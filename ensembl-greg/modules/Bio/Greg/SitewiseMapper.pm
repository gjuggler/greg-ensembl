package Bio::Greg::SitewiseMapper;

use strict;
use Getopt::Long;
use IO::File;
use File::Basename;
use File::Path;

use Time::HiRes qw(sleep);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;

use Bio::EnsEMBL::Hive;
use Bio::EnsEMBL::Hive::Process;

our @ISA = qw(Bio::EnsEMBL::Hive::Process);

my $dba;
my $pta;

my $params;

my $tree;
my $pos_values;
my $mapped_omegas;
my $gene_tags;

sub fetch_input {
  my $self = shift;

  # Load up the Compara DBAdaptor.
  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;

  $pta->local_mode(-1);

  ### DEFAULT PARAMETERS ###
  $params->{'sitewise_table'} = 'sitewise_aln';
  $params->{'do_mapping'} = 1;
  $params->{'collect_tags'} = 1;
  $params->{'collect_pfam'} = 1;
  $params->{'collect_go'} = 1;
  $params->{'create_plot'} = 0;
  $params->{'parameter_set_id'} = 1;
  #########################
  
  # Fetch parameters from the two possible locations. Input_id takes precedence!
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);

  #########################
  #
  # Load the tree.
  #
  my $node_id;
  $node_id = $params->{'protein_tree_id'};
  $node_id = $params->{'node_id'} if (!defined $node_id);

  my $param_set_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set($dba->dbc,$params->{'parameter_set_id'});
  my $new_params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($params,$param_set_params);

  $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($dba,$new_params);
  $tree = $tree->minimize_tree;
}

sub run {
  my $self = shift;

  $self->check_if_exit_cleanly;
  $self->{'start_time'} = time() * 1000;

  # Select all codons. 
  my $table = $params->{'sitewise_table'};
  my $node_id = $tree->node_id;

  if ($params->{'do_mapping'}) {
    print "Mapping sitewise to genome...\n";
    $self->do_mapping();
    print "  -> Finished mapping values!\n";
  }

  if ($params->{'collect_tags'}) {
    print "Collecting gene tags...\n";
    $self->collect_gene_tags();
    print "  -> Finished collecting tags!\n";
  }

  if ($params->{'collect_pfam'}) {
    print "Collecting Pfam annotations...\n";
    $self->collect_pfam();
    print "  -> Finished collecting Pfam!";
  }

  if ($params->{'collect_go'}) {
    print "Collecting GO annotations...\n";
    $self->collect_go();
    print "  -> Finished collecting GO terms!";
  }

  if ($params->{'create_plot'}) {
    $self->create_plot();
    print "  -> Finished plotting omegas!\n";
  }
}

sub do_mapping {
  my $self = shift;

  print "Mapping sitewise $node_id  to genome...\n";
    
  my $omega_cmd = qq^
      SELECT distinct(aln_position) aln_position FROM $table WHERE
      ncod >= 4
      AND note != "random"
      AND omega_upper > omega
      AND node_id=$node_id
      ;
    ^;
  my $sth = $dba->dbc->prepare($omega_cmd);
  $sth->execute();
  my @hashrefs = @{$sth->fetchall_arrayref({})};
  map {$pos_values->{$_->{'aln_position'}} = 1} @hashrefs;
  
  # Run the genomic mapping code.
  use Bio::Greg::EslrUtils;
  my $mapped_omegas = Bio::Greg::EslrUtils->mapSitewiseToGenome($tree,9606,$pos_values);

  #$dba->dbc->do("LOCK TABLE sitewise_genome WRITE;");
  my $insert_cmd = qq^
    REPLACE INTO sitewise_genome 
    (node_id,aln_position,member_id,chr_name,chr_start,chr_end,residue)
    VALUES (?,?,?,?,?,?,?);^;
  $sth = $dba->dbc->prepare($insert_cmd);

  foreach my $map (@{$mapped_omegas}) {
    $sth->execute($map->{'node_id'},
                  $map->{'aln_position'},
                  $map->{'member_id'},
                  $map->{'chr'},
                  $map->{'start'},
                  $map->{'end'},
		  $map->{'char'}
      );
    print $map->{'chr'}." ".$map->{'start'}."\n";
    sleep(0.05);
  }  
  $sth->finish();
}

sub collect_gene_tags {
  my $self = shift;

  #$gene_tags = Bio::Greg::EslrUtils->collectGeneTags($tree,$params);
  $gene_tags = Bio::Greg::EslrUtils->collectDuplicationTags($tree,$params);
}

sub collect_go_terms
{
    my $self = shift;

    my $tree = $pta->fetch_node_by_node_id($node_id);
    my @leaves = @{$tree->get_all_leaves};

    # GJ 2009-02-09 : priority for grabbing GO annotations. Human > mouse > zebrafish > drosophila > c.elegans
#    my @species_leaves = grep {$_->taxon_id == 9606 || $_->taxon_id == 10090 ||
#			       $_->taxon_id == 7955 || $_->taxon_id == 7227 ||
#			   $_->taxon_id == 6239} @leaves;
#    @leaves = @species_leaves if (scalar(@species_leaves) >= 1);
    my @leaves = grep {$_->taxon_id==9606} @leaves;

    my %go_terms;
    my $taxon_id;

    foreach my $leaf (@leaves)
    {
	$taxon_id = $leaf->taxon_id;
	my $ts = $leaf->transcript;
	next if (!defined $ts);

	my $db_entries = $ts->get_all_DBLinks;

	# Grep out all the GO xref entries.
	my @keepers = grep {$_->dbname eq "GO"} @{$db_entries};
	foreach my $db_e (@keepers)
	{
	    $self->insert_go_term($node_id,$leaf,$db_e);
	    sleep(0.1);
	}
    }
    
    my @gos = keys(%go_terms);
    return \@gos;
}

sub insert_go_term
{
    my $self = shift;
    my ($node_id,$leaf,$db_e) = @_;

    my $cmd = "INSERT IGNORE INTO go_terms (node_id,member_id,stable_id,go_term) values (?,?,?,?);";
    print $cmd."\n";
    my $sth = $dba->dbc->prepare($cmd);
    $sth->execute($node_id,$leaf->dbID,$leaf->stable_id,$db_e->display_id);
}


sub collect_pfam {
  my $self = shift;

  my $sa = $tree->get_SimpleAlign;
  my $pos_id_hash;
  foreach my $leaf ($tree->leaves) {
    print $leaf->stable_id."\n";
    my $name = $leaf->stable_id;
    my $tx = $leaf->get_Translation;
    my @features = @{$tx->get_all_ProteinFeatures('PFam')};
    foreach my $f (@features) {
      my $pf_id = $f->display_id;
      my $pf_lo = $f->hstart; # Start and end in the hit (Pfam domain) coordinates.
      my $pf_hi = $f->hend;
      my $lo = $f->start;
      my $hi = $f->end;
      foreach my $i (0 .. ($hi-$lo)) {
	my $pos = $lo + $i;
	my $aln_col = $sa->column_from_residue_number($name,$pos);
	my $obj = {
	  pf_pos => $pf_lo + $i,
	  score => $f->score
	};
	$pos_id_hash->{$pf_id."_".$aln_col} = $obj;
      }
    }
  }

  my $tree_node_id=0;
  $tree_node_id = $tree->subroot->node_id if ($tree->subroot);

  my $cmd = "REPLACE INTO sitewise_pfam (node_id,aln_position,pf_position,tree_node_id,pfam_id,score) VALUES (?,?,?,?,?,?)";
  my $sth = $tree->adaptor->prepare($cmd);

  use Time::HiRes qw(sleep);
  foreach my $key (sort keys %{$pos_id_hash}) {
    my $id;
    my $pos;
    ($id,$pos) = split("_",$key);
    my $obj = $pos_id_hash->{$key};
    my $pf_pos = $obj->{'pf_pos'};
    my $score = $obj->{'score'};
    printf("%s %s %s %s\n",
	   $id,
	   $pos,
	   $pf_pos,
	   $score);
    $sth->execute($tree->node_id,
		  $pos,
		  $pf_pos,
		  $tree_node_id,
		  $id,
		  $score);
    sleep(0.1);
  }

  $sth->finish;
}


sub create_plot {
  my $self = shift;

  my $base_dir = "/lustre/scratch103/ensembl/gj1/2xmammals_plots";

  my $node_id = $tree->node_id;

  use Digest::MD5 qw(md5_hex);
  my $digest = md5_hex($node_id);

  my $small_digest = substr($digest,0,10);

  my $subdir = substr($digest,0,1);
  my $sub_subdir = substr($digest,1,1);

  print "$digest\n";
  print "DIRS: $subdir $sub_subdir\n";
  my $output_dir = $base_dir."/".$subdir."/".$sub_subdir;
  use File::Path qw(mkpath);
  mkpath($output_dir,{mode => 0777});

  my @human_peps = grep {$_->taxon_id == 9606} $tree->leaves;
  foreach my $hum_pep (@human_peps) {
    my $stable_id = $hum_pep->stable_id;
    my $output_file = $output_dir."/${stable_id}.pdf";
    
    print "PLOT OUTPUT: $output_file\n";
    
    my $params = {
      parameter_set_id => 14,
      remove_blank_columns => 1,
      mask_outside_subtree => 1,
      remove_subtree => 0,
      mask_gblocks => 0,
      sitewise_table => 'sitewise_aln'
      };
    use Bio::Greg::EslrPlots;
    my $fresh_tree = $pta->fetch_node_by_node_id($node_id);
    Bio::Greg::EslrPlots->plotTreeWithOmegas($output_file,$params,$fresh_tree);
    $fresh_tree->release_tree;
  }

}

sub write_output {
  my $self = shift;

  $tree->release_tree;
}


1;
