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

  $tree = $pta->fetch_node_by_node_id($node_id);
  throw("No protein tree!") unless (defined $tree);
}

sub run {
  my $self = shift;

  $self->check_if_exit_cleanly;
  $self->{'start_time'} = time() * 1000;

    # Select all codons. 
  my $table = $params->{'sitewise_table'};
  my $node_id = $tree->node_id;

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
  $mapped_omegas = Bio::Greg::EslrUtils->mapSitewiseToGenome($tree,9606,$pos_values);
  print "  -> Finished mapping values!\n";

  # Collect various gene-centric tags.
  $gene_tags = Bio::Greg::EslrUtils->collectGeneTags($tree);
  print "  -> Finished collecting tags!\n";

  #$self->create_plot();
  #print "  -> Finished plotting omegas!\n";
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

  #$dba->dbc->do("LOCK TABLE sitewise_genome WRITE;");
  my $insert_cmd = qq^
REPLACE INTO sitewise_genome 
  (node_id,aln_position,member_id,chr_name,chr_start,chr_end,residue)
    VALUES (?,?,?,?,?,?,?);^;
  my $sth = $dba->dbc->prepare($insert_cmd);

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
  #$dba->dbc->do("UNLOCK TABLES;");

  #$dba->dbc->do("LOCK TABLE protein_tree_tag WRITE;");
  Bio::EnsEMBL::Compara::ComparaUtils->store_tags($tree,$gene_tags);
  #$dba->dbc->do("UNLOCK TABLES;");

  $tree->release_tree;
}


1;
