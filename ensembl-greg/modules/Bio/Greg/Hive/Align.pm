package Bio::Greg::Hive::Align;

use strict;
use Getopt::Long;
use IO::File;
use File::Basename;
use File::Path;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::BaseAlignFeature;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::PeptideAlignFeatureAdaptor;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;
#use Time::HiRes qw(time gettimeofday tv_interval sleep);
use POSIX qw(ceil floor);

use base ('Bio::Greg::Hive::Process');

#
# Some global-ish variables.
#

sub debug {1;}

sub fetch_input {
    my($self) = @_;
        
    ### DEFAULT PARAMETERS ###
    my $params = {
      alignment_table        => 'protein_tree_member',
      alignment_score_table  => 'protein_tree_member_score',
      alignment_method       => 'cmcoffee',
      alignment_max_gene_count         => 10000,
      alignment_use_exon_boundaries    => 0,
      alignment_executable             => '',
      t_coffee_executable => '',
      
      alignment_prank_codon_model => 0,
      alignment_prank_f           => 0
    };  
    
    #########################

    $self->load_all_params($params);

    # Load the tree.
    my $out_table = $self->param('alignment_table');
    $self->param('alignment_table','protein_tree_member');
    my $tree = $self->get_tree();
    $self->param('alignment_table',$out_table);
    
    my $method = $self->param('alignment_method');

    foreach my $leaf ($tree->leaves) {
      if ($leaf->distance_to_parent > 20) {
	$leaf->distance_to_parent(4);
      }
    }

    if ($method eq 'cmcoffee') {
      my $num_leaves = scalar(@{$tree->get_all_leaves});
      # Switch to fmcoffee if we have > 300 genes.
      if ($num_leaves > 300) {
	$self->param('alignment_method','fmcoffee');
      }
    }
    # Auto-switch to fmcoffee on single failure.

#    my @prank_order = qw(prank_f cmcoffee fmcoffee muscle);
#    my @mcoffee_order = qw(cmcoffee fmcoffee muscle);

    my $retry_count = $self->input_job->retry_count;
    #if ($retry_count > 0) {
    #  my $index = $retry_count;
    #  my @array;
    #  if ($self->param('alignment_method') =~ 'prank') {
#	@array = @prank_order;
#      } elsif ($self->param('alignment_method') =~ 'mcoffee') {
#	@array = @mcoffee_order;
#      }
#      $index = scalar(@array) -1 if ($index > scalar(@array) - 1);
#      if (defined $array[$index] && $array[$index] ne '') {
#        $self->param('alignment_method',$array[$index]);
#      }
#    }

    print "Alignment method: ".$self->param('alignment_method')."\n";

    # Fail if the gene count is too big.
    my $num_leaves = scalar(@{$tree->get_all_leaves});
    if ($num_leaves > $self->param('alignment_max_gene_count')) {
      $self->DESTROY;
      throw("Mcoffee job too big: try something else and FAIL it");
    }

    my $sa = $tree->get_SimpleAlign;

    # Export exon-cased if necessary.
    my $use_exons = $self->param('alignment_use_exon_boundaries');
    if ($use_exons) {
      my $sa_exons = Bio::EnsEMBL::Compara::ComparaUtils->get_ProteinTree_seqs($tree,1);
      $self->param('sa_exons',$sa_exons);
    }

    $self->param('tree',$tree);
    $self->param('aln',$sa);
}

sub run
{
    my $self = shift;

    $self->check_if_exit_cleanly;

    my $sa = $self->param('aln');
    my $tree = $self->param('tree');
    my $params = $self->params;

    my $sa_aligned;

    my $method = $self->param('alignment_method');
    if ($method =~ '(coffee|muscle|clustalw)') {
      $sa_aligned = $self->align_with_mcoffee($sa,$tree,$params);
    } elsif ($method =~ m/prank/) {
      if ($method =~ m/_f/) {
        $params->{alignment_prank_f} = 1;
      }
      if ($method =~ m/_codon/) {
        $params->{alignment_prank_codon_model} = 1;
        my $sa_cds = $tree->get_SimpleAlign(-cdna => 1);
        my $cds_aligned = $self->align_with_prank($sa_cds,$tree,$params);
        $sa_aligned = Bio::EnsEMBL::Compara::AlignUtils->translate($cds_aligned);
      } else {
        $sa_aligned = $self->align_with_prank($sa,$tree,$params);
      }
    } elsif ($method =~ m/papaya/) {
      $sa_aligned = $self->align_with_papaya($sa,$tree,$params);
    } elsif ($method =~ m/none/) {
      $sa_aligned = $sa;
    }

    $self->param('sa_aligned',$sa_aligned);
}


sub write_output {
    my $self = shift;
    $self->check_if_exit_cleanly;

    return if ($self->param('dont_write_output'));
    
    my $tree = $self->param('tree');
    my $sa_aligned = $self->param('sa_aligned');

    if (defined $tree) {
      $self->parse_and_store_alignment_into_proteintree($tree, $sa_aligned);
    }
    
    # Store our custom tags.
    #Bio::EnsEMBL::Compara::ComparaUtils->store_tags($tree,$tags);
}


sub DESTROY {
    my $self = shift;
    
    $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


##########################################
#
# internal methods
#
##########################################

sub align_with_prank {
  my $self = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;

  my $tmp = $self->worker_temp_directory;
  
  # Output alignment.
  my $aln_file = $tmp . "aln.fasta";
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$aln_file); # Write the alignment out to file.
  
  my $tree_file = $tmp . "tree.nh";
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  Bio::EnsEMBL::Compara::TreeUtils->to_file($treeI,$tree_file);

  my $output_file = $tmp . "output";
  
  my $executable = $params->{'alignment_executable'} || 'prank';
  my $extra_params = '';
  $extra_params .= ' -codon ' if ($params->{'alignment_prank_codon_model'});
  $extra_params .= ' +F ' if ($params->{'alignment_prank_f'});
  
  my $cmd = qq^$executable $extra_params -d=$aln_file -t=$tree_file -o=$output_file^;
  
  $output_file .= '.1.fas';

  # Run the command.
  $self->compara_dba->dbc->disconnect_when_inactive(1);
  my $rc = system($cmd);
  $self->compara_dba->dbc->disconnect_when_inactive(0);

  unless($rc == 0) {
    print "Prank error!\n";
    die;
  }
  
  use Bio::AlignIO;
  my $alignio = Bio::AlignIO->new(-file => $output_file,
                                  -format => "fasta");
  my $aln = $alignio->next_aln();
  return $aln;
}

sub align_with_papaya {
  my $self = shift;
  my $aln = shift;
  my $tree = shift;
  my $params = shift;

  my $tmp = $self->worker_temp_directory;
  print "$tmp\n";
  
  # Output alignment.
  my $aln_file = $tmp . "aln.fasta";
  Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$aln_file); # Write the alignment out to file.
  
  my $tree_file = $tmp . "tree.nh";
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
  Bio::EnsEMBL::Compara::TreeUtils->to_file($treeI,$tree_file);

  my $output_file = $tmp . "output";
  
  my $executable = $params->{'alignment_executable'} || 'papaya';
  my $extra_params = '';
  
  my $cmd = qq^$executable $extra_params --seqfile $aln_file --treefile $tree_file --outfile $output_file^;

  # Run the command.
  $self->compara_dba->dbc->disconnect_when_inactive(1);
  my $rc = system($cmd);
  $self->compara_dba->dbc->disconnect_when_inactive(0);

  unless($rc == 0) {
    print "Papaya error!\n";
    die;
  }
  
  $output_file .= '.fas';

  use Bio::AlignIO;
  my $alignio = Bio::AlignIO->new(-file => $output_file,
                                  -format => "fasta");
  my $aln = $alignio->next_aln();
  return $aln;
}


sub align_with_mcoffee
{
    my $self = shift;
    my $aln = shift;
    my $tree = shift;
    my $params = shift;

    my $default_params = {
      use_tree => 1,
      use_exons => 0,
      exon_aln => undef
    };

    $params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($default_params,$params);

    my $tmp = $self->worker_temp_directory;
    #my $tmp = "/tmp/mcoffee/";
    $tmp = $params->{'temp_dir'} if (defined $params->{'temp_dir'});
    mkdir($tmp);

    # Output alignment.
    my $input_fasta = $tmp . "input_seqs.fasta";
    Bio::EnsEMBL::Compara::AlignUtils->to_file($aln,$input_fasta); # Write the alignment out to file.

    my $input_tree = $tmp . "input_tree.nh";
    my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
    Bio::EnsEMBL::Compara::TreeUtils->to_file($treeI,$input_tree);
    
    my $output_file = $tmp . "output_aln.fasta";
    $output_file =~ s/\/\//\//g;  # converts any // in path to /

    my $tree_temp = $tmp . "tree_temp.dnd";
    $tree_temp =~ s/\/\//\//g;  # converts any // in path to /
    
    my $mcoffee_executable;
    $mcoffee_executable = $params->{'executable'};
    $mcoffee_executable = $params->{'t_coffee_executable'} if (!defined $mcoffee_executable);
    unless (-e $mcoffee_executable) {
	print "Using default T-Coffee executable!\n";
	$mcoffee_executable = "/nfs/users/nfs_g/gj1/bin/t_coffee";
    }
    #throw("can't find a M-Coffee executable to run\n") unless(-e $mcoffee_executable);
    
    # Output the params file.
    my $paramsfile = $tmp . "mcoffee.params";
    $paramsfile =~ s/\/\//\//g;  # converts any // in path to /
    open(OUTPARAMS, ">$paramsfile")
	or throw("Error opening $paramsfile for write");
    
    # GJ 2008-11-04: Variable args depending on method choice.
    my $method_string = '-method=';
    my $method = $params->{'alignment_method'};

    $method = "cmcoffee" unless (defined $method);
    if ($method eq 'cmcoffee' || $method eq 'mcoffee') {
	
      # CMCoffee, slow, comprehensive multiple alignments.
      $method_string .= "mafftgins_msa, muscle_msa, kalign_msa, probcons_msa"; #, t_coffee_msa";
    } elsif ($method eq 'fmcoffee') {
	
	# FMCoffee, fast but accurate alignments.
	$method_string .= "mafft_msa, muscle_msa, clustalw_msa, kalign_msa";
    } elsif ($method eq 'muscle') {
	
	# MUSCLE: very quick, kind of crappy alignments.
	$method_string .= "muscle_msa";
    } elsif ($method eq 'prank') {
	
	# PRANK: phylogeny-aware alignment.
	$method_string .= "prank_msa";
    }  elsif ($method eq 'clustalw') {
	$method_string .= "clustalw_msa";
    }
    # GJ 2008-11-17: Use exon boundaries if desired.
    if ($params->{'use_exons'}) {
	$method_string .= ", exon_pair";
	my $exons_fasta = $tmp . "exon_aln.fasta";
	my $exon_aln = $params->{'exon_aln'};
	Bio::EnsEMBL::Compara::AlignUtils->to_file($exon_aln,$exons_fasta);
	print OUTPARAMS "-template_file=$exons_fasta\n";
    }
    $method_string .= "\n";
    
    print OUTPARAMS $method_string;
    print OUTPARAMS "-mode=mcoffee\n";
    print OUTPARAMS "-output=fasta_aln,score_ascii\n";
    print OUTPARAMS "-outfile=$output_file\n";
    print OUTPARAMS "-newtree=$tree_temp\n";
    print OUTPARAMS "-usetree=$input_tree\n" if ($params->{'use_tree'});
    close(OUTPARAMS);
    
    my $cmd = $mcoffee_executable;
    $cmd .= " ".$input_fasta;
    $cmd .= " -parameters=$paramsfile";
    
    # Output some environment variables for tcoffee
#    $tmp = substr($tmp,0,-1);
    my $prefix = "export HOME_4_TCOFFEE=\"$tmp\";";
    $prefix .= "export DIR_4_TCOFFEE=\"$tmp\";";
    $prefix .= "export TMP_4_TCOFFEE=\"$tmp\";";
    $prefix .= "export CACHE_4_TCOFFEE=\"$tmp\";";
    $prefix .= "export NO_ERROR_REPORT_4_TCOFFEE=1;";
    $prefix .= "export MAFFT_BINARIES=/nfs/users/nfs_g/gj1/bin/mafft-bins/binaries;";  # GJ 2008-11-04. What a hack!
    
    # Run the command.
    $self->compara_dba->dbc->disconnect_when_inactive(1);
    my $rc = system($prefix.$cmd);
    $self->compara_dba->dbc->disconnect_when_inactive(0);
    
    unless($rc == 0) {
#	$self->DESTROY;
	print "MCoffee error!\n";
        die;
	#throw("MCoffee job, error running executable: $\n");
    }

    use Bio::AlignIO;
    my $alignio = Bio::AlignIO->new(-file => $output_file,
				    -format => "fasta");
    my $aln = $alignio->next_aln();
    return $aln;
}

########################################################
#
# ProteinTree input/output section
#
########################################################

sub update_single_peptide_tree
{
  my $self   = shift;
  my $tree   = shift;
  
  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;

  foreach my $member (@{$tree->get_all_leaves}) {
      next unless($member->isa('Bio::EnsEMBL::Compara::AlignedMember'));
      next unless($member->sequence);
      $member->cigar_line(length($member->sequence)."M");
      $pta->store($member);
      printf("single_pepide_tree %s : %s\n", $member->stable_id, $member->cigar_line) if($self->debug);
  }
}


sub parse_and_store_alignment_into_proteintree
{
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  
  my $pta = $self->compara_dba->get_ProteinTreeAdaptor;

  my %align_hash;
  foreach my $seq ($aln->each_seq) {
      my $id = $seq->display_id;
      $align_hash{$id} = $seq->seq;
  }

  # Convert alignment strings into cigar_lines
  my $alignment_length;
  foreach my $id (keys %align_hash) {
    next if ($id eq 'cons');
    my $alignment_string = $align_hash{$id};
    unless (defined $alignment_length) {
      $alignment_length = length($alignment_string);
    } else {
      if ($alignment_length != length($alignment_string)) {
        throw("While parsing the alignment, some id did not return the expected alignment length\n");
      }
    }
    # Call the method to do the actual conversion
    print "$id  $alignment_string\n";
    $align_hash{$id} = Bio::EnsEMBL::Compara::ComparaUtils->cigar_line($alignment_string);
  }
  
  my $table_name = $self->param('alignment_table');
  
  my $sth = $pta->prepare("INSERT INTO $table_name 
        (node_id,member_id,method_link_species_set_id,cigar_line)  VALUES (?,?,?,?)  ON DUPLICATE KEY UPDATE cigar_line=?");

  # Align cigar_lines to members and store
  foreach my $member (@{$tree->get_all_leaves}) {
      if ($align_hash{$member->stable_id} eq "") {
	  throw("mcoffee produced an empty cigar_line for ".$member->stable_id."\n");
      }
      $member->cigar_line($align_hash{$member->stable_id});
      ## Check that the cigar length (Ms) matches the sequence length
      my @cigar_match_lengths = map { if ($_ eq '') {$_ = 1} else {$_ = $_;} } map { $_ =~ /^(\d*)/ } ( $member->cigar_line =~ /(\d*[M])/g );
      my $seq_cigar_length; map { $seq_cigar_length += $_ } @cigar_match_lengths;
      my $member_sequence = $member->sequence;
      $member_sequence =~ s/\*//g;
      if ($seq_cigar_length != length($member_sequence)) {
	print "MC  ".$member->cigar_line."\n";
	print "MS  ".$member_sequence."\n";
	  throw("While storing the cigar line, the returned cigar length did not match the sequence length\n");
      }
            
      if ($table_name eq 'protein_tree_member') {
	  # We can use the default store method for the $member.
	  $pta->store($member);
      } else {
	  # Do a manual insert into the correct output table.
          my $cmd = "CREATE TABLE IF NOT EXISTS $table_name LIKE protein_tree_member;";
	  $pta->dbc->do($cmd);
	  printf("Updating $table_name %.10s : %.30s\n",$member->stable_id,$member->cigar_line) if ($self->debug);
	  $sth->execute($member->node_id,$member->member_id,$member->method_link_species_set_id,$member->cigar_line,$member->cigar_line);
      }
	  # Do a manual insert of the *scores* into the correct score output table.
      
      sleep(0.05);
  }

  $sth->finish;
}


1;
