#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 
    
=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::MCoffee

=cut

=head1 SYNOPSIS

my $db     = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $mcoffee = Bio::EnsEMBL::Compara::RunnableDB::Mcoffee->new ( 
                                                    -db      => $db,
                                                    -input_id   => $input_id,
                                                    -analysis   => $analysis );
$mcoffee->fetch_input(); #reads from DB
$mcoffee->run();
$mcoffee->output();
$mcoffee->write_output(); #writes to DB

=cut

=head1 DESCRIPTION

This Analysis/RunnableDB is designed to take a Family (or Homology) as input
Run a MCOFFEE multiple alignment on it, and store the resulting alignment
back into the family_member table.

input_id/parameters format eg: "{'family_id'=>1234,'options'=>'-maxiters 2'}"
    family_id       : use family_id to run multiple alignment on its members
    protein_tree_id : use 'id' to fetch a cluster from the ProteinTree
    options         : commandline options to pass to the 'mcoffee' program

=cut

=head1 CONTACT

  Contact Jessica Severin on module implemetation/design detail: jessica@ebi.ac.uk
  Contact Abel Ureta-Vidal on EnsEMBL/Compara: abel@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::RunnableDB::MCoffee;

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

use Bio::EnsEMBL::Hive;

use Bio::Greg::ProcessUtils;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);


#
# Some global-ish variables.
#

my $dba;
my $pta;

my $params;
my $tags;

my $tree;
my $sa;
my $sa_exons;

my $output_file;
my $sa_aligned;

my $dont_write_output;

sub debug {1;}

sub fetch_input {
    my($self) = @_;
    
    # Load up the Compara DBAdaptor.
    $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
    $pta = $dba->get_ProteinTreeAdaptor;
    
    ### DEFAULT PARAMETERS ###
    $params = {
      alignment_table        => 'protein_tree_member',
      alignment_score_table  => 'protein_tree_member_score',
      alignment_method       => 'cmcoffee',
      alignment_max_gene_count         => 10000,
      alignment_use_exon_boundaries    => 0,
      alignment_executable             => '',
      
      alignment_prank_codon_model => 0,
      alignment_prank_f           => 0
    };  
    
    #########################
    
    # Fetch parameters from the two possible locations. Input_id takes precedence!
    my $p_params = $self->get_params($self->parameters);
    my $i_params = $self->get_params($self->input_id);
    my $node_id = $i_params->{'protein_tree_id'};
    $node_id = $i_params->{'node_id'} if (!defined $node_id);
    my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags($dba,$node_id);

    $params = $self->replace_params($params,$p_params,$i_params,$t_params);
    Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);

    #########################

    $self->check_if_exit_cleanly;

    # Load the tree.
    if (defined $params->{'protein_tree_id'}) {
      $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($dba,$params);
    } elsif (defined $params->{'node_id'}) {
      $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis($dba,$params);
    } else {
	throw("No protein tree input ID!\n");
    }

    my $method = $params->{'alignment_method'};

    $dont_write_output = 0;
    if ($method eq 'none') {
      print "NO ALIGNMENT NECESSARY!! Skipping...\n";
      $dont_write_output = 1;
      return;
    }

    if ($method ne 'muscle') {
      if ($method eq 'cmcoffee') {
	my $num_leaves = scalar(@{$tree->get_all_leaves});
        # Switch to fmcoffee if we have > 300 genes.
	if ($num_leaves > 300) {
	  $params->{'alignment_method'} = 'fmcoffee';
	}
      }
      # Auto-switch to fmcoffee on single failure.
      if ($self->input_job->retry_count >= 1) {
	$params->{'method'} = 'fmcoffee';
      }
      # Auto-switch to muscle on a second failure.
      if ($self->input_job->retry_count >= 2) {
	$params->{'method'} = 'muscle';
      }
    }

    print "MCoffee alignment method: ".$params->{'alignment_method'}."\n";
    $tags->{'alignment_method'} = $params->{'alignment_method'};

    # Fail if the gene count is too big.
    my $num_leaves = scalar(@{$tree->get_all_leaves});
    if ($num_leaves > $params->{'alignment_max_gene_count'}) {
      #$self->dataflow_output_id($self->input_id, 2);
      $self->DESTROY;
      throw("Mcoffee job too big: try something else and FAIL it");
    }

    $sa = Bio::EnsEMBL::Compara::ComparaUtils->get_ProteinTree_seqs($tree,0);

    # Export exon-cased if necessary.
    my $use_exons = $params->{'alignment_use_exon_boundaries'};
    if ($use_exons) {
      $sa_exons = Bio::EnsEMBL::Compara::ComparaUtils->get_ProteinTree_seqs($tree,1);
    }
}

sub run
{
    my $self = shift;

    return if ($dont_write_output);    
    $self->check_if_exit_cleanly;

    my $method = $params->{'alignment_method'};
#    if ($method =~ 'coffee') {
#      $sa_aligned = $self->align_with_mcoffee($sa,$tree,$params);
#    } elsif ($method =~ 'prank') {
      $sa_aligned = $self->align_with_prank($sa,$tree,$params);
#    }
}


sub write_output {
    my $self = shift;
    $self->check_if_exit_cleanly;

    return if ($dont_write_output);
    
    if (defined $tree) {
      my $score_out = $output_file.".score_ascii";
      $self->parse_and_store_alignment_into_proteintree($tree, $sa_aligned, $score_out);
    }
    
    # Store our custom tags.
    Bio::EnsEMBL::Compara::ComparaUtils->store_tags($tree,$tags);
}


sub DESTROY {
    my $self = shift;
    
    if ($tree) {
	$tree->release_tree;
	$tree = undef;
    }

    unlink $output_file if (-e $output_file);
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

  $output_file = $tmp . "output";
  
  my $executable = $params->{'alignment_executable'} || 'prank';
  my $extra_params = '';
  $extra_params .= ' -codon ' if ($params->{'alignment_prank_codon_model'});
  $extra_params .= ' -F ' if ($params->{'alignment_prank_f'});
  
  my $cmd = qq^$executable -d=$aln_file -t=$tree_file -o=$output_file $extra_params^;
  
  $output_file .= '.1.fas';

  # Run the command.
  $dba->dbc->disconnect_when_inactive(1);
  my $rc = system($cmd);
  $dba->dbc->disconnect_when_inactive(0);

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
    
    $output_file = $tmp . "output_aln.fasta";
    $output_file =~ s/\/\//\//g;  # converts any // in path to /

    my $tree_temp = $tmp . "tree_temp.dnd";
    $tree_temp =~ s/\/\//\//g;  # converts any // in path to /
    
    my $mcoffee_executable;
    $mcoffee_executable = $params->{'executable'};
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
    }  else {
	
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
    $dba->dbc->disconnect_when_inactive(1) if ($dba);
    my $rc = system($prefix.$cmd);
    $dba->dbc->disconnect_when_inactive(0) if ($dba);
    
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
  my $mcoffee_scores_f = shift;
  
  my %align_hash;
  foreach my $seq ($aln->each_seq) {
      my $id = $seq->display_id;
      $align_hash{$id} = $seq->seq;
  }

  # Read in the scores file manually.
  my %score_hash;
  if ($mcoffee_scores_f && -e $mcoffee_scores_f) {
    my %overall_score_hash;
    my $FH = IO::File->new();
    $FH->open($mcoffee_scores_f) || throw("Could not open alignment scores file [$mcoffee_scores_f]");
    <$FH>; #skip header
    my $i=0;
    while(<$FH>) {
      $i++;
      next if ($i < 7); # skip first 7 lines.
      next if($_ =~ /^\s+/);  #skip lines that start with space
      if ($_ =~ /:/) {
        my ($id,$overall_score) = split(/:/,$_);
        $id =~ s/^\s+|\s+$//g;
        $overall_score =~ s/^\s+|\s+$//g;
        print "___".$id."___".$overall_score."___\n";
        next;
      }
      chomp;
      my ($id, $align) = split;
      $score_hash{$id} = '' if (!exists $score_hash{$id});
      $score_hash{$id} .= $align;
      print $id." ". $align."\n";
    }
    $FH->close;
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
  
  my $table_name = $params->{'alignment_table'};
  my $score_table = $params->{'alignment_score_table'};
  
  my $sth = $pta->prepare("INSERT INTO $table_name 
        (node_id,member_id,method_link_species_set_id,cigar_line)  VALUES (?,?,?,?)  ON DUPLICATE KEY UPDATE cigar_line=?");
  my $sth2 = $pta->prepare("INSERT INTO $score_table
        (node_id,member_id,method_link_species_set_id,cigar_line)  VALUES (?,?,?,?) ON DUPLICATE KEY UPDATE cigar_line=?");

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
      
      if (%score_hash) {
        my $score_string = $score_hash{$member->stable_id};
        $score_string =~ s/[^\d-]/9/g;   # Convert non-digits and non-dashes into 9s. This is necessary because t_coffee leaves some leftover letters.
        printf("Updating $score_table %.10s : %.30s ...\n",$member->stable_id,$score_string) if ($self->debug);
        $sth2->execute($member->node_id,$member->member_id,$member->method_link_species_set_id,$score_string,$score_string);
      }
      sleep(0.05);
  }

  $sth->finish;
  $sth2->finish;
}


1;
