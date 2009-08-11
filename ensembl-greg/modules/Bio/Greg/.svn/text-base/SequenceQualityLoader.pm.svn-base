package Bio::Greg::SequenceQualityLoader;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval sleep);
use Getopt::Long;
use IO::File;
use File::Basename;
use File::Path;

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
my $tags;

my $tree;
my $seq_id_quals_hash;

sub fetch_input {
  my $self = shift;

  # Load up the Compara DBAdaptor.
  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-DBCONN=>$self->db->dbc);
  $pta = $dba->get_ProteinTreeAdaptor;

  ### DEFAULT PARAMETERS ###
  my $p;
  
  $params->{'input_table_base'} = 'protein_tree';
  $params->{'output_table_base'} = 'protein_tree';
  #########################
  
  # Fetch parameters from the two possible locations. Input_id takes precedence!
  # (this utility method is from AlignmentProcess.pm)
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->parameters);
  $params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_string($params,$self->input_id);

  #########################
  #
  # Load the tree.
  #
  my $node_id;
  $node_id = $params->{'protein_tree_id'};
  $node_id = $params->{'node_id'} if (!defined $node_id);

  $pta->table_base($params->{'input_table_base'});
  $tree = $pta->fetch_node_by_node_id($node_id);

  throw("No protein tree!") unless (defined $tree);

  #
  # Some last-minute adjustments based on retry counts or somesuch.
  #

  #
  # Think of reasons why we want to fail the job.
  #

}

my @twox_ids = (9978,9371,9739,9478,42254,30538,10020,9365,59463,9358,9813,37347,9593,132908,30611,30608,9785,43179,9986,9685,9361);

sub run {
  my $self = shift;
  $self->check_if_exit_cleanly;
  $self->{'start_time'} = time()*1000;

  foreach my $leaf ($tree->leaves) {
    my $member = $leaf;
    
    # Only try with leaves that are 2x genomes.
    if (!grep {$_ eq $leaf->taxon_id} @twox_ids) {
      print $member->stable_id."\n";
      next;
    } else {
      print $member->stable_id." 2x!\n";
    }

    my $quals = Bio::EnsEMBL::Compara::ComparaUtils->get_quality_string_for_member($member);
    if (defined $quals) {
      my $pep = $member->sequence;
      #print $quals."\n";
      #print $pep."\n";
      #$quals = substr($quals,0,-1) if (length($quals) > length($pep));
      if (length($quals) == length($pep)+1) {
	$quals = substr($quals,0,-1);
      }
      die if (length($quals) != length($pep));
      
      printf "%s\n%s\n",$pep,$quals;
      $seq_id_quals_hash->{$member->sequence_id} = $quals;
    }
  }
  $self->{'end_time'} = time()*1000;
}


sub write_output {
  my $self = shift;

  foreach my $key (sort keys %{$seq_id_quals_hash}) {
      my $quals = $seq_id_quals_hash->{$key};
      
      # Now we can store it in the sequence_qual table!
      my $cmd = "REPLACE INTO sequence_quality (sequence_id,length,sequence) VALUES (?,?,?);";
      print $cmd."  ".$key."\n";
      my $sth = $dba->dbc->prepare($cmd);
      $sth->execute($key,length($quals),$quals);
      sleep(0.4);
  }
}


1;
