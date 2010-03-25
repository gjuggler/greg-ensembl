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

use Bio::Greg::ProcessUtils;

our @ISA = qw(Bio::EnsEMBL::Hive::Process Bio::Greg::ProcessUtils);

my $dba;
my $pta;

my $params;
my $tags;

my $tree;
my $seq_id_quals_hash;

sub fetch_input {
  my $self = shift;

  # Load up the Compara DBAdaptor.
  $dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new( -DBCONN => $self->db->dbc );
  $pta = $dba->get_ProteinTreeAdaptor;

  ### DEFAULT PARAMETERS ###
  $params = {
    alignment_table       => 'protein_tree_member',
    alignment_score_table => 'protein_tree_member_score'
  };

  #########################

  # Fetch parameters from the possible locations.
  my $p_params = $self->get_params( $self->parameters );
  my $i_params = $self->get_params( $self->input_id );
  my $node_id  = $i_params->{'protein_tree_id'};
  $node_id = $i_params->{'node_id'} if ( !defined $node_id );
  my $t_params = Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_tree_tags( $dba, $node_id );

  $params = $self->replace_params( $params, $p_params, $i_params, $t_params );
  Bio::EnsEMBL::Compara::ComparaUtils->hash_print($params);

  #########################

  # Load the tree.
  $params->{'alignment_table'} = 'protein_tree_member';
  if ( defined $params->{'node_id'} ) {
    $tree = Bio::EnsEMBL::Compara::ComparaUtils->get_tree_for_comparative_analysis( $dba, $params );
  } else {
    throw("No protein tree input ID!\n");
  }

  throw("No protein tree!") unless ( defined $tree );
}

sub run {
  my $self = shift;
  $self->check_if_exit_cleanly;
  $self->{'start_time'} = time() * 1000;

  foreach my $leaf ( $tree->leaves ) {
    my $member = $leaf;
    print $member->stable_id . "\n";

    #    next unless ($member->stable_id =~ m/(pvap|stop)/ig);

    my $gdb = $member->genome_db;

# This is finicky: we need to call the "db_adaptor" method to get the Bio::EnsEMBL::DBSQL::DBAdaptor object, and then the meta container.
    my $meta     = $gdb->db_adaptor->get_MetaContainer;
    my $coverage = @{ $meta->list_value_by_key('assembly.coverage_depth') }[0];
    if ( $coverage ne 'low' ) {
      print $member->stable_id . "\n";
      next;
    } else {
      print $member->stable_id . " 2x!\n";
    }

    my $quals = Bio::EnsEMBL::Compara::ComparaUtils->get_quality_string_for_member($member);
    if ( defined $quals ) {
      my $pep = $member->sequence;
      print $quals. "\n";
      print $pep. "\n";
      $quals = substr( $quals, 0, -1 ) if ( length($quals) > length($pep) );
      if ( length($quals) == length($pep) + 1 ) {
        $quals = substr( $quals, 0, -1 );
      }
      die if ( length($quals) != length($pep) );

      printf "%s\n%s\n", $pep, $quals;
      $seq_id_quals_hash->{ $member->sequence_id } = $quals;
    }
  }
  $self->{'end_time'} = time() * 1000;
}

sub write_output {
  my $self = shift;

  foreach my $key ( sort keys %{$seq_id_quals_hash} ) {
    my $quals = $seq_id_quals_hash->{$key};

    # Now we can store it in the sequence_qual table!
    my $cmd = "REPLACE INTO sequence_quality (sequence_id,length,sequence) VALUES (?,?,?);";
    print $cmd. "  " . $key . "\n";
    my $sth = $dba->dbc->prepare($cmd);
    $sth->execute( $key, length($quals), $quals );
    sleep(0.4);
  }
}

1;
