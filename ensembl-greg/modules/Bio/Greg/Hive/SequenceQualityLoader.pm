package Bio::Greg::Hive::SequenceQualityLoader;

use strict;
use Time::HiRes qw(time gettimeofday tv_interval sleep);
use Getopt::Long;
use IO::File;
use File::Basename;
use File::Path;

use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::ComparaUtils;

use base ('Bio::Greg::Hive::Process');

sub fetch_input {
  my $self = shift;

  ### DEFAULT PARAMETERS ###
  my $params = {
    alignment_table       => 'protein_tree_member',
    alignment_score_table => 'protein_tree_member_score'
  };

  $self->load_all_params($params);
}

sub run {
  my $self = shift;

  $self->{'start_time'} = time() * 1000;

  my $tree = $self->get_tree;
  my $seq_id_quals_hash;
  
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

  $self->param('seq_id_quals_hash',$seq_id_quals_hash);
}

sub write_output {
  my $self = shift;

  my $seq_id_quals_hash = $self->param('seq_id_quals_hash');

  foreach my $key ( sort keys %{$seq_id_quals_hash} ) {
    my $quals = $seq_id_quals_hash->{$key};

    # Now we can store it in the sequence_qual table!
    my $cmd = "REPLACE INTO sequence_quality (sequence_id,length,sequence) VALUES (?,?,?);";
    print $cmd. "  " . $key . "\n";
    my $sth = $self->compara_dba->dbc->prepare($cmd);
    $sth->execute( $key, length($quals), $quals );
    sleep(0.4);
  }
}

1;
