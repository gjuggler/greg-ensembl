package Bio::Greg::PrimateHIV::CollectPrimateHIVStats;

use strict;
use Time::HiRes qw(sleep);

use Cwd;
use Bio::AlignIO;

use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::Process;

use base ('Bio::Greg::Hive::CollectSitewiseStats', 'Bio::Greg::StatsCollectionUtils');

sub get_table_structure {
  my $self = shift;

  my $structure = {
    data_id   => 'int',
    gene_name => 'string',

    peptide_stable_id  => 'string',
    peptide_window_start => 'int',
    peptide_window_end => 'int',
    peptide_window_width => 'int',

    aln_window_start => 'int',
    aln_window_end => 'int',

    'hg19_chr_name'     => 'string',
    'hg19_window_start' => 'int',
    'hg19_window_end'   => 'int',

    'hg18_chr_name'     => 'string',
    'hg18_window_start' => 'int',
    'hg18_window_end'   => 'int',

    dnds_primates => 'float',
    dnds_hominids => 'float',

    pval_primates => 'float',
    pval_hominids => 'float',
    n_leaves_primates => 'int',
    n_sites_primates    => 'int',
    n_pos_sites_primates    => 'int',
    n_leaves_hominids => 'int',
    n_sites_hominids    => 'int',
    n_pos_sites_hominids    => 'int',
    unique_keys => 'data_id,peptide_window_start,peptide_window_end',
  };

  return $structure;
}

sub fetch_input {
  my ($self) = @_;

  my $params = {
    table             => 'stats_windows',
    reference_species => 9606,
    window_size => 10,
    window_step => 5
  };

  $self->load_all_params($params);

  # Create tables if necessary.
  $self->create_table_from_params( $self->hive_dba, $self->param('table'),
    $self->get_table_structure );
}

sub run {
  my $self = shift;

  foreach my $size (10, 30, 50, 100, 9999) {
    $self->run_with_windows($size,$size/2);
  }

}

sub gene_dnds_for_parameter_set {
  my $self = shift;
  my $parameter_set_id = shift;

  my $sth = $self->dbc->prepare("SELECT * FROM dnds_genes where parameter_set_id=? and node_id=?");
  $sth->execute($parameter_set_id,$self->param('node_id'));
  my $obj = $sth->fetchrow_hashref;
  $sth->finish;
  return $obj->{dnds};
}

sub run_with_windows {
  my $self = shift;
  my $w_size = shift;
  my $w_step = shift;

  my $params = $self->params;
  my $cur_params = $self->replace( $params, {} );

  my $ref_species = $self->param('reference_species');

  my $aln_aa = $self->get_aln;
  $self->pretty_print($aln_aa,{full => 1,length => 100});
  my $tree = $self->get_tree;

  my @members = $tree->leaves;
  my ($ref_member) = grep { $_->taxon_id == $self->param('reference_species') } @members;

  print "REF: ".$ref_member->stable_id."\n";

#  return unless (defined $ref_member);

  my @seqs = $aln_aa->each_seq;
  my ($ref_seq) = grep { $_->id eq $ref_member->stable_id } @seqs;
  my $ref_tx = $ref_member->get_Transcript;

  my $c_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(-url => 'mysql://ensro@ens-livemirror:3306/ensembl_compara_58');
  my ($cdna,$aa) = Bio::EnsEMBL::Compara::ComparaUtils->genomic_aln_for_member($c_dba,$ref_member);
  $self->pretty_print($aa,{full => 1, length => 100});

  my $len = $ref_member->seq_length;

  my $hominid_sites = $self->get_psc_hash($self->dbc,$self->params);
  my $hominid_tree = $self->get_tree($self->params);
  my $hominid_dnds = $self->gene_dnds_for_parameter_set(1);
  my $cur_params = $self->params;
  # Set to primate pset ID
  $cur_params->{parameter_set_id} = 2;
  $self->param('parameter_set_id',2);
  my $primate_sites = $self->get_psc_hash($self->dbc,$cur_params);
  my $primate_tree = $self->get_tree($cur_params);
  my $primate_dnds = $self->gene_dnds_for_parameter_set(2);

  # Back to normal pset ID
  $self->param('parameter_set_id',1);

  my $no_windows_yet = 1;
  for (my $i=1; $i < $len; $i += $w_step) {
    my $lo = $i;
    my $hi = $i + $w_size;
    $hi = $len if ($hi > $len);

    last if (!$no_windows_yet && $hi - $lo < $w_size - $w_step);
    $no_windows_yet = 0;
    printf ">>>> PEPTIDE WINDOW: %d %d\n",$lo,$hi;

    my $lo_coords = $self->get_coords_from_pep_position($ref_member,$lo);
    my $hi_coords = $self->get_coords_from_pep_position($ref_member,$hi);
    my $aln_coord_lo = $aln_aa->column_from_residue_number($ref_seq->id,$lo);
    my $aln_coord_hi = $aln_aa->column_from_residue_number($ref_seq->id,$hi);

    my $cur_params = $self->replace($params,{
      peptide_stable_id => $ref_tx->stable_id,
      peptide_window_start => $lo,
      peptide_window_end => $hi,
      peptide_window_width => $w_size,
      aln_window_start => $aln_coord_lo,
      aln_window_end => $aln_coord_hi,
      hg19_chr_name => $lo_coords->{hg19_name},
      hg18_chr_name => $lo_coords->{hg18_name},
      hg19_window_start => $lo_coords->{hg19_pos},
      hg18_window_start => $lo_coords->{hg18_pos},
      hg19_window_end => $hi_coords->{hg19_pos},
      hg18_window_end => $hi_coords->{hg18_pos},
      gene_name => $ref_member->get_Gene->external_name,
      dnds_primates => $primate_dnds,
      dnds_hominids => $hominid_dnds
                                    });

    foreach my $sites ($hominid_sites,$primate_sites) {      
      my @window_sites = map {$sites->{$_}} keys %$sites;
      @window_sites = grep {
        my $pos = $_->{aln_position};
        ($pos >= $aln_coord_lo && $pos <= $aln_coord_hi);
      } @window_sites;
      
      my $pval = $self->combined_pval(\@window_sites,'fisher');
      
      my @pos_sites = grep {$_->{omega} > 1} @window_sites;

      my $prefix;
      my $tree;
      $prefix = 'hominids' if ($sites == $hominid_sites);
      $prefix = 'primates' if ($sites == $primate_sites);
      $tree = $primate_tree if ($sites == $primate_sites);
      $tree = $hominid_tree if ($sites == $hominid_sites);
      my $added_params = {
        'n_leaves_'.$prefix => scalar($tree->leaves),
        'pval_'.$prefix => $pval,
        'n_sites_'.$prefix => scalar(@window_sites),
        'n_pos_sites_'.$prefix => scalar(@pos_sites),
      };
      $cur_params = $self->replace($cur_params,$added_params);
    }

    $cur_params->{data_id} = $cur_params->{node_id};
    $self->hash_print($cur_params);
    $self->store_params_in_table($self->dbc,$self->param('table'),$cur_params);

  }

}

sub get_coords_from_pep_position {
  my $self   = shift;
  my $member = shift;
  my $pos = shift;

  my $ref_tx = $member->get_Transcript;

  # Get the dba for the reference species.
  my $alias = $member->taxon->ensembl_alias;
  my $dba   = Bio::EnsEMBL::Registry->get_DBAdaptor( $alias, 'core' );
  my $asma  = $dba->get_AssemblyMapperAdaptor;
  my $csa   = $dba->get_CoordSystemAdaptor;

  my $new_cs = $csa->fetch_by_name('chromosome');
  my $old_cs = $csa->fetch_by_name( 'chromosome', 'NCBI36' );

  my $asm_mapper = $asma->fetch_by_CoordSystems( $new_cs, $old_cs );

  my @genomic1 = $ref_tx->pep2genomic( $pos, $pos );
  if (@genomic1) {
    my $coord1 = $genomic1[0];
    my $chr    = $ref_tx->seq_region_name;
    my $start  = $coord1->start;
    my $end    = $coord1->end;
    my $strand = $coord1->strand;

    # Map to old coordinates.
    my @coords = $asm_mapper->map( $chr, $start, $end, $strand, $new_cs );
    my ($coord) = @coords;
    if ($coord) {
      return {
        'hg19_name' => $chr,
        'hg18_name' => $chr,
        'hg19_pos'  => $start,
        'hg18_pos'  => $coord->start,
      };
    }

  }
}

1;
