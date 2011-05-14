package Bio::Greg::Hive::CountOrthologs;

use strict;
use warnings;

use base (
  'Bio::Greg::Hive::Process',
);

sub param_defaults {
    return {

 };
}

sub fetch_input {
  my $self = shift;

  $self->{_genes_structure} = $self->table_def;

  $self->store_param('job_id', $self->job_id);

}

sub run {
  my $self = shift;

  my $gene_id = $self->param('gene_id');
  my $mba = $self->compara_dba->get_MemberAdaptor;
  my $member = $mba->fetch_by_source_stable_id(undef, $gene_id);
  my $gene = $member->get_Gene;
  my $tx_member = $member->get_canonical_peptide_Member;
  if ($tx_member) {
    $self->store_param('protein_id', $tx_member->stable_id);
    $tx_member->add_tag('orth_type', 'ortholog_one2one');
    $self->store_param('o_'.9606, $member->stable_id);
  }
  $self->store_param('gene_id', $gene->stable_id);
  $self->store_param('gene_name', $gene->external_name);
  $self->store_param('data_id', $member->dbID);


  my $run_quickly = 1;
  if ($run_quickly == 1) {
    my $homology_adaptor = $self->compara_dba->get_HomologyAdaptor;
    my $homs = $homology_adaptor->fetch_all_by_Member($member);
    my $taxon_counts;
    foreach my $hom (@$homs) {
      my $desc = $hom->description;
      next unless ($desc =~ m/(one2one|one2many)/i );
      my @members = @{$hom->gene_list};
      foreach my $m (@members) {
        my $tx_id = $m->taxon_id;
        next if ($tx_id == $member->taxon_id);
        
        $taxon_counts->{$tx_id} = 0 if (!defined $taxon_counts->{$tx_id});
        $taxon_counts->{$tx_id}++;
        #print $tx_id."  ".$desc."\n";
        if ($desc =~ m/one2one/) {
          $self->store_param('o_'.$tx_id, $m->stable_id, 'string');
          delete $taxon_counts->{$tx_id};
        }
      }
    }

    foreach my $tx_id (keys %$taxon_counts) {
      my $count = $taxon_counts->{$tx_id};
      if ($count > 1) {
        $self->store_param('o_'.$tx_id, $count, 'string');
      } else {
        # We've got a one2many ortholog, but only one non-human member.
        # This means it's a human duplication.
        $self->store_param('o_'.$tx_id, 'human_dup', 'string');
      }
    }
  } else {    
    my $ortholog_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_one_to_one_ortholog_tree($self->compara_dba, $member, '(ortholog_one2one|ortholog_one2many)');  
    if (defined $ortholog_tree) {
      # Count how many we get from each taxa.
      my $taxon_counts;
      foreach my $leaf ($ortholog_tree->leaves) {
        next if ($leaf->taxon_id == $member->taxon_id);
        my $tx_id = $leaf->taxon_id;
        $taxon_counts->{$tx_id} = 0 if (!defined $taxon_counts->{$tx_id});
        $taxon_counts->{$tx_id}++;
      }
      
      foreach my $leaf ($ortholog_tree->leaves) {
        next if ($leaf->taxon_id == $member->taxon_id);
        my $tx_id = $leaf->taxon_id;
        if ($taxon_counts->{$tx_id} == 1) {
          $self->_store_member($leaf);
        } else {
          $self->store_param('o_'.$tx_id, $taxon_counts->{$tx_id}.'');
        }
      }
      $self->store_param('one2many_tree', $ortholog_tree->newick_format);
    }
    
    my $one2one_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_one_to_one_ortholog_tree($self->compara_dba, $member, 'ortholog_one2one');  
    if (defined $one2one_tree) {
      $self->store_param('one2one_full', $one2one_tree->newick_format);
      
      my $output_base = $self->get_output_folder;
      my $scratch_file = "$output_base/genome_tree.nh";
      if (!-e $scratch_file) {
        my $genome_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree($self->compara_dba);
        map {
          my $name = $_->name;
          if (defined $_->taxon_id) {
            $_->{_tags} = {
              name => $name,
              taxon_id => $_->ncbi_taxid
            };
          } else{
            $_->{_tags} = {
              name => $name
            };
          }
        } $genome_tree->nodes;
        print $genome_tree->nhx_format."\n";
        Bio::EnsEMBL::Compara::TreeUtils->to_file($genome_tree, $scratch_file, {nhx_format => 1});
      }
      my $genome_tree = Bio::EnsEMBL::Compara::TreeUtils->from_file($scratch_file);
      print $genome_tree->newick_format."\n";
      
      my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($one2one_tree);
      my $mammal_tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_clade($self->compara_dba, $tree_copy, 'Mammalia', $genome_tree);
      $self->store_param('one2one_mammals', $mammal_tree->newick_format);
      
      $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($one2one_tree);
      my $primate_tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_clade($self->compara_dba, $tree_copy, 'Primates', $genome_tree);
      $self->store_param('one2one_primates', $primate_tree->newick_format);
    }
  }
}

sub _store_member {
  my $self = shift;
  my $member = shift;

  my $taxon_id = $member->taxon_id;
  my $type = $member->get_tagvalue('orth_type');
  
  if ($type eq 'ortholog_one2many') {
    $self->store_param('o_'.$taxon_id, 'one2many');
  } elsif ($type eq 'ortholog_one2one') {
    $self->store_param('o_'.$taxon_id, $member->stable_id);
  }
}

sub write_output {
  my $self = shift;

  $self->create_table_from_params( $self->dbc, 'orthologs', $self->{_genes_structure});
  $self->store_params_in_table($self->dbc, 'orthologs', $self->params);
  
}

sub table_def {
  my $self = shift;

  return {
    job_id => 'int',
    unique_keys => 'data_id'
  };
}

1;
