package Bio::Greg::Gorilla::BranchTests;

use strict;
use File::Path;
use Bio::Greg::Codeml;
use FreezeThaw qw(freeze thaw cmpStr safeFreeze cmpStrHard);
use Data::Types qw(:all);

use base (
  'Bio::Greg::StatsCollectionUtils',
  'Bio::Greg::Hive::Process',
  'Bio::Greg::Hive::Align',
  'Bio::Greg::Hive::AlignmentScores',
  'Bio::Greg::Hive::CountSubstitutions'
);

sub param_defaults {
  return {
  };
}

sub fetch_input {
  my ($self) = @_;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  $self->{_genes_structure} = $self->genes_table;

  # Create tables if necessary.
#  $self->create_table_from_params( $self->compara_dba, 'subs',
#    $self->sites_table );
  $self->create_table_from_params( $self->compara_dba, 'subs',
    $self->codon_subs_table_def );
}

sub run {
  my $self = shift;

  $self->param('force_recalc', 0);

  my $gene_id = $self->param('gene_id');
  my $mba = $self->compara_dba->get_MemberAdaptor;
  my $member = $mba->fetch_by_source_stable_id(undef, $gene_id);
  my $gene = $member->get_Gene;
  $self->param('gene_name', $gene->external_name);
  $self->param('data_id', $member->dbID);

  my $ortholog_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_one_to_one_ortholog_tree($self->compara_dba, $member);

  my ($tree, $aln) = $self->_get_tree_and_alignment($ortholog_tree);
  $self->pretty_print($aln);

  $self->param('force_recalc', 0);

  my $params = {
    keep_species => join(',', map {$_->taxon_id} $tree->leaves)
  };
  my $genome_tree = Bio::EnsEMBL::Compara::ComparaUtils->get_genome_tree_subset($self->compara_dba, $params);
  # Match up genome tree names to gene tree IDs.
  my $id_hash;
  map {$id_hash->{$_->taxon_id} = $_->name} $ortholog_tree->leaves;
  foreach my $node ($genome_tree->nodes) {
    if ($node->is_leaf) {
      $node->name($id_hash->{$node->taxon_id});
    } else {
      $node->name('');
    }
  }

  $self->_run_tests($genome_tree, $aln);

  $self->_plot_subs($genome_tree, $aln);
}

sub write_output {
  my $self = shift;

  $self->create_table_from_params( $self->compara_dba, 'genes', $self->{_genes_structure});
  $self->store_params_in_table($self->compara_dba, 'genes', $self->params);
}

sub _tx_aln {
  my $self = shift;
  my $aln = shift;
  return Bio::EnsEMBL::Compara::AlignUtils->translate($aln);
}

sub _get_tree_and_alignment {
  my $self = shift;
  my $tree = shift;

  my @aln_types = ('compara', 'genomic_primates', 'genomic_mammals', 'genomic_all');

  my $keep_type = 'genomic_all';
  my $keep_file;

  foreach my $type (@aln_types) {
    my $f = $self->_save_file('aln_'.$type, 'fasta');
    my $file = $f->{full_file};
    $self->store_param('aln_'.$type.'_file', $f->{rel_file});

    if ($type eq $keep_type) {
      $keep_file = $f;
    }

    my $raw_f = $self->_save_file('aln_'.$type.'_raw', 'fasta');
    my $raw_file = $raw_f->{full_file};
    $self->store_param('aln_'.$type.'_raw_file', $raw_f->{rel_file});
    
    if (!-e $file || $self->param('force_recalc')) {
      my $ref_species = 9606;
      my @members     = $tree->leaves;
      my ($ref_member) = grep { $_->taxon_id == $ref_species } @members;
      die("No ref member!") unless (defined $ref_member);

      # Fetch the alignment.
      my $params = $self->params;
      $params->{quality_threshold} = 30;
      my $tree_aln_obj =
        Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( 
          $self->compara_dba, $tree, $ref_member, $params
        );
      my $aln = $tree_aln_obj->{aln};
      my $aln_tree = $tree_aln_obj->{tree};
      
      # Realign compara alignments w/ PRANK...
      if ($type eq 'compara') {
        $params->{alignment_prank_codon_model} = 1;
        $aln = $self->align_with_prank($aln, $aln_tree, $params);
      }

      # Get the raw (unfiltered, not-realigned) alignment.
      $params->{quality_threshold} = 0;
      my $tree_aln_obj =
        Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( 
          $self->compara_dba, $tree, $ref_member, $params
        );
      my $raw_aln = $tree_aln_obj->{aln};
    
      Bio::EnsEMBL::Compara::AlignUtils->to_file($aln, $file);
      Bio::EnsEMBL::Compara::AlignUtils->to_file($raw_aln, $raw_file);
    }
  }

  my $tree_f = $self->_save_file('tree', 'nh');
  my $tree_file = $tree_f->{full_file};
  $self->param('tree_file', $tree_f->{rel_file});  
  if (!-e $tree_file || $self->param('force_recalc')) {
    Bio::EnsEMBL::Compara::TreeUtils->to_file($tree, $tree_file);
  }
  my $tree_from_file = Bio::EnsEMBL::Compara::TreeUtils->from_file($tree_file);

  my $aln = Bio::EnsEMBL::Compara::AlignUtils->from_file($keep_file->{full_file});

  # Get sub-tree of primate sequences.
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_clade(
    $self->compara_dba, $tree, "Primates"
    );
  
  # Restrict the tree & alignment to each other.
  $tree = Bio::EnsEMBL::Compara::ComparaUtils->restrict_tree_to_aln($tree, $aln);
  $aln = Bio::EnsEMBL::Compara::ComparaUtils->restrict_aln_to_tree($aln, $tree);
  
  return ($tree, $aln);
}

sub _run_tests {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  my $branch_f = $self->_save_file('paml_branch_results', 'perlobj');
  my $branch_file = $branch_f->{full_file};
  $self->param('branch_results_file', $branch_f->{rel_file});  
  if (!-e $branch_file || $self->param('force_recalc')) {  

    my $tree_copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);

    my $config = {
      0 => [],
      1 => ['H'],
      2 => ['C'],
#      3 => ['H', 'CH'],
#      4 => ['C', 'CH'],
#      5 => ['G'],
#      6 => ['CGH'],
#      7 => ['CH', 'G', 'H', 'C'],
#      8 => ['H', 'C'],
#      9 => ['H', 'C'],
#      10 => ['C', 'G']
    };
    my $results = $self->run_branches($tree_copy, $aln, $config);
    my $obj = $results;
    open(OUT,">$branch_file");
    print OUT freeze($obj);
    close(OUT);
  }
  open(IN, $branch_file);
  my @lines = <IN>;
  close(IN);
  my ($branch_results) = thaw(join('',@lines));

  $self->parse_branch_results($branch_results, $tree, $aln);  
  return;

  my $sites_f = $self->_save_file('paml_branchsites_results', 'perlobj');
  my $sites_file = $sites_f->{full_file};
  $self->param('branchsites_results_file', $sites_f->{rel_file});  
  if (!-e $sites_file || $self->param('force_recalc')) {  
  
    my $config = {
      11 => ['CGH'],
      12 => ['H'],
      13 => ['C'],
      14 => ['G'],
      15 => ['O'],
      16 => ['CGHO']
    };

    my $results = $self->run_branchsites($tree, $aln, $config);
    my $obj = $results;
    open(OUT,">$sites_file");
    print OUT freeze($obj);
    close(OUT);
  }
  open(IN, $sites_file);
  my @lines = <IN>;
  close(IN);
  my ($branchsites_results) = thaw(join('',@lines));

  $self->parse_results($branchsites_results);  
}

sub run_branches {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $config = shift;

  my $lines_hash;

  foreach my $key (keys %$config) {
    my $fg_branches = $config->{$key};

    print "  running $key...\n";
    my $copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);
    my $lines = $self->run_branch($copy, $aln, $fg_branches, $key);
    $lines_hash->{$key} = $lines;
  }
  return $lines_hash;
}

sub parse_branch_results {
  my $self = shift;
  my $results = shift;
  my $tree = shift;
  my $aln = shift;

  foreach my $key (keys %$results) {
    my $lines = $results->{$key};

    my $lnL = Bio::Greg::Codeml->extract_lnL($lines);
    my @omegas = Bio::Greg::Codeml->extract_omegas($lines);

    $self->param("${key}_lnL", $lnL);
    for ( my $i = 0 ; $i < scalar @omegas ; $i++ ) {
      my $omega = $omegas[$i];
      $self->param("${key}_omega_${i}", $omega);
    }

    if ($key == 0) {
      $self->extract_subs($lines, $tree, $aln);
    }
  }
}

sub run_branch {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $fg_branches = shift;
  my $test_num = shift;

  my $key = "$test_num";

  my ($labeled_tree, $labeled_aln) = $self->label_branches($tree, $aln, $fg_branches);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($labeled_tree);
  print $treeI->ascii;
  my $res = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $labeled_aln, $self->worker_temp_directory, $self->params );

  return $res->{lines};
}

sub extract_subs {
  my $self = shift;
  my $lines = shift;
  my $tree = shift;
  my $aln = shift;

  my $pep_aln = Bio::EnsEMBL::Compara::AlignUtils->translate($aln);

  my $codeml_tree = Bio::Greg::Codeml->parse_codeml_results($lines);

  foreach my $node ( $codeml_tree->get_nodes ) {
    my $id = $node->id;
    
    # Grab the substitutions from the codeml tree node.
    my $subs_hash = $node->get_tag_values('substitutions');
    my $i         = 0;
    foreach my $sub_key ( sort { $a <=> $b } keys %$subs_hash ) {
      my $subst = $subs_hash->{$sub_key};
      
      #$self->hash_print($subst);
      my $obj = $self->extract_substitution_info(
        $node,
        $subst,
        $tree,
        $pep_aln,
        $aln
        );
      if ( defined $obj ) {
        $i++;
        my $params = $self->replace( $self->params, $obj );
        $self->store_params_in_table( $self->dbc, 'subs', $params );
      } else {
        
        #print "No results for subst:\n";
        #$self->hash_print($subst);
      }
    }  
    print "Stored $i substitutions for $id\n" if ( $self->debug );
  }
}

sub run_branchsite {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $fg_branches = shift;
  my $test_num = shift;

  my $key = "$test_num";

  my ($labeled_tree, $labeled_aln) = $self->label_branches($tree, $aln, $fg_branches);
  my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($labeled_tree);
  print $treeI->ascii;

  # Null model.
  my $null = $self->params;
  $null->{model} = 2;
  $null->{NSsites} = 2;
  $null->{fix_omega} = 1;
  $null->{omega} = 1;

  # Alternative model.
  my $alt = $self->params;
  $alt->{model} = 2;
  $alt->{NSsites} = 2;
  $alt->{fix_omega} = 0;
  $alt->{omega} = 0.3;

  my $res_null = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $labeled_aln, $self->worker_temp_directory, $null );
  my @null_omegas = @{ $res_null->{omegas} };
  for ( my $i = 0 ; $i < scalar @null_omegas ; $i++ ) {
    my $omega = $null_omegas[$i];
    $self->param( "omega_${key}_null_${i}", $omega );
  }
  $self->param("lnL_${key}_null", $res_null->{lnL});

  my $res_alt = Bio::Greg::Codeml->branch_model_likelihood( $treeI, $labeled_aln, $self->worker_temp_directory, $alt );
  my @alt_omegas = @{ $res_alt->{omegas} };
  for ( my $i = 0 ; $i < scalar @alt_omegas ; $i++ ) {
    my $omega = $alt_omegas[$i];
    $self->param( "omega_${key}_alt_${i}", $omega );
  }
  $self->param("lnL_${key}_alt", $res_alt->{lnL});
}


sub label_branches {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $fg_branches = shift;


  my $copy = Bio::EnsEMBL::Compara::TreeUtils->copy_tree($tree);

  my $char_to_taxon = {
    'H' => 9606,
    9606 => 'H',
    'C' => 9598,
    9598 => 'C',
    'G' => 9593,
    9593 => 'G',
    'O' => 9600,
    9600 => 'O'
  };

  foreach my $fg_branch (@$fg_branches) {
    foreach my $node ($copy->nodes) {
      my @leaves = $node->leaves;
      my @leaf_strings = map {
        if ($char_to_taxon->{$_->taxon_id}) {
          $char_to_taxon->{$_->taxon_id};
        }
      } @leaves;
      my $leaf_string = join('',sort {$a cmp $b} @leaf_strings);

      if ($leaf_string eq $fg_branch) {
        $node->name($node->name . '#1');
      }
    }
  }

  return ($copy, $aln);
}

sub _plot_subs {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;

  $self->param('force_recalc', 1);

  my $plot_aln_f = $self->_save_file('plot_aln', 'pdf');
  my $plot_aln_file = $plot_aln_f->{full_file};
  $self->param('plot_aln_file', $plot_aln_f->{rel_file});  

  my $plot_tree_f = $self->_save_file('plot_tree', 'pdf');
  my $plot_tree_file = $plot_tree_f->{full_file};
  $self->param('plot_tree_file', $plot_tree_f->{rel_file});  

  if (!-e $plot_aln_file || $self->param('force_recalc')) {  
    
    my @features;
    
    my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);
    
    my $id_to_species;
    map {
      $id_to_species->{$_->name} = $_->ensembl_alias_name;
    } $tree->leaves;
    
    map {
      $_->branch_length(1);
    } $treeI->nodes;
    
    my $data_id = $self->param('data_id');
    my $sth = $self->dbc->prepare("SELECT * FROM subs where data_id=$data_id");
    
    $sth->execute;
    my $muts;
    my @rows = @{ $sth->fetchall_arrayref({}) }; # Fetch all rows into an arrayref of hashrefs.
    $sth->finish;
    foreach my $row (@rows) {
      my $key = $row->{leaves_beneath};
      $muts->{$key} = [] if (!$muts->{$key});
      push @{$muts->{$key}}, $row;
    }

    my @keys = keys %$muts;
    @keys = sort {
      scalar(split(',', $b)) <=> scalar(split(',', $a))
    } @keys;
    
    foreach my $key (@keys) {
      # Take all mutations for a given branch, split the syn ones in half and
      # 'engulf' the nsyn's w/ the syn's.
      my @cur_branch_muts = @{$muts->{$key}};
      my @syn_rows = grep {$_->{mut_nsyn} == 0} @cur_branch_muts;
      my @nsyn_rows = grep {$_->{mut_nsyn} == 1} @cur_branch_muts;

      foreach my $row (@syn_rows, @nsyn_rows) {
        # Using the 'leaves_beneath' string from the table,
        # parse out the taxon IDs of the leaves beneath,
        # and match up to the node from our TreeI object
        my $leaves = $row->{leaves_beneath};
        
        my @ids = split(",", $leaves);
        
        my $org_string = join("/", map {$id_to_species->{$_}} @ids);
        
        my @id_nodes = map {$treeI->find($_)} @ids;
        my $tree_node;
        if (scalar(@id_nodes) > 1) {
          $tree_node = $treeI->lca(@id_nodes);
        } else {
          $tree_node = $id_nodes[0];
        }
        
        # Create a new internal node representing this mutation event in the tree.
        my $new_node = new $tree_node;
        $new_node->id('m_'.$row->{aln_pos});
        $new_node->set_tag_value('nsyn', $row->{mut_nsyn});
        $new_node->set_tag_value('from', $row->{nuc_from});
        $new_node->set_tag_value('to', $row->{nuc_to});
        
        $tree_node->split_branch_with_node($new_node, 1);
        $new_node->branch_length(1);
        $tree_node->branch_length(0);
        
        my $tax_name = $tree_node->id;
        if ($id_to_species->{$tree_node->id}) {
          $tax_name = $id_to_species->{$tree_node->id};
        }
        
        my $f = Bio::SeqFeature::Generic->new(
          -start => $row->{aln_pos},
          -end => $row->{aln_pos}+1,
          -score => $row->{mut_nsyn},
          -source => $org_string." substitutions"
          );
        push @features, $f;
      }
    }
    
    map {
      $_->id($id_to_species->{$_->id});
    } $treeI->leaves;
    
    map {
      $_->id($id_to_species->{$_->id});
    } $aln->each_seq;
    
    my $str = $treeI->root->as_text('nhx');
    my $gene_name = $self->param('gene_name');
    
    my $tmp = $self->worker_temp_directory;
    
    my $tree_subs_file = "$tmp/${gene_name}_subs.nhx";
    my $tree_file = "$tmp/${gene_name}.nh";
    my $aln_file = "$tmp/${gene_name}.fasta";
    my $gff_file = "$tmp/${gene_name}.gff";
    
    open(OUT, ">$tree_subs_file");
    print OUT $str."\n";
    close(OUT);
    
    Bio::EnsEMBL::Compara::TreeUtils->to_file($treeI, $tree_file);  
    
    Bio::EnsEMBL::Compara::AlignUtils->to_file($self->_tx_aln($aln), $aln_file);  
    
    my $gff = Bio::Tools::GFF->new(
      -file => ">$gff_file",
      -gff_version => 3
      );
    foreach my $f (@features) {
      $f->seq_id('chr1');
      $f->strand('+');
      $gff->write_feature($f);
    }
    $gff->close;
    
    my $rcmd = qq^
library(phylosim)
source("~/src/greg-ensembl/projects/phylosim/PhyloSimPlots.R")

sim <- PhyloSim()
readAlignment(sim, "${aln_file}")
pdf(file="${plot_aln_file}", width=30, height=6)
  plotAlignment(sim, aln.gff.file="${gff_file}", axis.text.size=5, aln.plot.chars=T)
dev.off()

sim <- PhyloSim()
phylo <- read.nhx("$str")
sim\$.phylo <- phylo

pdf(file="${plot_tree_file}")

xlab <- "Number of inferred substitutions since LCA"
obj <- plotTree(sim, tree.do.plot=F, color.by='nsyn', tree.xlab=xlab, tree.xlim.expand=c(0.05, 0.5))
p <- obj
p <- p + scale_colour_manual("Substitution Type", 
  values=c("gray", "red"), 
  breaks=c(0,1), 
  labels=c('Synonymous', 'Nonsynonymous')
)
p <- p + opts(title="${gene_name}")
print(p)
dev.off()
^;
    Bio::Greg::EslrUtils->run_r($rcmd);  
  }
}

sub store_param {
  my $self = shift;
  my $param = shift;
  my $value = shift;

  # First, store in the local param hashref.
  $self->param($param, $value);

  # Now, make sure this param is stored in our tables structure.
  my $str = $self->{_genes_structure};
  $str->{$param} = 'string';
  if (is_int($value)) {
    $str->{$param} = 'int';
  } elsif (is_float($value)) {
    $str->{$param} = 'float';
  }
}

sub _save_file {
  my $self = shift;
  my $filename_base = shift;
  my $ext = shift;

  my $id = $self->param('gene_name');

  my $filename = "${id}_${filename_base}";

  my $file_params = {
    id => $id,
    filename => $filename,
    extension => $ext,
    subfolder => 'data',
  };

  my $file_obj = $self->save_file($file_params);

  if ($self->param('data_tarball') && !-e $file_obj->{full_file}) {
    # Load the file from the tarball into our temp directory.
    #my $tarball = $self->param('data_tarball');
    #my $tmp = $self->worker_temp_directory;
    #my $rel_f = $file_obj->{rel_file};

    #my $cmd = qq^tar -zxvf $tarball $rel_f^;
    #print $cmd."\n";
    #$file_obj->{full_file} = $tmp . $filename;
  }

  return $file_obj;
}


sub fail {
  my $self = shift;
  my $reason_type = shift;
  my $reason_msg = shift;

  $self->param('fail_type', $reason_type);
  $self->param('fail_msg', $reason_msg);

  $self->fail_and_die("Job fail: $reason_msg");
}

sub genes_table {
  my $self = shift;
}

sub sites_table {
  my $self = shift;
}


1;
