package Bio::Greg::Mammals::Webpages;

use strict;
use Bio::Greg::Codeml;
use File::Path;
use JSON;

use base (
  'Bio::Greg::Hive::Process', 'Bio::Greg::StatsCollectionUtils'
);

my $TREE = 'Bio::EnsEMBL::Compara::TreeUtils';

sub param_defaults {
  return {
    parameter_set_id => 1,
    table => 'web_data',
    skip_file_output => 1
  };
}

sub genes_table_structure {
  my $self = shift;

  return {
    node_id => 'int',
    parameter_set_id => 'tinyint',
    stable_id => 'string',
    gene_name => 'string',

    sites_csv => 'string',
    sites_pdf => 'string',
    sites_png => 'string',
    
    aln_file => 'string',
    aln_pdf => 'string',
    aln_png => 'string',

    tree_file => 'string',
    params_file => 'string',
    data_file => 'string',
    
    ## Potentially useful data to sort by.
    # Comes from stats_genes table.
    chr_name => 'string',
    chr_start => 'string',
    chr_end => 'string',
    chr_strand => 'string',
    pval_fisher => 'float',
    seq_length_mean => 'float',
    slr_dnds => 'float',
    slr_kappa => 'float',
    gc_content_mean => 'float',
    tree_mean_path => 'float',
    leaf_count => 'float',
    sitewise_value_count => 'int',

    # Comes from on-the-fly calculations on sites data.
    f_neg => 'float',
    f_pos => 'float',
    f_neut => 'float',    

    unique_keys => 'node_id,parameter_set_id'
  };
}

sub fetch_input {
  my $self = shift;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  
  $self->create_table_from_params( $self->compara_dba, $self->param('table'),
                                   $self->genes_table_structure );
  
}

sub run {
  my $self = shift;

  # Gather data: Tree, align, sitewise info, params.
  my $params = $self->params;
  my $aln = $self->get_aln;
  my $tree = $self->get_tree;

  my $id = $self->stable_id;
  $self->param('stable_id',$self->stable_id);
  $self->param('gene_name',$self->gene_name);

  $self->pretty_print($aln);
  $self->pretty_print($tree);

  my $output_base = 'web';

  # Sites data.
  my $csv_f = $self->save_file('_sites','csv');
  if (!-e $csv_f->{full_file}) {
    $self->output_sites($csv_f->{full_file});
  }  
  if (-e $csv_f->{full_file}) {
    $self->param('sites_csv',$csv_f->{rel_file});
  }

  # Sites plots.
  eval {
    my $sites_pdf = $self->save_file('_sites','pdf');
    my $sites_png = $self->save_file('_sites','png');
    $self->plot_sites($csv_f->{full_file},$sites_pdf->{full_file},$sites_png->{full_file},$params);

    if (-e $sites_pdf->{full_file}) {
      $self->param('sites_pdf',$sites_pdf->{rel_file});
      $self->param('sites_png',$sites_png->{rel_file});
    }
  };

  # Alignment data.
  my $aln_f = $self->save_file('_aln','fasta');
  $self->save_aln($aln,$aln_f->{full_file});
    
  # Alignment plots.
  my $aln_pdf = $self->save_file('_aln','pdf');
  my $aln_png = $self->save_file('_aln','png');
  $self->plot_aln($aln_f->{full_file},$aln_pdf->{full_file},$aln_png->{full_file},$params);  

  # Tree data.
  my $tree_f = $self->save_file('_tree','nhx');
  $self->save_tree($tree,$tree_f->{full_file});

  # Params.
  my $params_f = $self->save_file('_params','txt');
  $self->save_params($params_f->{full_file});

  # Data.
  my $data_f = $self->save_file('_data','txt');
  $self->save_data($data_f->{full_file});

  # Load a few more calculations into the current params.
  $self->calculate_fractions($csv_f->{full_file});

  # Load params from the stats_genes table.
  $self->load_params;

  # Store refs in database.
  my $row = $self->replace($self->params,{
    aln_file => $aln_f->{rel_file},
    aln_pdf => $aln_pdf->{rel_file},
    aln_png => $aln_png->{rel_file},

    tree_file => $tree_f->{rel_file},
    params_file => $params_f->{rel_file},
    data_file => $data_f->{rel_file}
  });
  $self->store_params_in_table($self->dbc,$self->param('table'),$row);

}


sub save_file {
  my $self = shift;
  my $file_suffix = shift;
  my $file_ext = shift;

  #return if ($self->param('skip_file_output'));

  my $id = $self->stable_id;
  my $subfolder = 'web';
  
  my $params = {
    id => $id,
    filename => $id . $file_suffix,
    extension => $file_ext,
    subfolder => $subfolder
  };
  $self->SUPER::save_file($params);
  
}

sub save_aln {
  my $self = shift;
  my $aln = shift;
  my $aln_file = shift;

  return if ($self->param('skip_file_output'));

  Bio::EnsEMBL::Compara::AlignUtils->to_file( $aln, $aln_file );
}

sub save_data {
  my $self = shift;
  my $filename = shift;

  return if ($self->param('skip_file_output'));

  my $obj = $self->_get_gene_params;

  open(my $out,">${filename}");
  $self->hash_print($obj,$out);
  close($out);  
}

sub load_params {
  my $self = shift;

  my $obj = $self->_get_gene_params;

  $self->set_params($obj);
}

sub _get_gene_params {
  my $self = shift;

  my $node_id = $self->param('node_id');
  my $parameter_set_id = $self->param('parameter_set_id');

  my $sth = $self->dbc->prepare("SELECT * from stats_genes where node_id=${node_id} and parameter_set_id=${parameter_set_id};");
  $sth->execute;
  while ( my $obj = $sth->fetchrow_hashref ) {
    $sth->finish;
    return $obj;
  }
}

sub save_params {
  my $self = shift;
  my $filename = shift;

  return if ($self->param('skip_file_output'));

  open(my $out,">${filename}");
  $self->hash_print($self->params,$out);
  close($out);
}

sub save_tree {
  my $self = shift;
  my $tree = shift;
  my $filename = shift;

  return if ($self->param('skip_file_output'));

  Bio::EnsEMBL::Compara::TreeUtils->to_file( $tree, $filename );
}

sub plot_aln {
  my $self = shift;
  my $aln_file = shift;
  my $pdf_filename = shift;
  my $png_filename = shift;
  my $params = shift;

  return if ($self->param('skip_file_output'));

  my $cmd = qq^
  library(phylosim)
 
  sim <- PhyloSim();
  readAlignment(sim,file="$aln_file")

  pdf(file="$pdf_filename",width=15,height=15)
  plot(sim,plot.chars=T)
  dev.off()
  png(file="$png_filename",width=1500,height=300)
  plot(sim,plot.chars=F,num.pages=1)
  dev.off()

^;
  Bio::Greg::EslrUtils->run_r($cmd);
  


}

sub plot_sites {
  my $self = shift;
  my $sites_csv = shift;
  my $pdf_filename = shift;
  my $png_filename = shift;
  my $params = shift;

  return if ($self->param('skip_file_output'));

  my $folder = $self->get_output_folder;
  my $sites_file = "${folder}/sites_1.Rdata";

  my $node_id = $self->param('node_id');

  my $script = Bio::Greg::EslrUtils->baseDirectory . "/scripts/plot_sitewise.R";

  my $cmd = qq^
  source("$script")
  sites <- read.csv("$sites_csv",stringsAsFactors=F)

  if (nrow(sites) == 0) { q() }

  library(ggplot2)
  pdf("$pdf_filename",width=10,height=3)
  p <- ggplot.gene(sites,'aln_position')
  print(p)
  dev.off()
  png("$png_filename",width=500,height=150)
  p <- ggplot.gene(sites,'aln_position')
  print(p)
  dev.off()
^;
  print $cmd;

  Bio::Greg::EslrUtils->run_r($cmd);

}

sub output_sites {
  my $self = shift;
  my $filename = shift;
  my $params = shift;

  return if ($self->param('skip_file_output'));

  my $folder = $self->get_output_folder;
  my $sites_file = "${folder}/sites_1.Rdata";
  my $node_id = $self->param('node_id');
  my $script = Bio::Greg::EslrUtils->baseDirectory . "/scripts/plot_sitewise.R";

  my $cmd = qq^
  source("$script")
  load("$sites_file")
  sites <- subset(sites,node_id==${node_id})
  sites <- sign.lrt(sites)
  if (nrow(sites) > 0) {
    write.csv(sites,file="$filename",col.names=F,quote=F)
  }
^;
  Bio::Greg::EslrUtils->run_r($cmd);
}


sub write_output {

}
