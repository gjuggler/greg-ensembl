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
    skip_file_output => 0
  };
}

sub genes_table_structure {
  my $self = shift;

  return {
    job_id => 'int',
    node_id => 'int',
    parameter_set_id => 'tinyint',
    parameter_set_shortname => 'char8',

    # Reference gene info.
    gene_id => 'string',
    gene_name => 'string',
    gene_description => 'string',
    gene_taxon_id => 'int',
    protein_id => 'string',
    transcript_id => 'string',
    chr_name => 'string',
    chr_start => 'string',
    chr_end => 'string',
    chr_strand => 'string',  # Forgot to include this in SitewiseMammals.pm, so we need to collect this temporarily.

    # Sitewise data.
    sitewise_value_count => 'int',
    pval_fisher => 'float',
    f_neg => 'float',
    f_pos => 'float',
    f_neut => 'float',    
    n_neg => 'int',
    n_pos => 'int',
    n_neut => 'int',
    slr_dnds => 'float',
    slr_kappa => 'float',
    pfam_domains => 'string',  # We're collecting these here, but it could/should probably be moved to SitewiseMammals.pm

    # Alignment properties
    aln_length => 'int',
    seq_length => 'int',
    gc_cds => 'float',
    'gc_3' => 'float',
    gc_genomic => 'float',
    tree_mean_path => 'float',
    leaf_count => 'float',
    duplication_count => 'int',

    # Files we're generating.
    sites_file => 'string',    
    tree_file => 'string',
    aln_file => 'string',
    pep_aln_file => 'string',
    params_file => 'string',
    data_file => 'string',
    aln_pdf => 'string',
    aln_png => 'string',
    
    unique_keys => 'node_id,parameter_set_id'
  };
}

sub fetch_input {
  my $self = shift;

  # Fetch parameters from all possible locations.
  $self->load_all_params();

  $self->create_table_from_params( $self->compara_dba, $self->_get_output_table_name,
                                   $self->genes_table_structure );  
}

sub _get_output_table_name {
  my $self = shift;

  my $shortname = $self->param('parameter_set_shortname');
  my $table = "web_data_${shortname}";
  return $table;
}

sub _get_sitewise_output_dir {
  my $self = shift;

  my $output_dir = $self->get_output_folder;
  return "${output_dir}/data";
}

sub run {
  my $self = shift;

  $self->_load_params_from_genes_table;

  my $aln_file = $self->param('aln_file');
  # Turn the relative filename into something meaningful.
  $aln_file = $self->_get_sitewise_output_dir . "/$aln_file";

  if (!-e $aln_file) {
    die("Alignment file not found!");
  }
  # Set the existing_alignment_file flag to be used by the ComparaUtils method to fetch
  # the stored alignment file.
  $self->param('existing_alignment_file', $aln_file);
  my $tree = $self->get_tree;
  my $ref_member = $self->get_ref_member($tree); # Defined in StatsCollectionUtils.

  my $tree_aln_obj =
    Bio::EnsEMBL::Compara::ComparaUtils->get_compara_or_genomic_aln( $self->compara_dba, $tree, $ref_member,
                                                                     $self->params );  
  my $aln = $tree_aln_obj->{aln};
  my $pep_aln = $tree_aln_obj->{pep_aln};
  my $tree = $tree_aln_obj->{tree};
  $self->_output_files($tree, $aln, $pep_aln);

  $self->pretty_print($pep_aln);

  # Collect stuff we forgot in SitewiseMammals.pm.
  #   - first the strand.
  my ($ref_member) = grep { $_->stable_id == $self->param('protein_id') } $tree->leaves;
  my $strand = $ref_member->get_Gene->seq_region_strand;
  $self->param('chr_strand', $strand);
  $self->_collect_pfam_ids;

  # Sites data.
  my $sites_f = $self->_save_file('_sites','csv');
  if (!-e $sites_f->{full_file}) {
    $self->_output_sites($self->param('tree_file_full'), $self->param('pep_aln_file_full'), $sites_f->{full_file});
  }  
  if (-e $sites_f->{full_file}) {
    $self->param('sites_file',$sites_f->{rel_file});
    $self->param('sites_file_full',$sites_f->{full_file});
  }

  # Alignment plots.
  my $aln_pdf = $self->_save_file('_aln','pdf');
  my $aln_png = $self->_save_file('_aln','png');
  if (!-e $aln_pdf->{full_file}) {
    $self->_plot_aln($self->param('tree_file_full'), $self->param('pep_aln_file_full'),$self->param('sites_file_full'),$aln_pdf->{full_file},$aln_png->{full_file},$self->params);
  }
  $self->param('aln_pdf',$aln_pdf->{rel_file});
  $self->param('aln_png',$aln_png->{rel_file});

  # Params.
  my $params_f = $self->_save_file('_params','txt');
  $self->_save_params($params_f->{full_file});
  $self->param('params_file', $params_f->{rel_file});

  # Store refs in database.
  $self->store_params_in_table($self->dbc,$self->_get_output_table_name,$self->params);
}

sub _sites_table_name {
  my $self = shift;

  my $parameter_set_shortname = $self->param('parameter_set_shortname');
  my $table = "sites_c";
  $table = "sites_g" if ($parameter_set_shortname eq 'm_g');
  return $table;
}

sub _genes_table_name {
  my $self = shift;

  my $parameter_set_shortname = $self->param('parameter_set_shortname');
  my $table = "genes_c";
  $table = "genes_g" if ($parameter_set_shortname eq 'm_g');
  return $table;
}

sub _collect_pfam_ids {
  my $self = shift;

  my $node_id = $self->param('node_id');
  my $table = $self->_sites_table_name;

  my %pfam_ids;
  my $sth = $self->dbc->prepare("SELECT * from ${table} where node_id=${node_id} and pfam_domain IS NOT NULL;");
  $sth->execute;
  while ( my $obj = $sth->fetchrow_hashref ) {
    $pfam_ids{$obj->{pfam_domain}} = 1;
  }
  $sth->finish;

  my $pfam_string = join(',', sort keys %pfam_ids);
  $self->param('pfam_domains', $pfam_string);
}

sub _output_files {
  my $self = shift;
  my $tree = shift;
  my $aln = shift;
  my $pep_aln = shift;

  my $aln_f = $self->_save_file('_aln','fasta');
  my $pep_aln_f = $self->_save_file('_pep_aln','fasta');
  my $tree_f = $self->_save_file('_tree','nh');

  Bio::EnsEMBL::Compara::TreeUtils->to_file($tree,$tree_f->{full_file});
  Bio::EnsEMBL::Compara::AlignUtils->to_file( $aln, $aln_f->{full_file} );
  Bio::EnsEMBL::Compara::AlignUtils->to_file( $pep_aln, $pep_aln_f->{full_file} );

  $self->param('aln_file',$aln_f->{rel_file});
  $self->param('pep_aln_file',$pep_aln_f->{rel_file});
  $self->param('tree_file',$tree_f->{rel_file});

  $self->param('aln_file_full', $aln_f->{full_file});
  $self->param('pep_aln_file_full', $pep_aln_f->{full_file});
  $self->param('tree_file_full', $tree_f->{full_file});
}

sub _save_file {
  my $self = shift;
  my $file_suffix = shift;
  my $file_ext = shift;

  my $id = $self->param('protein_id');

  my $shortname = $self->param('parameter_set_shortname');
  my $subfolder = "web_data_${shortname}";
  
  my $params = {
    id => $id,
    filename => $id . '_' . $shortname . $file_suffix,
    extension => $file_ext,
    subfolder => $subfolder
  };
  return $self->save_file($params);
}

sub _load_params_from_genes_table {
  my $self = shift;

  my $node_id = $self->param('node_id');
  my $table = $self->_genes_table_name;
  print "GENES TABLE: $table\n";
  my $sth = $self->dbc->prepare("SELECT * from ${table} where node_id=${node_id};");
  $sth->execute;
  while ( my $obj = $sth->fetchrow_hashref ) {
    $sth->finish;
    $self->hash_print($obj);
    $self->set_params($obj);
    return;
  }

  $self->fail_and_die("No row found in the genes table -- must have failed!");

}

sub _save_params {
  my $self = shift;
  my $filename = shift;

  open(my $out,">${filename}");
  $self->hash_print($self->params,$out);
  close($out);
}

sub _plot_aln {
  my $self = shift;
  my $tree_file = shift;
  my $aln_file = shift;
  my $sites_file = shift;
  my $pdf_filename = shift;
  my $png_filename = shift;
  my $params = shift;

  my $cmd = qq^
  library(phylosim)

  # Create the phylosim object and read the alignment.
  sim <- PhyloSim();
  readAlignment(sim,file="$aln_file")
  readTree(sim,file="$tree_file")

  tracks <- NULL
  if (file.exists("${sites_file}")) {
    # Create the sites track to plot above the alignment.
    aln.row.count <- length(sim\$alignment[,1])
    aln.col.count <- length(sim\$alignment[1,])

    sites <- read.csv("$sites_file",stringsAsFactors=F)
    clade <- sites[,'clade']
    pos <- sites[,'aln_position']
    y.max <- log10(10)
    y.min <- log10(0.01)
    y_lo <- log10(sites[,'omega_lower'] + 0.01)
    y_hi <- log10(sites[,'omega_upper'] + 0.01)
    y_lo <- (y_lo - y.min) / (y.max-y.min)
    y_hi <- (y_hi - y.min) / (y.max-y.min)
    y_lo <- pmax(y_lo,0)
    y_hi <- pmin(y_hi,1)

    score <- sites[,'signed_lrt'] 
    max.lrt <- 20
    score <- (score + max.lrt) / (max.lrt*2)
    score <- pmax(score,0)
    score <- pmin(score,1)

    # Some predefined track heights.
    track.height <- as.integer(aln.row.count/2)
    small.track.height <- as.integer(aln.row.count/4)
    tiny.track.height <- as.integer(aln.row.count/8)

    # Create the main sites track.
    slr.track <- data.frame(
      id="Sitewise dN/dS (log scale)",height=as.integer(aln.row.count),
      color.gradient='blue,gray,gray,red',
      pos=pos,
      y_lo=y_lo,
      y_hi=y_hi,
      score=score,
      clade=clade,
      background='white',
      height=track.height
    )
    # Limit the SLR track to the max length of the alignment.
    slr.track <- subset(slr.track, pos <= aln.col.count)

    space.track <- data.frame(id='',height=tiny.track.height,pos=1)
    tracks <- list(space.track,slr.track)

    # Collect the nongap branch lengths and add a track for that.
    if (length(sites[,'nongap.bl']) > 0) {
      nongap.bl <- sites[,'nongap.bl']
      pos <- sites[,'aln_position']

      max.nongap.bl <- max(nongap.bl, na.rm=T)
      
      nongap.track <- data.frame(
        id="Non-gapped branch length",
        pos=pos,
        background='white',
        score=1,
        y_lo=0,
        y_hi=(nongap.bl/(max.nongap.bl*1.1)),
        col.gradient='gray60,gray60',
        layout='below'
      )
      tracks <- c(tracks, list(space.track, nongap.track))
    }

    # Collect any Pfam domains and add tracks for those.
    pfam.sites <- subset(sites, !is.na(pfam_domain))
    if (nrow(pfam.sites) > 0) {
      pfam.tracks <- list()
      pfam.ids <- sort(unique(pfam.sites[, 'pfam_domain']))
      for (pfam in pfam.ids) {
        cur.pfam.sites <- subset(pfam.sites, pfam_domain==pfam)
        cur.pfam.track <- data.frame(
          pos=cur.pfam.sites[, 'aln_position'],
          background='white',
          height=small.track.height,
          id=pfam,
          score=1,
          y_lo=0,
          y_hi=1,
          col.gradient='black,blue'
        )
        pfam.tracks <- c(pfam.tracks, list(cur.pfam.track))
      }
      tracks <- c(list(space.track), pfam.tracks,list(space.track), tracks)
    }
  }

  pdf(file="$pdf_filename", width=15, height=15)
  plot(sim, plot.chars=T, tracks=tracks, tree.xlim=c(0,1))
  dev.off()
  png(file="$png_filename", width=1500, height=300)
  plot(sim, plot.chars=F, plot.labels=F, num.pages=1, tracks=tracks, tree.xlim=c(0,1))
  dev.off()
^;
  Bio::Greg::EslrUtils->run_r($cmd);
}

sub _output_sites {
  my $self = shift;
  my $tree_file = shift;
  my $pep_aln_file = shift;
  my $filename = shift;
  my $params = shift;

  my $node_id = $self->param('node_id');
  my $sites_table = $self->_sites_table_name;

  my $collect_script = Bio::Greg::EslrUtils->baseDirectory . "/scripts/collect_sitewise.R";
  my $plot_script = Bio::Greg::EslrUtils->baseDirectory . "/scripts/plot_sitewise.R";

  my $dbname = $self->hive_dba->dbc->dbname;

  my $clade = $self->param('parameter_set_name');

  my $cmd = qq^
  dbname <- "${dbname}"
  source("${collect_script}")
  source("${plot_script}")
  sites <- get.vector(con, "SELECT * from ${sites_table} where node_id=${node_id}")

  if (nrow(sites) > 0) {
    # These methods are defined in plot_sitewise.R
    sites <- sign.lrt(sites)
    sites <- add.pval(sites)

    # Function 'add.nongap.bl' is also defined in plot_sitewise.R
    library(phylosim)
    sim <- PhyloSim()
    readTree(sim, "${tree_file}")
    readAlignment(sim, "${pep_aln_file}")
    sites <- add.nongap.bl(sites, sim\$.phylo, sim\$.alignment)

    sites[,'clade'] <- "${clade}"
    write.csv(sites,file="$filename",quote=F)
  }
^;
  Bio::Greg::EslrUtils->run_r($cmd);
}


sub write_output {

}
