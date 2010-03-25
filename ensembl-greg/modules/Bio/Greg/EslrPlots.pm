package Bio::Greg::EslrPlots;

use strict;
use File::Path;
use Bio::EnsEMBL::Compara::ComparaUtils;
use Bio::EnsEMBL::Compara::AlignUtils;
use Bio::EnsEMBL::Compara::TreeUtils;

sub plotTrees {
  my $class    = shift;
  my $file_out = shift;
  my $params   = shift;
  my $tree_in  = shift;

  my @node_ids = @{ $params->{'node_ids'} };

  my $height = $params->{'height'} || 30;
  my $width  = $params->{'width'}  || 20;

  my $n = scalar(@node_ids);

  my $r_start = qq{
library(ape);
source("aln_tools/aln.tools.R");
source("aln_tools/slr.tools.R");
source("aln_tools/phylo.tools.R");
source("aln_tools/plot.phylo.greg.R");

width = $width;
height = $height;

mult = 50 / width;
width = width*mult;
height = height*mult;

dpar = par(no.readonly=T)
print(paste("Width:",width," Height:",height));
pdf(file="${file_out}",width=$width,height=$height,pointsize=12);

par(mfrow=c($n,1));

};

  my $r_end = qq{
dev.off();
q();
};

  my @cmds;
  foreach my $node_id (@node_ids) {
    $params->{'node_id'} = $node_id;
    push @cmds, $class->plotTree( $file_out, $params, $tree_in, 1 );
  }

  my $full_cmd = $r_start . join( "\n", @cmds ) . $r_end;

  my $dir = "/tmp/eslrplots";
  mkpath( [$dir] );
  my $temp = "$dir/rcmd.txt";
  open( OUT, ">$temp" );
  print OUT $full_cmd . "\n";
  close(OUT);

  # Create the PDF, then use ImageMagick to create a PNG and thumbnail.
  # cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  my $rc      = 0;
  my $nullify = "";
  if ( $params->{'quiet'} ) {
    $nullify = ">/dev/null";
  }
  if ( $ENV{USER} =~ /greg/ ) {
    $rc = system("R --vanilla < $temp $nullify");

    #$rc = system("/ebi/research/software/Linux_x86_64/bin/R-2.7.0 --vanilla < $temp");
  } else {
    $rc = system("/software/bin/R-2.9.0 --vanilla < $temp $nullify");
  }

  if ( $params->{'cleanup_temp'} ) {
    use File::Path qw(rmtree);
    rmtree($dir);
  }
}

sub plotTree {
  my $class          = shift;
  my $file_out       = shift;
  my $params         = shift;
  my $tree_in        = shift;
  my $return_mid_cmd = shift;

  my $dba;
  if ($tree_in) {
    $dba = $tree_in->adaptor->db;
  } else {
    $dba = $params->{'dba'};
  }

  my $defaults = {
    max_tree_length                  => 0,
    alignment_quality_mask_character => 'O',                # Goes to dark gray.
    sequence_quality_mask_character  => 'B',                # Goes to lighter gray.
    gblocks_mask_character           => 'U',
    subtree_mask_character           => 'Z',
    remove_subtree                   => 0,
    remove_blank_columns             => 1,
    output_gene_info                 => 0,
    sitewise_table                   => 'sitewise_omega',
    exon_cased                       => 0,
    node_id                          => '',
    subtrees                         => '',
    cleanup_temp                     => 1
  };
  $params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params( $defaults, $params );

  my $param_set = $params->{'parameter_set'};
  if ($param_set) {
    print "Loading parameters from param_set $param_set...\n";
    my $param_set_params =
      Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set( $dba->dbc, $param_set );
    $params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params( $params, $param_set_params );
  }

  # Prune the tree.
  $params->{'node_id'} = $tree_in->node_id if ( defined $tree_in );
  my $node_id = $params->{'node_id'};

  # Define paths.
  my $dir = "/tmp/eslrplots";
  mkpath( [$dir] );
  my $tree_f = "$dir/tree_${node_id}.nh";
  my $aln_f  = "$dir/aln_${node_id}.fasta";

  my $url         = "";
  my $spf         = "";
  my $skip_output = 0;
  unless ($skip_output) {
    my ( $tree, $sa ) =
      Bio::EnsEMBL::Compara::ComparaUtils->get_tree_and_alignment( $dba, $params );

    Bio::EnsEMBL::Compara::AlignUtils->pretty_print( $sa, { length => 200 } );

    # Output the alignment.
    #my $sa = $tree->get_SimpleAlign;
    my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

    #$sa = Bio::EnsEMBL::Compara::ComparaUtils->fetch_masked_alignment($sa,$tree,$params);

    my $old_to_new;
    if ( $params->{'remove_blank_columns'} ) {
      my ( $flat_aln, $new_to_old, $o_t_n ) =
        @{ Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($sa) };
      $sa         = $flat_aln;
      $old_to_new = $o_t_n;
    }

    open( OUT, ">$aln_f" );
    foreach my $seq ( $sa->each_seq ) {
      print OUT ">" . $seq->id . "\n";
      print OUT $seq->seq . "\n";
    }
    close(OUT);

    # Output the tree.
    open( OUT, ">$tree_f" );
    my $tree2 = $tree->copy;
    $tree2 = $tree2->minimize_tree;
    print OUT $tree2->newick_format() . "\n";
    print $tree2->newick_format() . "\n";

    close(OUT);
  }

  my $max_tree_length = $params->{'max_tree_length'};

  my $out_cmd = qq^
    pdf(file="$file_out",width=width*10,height=height*10);
  ^;
  $out_cmd = qq^
    png(file="$file_out",width=width*10*72,height=height*10*72);
  ^ if ( $file_out =~ /\.png/i );

  my $aln_file_list  = "c(" . join( ", ", ("\"$aln_f\"") ) . ")";
  my $tree_file_list = "c(" . join( ", ", ("\"$tree_f\"") ) . ")";

  my $r_start = qq{
library(ape);
source("aln-tools/aln.tools.R");
source("aln-tools/slr.tools.R");
source("aln-tools/phylo.tools.R");
source("aln-tools/plot.phylo.greg.R");

print("$aln_f");
aln = read.aln("$aln_f");
aln\$tree = read.tree("$tree_f");

height = max(1,aln\$num_seqs / 10);
width = min(48,aln\$length/10);

mult = 50 / width;
width = width*mult;
height = height*mult;

dpar = par(no.readonly=T)
print(paste("Width:",width," Height:",height));
pdf(file="${file_out}",width=width,height=height,pointsize=12);
};

  my $r_mid = qq{
# Use layout to arrange the regions for SLR values, tracks, tree, and alignment.
aln = read.aln("$aln_f");
aln\$tree = read.tree("$tree_f");

num_tracks = 1
track_heights = c(.2)
aln_height = 1
heights=c(track_heights,aln_height)
tree_width=0.2
aln_width=1
pad_width = 1/50
widths=c(tree_width,aln_width,pad_width)

# Plotting order: aln, tree, track 1,2,3...n
num_rows = num_tracks+1
num_cols = 3
numbers = c()
for (i in 1:num_tracks) {
  numbers = c(numbers,3,i+3,0)
}
numbers = c(numbers,2,1,0)
print(numbers)
layout(matrix(numbers,nrow=num_rows,ncol=num_cols,byrow=T),heights=heights,widths=widths)
par(xaxs='i',mai=rep(0,4))

tree = aln\$tree
a = plot.aln(aln,overlay=F,axes=F,draw.chars=F,plot.tree=F)

xl = tree_length(tree);
if ($max_tree_length > 0) {
  low.lim = -$max_tree_length;
  hi.lim = xl*1.1;
} else {
  low.lim = -xl;
  hi.lim = xl*1.1;
}
plot.phylo.greg2(tree,y.lim=a\$ylim,x.lim=c(low.lim,hi.lim),show.tip.label=T)
axis.phylo(tree,side=3)
plot(c(1),type='n',axes=F) # The blank upper-left hand corner.

# Now plot the tracks.
plot.aln.bars(aln)

};

  my $r_end = qq{
dev.off()
q();
};

  if ($return_mid_cmd) {
    return $r_mid;
  }

  my $rcmd = "$r_start \n $r_mid \n $r_end \n";

  my $temp = "$dir/rcmd.txt";
  open( OUT, ">$temp" );
  print OUT $rcmd . "\n";
  close(OUT);

  # Create the PDF, then use ImageMagick to create a PNG and thumbnail.
  # cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  my $rc      = 0;
  my $nullify = "";
  if ( $params->{'quiet'} ) {
    $nullify = ">/dev/null";
  }
  if ( $ENV{USER} =~ /greg/ ) {
    $rc = system("R --vanilla < $temp $nullify");

    #$rc = system("/ebi/research/software/Linux_x86_64/bin/R-2.7.0 --vanilla < $temp");
  } else {
    $rc = system("/software/bin/R-2.9.0 --vanilla < $temp $nullify");
  }

  if ( $params->{'cleanup_temp'} ) {
    use File::Path qw(rmtree);
    rmtree($dir);
  }
}

sub plotTreeWithOmegas {
  my $class    = shift;
  my $file_out = shift;
  my $params   = shift;
  my $tree_in  = shift;

  my $skip_output = 0;
  my $dba         = $tree_in->adaptor->db;

  my $defaults = {
    alignment_quality_mask_character => 'O',                # Goes to dark gray.
    sequence_quality_mask_character  => 'B',                # Goes to lighter gray.
    gblocks_mask_character           => 'U',
    subtree_mask_character           => 'Z',
    remove_subtree                   => 0,
    remove_blank_columns             => 0,
    output_gene_info                 => 0,
    sitewise_table                   => 'sitewise_omega',
    exon_cased                       => 0,
    node_id                          => '',
    subtrees                         => '',
    cleanup_temp                     => 1
  };
  $params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params( $defaults, $params );
  my $param_set_params =
    Bio::EnsEMBL::Compara::ComparaUtils->load_params_from_param_set( $dba->dbc,
    $params->{'parameter_set_id'} );
  $params = Bio::EnsEMBL::Compara::ComparaUtils->replace_params( $params, $param_set_params );

  # Define paths.
  my $dir = "/tmp/eslrplots";
  mkpath( [$dir] );
  my $tree_f  = "$dir/tree.nh";
  my $aln_f   = "$dir/aln.fasta";
  my $omega_f = "$dir/omegas.txt";

  #unlink($omega_f) if (-e $omega_f);
  #unlink($tree_f) if (-e $tree_f);
  #unlink($aln_f) if (-e $aln_f);

  # Prune the tree.
  $params->{'node_id'} = $tree_in->node_id if ( defined $tree_in );

  my @subtrees = ();

  #  if ($params->{'subtrees'}) {
  #    foreach my $st (@{$params->{'subtrees'}}) {
  #      push @subtrees, $st;
  #    }
  #  } else {
  push @subtrees, {
    name  => "",
    table => "sitewise_aln",
    pset  => 1
    };

  #  }

  my $url = "";
  my $spf = "";
  unless ($skip_output) {
    my ( $tree, $sa ) =
      Bio::EnsEMBL::Compara::ComparaUtils->get_tree_and_alignment( $dba, $params );

    # Output the alignment.
    my $treeI = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($tree);

    my $old_to_new;
    if ( $params->{'remove_blank_columns'} ) {
      my ( $flat_aln, $new_to_old, $o_t_n ) =
        @{ Bio::EnsEMBL::Compara::AlignUtils->remove_blank_columns($sa) };
      $sa         = $flat_aln;
      $old_to_new = $o_t_n;
    }

    open( OUT, ">$aln_f" );
    foreach my $seq ( $sa->each_seq ) {
      print OUT ">" . $seq->id . "\n";
      print OUT $seq->seq . "\n";
    }
    close(OUT);

    # Output the tree.
    open( OUT, ">$tree_f" );
    my $tree2 = $tree->copy;
    $tree2 = $tree2->minimize_tree;
    print OUT $tree2->newick_format() . "\n";
    close(OUT);

    # Gather gene / protein info.
    my @human_leaves = grep { $_->taxon_id == 9606 } $tree->leaves;
    my $gene_obj = {};
    if ( scalar @human_leaves > 0 ) {
      my $leaf = $human_leaves[0];

      #my $gene = $leaf->get_Gene;
      #my $tscript = $leaf->get_Transcript;
      #my $tslation = $leaf->get_Translation;

      $gene_obj->{'Species'} = $leaf->taxon->genus . "_" . $leaf->taxon->species;

      #$gene_obj->{'Ensembl Gene'} = $gene->stable_id;
      #$gene_obj->{'Ensembl Peptide'} = $tslation->stable_id;
      #$gene_obj->{'Gene Name'} = $gene->external_name || $gene->description;
      #$url = "\"http://www.ensembl.org/Homo_sapiens/Gene/Compara_Tree?g=".$gene->stable_id."\"";
    }
    my $f = "";
    map { $f .= "%-20s: %s\\n" } keys( %{$gene_obj} );
    my $vals = "";
    $vals = join( ",", map { "\"$_\", \"" . $gene_obj->{$_} . "\"" } sort( keys( %{$gene_obj} ) ) );
    $vals = "\"\"" if ( $vals eq '' );
    $spf = "sprintf(\"$f\",$vals)";

    # Output omega values.

    foreach my $subtree (@subtrees) {
      my $table = $subtree->{'table'};
      my $pset  = $subtree->{'pset'};

      my $id = "${table}.${pset}";
      $omega_f = "$dir/$id.txt";
      open( OUT, ">$omega_f" );
      print OUT join( "\t", "site", "omega", "lower", "upper", "lrt_stat", "ncod", "note", "type" )
        . "\n";
      my $param_set_str = "";
      if ( exists $params->{'parameter_set_id'} ) {
        $param_set_str = "AND parameter_set_id=" . $params->{'parameter_set_id'};
      }
      my $node_id        = $tree->node_id;
      my $sitewise_table = $params->{'sitewise_table'};
      my $cmd            = qq^
   SELECT aln_position,omega,omega_lower,
          omega_upper,lrt_stat,ncod,note,type 
   FROM ${sitewise_table} WHERE 
   node_id=${node_id} ${param_set_str} order by aln_position;^;
      print $cmd. "\n";
      my $sth = $dba->dbc->prepare($cmd);
      $sth->execute();

      while ( my $arr = $sth->fetchrow_hashref ) {
        my $pos = $arr->{'aln_position'};
        $pos = $old_to_new->{$pos} if ( defined $old_to_new->{$pos} );
        my $omega    = $arr->{'omega'};
        my $lower    = $arr->{'omega_lower'};
        my $upper    = $arr->{'omega_upper'};
        my $lrt_stat = $arr->{'lrt_stat'};

        $lrt_stat *= -1 if ( $omega < 1 );

        my $ncod = $arr->{'ncod'} || "NA";
        my $note = $arr->{'note'} || "NA";
        my $type = $arr->{'type'} || "NA";
        print OUT join( "\t", ( $pos, $omega, $lower, $upper, $lrt_stat, $ncod, $note, $type ) )
          . "\n";

        #print join("\t",($pos,$omega,$lower,$upper,$ncod,$note,$type))."\n";
      }
    }

    close(OUT);
  }

  my $out_cmd = qq^
    pdf(file="$file_out",width=width*10,height=height*10);
  ^;
  $out_cmd = qq^
    png(file="$file_out",width=width*10*72,height=height*10*72);
  ^ if ( $file_out =~ /\.png/i );

  my $aln_file_list  = "c(" . join( ", ", ("\"$aln_f\"") ) . ")";
  my $tree_file_list = "c(" . join( ", ", ("\"$tree_f\"") ) . ")";
  my $slr_file_list  = "c(" . join( ", ", ("\"$omega_f\"") ) . ")";

  my $rcmd = qq{
    library(ape);
    source("aln-tools/aln.tools.R");
    source("aln-tools/phylo.tools.R");
    source("aln-tools/plot.phylo.greg.R");
    source("aln-tools/slr.tools.R");

    files = $aln_file_list;
    trees = $tree_file_list;
    slrs = $slr_file_list;
    alns = list();
    for (i in 1:length(files)) {
	aln = read.aln(files[i]);
	aln\$tree = read.tree(trees[i]);
	aln\$omegas = slrs[i];
	alns = append(alns,list(aln));
    }

    aln1 = alns[[1]];

    height = 5;
    width = 25;
    #height = max(1,aln1\$num_seqs / 10);
    #width = min(24,aln1\$length/10);

    mult = 50 / width;
    width = width*mult;
    height = height*mult;


dpar = par(no.readonly=T)

#par(cex=1/width);
#par(ps=width/10);

print(paste("Width:",width," Height:",height));
#library(Cairo);
#Cairo(file="${file_out}",width=width,height=height,type="svg",units="in");
pdf(file="${file_out}",width=width,height=height);

# Use layout to arrange the regions for SLR values, tracks, tree, and alignment.
track_heights = c(1,.1,.1,.1)
aln_height = 1
num_tracks = length(track_heights)
heights=c(track_heights,aln_height)
tree_width=0.2
aln_width=1
pad_width = aln_width/10
widths=c(tree_width,aln_width,pad_width)

# Plotting order: aln, tree, track 1,2,3...n
num_rows = num_tracks+1
num_cols = 3
numbers = c()
for (i in 1:num_tracks) {
  numbers = c(numbers,3,i+3,0)
}
numbers = c(numbers,2,1,0)

layout(matrix(numbers,nrow=num_rows,ncol=num_cols,byrow=T),heights=heights,widths=widths)
#layout.show(num_tracks+3)
par(xaxs='i',mai=rep(0,4))
source("aln-tools/aln.tools.R")
source("aln-tools/slr.tools.R")
source("aln-tools/plot.phylo.greg.R")
source("aln-tools/phylo.tools.R")
aln = read.aln(files[1])
tree = read.tree(trees[1])
slr = read.slr(slrs[1])

# Plot aln.
a = plot.aln(aln,overlay=F,axes=F,draw.chars=F)

# Plot tree.
xl = tree_length(tree)
plot.phylo.greg2(tree,y.lim=a\$ylim,x.lim=c(xl-2,xl*1.1),show.tip.label=T,edge.width=2)
axis.phylo(tree,side=3)
plot(c(1),type='n',axes=F,lwd=4)

# Plot SLR.
plot.slr.blocks(slr,xlim=a\$xlim,log=T)
axis(side=2,at=c(0.01,1,10))

# Plot tracks.
plot.slr.type(slr,xlim=a\$xlim,type='positive')
plot.slr.type(slr,xlim=a\$xlim,type='negative')
plot.aln.bars(aln)
dev.off()
q();
};

  my $temp = "$dir/rcmd.txt";
  open( OUT, ">$temp" );
  print OUT $rcmd . "\n";
  close(OUT);

  # Create the PDF, then use ImageMagick to create a PNG and thumbnail.
  # cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  my $rc      = 0;
  my $nullify = "";

  #my $nullify = ">/dev/null";
  if ( $ENV{USER} =~ /greg/ ) {
    $rc = system("R --vanilla < $temp");

    #$rc = system("/ebi/research/software/Linux_x86_64/bin/R-2.7.0 --vanilla < $temp");
  } else {
    $rc = system("/software/bin/R-2.9.0 --vanilla < $temp >/dev/null");
  }

  if ( $params->{'cleanup_temp'} ) {
    use File::Path qw(rmtree);
    rmtree($dir);
  }

}

sub plotDuplicationComparison {

  #     my $node_id = $tree->node_id;
  #     print "Creating alignment plots for $node_id...\n";

  #     my $dir = "$html_dir/node_id/$node_id";

  #     my @files = ();
  #     my @tree_files = ();
  #     my @slr_files = ();
  #     foreach my $f ("tree_a","tree_b","tree_full")
  #     {
  # 	my $aln_file = "$dir/${f}.fasta";
  # 	push @files, "\"$aln_file\"";
  # 	my $tree_file = "$dir/${f}.nh";
  # 	push @tree_files, "\"$tree_file\"";
  # 	my $slr_file = "$dir/${f}_omegas.txt";
  # 	push @slr_files, "\"$slr_file\"";

  #     }
  #     my $aln_file_list = "c(" . join(", ",@files) . ")";
  #     my $tree_file_list = "c(" . join(", ",@tree_files) . ")";
  #     my $slr_file_list = "c(" . join(", ",@slr_files) . ")";

  #     my $pdf_output = "$dir/plot.pdf";
  #     my $pdf_output_2 = "$dir/plot2.pdf";
  #     my $pdf_output_3 = "$dir/plot3.pdf";
  #     my $pdf_output_4 = "$dir/plot4.pdf";

  #     my $rcmd = <<"RCMD";
  #     library(ape);
  #     source("aln-tools/aln.tools.R");
  #     source("aln-tools/phylo.tools.R");
  #     source("aln-tools/plot.phylo.greg.R");
  #     source("aln-tools/slr.tools.R");

  #     files = $aln_file_list;
  #     trees = $tree_file_list;
  #     slrs = $slr_file_list;
  #     alns = list();
  #     for (i in 1:length(files)) {
  # 	aln = read.aln(files[i]);
  # 	aln\$tree = read.tree(trees[i]);
  # 	aln\$omegas = slrs[i];
  # 	alns = append(alns,list(aln));
  #     }

  #     aln1 = alns[[1]];
  #     aln2 = alns[[2]];
  #     aln_full = alns[[3]];

  #     height = max(2.5,aln1\$num_seqs / 50);
  #     width = min(16,aln1\$length/10);
  #     rowHeight = height/aln1\$num_seqs
  #     cex = rowHeight*5
  #     pdf(file="$pdf_output",width=width*10,height=height*10);

  #     omd = par("omd");
  #     pad = 0.01;
  #     alignment_height = 0.3;

  #     print(aln1\$length);
  #     print(aln2\$length);

  #     fit_to_track = function(
  # 	track_index,  # 1-based index.
  # 	num_tracks,   # total # of tracks.
  # 	height=1,     # height multiplier.
  # 	padding=0.1
  # 	) {
  # 	# Find the appropriate lower and upper bounds for this track.
  # 	track_height = 1/num_tracks;
  # 	# My actual height is decreased by the padding.
  # 	my_height = track_height - 2*(track_height*padding);

  # 	# Decrease height again by the height multiplier.
  # 	if (height < 1) my_height = track_height * height;

  # 	# Center the track in the "div" if the given height is less than the div height.
  # 	track_center_padding = (track_height - my_height) / 2;

  # 	lower_bound = (track_index-1)*(track_height);
  # 	upper_bound = track_index*(track_height);
  # 	low_y = lower_bound + track_center_padding;
  # 	hi_y = upper_bound - track_center_padding;

  # 	par(fig=c(0,1,low_y,hi_y),new=TRUE);
  #     }

  #     max_num_seqs = max(aln2\$num_seqs,aln1\$num_seqs);
  #     tree_width = 20;

  #     # Top alignment.
  #     par(new=FALSE);
  #     plot.new();
  #     fit_to_track(3,3,height=aln1\$num_seqs/max_num_seqs);
  #     plot.aln(aln1,square=FALSE,overlay=FALSE,tree.width=tree_width);

  #     # Bottom alignment.
  #     fit_to_track(1,3,height=aln2\$num_seqs/max_num_seqs);
  #     plot.aln(aln2,square=FALSE,overlay=FALSE,tree.width=tree_width);

  #     # Omegas.
  #     par(mar=c(0,0,0,0));
  #     fit_to_track(2,3,padding=0);
  #     plot.new();
  #     plot.window(xlim=c(0+0.5-tree_width,aln1\$length+0.5),ylim=c(0,2));
  #     plot.slr.blocks(aln1\$omegas,overlay=TRUE,block.xlim=c(.2,0.4));
  #     plot.slr.blocks(aln_full\$omegas,overlay=TRUE,block.xlim=c(0.4,0.6));
  #     plot.slr.blocks(aln2\$omegas,overlay=TRUE,block.xlim=c(0.6,.8));

#     o1 = read.table(aln1\$omegas,header=TRUE,sep="\t");
#     o2 = read.table(aln2\$omegas,header=TRUE,sep="\t");
#     o_full = read.table(aln_full\$omegas,header=TRUE,sep="\t");
#     for (i in 1:nrow(o_full)) {
#        lines(x=(i-.5)+c(0.3,0.5,0.7),y=c(o1[i,]\$omega,o_full[i,]\$omega,o2[i,]\$omega),col='black',lwd=5);
#     }

  #     dev.off();

#     pdf(file="$pdf_output_2");
#     par(new=FALSE);
#     good_rows = o2\$ncod > 4 & o1\$ncod > 4;
#     plot(x=o2[good_rows,]\$omega,y=o1[good_rows,]\$omega,xlim=c(0,2),ylim=c(0,2),main="",xlab="",ylab="");
#     title(main="Sitewise scatterplot of omegas",xlab="dup_b",ylab="dup_a");
#     abline(coef=c(0,1));
#     dev.off();

  #     pdf(file="$pdf_output_3");
  #     par(new=FALSE);
  #     diff = o2[good_rows,]\$omega - o1[good_rows,]\$omega;
  #     diff = diff[abs(diff) < 2]
  #     hist(diff,breaks=seq(from=-2,to=2,by=0.1),
  # 	 xlim=c(-2,2),
  # 	main="Distribution of sitewise omega differences",
  # 	xlab="(omega_b - omega_a)");
  #     abline(v=0,col='red',lwd=2);
  #     dev.off();

  #     pdf(file="$pdf_output_4");
  #     plot.phylo.greg(aln_full\$tree,plot=TRUE,max.pointsize=4);
  #     axis(1);
  #     dev.off();

  # RCMD
  #     ;

  #     my $temp = "temp.txt";
  #     open(OUT,">$temp");
  #     print OUT $rcmd."\n";
  # #    print $rcmd."\n";
  #     close(OUT);

  #     # Create the PDF, then use ImageMagick to create a PNG and thumbnail.
  #     # cmd to run: /software/R-2.7.1/bin/R CMD BATCH $filename
  #     my $rc = 0;
  #     if ($url =~ /127/) {
  # 	my $nullify = "";
  # 	#my $nullify = ">/dev/null";
  # 	$rc = system("/ebi/research/software/Linux_x86_64/bin/R-2.7.0 --vanilla < $temp $nullify");
  #     } else {
  # 	$rc = system("R --vanilla < $temp >/dev/null");
  #     }

  #     die "R returned an error!" if ($rc);
  #     unlink("$temp");

  #     print "Rasterizing PDFs ...\n";

  #     # Use ImageMagick to convert from PDF to PNG.
  #     my @pdf_files = ($pdf_output,$pdf_output_2,$pdf_output_3,$pdf_output_4);

  #     foreach my $file (@pdf_files) {
  # 	my ($name,$path,$suffix) = fileparse($file,(".pdf"));

  # 	my $cmd = "convert $file -thumbnail x300 ${path}/${name}_thumb.jpg";
  # 	$rc = system($cmd);

  # 	my $sa = $tree->get_SimpleAlign;
  # 	my $full_width = $img_full_width * ($sa->length + 150);
  # 	$full_width = 3000 if ($full_width < 3000);

  # #	$cmd = "convert $file -quality 80 -resize ${full_width}x5000 ${path}/${name}_full.jpg";
  # #	$rc = system($cmd);
  #     }

}

1;
