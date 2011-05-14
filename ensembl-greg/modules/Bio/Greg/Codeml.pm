# $Id: Codeml.pm 14743 2008-06-19 23:16:13Z jason $
#
# BioPerl module for Bio::Tools::Run::Phylo::PAML::Codeml
#
# Cared for by Jason Stajich <jason-at-bioperl-dot-org>
#
# Copyright Jason Stajich
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Tools::Run::Phylo::PAML::Codeml - Wrapper aroud the PAML program codeml

=head1 SYNOPSIS

  use Bio::Tools::Run::Phylo::PAML::Codeml;
  use Bio::AlignIO;

  my $alignio = Bio::AlignIO->new(-format => 'phylip',
  			         -file   => 't/data/gf-s85.phylip');

  my $aln = $alignio->next_aln;

  my $codeml = Bio::Tools::Run::Phylo::PAML::Codeml->new();
  $codeml->alignment($aln);
  my ($rc,$parser) = $codeml->run();
  my $result = $parser->next_result;
  my $MLmatrix = $result->get_MLmatrix();
  print "Ka = ", $MLmatrix->[0]->[1]->{'dN'},"\n";
  print "Ks = ", $MLmatrix->[0]->[1]->{'dS'},"\n";
  print "Ka/Ks = ", $MLmatrix->[0]->[1]->{'omega'},"\n";

=head1 DESCRIPTION

This is a wrapper around the codeml program of PAML (Phylogenetic
Analysis by Maximum Likelihood) package of Ziheng Yang.  See
http://abacus.gene.ucl.ac.uk/software/paml.html for more information.

This module is more about generating the properl codeml.ctl file and
will run the program in a separate temporary directory to avoid
creating temp files all over the place.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Jason Stajich

Email jason-at-bioperl-dot-org

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::Greg::Codeml;
use vars qw(@ISA %VALIDVALUES $MINNAMELEN $PROGRAMNAME $PROGRAM);
use strict;
use Bio::Root::Root;
use Bio::AlignIO;
use Bio::TreeIO;
use Bio::Tools::Run::WrapperBase;
use Bio::Tools::Phylo::PAML;
use Cwd;

@ISA = qw(Bio::Root::Root Bio::Tools::Run::WrapperBase);

=head2 Default Values

Valid and default values for codeml programs are listed below.  The
default values are always the first one listed.  These descriptions
are essentially lifted from the example codeml.ctl file and pamlDOC
documentation provided by the author.

B<CodonFreq> specifies the equilibrium codon frequencies in codon
substitution model. These frequencies can be assumed to be equal (1/61
each for the standard genetic code, B<CodonFreq> = 0), calculated from
the average nucleotide frequencies (B<CodonFreq> = 1), from the average
nucleotide frequencies at the three codon positions (B<CodonFreq> = 2),
or used as free parameters (B<CodonFreq> = 3). The number of parameters
involved in those models of codon frequencies is 0, 3, 9, and 60
(under the universal code), for B<CodonFreq> = 0, 1, 2, and 3
respectively.

B<aaDist> specifies whether equal amino acid distances are assumed (=
0) or Grantham's matrix is used (= 1) (Yang et al. 1998).

B<runmode> = -2 performs ML estimation of dS and dN in pairwise
comparisons. The program will collect estimates of dS and dN into the
files 2ML.dS and 2ML.dN. Since many users seem interested in looking
at dN /dS ratios among lineages, examination of the tree shapes
indicated by branch lengths calculated from the two rates may be
interesting although the analysis is ad hoc. If your species names
have no more than 10 characters, you can use the output distance
matrices as input to Phylip programs such as neighbor without
change. Otherwise you need to edit the files to cut the names short.

B<model> concerns assumptions about the dN/dS rate ratios among
branches (Yang 1998; Yang and Nielsen 1998). B<model> =0 means a single
dN/dS ratio for all lineages (branches), 1 means one ratio for each
branch (free ratio model), and 2 means arbitrary number of rations
(such as the 2-ratios or 3-ratios models. with B<model> =2, you may
specify the omega ratios for the branches using branch labels (read
about the tree structure file in the document).  This option seems
rather easy to use. Otherwise, the program will ask the user to input
a branch mark for the dN/dS ratio assumed for each branch. This should
be an integral number between 0 to k - 1 if k different dN/dS ratios
(omega_0 - omega_k - 1) are assumed for the branches of the
tree. B<Bioperl> note basically, doing this interactively is not going
to work very well, so this module is really focused around using the 0
or 1 parameters.  Read the program documentation if you'd like some more 
detailed instructions.

B<NSsites> specifies models that allow the dN/dS ratio (omega) to vary
among sites (Nielsen and Yang 1998, Yang et al. 2000) B<Nssites> = m
corresponds to model Mm in Yang et al (2000).  The variable B<ncatG>
is used to specify the number of categories in the omega distribution
under some models.  The values of ncatG() used to perform our
analyses are 3 for M3 (discrete), 5 for M4 (freq), 10 for the
continuous distributions (M5: gamma, M6: 2gamma, M7: beta, M8:beta&amp;w,
M9:beta&amp;gamma, M10: beta&gamma+1, M11:beta&amp;normal&gt;1, and
M12:0&amp;2normal&gt;1, M13:3normal&gt;0). This means M8 will have 11 site
classes (10 from the beta distribution plus 1 additional class). The
posterior probabilities for site classes as well as the expected omega
values for sites are listed in the file rst, which may be useful to
pinpoint sites under positive selection, if they exist. 

To make it easy to run several B<Nssites> models in one go, the
executable L<Bio::Tools::Run::Phylo::PAML::Codemlsites> can be used,
which asks you how many and which models to run at the start of the
program. The number of categories used will then match those used in
Yang et al(2000).

As noted in that paper, some of the models are hard to use, in
particular, M12 and M13. Recommended models are 0 (one-ratio), 1
(neutral), 2 (selection), 3 (discrete), 7 (beta), and 8
(beta&amp;omega ). Some of the models like M2 and M8 are noted to be
prone to the problem of multiple local optima. You are advised to run
the program at least twice, once with a starting omega value E<lt>1 and a
second time with a value E<gt>1, and use the results corresponding to the
highest likelihood. The continuous neutral and selection models of
Nielsen and Yang (1998) are not implemented in the program.


B<icode> for genetic code and these correspond to 1-11 in the genbank 
transl table.
  0:universal code
  1:mamalian mt
  2:yeast mt
  3:mold mt,
  4:invertebrate mt
  5:ciliate nuclear
  6:echinoderm mt
  7:euplotid mt
  8:alternative yeast nu.
  9:ascidian mt
  10:blepharisma nu

B<RateAncestor> For codon sequences, ancestral reconstruction is not
implemented for the models of variable dN/dS ratios among sites. The
output under codon-based models usually shows the encoded amino acid
for each codon. The output under "Prob of best character at each node,
listed by site" has two posterior probabilities for each node at each
codon (amino acid) site. The first is for the best codon. The second,
in parentheses, is for the most likely amino acid under the codon
substitution model. This is a sum of posterior probabilities across
synonymous codons. In theory it is possible although rare for the most
likely amino acid not to match the most likely codon.

B<Output> for codon sequences (seqtype = 1): The codon frequencies in
each sequence are counted and listed in a genetic code table, together
with their sums across species. Each table contains six or fewer
species. For data of multiple genes (option G in the sequence file),
codon frequencies in each gene (summed over species) are also
listed. The nucleotide distributions at the three codon positions are
also listed. The method of Nei and Gojobori (1986) is used to
calculate the number of synonymous substitutions per synonymous site
(dS ) and the number of nonsynonymous substitutions per nonsynonymous
site (dN ) and their ratio (dN /dS ). These are used to construct
initial estimates of branch lengths for the likelihood analysis but
are not MLEs themselves. Note that the estimates of these quantities
for the a- and b-globin genes shown in Table 2 of Goldman and Yang
(1994), calculated using the MEGA package (Kumar et al., 1993), are
not accurate.  

Results of ancestral reconstructions (B<RateAncestor> = 1) are collected 
in the file rst. Under models of variable dN/dS ratios among sites (NSsites models), 
the posterior probabilities for site classes as well as positively 
selected sites are listed in rst.

INCOMPLETE DOCUMENTATION OF ALL METHODS

=cut

BEGIN {

  $MINNAMELEN = 25;
  $PROGRAMNAME = 'codeml' . ( $^O =~ /mswin/i ? '.exe' : '' );
  if ( defined $ENV{'PAMLDIR'} ) {
    $PROGRAM =
      Bio::Root::IO->catfile( $ENV{'PAMLDIR'}, $PROGRAMNAME ) . ( $^O =~ /mswin/i ? '.exe' : '' );
  }

  # valid values for parameters, the default one is always
  # the first one in the array
  # much of the documentation here is lifted directly from the codeml.ctl
  # example file provided with the package
  %VALIDVALUES = (
    'outfile' => 'mlc',
    'noisy'   => [ 0 .. 3, 9 ],
    'verbose' => [ 1, 0, 2 ],     # 0:concise, 1:detailed, 2:too much

    # (runmode) 0:user tree, 1:semi-autmatic, 2:automatic
    #           3:stepwise addition, 4,5:PerturbationNNI
    #           -2:pairwise
    'runmode' => [ 0 .. 5, -2 ],

    'seqtype' => [ 1 .. 3 ],      # 1:codons, 2:AAs, 3:codons->AAs

    'CodonFreq' => [ 2, 0, 1, 3, 4, 5, 6, 7 ],    # 0:1/61 each, 1:F1X4,
                                                  # 2:F3X4, 3:codon table

    # (aaDist) 0:equal, +:geometric, -:linear,
    #          1-6:G1974,Miyata, c,p,v,a
    'aaDist' => [ 0, '+', '-', 1 .. 6 ],

    # (aaRatefile) only used for aa seqs
    # with model=empirical(_F)
    # default is usually 'wag.dat', also
    # dayhoff.dat, jones.dat, mtmam.dat, or your own
    'aaRatefile' => 'wag.dat',

    # (model) models for codons
    # 0: one, 1:b, 2:2 or more dN/dS ratios for branches
    'model' => [ 0 .. 2, 7 ],

    # (NSsites) number of S sites
    # 0: one w;1:neutral;2:selection; 3:discrete;4:freqs;
    # 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
    # 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
    # 13:3normal>0
    'NSsites' => [ 0 .. 13 ],

    # (icode) genetic code
    # 0:universal code
    # 1:mamalian mt
    # 2:yeast mt
    # 3:mold mt,
    # 4:invertebrate mt
    # 5:ciliate nuclear
    # 6:echinoderm mt
    # 7:euplotid mt
    # 8:alternative yeast nu.
    # 9:ascidian mt
    #10:blepharisma nu
    # these correspond to 1-11 in the genbank transl table

    'icode' => [ 0 .. 10 ],

    'Mgene' => [ 0, 1, 2, 3, 4 ]
    ,    # 0:rates, 1:separate analysis, 2: same (k,w) but diff pi and b.l.s
         # 3: same pi, diff. (k,w) and b.l.s, 4: different (k,w,), pis, and b.l.s
    'gene_codon_counts' => '',
    'fix_kappa'         => [ 0, 1 ],                # 0:estimate kappa, 1:fix kappa
    'kappa'             => '2',                     # initial or fixed kappa
    'fix_omega'         => [ 0, 1 ],                # 0: estimate omega, 1: fix omega
    'omega'             => '.3',                    # initial or fixed omega for
                                                    # codons or codon-base AAs
    'fix_alpha'         => [ 1, 0 ],                # 0: estimate gamma shape param
                                                    # 1: fix it at alpha
    'alpha'             => '0.',                    # initial or fixed alpha
                                                    # 0: infinity (constant rate)
    'Malpha'            => [ 0, 1 ],                # different alphas for genes
    'ncatG'             => [ 3, 1, 2, 4 .. 10 ],    # number of categories in
                                                    # dG of NSsites models

    # (clock)
    # 0: no clock, 1: global clock, 2: local clock
    # 3: TipDate
    'clock' => [ 0 .. 3 ],

    # (getSE) Standard Error:
    # 0:don't want them, 1: want S.E.
    'getSE' => [ 0, 1 ],

    # (RateAncestor)
    # 0,1,2 rates (alpha>0) or
    # ancestral states (1 or 2)
    'RateAncestor' => [ 1, 0, 2 ],
    'Small_Diff'   => '.5e-6',

    # (cleandata) remove sites with ambiguity data
    # 1: yes, 0:no
    'cleandata' => [ 0, 1 ],

    # this is the number of datasets in
    # the file - we would need to change
    # our api to allow >1 alignment object
    # to be referenced at time
    'ndata' => 1,

    # (method)
    # 0: simultaneous,1: 1 branch at a time
    'method' => [ 0, 1 ],

    # allow branch lengths to be fixed
    # 0 ignore
    # -1 use random starting points
    # 1 use the branch lengths in initial ML iteration
    # 2 branch lengths are fixed
    'fix_blength' => [ 0, -1, 1, 2 ],
  );
}

=head2 program_name

 Title   : program_name
 Usage   : $factory->program_name()
 Function: holds the program name
 Returns:  string
 Args    : None

=cut

sub program_name {
  return 'codeml';
}

=head2 program_dir

 Title   : program_dir
 Usage   : ->program_dir()
 Function: returns the program directory, obtiained from ENV variable.
 Returns:  string
 Args    :

=cut

sub program_dir {
  return Bio::Root::IO->catfile( $ENV{PAMLDIR} ) if $ENV{PAMLDIR};
}

=head2 new

 Title   : new
 Usage   : my $obj = Bio::Tools::Run::Phylo::PAML::Codeml->new();
 Function: Builds a new Bio::Tools::Run::Phylo::PAML::Codeml object 
 Returns : Bio::Tools::Run::Phylo::PAML::Codeml
 Args    : -alignment => the Bio::Align::AlignI object
           -save_tempfiles => boolean to save the generated tempfiles and
                              NOT cleanup after onesself (default FALSE)
           -tree => the Bio::Tree::TreeI object
           -branchlengths => 0: ignore any branch lengths found on the tree
                             1: use as initial values
                             2: fix branch lengths
           -params => a hashref of PAML parameters (all passed to set_parameter)
           -executable => where the codeml executable resides

See also: L<Bio::Tree::TreeI>, L<Bio::Align::AlignI>

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  $self->{'_branchLengths'} = 0;
  my ( $aln, $tree, $st, $params, $exe, $tempdir, $ubl ) = $self->_rearrange( [
      qw(ALIGNMENT TREE SAVE_TEMPFILES
        PARAMS EXECUTABLE TEMPDIR BRANCHLENGTHS)
    ],
    @args
  );
  defined $aln     && $self->alignment($aln);
  defined $tree    && $self->tree( $tree, branchLengths => ( $ubl || 0 ) );
  defined $st      && $self->save_tempfiles($st);
  defined $exe     && $self->executable($exe);
  defined $tempdir && $self->tempdir($tempdir);

  $self->set_default_parameters();
  if ( defined $params ) {
    if ( ref($params) !~ /HASH/i ) {
      $self->warn("Must provide a valid hash ref for parameter -FLAGS");
    } else {
      map { $self->set_parameter( $_, $$params{$_} ) } keys %$params;
    }
  }
  return $self;
}

=head2 prepare

 Title   : prepare
 Usage   : my $rundir = $codeml->prepare($aln);
 Function: prepare the codeml analysis using the default or updated parameters
           the alignment parameter must have been set
 Returns : value of rundir
 Args    : L<Bio::Align::AlignI> object,
	   L<Bio::Tree::TreeI> object [optional]

=cut

sub prepare {
  my ( $self, $aln, $tree ) = @_;

  unless ( $self->save_tempfiles ) {

    # brush so we don't get plaque buildup ;)
    $self->cleanup();
  }

  $tree = $self->tree      unless $tree;
  $aln  = $self->alignment unless $aln;
  my $tempdir = $self->tempdir();

  if ( !$aln ) {
    $self->warn("must have supplied a valid alignment file in order to run codeml");
    return 0;
  }
  my ( $tempseqFH, $tempseqfile );
  if ( !ref($aln) && -e $aln ) {
    $tempseqfile = $aln;
  } else {
    ( $tempseqFH, $tempseqfile ) = $self->io->tempfile(
      '-dir' => $tempdir,
      UNLINK => ( $self->save_tempfiles ? 0 : 1 )
    );
    #print "TEMP DIR: $tempdir\n";
    my $gene_codon_counts = $self->{_codemlparams}->{gene_codon_counts};
    delete $self->{_codemlparams}->{gene_codon_counts};
    my $alnout = Bio::AlignIO->new(
      '-format'      => 'phylip',
      '-fh'          => $tempseqFH,
      '-interleaved' => 0,
      '-idlength'    => $MINNAMELEN > $aln->maxdisplayname_length()
      ? $MINNAMELEN
      : $aln->maxdisplayname_length() + 1,
      '-paml_mgenes' => $gene_codon_counts
    );

    $alnout->write_aln($aln);
    $alnout->close();
    undef $alnout;
    close($tempseqFH);
  }

  # now let's print the codeml.ctl file.
  # many of the these programs are finicky about what the filename is
  # and won't even run without the properly named file.  Ack

  my $codeml_ctl = "$tempdir/codeml.ctl";
  open( CODEML, ">$codeml_ctl" ) or $self->throw("cannot open $codeml_ctl for writing");
  print CODEML "seqfile = $tempseqfile\n";
  my $outfile = $self->outfile_name;
  print CODEML "outfile = $outfile\n";

  if ($tree) {
    my ( $temptreeFH, $temptreefile );
    if ( !ref($tree) && -e $tree ) {
      $temptreefile = $tree;
    } else {
      ( $temptreeFH, $temptreefile ) = $self->io->tempfile(
        '-dir' => $tempdir,
        UNLINK => ( $self->save_tempfiles ? 0 : 1 )
      );

      my $treeout = Bio::TreeIO->new(
        '-format' => 'newick',
        '-fh'     => $temptreeFH
      );
      $treeout->write_tree($tree);
      $treeout->close();
      close($temptreeFH);
    }
    print CODEML "treefile = $temptreefile\n";
  }

  my %params = $self->get_parameters;
  while ( my ( $param, $val ) = each %params ) {
    next if $param eq 'outfile';
    print CODEML "$param = $val\n";
  }
  close(CODEML);

  # Print a copy of the grantham.dat if needed.
  if ($params{'aaDist'} == 1) {
    print STDERR "PRINTING OUTPUT TABLE!!!\n";
    $self->output_grantham_table;
  }

  return $tempdir;
}

=head2 run

 Title   : run
 Usage   : my ($rc,$parser) = $codeml->run($aln,$tree);
 Function: run the codeml analysis using the default or updated parameters
           the alignment parameter must have been set
 Returns : Return code, L<Bio::Tools::Phylo::PAML>
 Args    : L<Bio::Align::AlignI> object,
	   L<Bio::Tree::TreeI> object [optional]


=cut

sub run {
  my ($self)  = shift;
  my $tmpdir  = $self->prepare(@_);
  my $outfile = $self->outfile_name;

  my ( $rc, $parser ) = (1);
  {
    my $cwd = cwd();
    my $exit_status;
    chdir($tmpdir);
    my $codemlexe = $self->executable();
    $self->throw("unable to find or run executable for 'codeml'")
      unless $codemlexe && -e $codemlexe && -x _;
    my $run;
    if ( $self->{'_branchLengths'} ) {
      open( $run, "echo $self->{'_branchLengths'} | $codemlexe |" )
        or $self->throw("Cannot open exe $codemlexe");
    } else {
      open( $run, "$codemlexe |" ) or $self->throw("Cannot open exe $codemlexe");
    }
    my @output = <$run>;

    $exit_status = close($run);
    $self->error_string( join( '', @output ) );
    if ( ( grep { /error/io } @output ) || !$exit_status ) {

      if (!$exit_status) {
        $self->throw("Bad codeml exist status!\n". $self->error_string);
      } else {
        $self->throw( "ERROR RUNNING CODEML:\n" . $self->error_string );
        $rc = 0;
      }
    }

    # GJ 2009-01-08: Put the main results into a string and store it.
    open( IN, "$tmpdir/$outfile" );
    my @main_results_lines = <IN>;
    $self->main_results( \@main_results_lines );
    close(IN);

    eval {
      $parser = Bio::Tools::Phylo::PAML->new(
        -file    => "$tmpdir/$outfile",
        -verbose => $self->verbose,
        -dir     => "$tmpdir"
      );
    };
    warn() if $@;
    eval {

      #
      # GJ 2009-01-08 : Parse the supplementary results.
      #
      open( IN, "$tmpdir/rst" );
      my @supps = <IN>;
      $self->supplementary_results( \@supps );
      close(IN);
    };
    warn() if $@;

    chdir($cwd);
  }
  return ( $rc, $parser );
}

# GJ 2009-01-08 : store the supplementary results for later retrieval.
sub supplementary_results {
  my $self           = shift;
  my $suppl_arrayref = shift;

  $self->{'_suppl'} = $suppl_arrayref if ( defined $suppl_arrayref );
  return $self->{'_suppl'};
}

sub main_results {
  my $self          = shift;
  my $main_arrayref = shift;

  $self->{'_main'} = $main_arrayref if ( defined $main_arrayref );
  return $self->{'_main'};
}

sub all_lines {
  my $self = shift;
  
  my @main_lines = @{$self->main_results};
  my @supp_lines = @{$self->supplementary_results};
  
  my @all_lines = (@main_lines,@supp_lines);
  return \@all_lines;
}

sub parse_results {
  my $self = shift;
  my $aln = shift;
  my $pep_aln = shift;
  my $lines_arrayref = shift;

  my $results;

  my $bayes_results = $self->extract_empirical_bayes($lines_arrayref);
  $results = Bio::EnsEMBL::Compara::ComparaUtils->replace_params($results, $bayes_results);

  my $lnl = $self->extract_lnL($lines_arrayref);
  print "LNL: $lnl\n";

  return $results;
}

# GJ 2009-01-08 : extracts the naive bayes predictions from the supplementary string.
sub extract_empirical_bayes {
  my $self           = shift;
  my $suppl_arrayref = shift;

  $suppl_arrayref = $self->supplementary_results() unless ( defined $suppl_arrayref );
  my @suppl = @{$suppl_arrayref};

# Naive:
#Naive Empirical Bayes (NEB) probabilities for 10 classes& postmean_w
# 1 M   0.00320 0.01289 0.02424 0.03656 0.04965 0.06350 0.07828 0.09450 0.11369 0.52349 (10)  0.642
# Bayesian:
#Bayes Empirical Bayes (BEB) probabilities for 11 classes (class)& postmean_w
# 1 M   0.03403 0.08222 0.10935 0.11735 0.11455 0.10627 0.09566 0.08451 0.07399 0.07080 0.11128 ( 4)  0.745 +-  1.067
# NOTE FROM PAML DOCUMENTATION, p.29: "We suggest that you ignore the NEB output and use the BEB results only."

  my $results;
  my $bayes_results;
  my $pos_sites;

  my $site_class_omegas;

  my $has_bayes_section = 0;
  foreach my $line (@suppl) {
    #print "$line";
    $has_bayes_section = 1 if ( $line =~ m/Bayes Empirical Bayes/i );
  }

  #print "Has bayes: $has_bayes_section\n";
  my $in_supplement = 0;
  my $seen_bayes         = 0;
  my $number_of_bayes_classes = 0;
  my $seen_pos_sel_sites = 0;
  my $seen_naive         = 0;
  my $naive_output_has_p_gt_one_at_end = 0;

  foreach my $line (@suppl) {
    chomp($line);
    $line = strip($line);
    #print "$line\n";
    next if ( length($line) == 0 );
    next
      if ( $line =~ /amino acids/i )
      ;    # Skip lines like: (amino acids refer to 1st sequence: ENSDARP00000087283)

    # Wait until we get to the supplemental part
    $in_supplement = 1 if ( $line =~ m/supplemental results/i);
    next unless ($in_supplement);

    

    if ( $line =~ /w:\s+(.*)/ ) {
      my $stuff = $1;
      chomp $stuff;
      my @toks = split( /\s+/, $stuff );
      for ( my $i = 0 ; $i < scalar(@toks) ; $i++ ) {
        $site_class_omegas->{$i} = $toks[$i];
        #print "Site class:$i w:" . $toks[$i] . "\n";
      }
    }

    $seen_bayes++ if ( $line =~ m/Bayes Empirical Bayes/i );
    $seen_naive++ if ( $line =~ m/Naive Empirical Bayes/i );

    if ($seen_bayes && $line =~ m/probabilities for (\d+) classes/i) {
      $number_of_bayes_classes = $1;
      #print "  NUMBER OF CLASSES: $number_of_bayes_classes\n";
    }

    # If we have a bayes section, wait for that one.
    next if ($has_bayes_section && $seen_naive && !$seen_bayes);

    $seen_pos_sel_sites++ if ( $line =~ m/positively selected sites/i );

    #print "[$seen_bayes $seen_naive $seen_pos_sel_sites]";
    #print "$line\n";
    
    if ( !$has_bayes_section && $seen_naive && !$seen_pos_sel_sites ) {

      if ( $line =~ /empirical/i ) {
        if ( $line =~ /w>1/i ) {
          $naive_output_has_p_gt_one_at_end = 1;
        }
      }
      
      #Naive Empirical Bayes (NEB) probabilities for 10 classes& postmean_w
      # 1 M   0.00320 0.01289 0.02424 0.03656 0.04965 0.06350 0.07828 0.09450 0.11369 0.52349 (10)  0.642
      chomp $line;
      my @bits = split( /\s+/, $line );
      my $pos = shift @bits;
      my $residue = shift @bits;

      # Gets rid of the header line for the NEB output.
      next if ($line =~ m/w>1/i);
      
      my $prob_gt_one = 0;
      if ($naive_output_has_p_gt_one_at_end) {
        $prob_gt_one = pop @bits;
      } else {
        foreach my $site_class ( sort keys %{$site_class_omegas} ) {
          my $post_prob = $bits[ $site_class + 3 ];
          
          print "  Post prob: ${site_class} $post_prob\n";
          if ( $site_class_omegas->{$site_class} > 1 ) {
            $prob_gt_one += $post_prob;
          }
        }
      }

      my $omega = pop @bits;
      
      my $obj = {
        omega => $omega,
        prob_gt_one => $prob_gt_one,
        se => 0,
        is_positive_site => 0
      };
      $bayes_results->{$pos} = $obj;
    }

    if ( $has_bayes_section && $seen_bayes && !$seen_pos_sel_sites) {
      next if ( $line =~ /empirical/i );

      #Bayes Empirical Bayes (BEB) probabilities for 11 classes (class)& postmean_w
      # 1 M   0.03403 0.08222 0.10935 0.11735 0.11455 0.10627 0.09566 0.08451 0.07399 0.07080 0.11128 ( 4)  0.745 +-  1.067

      # GJ 2010-01-25 - PAML docs say nothing, but the the first 10 classes are probably 
      # equivalent to a 10-site class model, with the 11th position position being the
      # probability of w>1 integrated over all possible parameter values. So we'll take 
      # the last category as the p(w>1) and ignore the rest of the site classes.

      my @bits            = split( /\s+/, $line );

      # Shave off the beginning.
      my $pos = shift @bits;
      my $residue = shift @bits;
      # Shave off the end.
      my $se              = pop @bits;
      my $plus_minus_sign = pop @bits;
      my $omega           = pop @bits;

      my $prob_gt_one = 0;

      #print "$residue $pos $omega\n";

      # The rest of @bits is the probability of being in each site class.
      # NOTE: these are NOT necessarily the same number of site classes as chosen in the NSSites
      # config parameter. Apparently, PAML has a hard-coded 5 classes for M2 and 11 classes for M8.
      # With both situations, the last class is the only one with p(w)>1, so we'll use the
      # posterior of being in that class as the posterior of being under pos-sel.
      my $last_class_posterior = $bits[ $number_of_bayes_classes - 1 ];
      my $prob_gt_one = $last_class_posterior;

      my $obj = {
        omega => $omega,
        prob_gt_one => $prob_gt_one,
        se => $se,
        is_positive_site => 0
      };
      $bayes_results->{$pos} = $obj;
    }

    if ($seen_pos_sel_sites) {
      next if ( $line =~ m/positively/i );
      next if ( $line =~ m/Prob\(w>1\)/i );
      last if ( $line =~ m/lnL/i );
      last if ( $line =~ m/reconstruction/i );

      my @bits = split( /\s+/, $line );

      #Positively selected sites
      #     7 T      0.705         3.979 +- 3.258
      #    41 S      0.830         4.654 +- 3.213
      my $site = shift @bits;
      #print "$site\n";
      $bayes_results->{$site}->{is_positive_site} = 1;
    }
  }

  my $results = $bayes_results;
  return $results;
}

# GJ 2009-01-08 : Parse the tree from a block of CODEML-formatted text.
sub extract_tree {
  my $self          = shift;
  my $main_arrayref = shift;
  $main_arrayref = $self->main_results unless ( defined $main_arrayref );

  my @main = @{$main_arrayref};

  #print "@main\n";

  my $treeI;
  my $spotted_tree_count = 0;
  foreach my $line (@main) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blanks.

    # GJ 2009-01-08 : we should see 3 trees, and the third one will have the proper labels in them.

    my $line_contains_tree = ( $line =~ /\(.*?\)\;/ );
    $spotted_tree_count++ if ($line_contains_tree);

    if ( $spotted_tree_count == 3 && $line_contains_tree ) {
      $line =~ s/: /:/g;               # Remove spaces in the string.
      $line =~ s/, /,/g;

      #print "$line\n";
      open( my $fake_fh, "<", \$line );
      my $treein = new Bio::TreeIO( -fh => $fake_fh, -format => 'newick' );
      $treeI = $treein->next_tree;
      $treein->close();
    }
  }
  return $treeI;
}

# GJ 2009-01-09 Extracting log-likelihoods from runs.
sub extract_lnL {
  my $self          = shift;
  my $main_arrayref = shift;

  my @main = @{$main_arrayref};
  foreach my $line (@main) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blank lines.

    # Example line:
    # lnL(ntime:  8  np: 13):  -4820.123143     +0.000000
    if ( $line =~ /lnL(.*):\s+(\S*?)\s+(\S*?)/ ) {

      #print "$1 $2\n";
      my $lnL = $2;
      return $lnL;
    }
  }
  return undef;
}

sub parse_params {
  my $self = shift;
  my $main_arrayref = shift;

  my $seen_lnl = 0;
  my $n_past_lnl = 0;

  my @branches;
  my @branch_lengths;
  my @branch_standard_errors;
  my @params;
  my @standard_errors;

  my @main = @{$main_arrayref};
  foreach my $line (@main) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blank lines.
 
    if ($seen_lnl) {
      $n_past_lnl++;
    }
    if ($n_past_lnl == 1) {
      my $str = strip($line);
      @branches = split(/\s+/,$str);
    }
    if ($n_past_lnl == 2) {
      # We're in the parameters line.
      my @tokens = split(/\s+/,strip($line));
      my $divider = scalar(@branches) - 1;
      @branch_lengths = @tokens[0..$divider];
      my @param_tokens = @tokens[$divider+1..scalar(@tokens)-1];
      @params = @param_tokens;
    }
    if ($n_past_lnl == 4) {
      # We're in the parameters line.
      my @tokens = split(/\s+/,strip($line));
      my $divider = scalar(@branches) - 1;
      @branch_standard_errors = @tokens[0..$divider];
      my @param_tokens = @tokens[$divider+1..scalar(@tokens)-1];
      @standard_errors = @param_tokens;
    }
    
    $seen_lnl = 1 if ( $line =~ /lnL(.*):\s+(\S*?)\s+(\S*?)/ );
  }

  my $branch_map;
  for (my $i=0; $i < scalar(@branches); $i++) {
    my $key = $branches[$i];
    my $length = $branch_lengths[$i];
    my $se = $branch_standard_errors[$i];
    $branch_map->{$key} = [$length,$se];
  }

  return {
    branches => $branch_map,
    params => \@params,
    params_se => \@standard_errors
  };
}

# GJ 2009-01-09 Extracting log-likelihoods from runs.
sub extract_kappa {
  my $self          = shift;
  my $main_arrayref = shift;

  my @main = @{$main_arrayref};
  foreach my $line (@main) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blank lines.

    # Example line:
    # lnL(ntime:  8  np: 13):  -4820.123143     +0.000000
    if ( $line =~ /kappa.*=\s*(\S+?)\s*/ ) {
      my $kappa = $1;
      return $1;
    }
  }
  return undef;
}

# GJ 2009-01-09 Extracting log-likelihoods from runs.
sub extract_omegas {
  my $self          = shift;
  my $main_arrayref = shift;

  my @main = @{$main_arrayref};

  foreach my $line (@main) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blank lines.

    #print "$line\n";
    if ( $line =~ m/omega/ ) {
      my @tokens = split( "=", $line );
      my $value_string = $tokens[1];
      $value_string = strip($value_string);
      return ($value_string);
    } elsif ( $line =~ m/w \(dN\/dS\) for branches/ ) {

      #print "$line\n";
      my @tokens = split( ":", $line );

      #print "Tokens: @tokens\n";
      my $values_string = $tokens[1];
      $values_string = strip($values_string);
      my @values = split( /\s/, $values_string );

      #print "Omegas: @values\n";
      return @values;
    }
  }
  return undef;
}

sub each_line {
  my $self          = shift;
  my $main_arrayref = shift;
  $main_arrayref = $self->main_results unless ( defined $main_arrayref );
  my @main = @{$main_arrayref};
  return @main;
}

sub get_t_tree {
  my $self = shift;
  my $tree = shift;

  return $self->get_tree( $tree, 't' );
}

sub get_ds_tree {
  my $self = shift;
  my $tree = shift;

  return $self->get_tree( $tree, 'dS' );
}

sub get_dnds_tree {
  my $self = shift;
  my $tree = shift;

  return $self->get_tree( $tree, 'dN/dS' );
}

sub get_tree {
  my $self  = shift;
  my $tree  = shift;
  my $value = shift || 't';

  # Create a copy of the tree.
  my $newick   = Bio::EnsEMBL::Compara::TreeUtils->to_newick($tree);
  my $new_tree = Bio::EnsEMBL::Compara::TreeUtils->to_treeI($newick);

  my $hash = $self->parse_branch_params();

  #print join(",",keys(%$hash))."\n";

  foreach my $node ( $new_tree->get_leaf_nodes ) {
    my $node_id = $node->id;
    #print "id:$node_id\n";
    warn "No dnds id for " . $node_id . "\n" unless ( defined $hash->{ $node_id } );
    
    if (!defined $hash->{$node_id}) {
      # Try removing the branch model specs from the node ID.
      $node_id =~ s/#\d+//g;
      $node_id =~ s/\$\d+//g;
      $node->id($node_id);
    }

    return undef if ( !defined $hash->{ $node_id } );
    while ( defined $hash->{ $node->id } ) {
      #print "Node: ".$node->id."\n";
      $self->_set_branch_length( $node, $hash, $value );

      my $obj         = $hash->{ $node->id };
      my $parent_id   = $obj->{'parent'};
      my $parent_node = $node->ancestor;
      if ( defined $parent_node ) {
        $parent_node->id($parent_id);
        $node = $parent_node;
      } else {
        return $new_tree;
      }
    }
  }
  return $new_tree;
}

sub _set_branch_length {
  my $self         = shift;
  my $node         = shift;
  my $dnds_hash    = shift;
  my $value_to_set = shift;

  die "No branch length ID!" unless ( defined $dnds_hash->{ $node->id } );
  my $obj  = $dnds_hash->{ $node->id };
  my $dnds = $obj->{$value_to_set};
  $node->branch_length($dnds);
}

sub parse_branch_params {
  my $class = shift;
  my $line_arrayref = shift;

  my $param_objs = $class->parse_params($line_arrayref);
  my $branches = $param_objs->{branches};

  my @params = @{$param_objs->{params}};
  my @params_se = @{$param_objs->{params_se}};

  my $kappa = shift @params;
  my $kappa_se = shift @params_se;

  my $values;
  my $looking = 0;
  my $free_model = 0;
  foreach my $line ( @$line_arrayref ) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blank lines.

    if ( $line =~ m/dN & dS for each branch/i ) {
      $looking = 1;
      next;
    }
    if ($line =~ m/tree length for/i) {
      $looking = 0;
    }

    if ($line =~ m/Model: free dN\/dS Ratios for branches for branches/i) {
      $free_model = 1;
    }

    if ($looking) {
      if ( $line =~ m/\d\.\.\d/ ) {
        $line = strip($line);

        #	print $line."\n";
        # branch           t        N        S    dN/dS       dN       dS   N*dN   S*dS
        #   5..6       0.101   2374.1    973.9   0.2223   0.0166   0.0747   39.4   72.8
        my @tokens   = split( /\s+/,  $line );
        my @branches = split( /\.\./, $tokens[0] );
        my $child    = $branches[1];
        my $parent = $branches[0];
        my $id = $child;
        my $parent_id = $parent;

        my $branch_label = $tokens[0];
        my $branch_bl_arrayref = $branches->{$branch_label};
        my ($bl,$se) = @{$branch_bl_arrayref};

        my $dnds_se = -1;
        $dnds_se = shift @params_se if ($free_model);

        #	print "$child!!\n";
        my $obj = {
          't'      => $tokens[1],
          't_se'   => $se,
          'N'      => $tokens[2],
          'S'      => $tokens[3],
          'dN/dS'  => $tokens[4],
          'dN/dS_se' => $dnds_se,
          'dN'     => $tokens[5],
          'dS'     => $tokens[6],
          'node_a' => $parent,
          'node_b' => $child,
          'id' => $id,
          'branch_label' => $branch_label,
          'parent' => $parent_id
        };
        $values->{$id} = $obj;
      }
    }
  }
  return $values;
}

sub parse_number_tree {
  my $class = shift;
  my $input_lines_arrayref = shift;

  return $class->parse_tree($input_lines_arrayref,0);
}

sub parse_id_tree {
  my $class = shift;
  my $input_lines_arrayref = shift;

  return $class->parse_tree($input_lines_arrayref,1);
}

sub parse_tree {
  my $class = shift;
  my $input_lines_arrayref = shift;
  my $return_id_tree = shift; # 0 = return numbers, 1 = return IDs

  my $look_for_tree_one = 0;
  my $look_for_tree_two = 0;
  my $tree_one;
  my $tree_two;

  foreach my $line ( @$input_lines_arrayref ) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blank lines.
    if ( $line =~ m/^tree length =/i ) {
      $look_for_tree_one = 1;
      next;
    }
    if ($look_for_tree_one) {
      $line =~ s/ //g; # Remove paml's extra spaces.
      $tree_one          = Bio::TreeIO->new( -string => $line )->next_tree;
      $look_for_tree_one = 0;
      $look_for_tree_two = 1;
      next;
    }
    if ($look_for_tree_two) {
      $line =~ s/ //g; # Remove paml's extra spaces.
      $tree_two = Bio::TreeIO->new( -string => $line )->next_tree;
      $look_for_tree_two = 0;
    }
  }

  return unless ( defined $tree_one );

  # Get the internal node IDs from PAML's branch labels.
  my $child_to_parent_id;
  my $looking = 0;
  foreach my $line ( @$input_lines_arrayref) {
    chomp $line;
    next if ( length($line) == 0 );    # skip blank lines.
    if ( $line =~ m/dN & dS for each branch/i ) {
      $looking = 1;
      next;
    }
    if ($looking) {
      if ( $line =~ m/\d\.\.\d/ ) {
        $line = strip($line);
        # branch           t        N        S    dN/dS       dN       dS   N*dN   S*dS
        #   5..6       0.101   2374.1    973.9   0.2223   0.0166   0.0747   39.4   72.8
        my @tokens   = split( /\s+/,  $line );
        my @branches = split( /\.\./, $tokens[0] );
        my $child    = $branches[1];
        my $parent = $branches[0];
        $child_to_parent_id->{$child} = $parent;
      }
      if ($line =~ m/tree length/i) {
        $looking = 0; # Stop looking, we're done!
      }
    }
  }
  # Assign internal node IDs.
  foreach my $leaf ($tree_one->get_leaf_nodes) {
    my $node = $leaf;
    while (my $parent = $node->ancestor) {
      my $child_id = $node->id;
      if (defined $child_to_parent_id->{$child_id}) {
        $parent->id($child_to_parent_id->{$child_id});
      }
      $node = $parent;
    }
  }

  if ($return_id_tree) {
    return $tree_two;
  } else {
    return $tree_one;
  }
}

sub parse_leaf_number_map {
  my $class = shift;
  my $input_lines_arrayref = shift;

  my $map;

  my $number_tree = $class->parse_number_tree($input_lines_arrayref);
  my $id_tree = $class->parse_id_tree($input_lines_arrayref);
  
  # Tree one contains the PAML node numbers.
  my @number_nodes = $number_tree->get_nodes;
  # Tree two contains the leaf node names.
  my @id_nodes = $id_tree->get_nodes;
  for ( my $i = 0 ; $i < scalar(@number_nodes) ; $i++ ) {
    my $number = $number_nodes[$i];
    my $id = $id_nodes[$i];
    #print $number->id."  ".$id->id."\n";
    if ( $number->id && $id->id ) {
      $map->{ strip( $number->id ) } = strip( $id->id );
    }
  }
  return $map;
}

sub parse_branch_substitution_map {
  my $class = shift;
  my $lines_arrayref = shift;

  my $map = {};

  my $within_changes = 0;
  my $node_id;
  my $parent_id;
  my $branch_line;
  foreach my $line (@$lines_arrayref) {
    $within_changes = 1 if ($line =~ m/Summary of changes along branches/i);
    $within_changes = 0 if ($line =~ m/List of extant and reconstructed/i);

    if ($within_changes) {
      #print $line;
      # Branch 5:    9..5  (ENSGGOP00000012863)  (n= 2.0 s=10.0)
      if ($line =~ m/Branch (\d+):\s+(\d+)\.\.(\d+)/gi) {
        $branch_line = $line;
        $parent_id = $2;
        $node_id = $3;
        $map->{$node_id} = {};
      }
      if ($line =~ m/\s*(\d+) (\S+) \((\S)\) (.*) -> (\S+) \((\S)\)/g) {
        my $obj = {
          id => $node_id,
          parent_id => $parent_id,
          pos => $1,
          codon_a => $2,
          codon_b => $5,
          aa_a => $3,
          aa_b => $6,
          confidence => $4,
          line => $line,
          branch_line => $branch_line
        };
        next if ($obj->{codon_b} eq '---');

        $map->{$node_id}->{$1} = $obj;
      }
    }
  }
  return $map;
}

sub strip {
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
}

# GJ 2009-01-08 : Convenience method for re-estimating branch lengths using the M0 model.
sub get_m0_tree {
  my $class     = shift;
  my $tree      = shift;
  my $codon_aln = shift;

  # M0 parameters:
  my $params = {
    NSsites     => 2,
    fix_blength => 1
  };

  my $codeml = $class->new( -params => $params );
  $codeml->tree($tree);
  $codeml->alignment($codon_aln);
  $codeml->run();

  my $new_tree = $codeml->extract_tree();

  return $new_tree;
}

sub parse_codeml_results {
  my $class = shift;
  my $line_ref = shift; # Arrayref with ALL output lines (main + suppl)

  my $tree = $class->parse_number_tree($line_ref);

  my $subs = $class->parse_branch_substitution_map($line_ref);
  my $branch_params = $class->parse_branch_params($line_ref);
  my $map = $class->parse_leaf_number_map($line_ref);

  foreach my $node ($tree->get_nodes) {
    my $id = $node->id;
    my $subst = $subs->{$id};
    $node->set_tag_value('substitutions',$subst);

    my $params = $branch_params->{$id};
    foreach my $key (keys %$params) {
      $node->set_tag_value($key,$params->{$key});
    }

    if (defined $map->{$id}) {
      $node->id($map->{$id});
    }
  }
  
  return $tree;
}

sub codon_model_likelihood {
  my $class     = shift;
  my $tree      = shift;
  my $codon_aln = shift;
  my $tempdir   = shift;
  my $params    = shift;

  my $default_params = {
    NSsites     => 0,    # No sites categories.
    model       => 0,    # One rate throughout tree.
    fix_blength => 0,    # Ignore input branch lengths
  };

  my $final_params = $default_params;

  my $codeml = $class->new(
    -params    => $final_params,
    -tree      => $tree,
    -alignment => $codon_aln,
    -tempdir   => $tempdir
  );

  #$codeml->save_tempfiles(1);
  my ( $rs, $parser ) = $codeml->run();

  my $lines_ref = $codeml->all_lines;

  my $lnL    = $codeml->extract_lnL($lines_ref);
  my @omegas = $codeml->extract_omegas($lines_ref);

  return {
    lnL    => $lnL,
    omegas => \@omegas
  };
}

sub branch_model_likelihood {
  my $class     = shift;
  my $tree      = shift;
  my $codon_aln = shift;
  my $tempdir   = shift;
  my $params    = shift;

  my $default_params = {
    fix_blength => 1,    # Use initial branch lengths as estimates.
    model       => 0,
    cleandata   => 0,
    getSE => 0,
    fix_omega => 0,
    omega => 0.3,
    NSsites => 0,
  };

  my $has_any_fg_branches = 0;
  foreach my $node ($tree->nodes) {
    my $id = $node->id;
    if ($id =~ m/[\#\$]/) {
      $has_any_fg_branches = 1;
    }
  }
  if ($has_any_fg_branches) {
    $default_params->{model} = 2;
  }


  my @variable_params = qw(model omega Mgene gene_codon_counts cleandata aaDist NSsites fix_omega);

  # If certain parameters are given, apply them to the parameter object
  # passed to the new Codeml instance.
  my $final_params = {%$default_params};
  foreach my $param (@variable_params) {
    $final_params->{$param} = $params->{$param} if ( defined $params->{$param} );
  }
  $final_params->{verbose} = 1;

  # Create the new Codeml object and run it.
  my $codeml = $class->new(
    -params         => $final_params,
    -tree           => $tree,
    -alignment      => $codon_aln,
    -tempdir        => $tempdir,
    no_param_checks => 1
  );
  $codeml->save_tempfiles(1);

  # Note: ignore the $parser object here because it doesn't seem to work...
  my ( $rs, $parser ) = $codeml->run();
  my $lines_ref = $codeml->all_lines;
  my @lines = @$lines_ref;

  # Manually extract likelihood and omega values from the results files.
  my $lnL    = $class->extract_lnL($lines_ref);
  my $kappa = $class->extract_kappa($lines_ref);
  my @omegas = $class->extract_omegas($lines_ref);
  #print "Omegas: @omegas\n";
  my $codeml_tree = $class->parse_codeml_results($lines_ref);

  # Return an object with our results of interest.
  return {
    lnL       => $lnL,
    omegas    => \@omegas,
    tree => $codeml_tree,
    lines => $lines_ref
  };
}

sub file_array {
  my $file = shift;

  open( IN, $file );
  my @lines = <IN>;
  close(IN);
  return @lines;
}

# GJ 2009-01-09 : Ratio tests for pos-sel.
sub NSsites_ratio_test {
  my $class     = shift;
  my $tree      = shift;
  my $codon_aln = shift;
  my $model_a   = shift;
  my $model_b   = shift;
  my $tempdir   = shift;

  # Model A should be a model nested within B, i.e. A=m2 b=m3, A=m7 b=m8
  $model_a = 7 unless ( defined $model_a );
  $model_b = 8 unless ( defined $model_b );

  my $fix_omega = 0;
  my $omega     = 1.3;
  if ( $model_a eq '8a' ) {
    $model_a   = '8';
    $fix_omega = 1;
    $omega     = 1;
  }

  my $params = {
    NSsites     => $model_a,
    fix_blength => 1,
    fix_omega   => $fix_omega,
    omega       => $omega
  };
  my $codemla = $class->new(
    -params    => $params,
    -tree      => $tree,
    -alignment => $codon_aln,
    -tempdir   => $tempdir
  );
  $codemla->run();
  my $a_lines = $codemla->all_lines;
  my $ma_lnL = $class->extract_lnL($a_lines);

  $params = {
    NSsites     => $model_b,
    fix_blength => 1,
    fix_omega   => 0,
    omega       => 1.3
  };
  my $codemlb = $class->new(
    -params    => $params,
    -tree      => $tree,
    -alignment => $codon_aln,
    -tempdir   => $tempdir
  );
  $codemlb->run();
  my $b_lines = $codemlb->all_lines;
  my $mb_lnL = $class->extract_lnL($b_lines);

  if ( $model_b == 8 ) {

    # We should probably run M8 again here.
    # From Anisimova et al. 2002 (http://www.citeulike.org/user/gjuggler/article/1597691):
    #   "To avoid being trapped at a local optimum, it is important to run
    #   M8 at least twice, once with initial omega > 1 and once with
    #   omega < 1, and results corresponding to the highest likelihood
    #   value should be used."

    $params = {
      NSsites     => 8,
      fix_blength => 1,
      omega       => 0.3
    };
    my $codemlb2 = $class->new(
      -params    => $params,
      -tree      => $tree,
      -alignment => $codon_aln,
      -tempdir   => $tempdir
    );
    $codemlb2->run();
    my $b2_lines = $codemlb2->all_lines;
    my $mb_lnL2 = $codemlb2->extract_lnL($b2_lines);

    # Use the 2nd run's results if they have a higher likelihood.
    if ( $mb_lnL2 > $mb_lnL ) {
      $codemlb = $codemlb2;
      $mb_lnL  = $mb_lnL2;
    }
  }

  #print "model A: $ma_lnL  MB: $mb_lnL\n";

  my $twice_lnL = ( $mb_lnL - $ma_lnL ) * 2;

  my $chi_nf = 5.991;    # Chi2 value at 0.05 with 2df
  my $chi_nn = 9.210;    # Chi2 value at 0.01 with 2df

  return ( $twice_lnL, $codemla, $codemlb );
}

# GJ 2009-01-08 : Added unrooted tree method to Codeml.
sub _ensure_unrooted_tree {
  my $self = shift;
  my $tree = shift;

# GJ 2008-12-15 : need to reroot the tree at an internal node so that it has three branches from the root --
# PAML requires this for "unrooted" input.

  # GJ 2009-01-08 : TODO: Actually make this work by manually forcing a trifurcation.

  my $root             = $tree->get_root_node;
  my @children         = $root->each_Descendent;
  my $root_child_count = scalar(@children);
  if ( $root_child_count == 2 ) {

    #      print "CHILD COUNT: $root_child_count\n";
    my $child = $children[0];

    #	$tree->reroot($child);
  }
  return $tree;
}

sub output_grantham_table {
  my $self = shift;

  my $filename = $self->tempdir."grantham.dat";
  my $grantham = $self->get_grantham_table_string;
  open(OUT,">$filename");
  print OUT $grantham;
  close(OUT);
}

sub get_grantham_score_hash {
  my $self = shift;

  my $g = $self->get_grantham_table_string;

  my @number_rows = ();
  my @lines = split("\n",$g);
  foreach my $line (@lines) {
    last if (strip($line) eq '');
    push @number_rows,$line;
  }
  #print "number_rows:[@number_rows]\n";
  
  my $scores;

  my ($order_line) = grep {$_ =~ m/A\s+R\s+N/gi} @lines;
  my @residues = split(/\s+/,$order_line);
  @residues = grep {strip($_) ne ''} @residues;
  #print "residues: [@residues]\n";
  #return;
  
  my $i=0;
  my @residues_minus_one = @residues;
  shift @residues_minus_one;
  foreach my $residue (@residues_minus_one) {
    my $i_residue = $residue;
    $scores->{$i_residue.$i_residue} = 0;
    my $row = $number_rows[$i];
    #print "row:[$row]\n";
    my @numbers = split(/\s+/,$row);
    @numbers = grep {strip($_) ne ''} @numbers;
    #print "numbers: [@numbers]\n";
    my $j=0;
    foreach my $num (@numbers) {
      my $j_residue = $residues[$j];
      #print "i[$i_residue] j[$j_residue] score[$num] \n";
      $scores->{$j_residue.$j_residue} = 0;
      $scores->{$i_residue.$j_residue} = strip($num) + 0.0;
      $scores->{$j_residue.$i_residue} = strip($num) + 0.0;
      $j++;
    }
    $i++;
  }
  return $scores;
}

sub get_grantham_table_string {
  my $self = shift;

  my $grantham = qq^111  
111  85  
126  96  23  
195 180 139 154  
 91  43  46  61 154  
107  54  41  45 170  29  
 60 125  79  94 158  87  98  
 86  29  68  81 174  24  41  98  
 94  98 149 168 197 109 134 136  94  
 96 102 153 172 198 113 139 138  99   5  
106  26  94 102 202  53  57 127  32 102 107  
 85  92 141 160 196 101 126 127  86  10  14  95  
113  97 158 177 205 116 140 153 100  21  22 102  29  
 27 103  90 108 169  75  94  42  76  96  98 103  87 114  
 99 109  46  66 112  68  80  55  89 142 144 121 135 155  73  
 58  71  65  85 149  41  66  59  47  89  92  78  81 103  38  58  
148 101 174 191 215 130 152 184 115  61  61 110  67  40 147 177 129  
112  77 142 160 194  99 123 147  83  33  36  85  35  22 110 143  92  37  
 65  96 133 152 191  96 121 109  84  30  32  97  22  50  68 123  70  88  55  

 A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val


Table 1 in Grantham (1974).

	c	p	v      aromaticity (from Xuhua Xia)

Ala	0	8.1	31      -0.11    
Arg	.65	10.5	124     0.079    
Asn	1.33	11.6	56      -0.136   
Asp	1.38	13.0	54      -0.285   
Cys	2.75	5.5	55      -0.184   
Gln	.89	10.5	85      -0.067   
Glu	.92	12.3	83      -0.246   
Gly	.74	9.0	3       -0.073   
His	.58	10.4	96      0.32     
Ile	0	5.2	111     0.001    
Leu	0	4.9	111     -0.008   
Lys	.33	11.3	119     0.049    
Met	0	5.7	105     -0.041   
Phe	0	5.2	132     0.438    
Pro	.39	8.0	32.5    -0.016   
Ser	1.42	9.2	32      -0.153   
Thr	.71	8.6	61      -0.208   
Trp	.13	5.4	170     0.493    
Tyr	.20	6.2	136     0.381    
Val	0	5.9	84      -0.155   
^;
  return $grantham;
}


=head2 error_string

 Title   : error_string
 Usage   : $obj->error_string($newval)
 Function: Where the output from the last analysus run is stored.
 Returns : value of error_string
 Args    : newvalue (optional)


=cut

sub error_string {
  my ( $self, $value ) = @_;
  if ( defined $value ) {
    $self->{'error_string'} = $value;
  }
  return $self->{'error_string'};

}

=head2 alignment

 Title   : alignment
 Usage   : $codeml->align($aln);
 Function: Get/Set the L<Bio::Align::AlignI> object
 Returns : L<Bio::Align::AlignI> object
 Args    : [optional] L<Bio::Align::AlignI>
 Comment : We could potentially add support for running directly on a file
           but we shall keep it simple
 See also: L<Bio::SimpleAlign>

=cut

sub alignment {
  my ( $self, $aln ) = @_;

  if ( defined $aln ) {
    if ( -e $aln ) {
      $self->{'_alignment'} = $aln;
    } elsif ( !ref($aln) || !$aln->isa('Bio::Align::AlignI') ) {
      $self->warn(
        "Must specify a valid Bio::Align::AlignI object to the alignment function not $aln");
      return undef;
    } else {
      $self->{'_alignment'} = $aln;
    }
  }
  return $self->{'_alignment'};
}

=head2 tree

 Title   : tree
 Usage   : $codeml->tree($tree, %params);
 Function: Get/Set the L<Bio::Tree::TreeI> object
 Returns : L<Bio::Tree::TreeI> 
 Args    : [optional] $tree => L<Bio::Tree::TreeI>,
           [optional] %parameters => hash of tree-specific parameters:
                  branchLengths: 0, 1 or 2
                  out

 Comment : We could potentially add support for running directly on a file
           but we shall keep it simple
 See also: L<Bio::Tree::Tree>

=cut

sub tree {
  my ( $self, $tree, %params ) = @_;
  if ( defined $tree ) {
    if ( !ref($tree) || !$tree->isa('Bio::Tree::TreeI') ) {
      $self->warn("Must specify a valid Bio::Tree::TreeI object to the alignment function");
    }

    # GJ 2009-01-08 : Ensure the tree is unrooted.
    $tree = $self->_ensure_unrooted_tree($tree);

    $self->{'_tree'} = $tree;
    if ( defined $params{'_branchLengths'} ) {
      my $ubl = $params{'_branchLengths'};
      if ( $ubl !~ m/^(0|1|2)$/ ) {
        $self->throw(
          "The branchLengths parameter to tree() must be 0 (ignore), 1 (initial values) or 2 (fixed values) only"
        );
      }
      $self->{'_branchLengths'} = $ubl;
    }
  }
  return $self->{'_tree'};
}

=head2 get_parameters

 Title   : get_parameters
 Usage   : my %params = $self->get_parameters();
 Function: returns the list of parameters as a hash
 Returns : associative array keyed on parameter names
 Args    : none


=cut

sub get_parameters {
  my ($self) = @_;

  # we're returning a copy of this
  return %{ $self->{'_codemlparams'} };
}

=head2 set_parameter

 Title   : set_parameter
 Usage   : $codeml->set_parameter($param,$val);
 Function: Sets a codeml parameter, will be validated against
           the valid values as set in the %VALIDVALUES class variable.  
           The checks can be ignored if one turns off param checks like this:
             $codeml->no_param_checks(1)
 Returns : boolean if set was success, if verbose is set to -1
           then no warning will be reported
 Args    : $param => name of the parameter
           $value => value to set the parameter to
 See also: L<no_param_checks()>

=cut

sub set_parameter {
  my ( $self, $param, $value ) = @_;
  unless ( defined $self->{'no_param_checks'} && $self->{'no_param_checks'} == 1 ) {
    if ( !defined $VALIDVALUES{$param} ) {
      $self->warn(
        "unknown parameter $param will not be set unless you force by setting no_param_checks to true"
      );
      return 0;
    }
    if ( ref( $VALIDVALUES{$param} ) =~ /ARRAY/i
      && scalar @{ $VALIDVALUES{$param} } > 0 ) {

      unless ( grep { $value eq $_ } @{ $VALIDVALUES{$param} } ) {
        $self->warn(
          "parameter $param specified value $value is not recognized, please see the documentation and the code for this module or set the no_param_checks to a true value"
        );
        return 0;
      }
    }
  }
  $self->{'_codemlparams'}->{$param} = $value;
  return 1;
}

=head2 set_default_parameters

 Title   : set_default_parameters
 Usage   : $codeml->set_default_parameters(0);
 Function: (Re)set the default parameters from the defaults
           (the first value in each array in the 
	    %VALIDVALUES class variable)
 Returns : none
 Args    : boolean: keep existing parameter values


=cut

sub set_default_parameters {
  my ( $self, $keepold ) = @_;
  $keepold = 0 unless defined $keepold;

  while ( my ( $param, $val ) = each %VALIDVALUES ) {

    # skip if we want to keep old values and it is already set
    next if ( defined $self->{'_codemlparams'}->{$param} && $keepold );
    if ( ref($val) =~ /ARRAY/i ) {
      $self->{'_codemlparams'}->{$param} = $val->[0];
    } else {
      $self->{'_codemlparams'}->{$param} = $val;
    }
  }
}

=head1 Bio::Tools::Run::WrapperBase methods

=cut

=head2 no_param_checks

 Title   : no_param_checks
 Usage   : $obj->no_param_checks($newval)
 Function: Boolean flag as to whether or not we should
           trust the sanity checks for parameter values  
 Returns : value of no_param_checks
 Args    : newvalue (optional)


=cut

sub no_param_checks {
  my ( $self, $value ) = @_;
  if ( defined $value ) {
    $self->{'no_param_checks'} = $value;
  }
  return $self->{'no_param_checks'};
}

=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $obj->save_tempfiles($newval)
 Function: 
 Returns : value of save_tempfiles
 Args    : newvalue (optional)


=cut

=head2 outfile_name

 Title   : outfile_name
 Usage   : my $outfile = $codeml->outfile_name();
 Function: Get/Set the name of the output file for this run
           (if you wanted to do something special)
 Returns : string
 Args    : [optional] string to set value to


=cut

sub outfile_name {
  my $self = shift;
  if (@_) {
    return $self->{'_codemlparams'}->{'outfile'} = shift @_;
  }
  unless ( defined $self->{'_codemlparams'}->{'outfile'} ) {
    $self->{'_codemlparams'}->{'outfile'} = 'mlc';
  }
  return $self->{'_codemlparams'}->{'outfile'};
}

=head2 tempdir

 Title   : tempdir
 Usage   : my $tmpdir = $self->tempdir();
 Function: Retrieve a temporary directory name (which is created)
 Returns : string which is the name of the temporary directory
 Args    : none


=cut

=head2 cleanup

 Title   : cleanup
 Usage   : $codeml->cleanup();
 Function: Will cleanup the tempdir directory after a PAML run
 Returns : none
 Args    : none


=cut

=head2 io

 Title   : io
 Usage   : $obj->io($newval)
 Function:  Gets a L<Bio::Root::IO> object
 Returns : L<Bio::Root::IO>
 Args    : none


=cut

sub DESTROY {
  my $self = shift;
  unless ( $self->save_tempfiles ) {
    $self->cleanup();
  }
  $self->SUPER::DESTROY();
}

sub tempdir {
  my $self         = shift;
  my $new_temp_dir = shift;

  if ( defined $new_temp_dir ) {
    $self->{'_tempdir'} = $new_temp_dir;
  }
  if ( defined $self->{'_tempdir'} ) {
    return $self->{'_tempdir'};
  }

  my $tempdir = $self->SUPER::tempdir;
  $self->SUPER::tempdir( $tempdir . "_codeml" );

  return $self->SUPER::tempdir();
}

1;
