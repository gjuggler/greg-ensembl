use Bio::TreeIO;

my $tree_file = "/nfs/users/nfs_g/gj1/src/greg-ensembl/projects/slrsim/trees/mammals44tree.nh";
my $treeio = Bio::TreeIO->new(-file => $tree_file);
my $tree = $treeio->next_tree;
map {$_->id(lc($_->id))} $tree->leaves;

my @eutherian = qw(alpaca armadillo bushbaby cat chimp cow dog dolphin elephant guinea_pig hedgehog horse human kangaroo_rat megabat microbat mouse mouse_lemur pika rabbit rat rhesus rock_hyrax shrew sloth squirrel tarsier tenrec treeshrew);  
my @full = qw(chicken chimp cow dog fugu guinea_pig horse human lizard monodelphis mouse platypus rat rhesus stickleback tetraodon zebrafish zebra_finch);
my @full_mammal = qw(chimp cow dog guinea_pig horse human mouse rat rhesus);
my @hmrd = qw(human mouse rat dog);

my @primates = qw(human chimp gorilla orangutan rhesus marmoset tarsier mouse_lemur bushbaby);
my @glires = qw(guinea_pig mouse rat kangaroo_rat squirrel rabbit pika);
my @laur = qw(alpaca dolphin cow horse cat dog microbat megabat hedgehog shrew);

my @glires_short = qw(mouse rat kangaroo_rat squirrel guinea_pig);
my @laur_short = qw(cow cat microbat megabat);

print "Primates:".$tree->slice_by_ids(@primates)->total_branch_length."\n";
print "Glires:".$tree->slice_by_ids(@glires)->total_branch_length."\n";
print "Glires short:".$tree->slice_by_ids(@glires_short)->total_branch_length."\n";
print "Laur:".$tree->slice_by_ids(@laur)->total_branch_length."\n";
print "Laur short:".$tree->slice_by_ids(@laur_short)->total_branch_length."\n";


my $tree_file = "/nfs/users/nfs_g/gj1/src/greg-ensembl/ensembl-compara-63/scripts/pipeline/species_tree_blength.nh";
my $treeio = Bio::TreeIO->new(-file => $tree_file);
my $tree = $treeio->next_tree;

print $tree->ascii(1,1,1)."\n";
