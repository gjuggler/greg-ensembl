package Bio::EnsEMBL::Compara::TreeAlignUtils;

sub sort_aln_by_tree {
    my ($aln, $tree) = @_;

    my (@seq,@ids,%order);
    foreach my $seq ( $aln->each_seq() ) {
        push @seq, $seq;
        push @ids, $seq->display_id;
    }

    if (! ref $tree) {
	# It's not a ref, so it should be a filename.
        # Check that the file exists, and try to load it.
	my $infile = $tree;
	$aln->throw("Tree file not found: $infile") unless (-e $infile);
	my $temp = Bio::TreeIO->new('-file' => $infile,
				 '-format' => 'newick');
	$tree = $temp->next_tree;
	$temp->close();
    } else {
	# It's an object ref, so we assume it's a TreeI object.
	if (! $tree->isa('Bio::Tree::TreeI')) {
	    warn "sort_by_tree: given a hashref that is not a TreeI object!\n";
	}
    }

    my $ct=1;
    foreach my $leaf ($tree->get_leaf_nodes) {
	my $name = $leaf->id;
	#$self->throw("Not found in alignment: $name") unless $aln->_in_aln($name, \@ids);
	$order{$name}=$ct++;
    }
    
    # use the map-sort-map idiom:
    my @sorted= map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [$order{$_->id()}, $_] } @seq;
    my $aln = $aln->new;
    foreach (@sorted) { $aln->add_seq($_) }
    return $aln;
}

sub load_tree_from_string {
    my $tree_string = shift;
    open(my $fake_fh, "+<", \$tree_string);
    my $treein = new Bio::TreeIO
	(-fh => $fake_fh,
	 -format => 'newick');
    my $tree = $treein->next_tree;
    $treein->close;
    return $tree;
}

1;
