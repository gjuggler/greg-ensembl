#
# Ensembl module for Bio::EnsEMBL::Compara::DBSQL::Sitewise_dNdSAdaptor
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Compara::DBSQL::Sitewise_dNdSAdaptor - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Albert Vilella

This modules is part of the Ensembl project http://www.ensembl.org

Email avilella@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::Compara::DBSQL::Sitewise_dNdSAdaptor;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Compara::Sitewise_dNdS;
use Bio::EnsEMBL::Utils::Exception;


# my $conservation_scores = $self->_fetch_all_by_GenomicAlignBlockId_WindowSize($genomic_align_block->dbID, $window_size, $PACKED);

sub fetch_all_by_ProteinTreeId {
    my ($self, $protein_tree_id) = @_;
    my $sitewise_dnds_values = [];

    my $sql = qq{
  	SELECT
	    sitewise_id,
	    aln_position,
	    node_id,
	    omega,
	    omega_lower,
	    omega_upper,
	    type
	FROM
	    sitewise_aln
	WHERE
	    node_id = ?
	};

    my $sth = $self->prepare($sql);
    $sth->execute($protein_tree_id);

    my $sitewise_dnds;
    while (my ($sitewise_id,$aln_position,$node_id,
               $omega,$omega_lower,$omega_upper,
               $type) = $sth->fetchrow_array()) {
	$sitewise_dnds = Bio::EnsEMBL::Compara::Sitewise_dNdS->new_fast(
				       {'adaptor' => $self,
					'_dbID' => $sitewise_id,
					'aln_position' => $aln_position,
					'node_id' => $node_id,
					'omega' => $omega,
					'omega_lower' => $omega_lower,
					'omega_upper' => $omega_upper,
					'type' => $type});
	push(@$sitewise_dnds_values, $sitewise_dnds);
    }

    #sort into numerical order based on aln_position
    my @sorted_values = sort {$a->{aln_position} <=> $b->{aln_position}} @$sitewise_dnds_values;
    return \@sorted_values;
}
