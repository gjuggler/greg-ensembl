# Author: jkh1
# 2006-08-18
#

=head1 NAME

 Treefam::FamilyHandle

=head1 SYNOPSIS

 use Treefam::DBConnection;

 my $dbc = new Treefam::DBConnection ();

 my $famh = $dbc->get_FamilyHandle();

 my $family = $famh->get_by_id('TF300011');

=head1 DESCRIPTION

 Enables retrieval of Treefam family objects.

=head1 CONTACT

 jkh1@sanger.ac.uk


=cut


package Treefam::FamilyHandle;


use strict;
use Carp;
use Treefam::Family;
use Treefam::DBConnection;
use Scalar::Util qw(weaken);


=head2 new

 Arg: Treefam::DBConnection
 Description: Creates a new family object handle
 Returntype: Treefam::FamilyHandle

=cut

sub new {

  my $class = shift;
  my $self = {};
  $self->{'DBConnection'} = shift;
  weaken($self->{'DBConnection'});
  bless ($self, $class);

  return $self;

}

=head2 get_by_id

 Arg: Treefam family database ID
 Description: Gets family with given Treefam database ID
 Returntype: Treefam::Family object

=cut

sub get_by_id {

  my ($self,$familyID) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $type;
  if ($dbc->get_Database eq 'treefam_3' || $dbc->get_Database eq 'treefam_4') {
    if    ($familyID =~ /^TF1/) { $type = "familyA";}
    elsif ($familyID =~ /^TF3/ || $familyID =~ /^TF5/) { $type = "familyB";}
    else { croak "\nERROR: unknown type for $familyID.\n";}
  }
  else { # treefam-5, treefam-6, etc.
    if    ($familyID =~ /^TF1/) { $type = "familyA";}
    elsif ($familyID =~ /^TF3/) { $type = "familyB";}
    elsif ($familyID =~ /^TF5/) { $type = "familyC";}
    else { croak "\nERROR: unknown type for $familyID.\n";}
  }

  my $query = qq( SELECT COUNT(*) FROM $type WHERE AC= ? );
  my $sth= $dbh->prepare ($query);
  # check if family exists in database
  $sth->execute($familyID);
  my ($count) = $sth->fetchrow_array();
  unless ($count) {
    return undef;
  }
  return new Treefam::Family($dbc,$familyID,$type);
}

=head2 get_by_ac

 Arg: Treefam family database ID
 Description: A synonym for get_by_id. Gets family with
              given Treefam database ID
 Returntype: Treefam::Family object

=cut

sub get_by_ac {

  my ($self,$familyID) = @_;
  return $self->get_by_id($familyID);

}


=head2 get_by_symbol

 Arg: family symbol
 Description: Gets family with given symbol. Only curated
              families have symbols.
 Returntype: Treefam::Family object

=cut

sub get_by_symbol {

  my ($self,$symbol) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $sth = $dbh->prepare("SELECT AC FROM familyA
                           WHERE symbol=?");
  $sth->execute($symbol);
  my ($familyID) = $sth->fetchrow_array();
  $sth->finish();

  return undef if( !defined $familyID );

  return $self->get_by_id($familyID);

}

=head2 get_by_gene

 Arg: Treefam::Gene object
 Description: Gets family containing given gene.
 Returntype: Treefam::Family object

=cut

sub get_by_gene {

  my ($self,$gene) = @_;
  my $geneID = ref($gene) ? $gene->ID() : $gene;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  if ($dbc->get_Database eq 'treefam_3' || $dbc->get_Database eq 'treefam_4') {
    foreach my $type('famB_gene','famA_gene') { # treefam_3, treefam_4
      # get family with the highest HMMER score
      my $query = qq( SELECT DISTINCT t.AC
                      FROM hmmer h, $type t
                      LEFT JOIN genes g ON t.ID=g.ID
                      WHERE g.GID= ?
                      AND h.ID=g.ID
                      AND t.AC<'TF500000'
                      ORDER BY h.SCORE DESC);
      my $sth = $dbh->prepare($query);
      $sth->execute($geneID);
      my ($familyID) = $sth->fetchrow_array();
      $sth->finish();
      return $self->get_by_id($familyID) if ($familyID);
    }
  }
  else { # treefam-5, treefam-6, etc.
    # get family with the highest HMMER score
    my $query = qq( SELECT DISTINCT t.AC
                    FROM hmmer_matches h, fam_genes t
                    LEFT JOIN genes g ON t.ID=g.ID
                    WHERE g.GID= ?
                    AND h.ID=g.ID
                    AND t.FAM_TYPE IN ('A','B')
                    ORDER BY h.SCORE DESC );
    my $sth = $dbh->prepare($query);
    $sth->execute($geneID);
    my ($familyID) = $sth->fetchrow_array();
    $sth->finish();
    return $self->get_by_id($familyID) if ($familyID);
  }
  return undef;
}


=head2 get_all_by_type

 Arg: string, type of families: A or B
 Description: Gets all families of given type
 Returntype: list of Treefam::Family objects

=cut


sub get_all_by_type {

  my ($self,$type) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  $type = ($type eq 'A') ? 'familyA' : 'familyB';
  my $query = qq( SELECT DISTINCT AC FROM $type );
  my $sth = $dbh->prepare($query);
  $sth->execute();
  my @families;
  while ( my ($familyID) = $sth->fetchrow_array()) {
    push @families,$self->get_by_id($familyID);
  }
  $sth->finish();
  return @families;

}

=head2 get_all_by_species

 Arg1: string, species (Swissprot code or latin name)
 Arg2: optional, string, type of families (A or B)
 Description: Gets all families that have genes of the
              given species
 Returntype: list of Treefam::Family objects

=cut

sub get_all_by_species {

  my ($self,$species,$type) = @_;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $query;
  my @families;
  if ($dbc->get_Database eq 'treefam_3' || $dbc->get_Database eq 'treefam_4') {
    if ($type) {
      if ($type eq 'A') {
	$query = qq(SELECT DISTINCT t.AC FROM famA_gene t, genes g, species sp
                    WHERE t.ID = g.ID
                    AND  g.TAX_ID = sp.TAX_ID
                    AND (sp.SWCODE = ? OR sp.TAXNAME = ?) );
      }
      elsif ($type eq 'B') {
	$query = qq(SELECT DISTINCT t.AC FROM famB_gene t, genes g, species sp
                    WHERE t.ID = g.ID
                    AND  g.TAX_ID = sp.TAX_ID
                    AND (sp.SWCODE = ? OR sp.TAXNAME = ?) );
      }
      my $sth = $dbh->prepare($query);
      $sth->execute($species,$species);
      while ( my ($familyID) = $sth->fetchrow_array()) {
	push @families,$self->get_by_id($familyID);
      }
      $sth->finish();
    }
    else {
      foreach my $type('famB_gene','famA_gene') {
	$query = qq( SELECT DISTINCT t1.AC FROM $type t1, genes g, species s
	      	     WHERE t1.ID = g.ID
                     AND  g.TAX_ID = s.TAX_ID
                     AND (s.SWCODE = ? OR s.TAXNAME = ?) );
	my $sth = $dbh->prepare($query);
	$sth->execute($species,$species);
	while ( my ($familyID) = $sth->fetchrow_array()) {
	  push @families,$self->get_by_id($familyID);
	}
	$sth->finish();
      }
    }
  }
  else {
    if ($type) {
      $query = qq( SELECT DISTINCT t.AC FROM fam_genes t, genes g, species sp
                   WHERE t.ID = g.ID
                   AND  g.TAX_ID = sp.TAX_ID
                   AND t.FAM_TYPE = '$type'
                   AND (sp.SWCODE = ? OR sp.TAXNAME = ?) );
    }
    else {
      $query = qq( SELECT DISTINCT t1.AC FROM fam_genes t1, genes g, species s
                   WHERE t1.ID = g.ID
                   AND  g.TAX_ID = s.TAX_ID
                   AND t1.FAM_TYPE IN ('A','B')
                   AND (s.SWCODE = ? OR s. TAXNAME = ?) );
    }
    my $sth = $dbh->prepare($query);
    $sth->execute($species,$species);
    while ( my ($familyID) = $sth->fetchrow_array()) {
      push @families,$self->get_by_id($familyID);
    }
    $sth->finish();
    return @families;
  }
}

=head2 get_all_by_domain

 Arg1: string, pfam id of domain
 Arg2: (optional) e-value cut-off
 Description: Gets all families that have genes with the
              given domain with e-value below given cut-off
              (default is 1e-2).
 Returntype: list of Treefam::Family objects

=cut

sub get_all_by_domain {

  my $self = shift;
  my $pfamid = shift;
  unless ($pfamid) {
    croak "Pfam domain id required";
  }
  my $cutoff = shift if @_;
  if (!$cutoff) {
    $cutoff = 1e-2;
  }
  my @families;
  my $dbc = $self->{'DBConnection'};
  my $dbh = $dbc->{'database_handle'};
  my $query = qq(SELECT DISTINCT t1.AC FROM fam_genes t1, pfam p
                 WHERE t1.ID = p.ID
                 AND p.PFAM_ID = ?
                 AND t1.FAM_TYPE IN ('A','B') );
  my $sth = $dbh->prepare($query);
  $sth->execute($pfamid);
  while (my ($familyID)=$sth->fetchrow_array()) {
    push (@families,$self->get_by_id($familyID));
  }
  $sth->finish();

  return @families;

}

=head2 get_all_by_domain_list

 Arg: list of strings, pfam ids of domains
 Description: Gets all families with genes that have
              all the given domains. No e-value cut-off.
 Returntype: list of Treefam::Gene objects

=cut

sub get_all_by_domain_list {

  my $self = shift;
  my @pfamids = @_;
  unless ($pfamids[0]) {
    croak "Pfam domain ids required";
  }
  my @genes;
  my %seen;
  # remove duplicates
  my @uniq_pfamids = grep { !$seen{$_}++ } @pfamids;
  %seen = ();
  foreach my $pfamid(@uniq_pfamids) {
    map { $seen{$_->ID}++ if $_ } $self->get_all_by_domain($pfamid,10000);
  }
  my $n = scalar(@uniq_pfamids);
  my @IDs = grep { $seen{$_}==$n } keys %seen;
  while (my $ID=shift @IDs) {
    push (@genes,$self->get_by_id($ID));
  }

  return @genes;

}


1;
