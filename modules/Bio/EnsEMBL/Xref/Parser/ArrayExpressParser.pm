=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Xref::Parser::ArrayExpressParser

=head1 DESCRIPTION

A parser class to parse the ArrayExpress source. Fetches the stable ids from core database
and populates the xref database's accession and label columns as DIRECT xref,
provided the species is active and declared in ArrayExpress.

-species = MULTI (Used by all ensembl species)

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::ArrayExpressParser->new(
    source_id  => 1,
    species_id => 9606,
    files      => ["project=>ensembl"],
    xref_dba   => $xref_dba,  # xref db adaptor
    dba        => $dba   # core db adaptor
  );

  $parser->run();
=cut


package Bio::EnsEMBL::Xref::Parser::ArrayExpressParser;

use strict;
use warnings;
use Carp;
use Bio::EnsEMBL::Registry;
use URI::ftp;

use parent qw( Bio::EnsEMBL::Xref::Parser );

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $species_name = $self->{species};
  my $files        = $self->{files};
  my $verbose      = $self->{verbose} // 0;
  my $xref_dba     = $self->{xref_dba};
  my $dba          = $self->{dba};

  defined $species_name or confess "Species name is required";

  my $species_lookup      = $self->_get_species($files, $verbose);
  my $active = $self->_is_active_species( $species_lookup, $species_name, $verbose );

  if ( !$active ) {
    return 0;
  }

  #get stable_ids from core and create xrefs
  my $gene_adaptor = $dba->get_GeneAdaptor();

  print "Finished loading the gene_adaptor\n" if $verbose;

  my @stable_ids = map { $_->stable_id } @{ $gene_adaptor->fetch_all() };

  my $xref_count = 0;
  foreach my $gene_stable_id (@stable_ids) {

    my $xref_id = $xref_dba->add_xref(
      {
        acc        => $gene_stable_id,
        label      => $gene_stable_id,
        source_id  => $source_id,
        species_id => $species_id,
        info_type  => "DIRECT"
      }
    );

    $xref_dba->add_direct_xref( $xref_id, $gene_stable_id, 'gene', '' );
    if ($xref_id) {
      $xref_count++;
    }
  }

  print "Added $xref_count DIRECT xrefs\n" if ($verbose);
  if ( !$xref_count ) {
   croak "No arrayexpress xref added\n";
  }

  return 0;      # successful

}

sub _get_species {
  my ( $self, $files, $verbose ) = @_;
  $verbose = ( defined $verbose ) ? $verbose : 0;

  my %species_lookup;
  foreach my $file (@$files) {
    my ($species) = split( /\./, $file );
    $species =~ s/^.*\///g;
    $species_lookup{$species} = 1;
  }
  return \%species_lookup;
}

# checks if the species is still active in ArrayExress
sub _is_active_species {
  my ( $self, $species_lookup, $species_name, $verbose ) = @_;

  if ( $species_lookup->{$species_name} ) {
    printf( 'Found ArrayExpress has declared the name "%s". This was an alias' . "\n", $species_name ) if $verbose;
    return 1;
  }
  return 0;
}

1;
