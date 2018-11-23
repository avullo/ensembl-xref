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
use Bio::EnsEMBL::Xref::FetchFiles;
use URI::ftp;
use List::Compare;

use parent qw( Bio::EnsEMBL::Xref::Parser );

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id             = $self->{source_id};
  my $species_id            = $self->{species_id};
  my $override_species_name = $self->{species};
  my $files                 = $self->{files};
  my $verbose               = $self->{verbose} // 0;
  my $xref_dba              = $self->{xref_dba};
  # my $core_dba              = $self->{dba};

  my $file = shift @{$files};

  # project could be ensembl or ensemblgenomes
  my $project;
  if ( $file =~ /project[=][>](\S+?)[,]/ ) {
    $project = $1;
  }

  # If a species is not listed on the ArrayExpress FTP site by any of
  # our aliases, we cannot load xrefs for this species
  my $ref_species_name = $self->is_active_species($species_id, $xref_dba, $override_species_name);
  if (! $ref_species_name ) {
    return;
  }
  
  #get stable_ids from core and create xrefs
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($ref_species_name,'core','Gene');

  print "Finished loading the gene_adaptor\n" if $verbose;

  my @stable_ids = map { $_->stable_id } @{ $gene_adaptor->fetch_all() };
  # FIXME: this would be a lot faster if we extended BaseFeatureAdaptor to
  # fetchall_stable_ids

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

=head2 is_active_species

Arg 1      :  species_id - the ID from species table in the Xref schema
Arg 2      :  xref_dba - DBAdaptor to the Xref database
Arg 3      :  species_name - An alternate name for the species to match the remote site

Description:  Determine whether the species (species_id) is available for this source 
              by consulting its FTP site for species-named files. Can supply a 
              custom species name as override to allow flexibility

Returntype :  String - The most production species name for Ensembl

=cut

sub is_active_species {
  my ($self, $species_id, $xref_dba, $species_name) = @_;
  
  # collect Ensembl names and species aliases for our current species
  my $species_record = $xref_dba->get_species_particulars($species_id);
  my $aliases = [$species_record->{aliases}];
  push @$aliases, $species_name if defined $species_name;

  my $client = Bio::EnsEMBL::Xref::FetchFiles->new();
  # File names end .ensgene.tsv or similar. Thes rest is a latin species name
  my $file_names = $client->list_ftp_files(
    'ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/bioentity_properties/ensembl/',
    qr/\..+/
  );
  
  my $comparison = List::Compare->new( $file_names, $aliases );
  my @match = $comparison->get_intersection;
  
  if (scalar @match > 0) {
    return $species_record->{name};
  } else {
    return 0;
  }
}

1;
