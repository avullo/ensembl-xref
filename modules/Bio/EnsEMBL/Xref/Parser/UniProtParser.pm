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



package Bio::EnsEMBL::Xref::Parser::UniProtParser;

use strict;
use warnings;

use Carp;
use Module::Load;
use Readonly;

use parent qw( Bio::EnsEMBL::Xref::Parser );


Readonly my $DEFAULT_LOADER_BATCH_SIZE         => 1000;
Readonly my $DEFAULT_LOADER_CHECKPOINT_SECONDS => 300;

Readonly my %source_name_for_section => (
  'Swiss-Prot' => 'Uniprot/SWISSPROT',
  'TrEMBL'     => 'Uniprot/SPTREMBL',
);


=head2 run
The run method does the actual parsing and creation of direct xrefs.
Parser gets initialized as noted above and run is called from
Bio::EnsEMBL::Production::Pipeline::Xrefs::ParseSource

my $parser = Bio::EnsEMBL::Xref::Parser::UniProtParser-E<gt>new(..)
$parser-E<gt>run();

=cut

sub run {
  my ( $self ) = @_;

  my $general_source_id = $self->{source_id};
  my $species_id        = $self->{species_id};
  my $species_name      = $self->{species};
  my $files             = $self->{files};
  my $release_file      = $self->{rel_file};
  my $verbose           = $self->{verbose} // 0;
  my $xref_dba          = $self->{xref_dba};

  my $loader_batch_size
    = $self->{loader_batch_size} // $DEFAULT_LOADER_BATCH_SIZE;
  my $loader_checkpoint_seconds
    = $self->{loader_checkpoint_seconds} // $DEFAULT_LOADER_CHECKPOINT_SECONDS;

  # Try to control where ETL modules can come from, just in case
  # someone does something really weird with configuration options.
  my $extractor_class   = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';
  my $transformer_class = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';
  my $loader_class      = 'Bio::EnsEMBL::Xref::Parser::UniProtParser::';

  $extractor_class   .= $self->{extractor_class}   // 'Extractor';
  $transformer_class .= $self->{transformer_class} // 'Transformer';
  $loader_class      .= $self->{loader_class}      // 'Loader';

  load $extractor_class;
  load $transformer_class;
  load $loader_class;

  my $extractor = $extractor_class->new({
    'file_names' => $files,
    'species_id' => $species_id,
    'xref_dba'   => $xref_dba,
  });
  my $transformer = $transformer_class->new({
    'species_id' => $species_id,
    'xref_dba'   => $xref_dba,
  });
  my $loader = $loader_class->new({
    'batch_size'         => $loader_batch_size,
    'checkpoint_seconds' => $loader_checkpoint_seconds,
    'xref_dba'           => $xref_dba,
  });

  # Generate a map of existing dependent-xref links, which will be
  # used to prevent insertion of duplicates
  my $source_id_map = $transformer->get_source_id_map();
  while ( my ( $section, $pri_ref ) = each %{ $source_id_map } ) {
    while ( my ( $priority, $source_id ) = each %{ $pri_ref } ) {

      if ( $priority ne 'direct' ) {
        $loader->prepare_source_for_dependent_xrefs( $source_id );
      }

      # Seeing as we are already looping over all sources.
      # FIXME: do we really care about the file name here?
      if ( $verbose ) {
        print "$section $priority source id for $files->[0]: $source_id\n";
      }

    }
  }

 RECORD:
  while ( $extractor->get_uniprot_record() ) {

    my $extracted_record = $extractor->extract();
    if ( $extracted_record eq 'SKIP' ) {
      next RECORD;
    }

    my $transformed_data
      = $transformer->transform( $extracted_record );

    $loader->load( $transformed_data );

  }

  $extractor->finish();
  $transformer->finish();
  $loader->finish();

  # Extract release numbers from the release file, if provided
  if ( defined $release_file ) {
    my $release_numbers = $self->_get_release_numbers_from_file( $release_file,
                                                                 $verbose );
    $self->_set_release_numbers_on_uniprot_sources( $source_id_map,
                                                    $release_numbers );
  }

  return 0;
} ## end sub run


=head2 _get_release_numbers_from_file
  Arg [1]    : release_file_name
  Arg [2]    : verbose
  Description: Extract Swiss-Prot and TrEMBL release info from the release file

=cut

sub _get_release_numbers_from_file {
  my ( $self, $release_file_name, $verbose ) = @_;

  my $xref_dba = $self->{xref_dba};

  my $release_io = $xref_dba->get_filehandle( $release_file_name );
  my $release_numbers = {};

  while ( my $line = $release_io->getline() ) {
    my ( $section, $release )
      = ( $line =~ m{
                      \A
                      UniProtKB/
                      (
                        Swiss-Prot
                      |
                        TrEMBL
                      )
                      \s+
                      Release
                      \s+
                      ( [^\n]+ )
                  }msx );
    if ( defined $section ) {
      $release_numbers->{ $source_name_for_section{$section} } = $release;
      if ( $verbose ) {
        print "$section release is '$release'\n";
      }
    }
  }

  $release_io->close();

  return $release_numbers;
} ## end sub _get_release_numbers_from_file


=head2
  Arg [1]    : Source ID Map
  Arg [2]    : Release number
  Description: Set the release number for the UniProt Source

=cut

sub _set_release_numbers_on_uniprot_sources {
  my ( $self, $source_id_map, $release_numbers ) = @_;

  my $xref_dba = $self->{xref_dba};

  foreach my $source ( keys %{ $source_id_map } ) {
    # Priority names are not important here, we only need source IDs
    foreach my $source_id ( values %{ $source_id_map->{$source} } ) {
      $xref_dba->set_release( $source_id, $release_numbers->{$source} );
    }
  }

  return;
} ## end sub _set_release_numbers_on_uniprot_sources

1;
