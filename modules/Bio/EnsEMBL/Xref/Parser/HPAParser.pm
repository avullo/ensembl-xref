
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

=head1 NAME

Bio::EnsEMBL::Xref::Parser::HPAParser

=head1 DESCRIPTION

This parser will read and creates direct xrefs from a simple comma-delimited file downloaded from the Human Protein Atlas (HPA) database.
The database contains two types of antibody, their own HPA antibodies and Collaborator antibody (CAB) commercial antibodies.

 data_uri        = http://www.proteinatlas.org/download/xref.php

The columns of the file should be the following:

 1)  Antibody
 2)  Antibody ID
 3)  Ensembl Peptide ID
 4)  Link (URL)

 Antibody,antibody_id,ensembl_peptide_id,link
 CAB000001,1,ENSP00000363822,http://www.proteinatlas.org/ENSG00000169083-AR
 CAB000001,1,ENSP00000379358,http://www.proteinatlas.org/ENSG00000169083-AR

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Xref::Parser::HPAParser->new(
    source_id  => 11,
    species_id => 9606,
    files      => ["hpa.txt"],
    xref_dba   => $xref_dba
  );

  $parser->run();

=cut

package Bio::EnsEMBL::Xref::Parser::HPAParser;

use strict;
use warnings;
use Carp;
use Text::CSV;


use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2
The run method does the actual parsing and creation of direct xrefs.
Parser gets initialized as noted above and run is called from
Bio::EnsEMBL::Production::Pipeline::Xrefs::ParseSource

my $parser = Bio::EnsEMBL::Xref::Parser::HPAParser-E<gt>new(..)
$parser-E<gt>run();

=cut

sub run {

    my ( $self, $ref_arg ) = @_;
    my $source_id  = $self->{source_id};
    my $species_id = $self->{species_id};
    my $files      = $self->{files};
    my $verbose    = $self->{verbose} // 0;
    my $xref_dba   = $self->{xref_dba};

    if (   ( !defined $source_id )
        or ( !defined $species_id )
        or ( !defined $files ) )
    {
        confess 'Need to pass source_id, species_id, and files';
    }

    my $file = shift @{$files};

    my $file_io = $xref_dba->get_filehandle($file);

    if ( !defined $file_io ) {
        confess "Could not open $file\n";
    }

    my $input_file = Text::CSV->new(
        {
            sep_char       => q{,},
            empty_is_undef => 1,
        }
    ) or confess "Cannot use file $file: " . Text::CSV->error_diag();

    my @expected_columns = qw(Antibody antibody_id ensembl_peptide_id link);
    my $header           = $input_file->getline($file_io);

    if ( scalar @{$header} != scalar @expected_columns ) {
        confess "input file $file has an incorrect number of columns";
    }

    $input_file->column_names(@expected_columns);
    my $parsed_count = 0;
    while ( my $data = $input_file->getline_hr($file_io) ) {

        $xref_dba->add_to_direct_xrefs(
            {
                acc        => $data->{'antibody_id'},
                version    => '1',
                label      => $data->{'Antibody'},
                stable_id  => $data->{'ensembl_peptide_id'},
                type       => 'translation',
                source_id  => $source_id,
                species_id => $species_id,
                info_type  => "DIRECT"
            }
        );

        $parsed_count++;
    }    #end while

    printf( "%d direct xrefs succesfully parsed\n", $parsed_count ) if ($verbose);

    $file_io->close();

    return 0;
} ## end sub run

1;
