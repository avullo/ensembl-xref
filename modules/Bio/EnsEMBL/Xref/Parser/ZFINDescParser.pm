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


package Bio::EnsEMBL::Xref::Parser::ZFINDescParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( XrefParser::BaseParser );


sub run {
  my ( $self ) = @_;

  my $source_id  = $self->{source_id};
  my $species_id = $self->{species_id};
  my $files      = $self->{files};
  my $xref_dba   = $self->{xref_dba};
  my $verbose    = $self->{verbose} // 0;

  if ( (!defined $source_id) or (!defined $species_id) or (!defined $files) ) {
    croak "Need to pass source_id, species_id and files as pairs";
  }

  my $file = shift @{$files};

#e.g.
#ZDB-GENE-050102-6       WITHDRAWN:zgc:92147     WITHDRAWN:zgc:92147     0
#ZDB-GENE-060824-3       apobec1 complementation factor  a1cf    0
#ZDB-GENE-090212-1       alpha-2-macroglobulin-like      a2ml    15      ZDB-PUB-030703-1


  my $count = 0;
  my $withdrawn = 0;

  my $file_io = $xref_dba->get_filehandle($file);

  if ( !defined $file_io ) {
    croak "Can't open ZFIN file $file\n";
  }

  my $input_file = Text::CSV->new({
    sep_char       => "\t",
    empty_is_undef => 1,
    binary         => 1
  }) or croak "Cannot use file $file: " . Text::CSV->error_diag ();


  # 2 extra columns are ignored
  $input_file->column_names( [ 'zfin', 'desc', 'label'] );

  while ( my $data = $input_file->getline_hr( $file_io ) ) {
    # skip if WITHDRAWN: this precedes both desc and label
    if ( $data->{'label'} =~ /\A WITHDRAWN:/xms ) {
      $withdrawn++;
    }
    else {
      $xref_dba->add_xref({
        acc        => $data->{'zfin'},
        label      => $data->{'label'},
        desc       => $data->{'desc'},
        source_id  => $source_id,
        species_id => $species_id,
        info_type  => "MISC"
      });
      $count++;
    }
  }

  $input_file->eof or croak "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if($verbose){
    print "$count ZFIN xrefs added, $withdrawn withdrawn entries ignored\n";
  }

  return 0;
}

1;
