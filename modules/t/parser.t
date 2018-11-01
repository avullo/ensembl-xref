=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

use Test::More;
use Test::Exception;

use_ok 'Bio::EnsEMBL::Xref::Parser';

throws_ok { Bio::EnsEMBL::Xref::Parser->new() } qr/Need to pass/, 'Throws with no arguments';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1) } qr/Need to pass/, 'Throws with not enough arguments';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2) } qr/Need to pass/, 'Throws with not enough arguments';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => 'dummy') }
  qr/check it is a reference/, 'Throws with scalar files arg';
throws_ok { Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => {'dummy'}) }
  qr/was expected/, 'Throws with non arrayref files arg';

my $parser = Bio::EnsEMBL::Xref::Parser->new(source_id => 1, species_id => 2, files => ['dummy']);
isa_ok($parser, 'Bio::EnsEMBL::Xref::Parser');

throws_ok { $parser->run() } qr/abstract method/, 'Throws calling abstract method';

done_testing();

