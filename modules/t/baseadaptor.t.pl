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

use strict;
use warnings;

use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Xref::Test::TestDB;
use Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor;

use_ok 'Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor';

my $db = Bio::EnsEMBL::Xref::Test::TestDB->new();
my %config = %{ $db->config };

my $xref_dba = Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor->new(
  host   => $config{host},
  dbname => $config{db},
  user   => $config{user},
  pass   => $config{pass},
  port   => $config{port}
);

ok( defined $xref_dba, 'BaseAdaptor handle returned');

# Test get_filehandle
throws_ok { $xref_dba->get_filehandle() } qr/No file name/, 'Throws with no arguments';
throws_ok { $xref_dba->get_filehandle('fake_file.tsv') } qr/Can not open file/, 'Throws fake file name';


# get_source_id_for_source_name
throws_ok { $xref_dba->get_source_id_for_source_name() } qr/source_name undefined/, 'Throws with no arguments';
throws_ok { $xref_dba->get_source_id_for_source_name('fake_name') } qr/No source_id/, 'Throws with no matching source_id';


# get_source_ids_for_source_name_pattern
throws_ok { $xref_dba->get_source_ids_for_source_name_pattern() } qr/source_name undefined/, 'Throws with no arguments';
ok( defined $xref_dba->get_source_ids_for_source_name_pattern('fake_name'), 'Array returned' );

# get_source_name_for_source_id
throws_ok { $xref_dba->get_source_name_for_source_id() } qr/source_id undefined/, 'Throws with no arguments';
is( $xref_dba->get_source_name_for_source_id('fake_name'), -1, 'source not present in db' );

# get_valid_xrefs_for_dependencies


done_testing();
