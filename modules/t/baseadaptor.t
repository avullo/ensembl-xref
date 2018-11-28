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

# Check that the module loaded correctly
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

# Ensure that the BaseAdaptor handle is returned
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



# get_valid_xrefs_for_direct_xrefs



# label_to_acc



# get_valid_codes



# upload_xref_object_graphs
my $fake_xref = {
  accession => 'NM01234',
  version => 1,
  label => 'NM01234.1',
  description => 'Fake RefSeq transcript',
  species_id => '9606',
  info_type => 'DIRECT',
  info_text => 'These are normally aligned',
  dumped => 'NO_DUMP_ANOTHER_PRIORITY'
};

my @xref_array_00 = ($fake_xref);

throws_ok { $xref_dba->upload_xref_object_graphs() } qr/Please give me some xrefs to load/, 'Throws with no arguments';
throws_ok { $xref_dba->upload_xref_object_graphs( \@xref_array_00 ) } qr/Your xref does not have an accession-number/, 'Throws with bad accession';

# This auto-deploys the schema
# $db = Bio::EnsEMBL::Xref::Test::TestDB->new(
#   config_file => 'testdb.conf'
#   # config => {
#   #   driver => 'SQLite',
#   #   file => 'test.db',
#   #   create => 1
#   # }
# );

ok($db, 'TestDB ready to go');


# Add an example source to the db
my $source = $db->schema->resultset('Source')->create({
  name => 'RefSeq',
  status => 'KNOWN',
  source_release => '38',
  download => 'Y',
  priority => 1,
  priority_description => 'Like a boss',
  ordered => 10
});

ok(defined $source->source_id, 'Was the source created in the DB?');

my $xref = $source->create_related('xrefs', {
  accession => 'NM01234',
  version => 1,
  label => 'NM01234.1',
  description => 'Fake RefSeq transcript',
  species_id => '9606',
  info_type => 'DIRECT',
  info_text => 'These are normally aligned',
  dumped => 'NO_DUMP_ANOTHER_PRIORITY'
});

my $rs = $db->schema->resultset('Xref')->search(
  { accession => 'NM01234'}
);

my $matching_xref = $rs->next;
ok(defined $matching_xref,'A result was pulled from the DB');


# Test to ensure that the code fails if there is no SOURCE_ID defined
my $new_xref_00 = {
  ACCESSION => 'NM01235',
  VERSION => 1,
  LABEL => 'NM01235.1',
  DESCRIPTION => 'Fake RefSeq transcript',
  SPECIES_ID => '9606',
  INFO_TYPE => 'DIRECT',
  INFO_TEXT => 'These are normally aligned'
};

my @xref_array_01 = ($new_xref_00);

throws_ok { $xref_dba->upload_xref_object_graphs( \@xref_array_01 ) } qr/Your xref: NM01235 does not have a source-id/, 'Throws with no source ID';


# Test that a valid xref hashref gets added to the db
my $new_xref_01 = {
  ACCESSION => 'NM01235',
  VERSION => 1,
  LABEL => 'NM01235.1',
  DESCRIPTION => 'Fake RefSeq transcript',
  SPECIES_ID => '9606',
  SOURCE_ID => $source->source_id,
  INFO_TYPE => 'DIRECT',
  INFO_TEXT => 'These are normally aligned'
};

my @xref_array_02 = ( $new_xref_01 );

ok( !defined $xref_dba->upload_xref_object_graphs( \@xref_array_02 ), "NM01235 was added to the xref table" );

# get_xref
my $xref_id = $xref_dba->get_xref('NM01235', $source->source_id, 9606);
ok( defined $xref_id, "NM01235 (xref_id: $xref_id) was added to the xref table" );


# # upload_direct_xrefs
# throws_ok { $xref_dba->upload_direct_xrefs( \@xref_array_02 ) } qr/Problem Could not find accession/, 'Throws with bad accession';


# add_meta_pair



# get_xref_sources
my $xref_sources = $xref_dba->get_xref_sources();
ok( defined $xref_sources, 'There are sources in the db');


# species_id2taxonomy



# species_id2name



# get_xref_id



# get_xref - Already tested for other loading functions


# get_object_xref



# add_xref
# Test for xref that already exists
my $new_xref_02 = {
  acc => 'NM01235',
  version => 1,
  label => 'NM01235.1',
  desc => 'Fake RefSeq transcript',
  species_id => '9606',
  source_id => $source->source_id,
  info_type => 'DIRECT',
  info_text => 'These are normally aligned',
  update_label => 1,
  update_desc => 1
};
my $xref_id_old = $xref_dba->add_xref($new_xref_02);
ok( defined $xref_id_old, "NM01235 (xref_id: $xref_id_old) was loaded" );
ok($xref_id_old == $xref_id, "NM01235 already existed in the db");

# Test for a new xref
my $new_xref_03 = {
  acc => 'NM01236',
  version => 1,
  label => 'NM01236.1',
  desc => 'Fake RefSeq transcript',
  species_id => '9606',
  source_id => $source->source_id,
  info_type => 'DIRECT',
  info_text => 'These are normally aligned',
  update_label => 1,
  update_desc => 1
};
my $xref_id_new = $xref_dba->add_xref($new_xref_03);
ok( defined $xref_id_new, "NM01236 (xref_id: $xref_id_new) was loaded to the xref table" );


# add_object_xref
throws_ok { $xref_dba->add_object_xref() } qr/add_object_xref needs an xref_id/, 'Throws with no arguments';
throws_ok { $xref_dba->add_object_xref({xref_id => $xref_id}) } qr/add_object_xref needs an ensembl_id/, 'Throws with no arguments';
throws_ok { $xref_dba->add_object_xref({xref_id => $xref_id, ensembl_id => 1}) } qr/add_object_xref needs an object_type/, 'Throws with no arguments';

my $object_xref_id = $xref_dba->add_object_xref({xref_id => $xref_id, ensembl_id => 1, object_type => 'Gene'});
ok( defined $object_xref_id, "Object_xref entry inserted - $object_xref_id" );

done_testing();
