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

use utf8;
package Bio::EnsEMBL::Xref::Schema::Result::Source;

=head1 NAME

Bio::EnsEMBL::Xref::Schema::Result::Source

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<source>

=cut

__PACKAGE__->table("source");

=head1 ACCESSORS

=head2 source_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_auto_increment: 1
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 status

  data_type: 'enum'
  default_value: 'NOIDEA'
  extra: {list => ["KNOWN","XREF","PRED","ORTH","PSEUDO","LOWEVIDENCE","NOIDEA"]}
  is_nullable: 0

=head2 source_release

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=head2 download

  data_type: 'enum'
  default_value: 'Y'
  extra: {list => ["Y","N"]}
  is_nullable: 1

=head2 ordered

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 priority

  data_type: 'integer'
  default_value: 1
  extra: {unsigned => 1}
  is_nullable: 1

=head2 priority_description

  data_type: 'varchar'
  default_value: (empty string)
  is_nullable: 1
  size: 40

=cut

__PACKAGE__->add_columns(
  "source_id",
  {
    data_type => "integer",
    extra => { unsigned => 1 },
    is_auto_increment => 1,
    is_nullable => 0,
  },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "status",
  {
    data_type => "enum",
    default_value => "NOIDEA",
    extra => {
      list => ["KNOWN", "XREF", "PRED", "ORTH", "PSEUDO", "LOWEVIDENCE", "NOIDEA"],
    },
    is_nullable => 0,
  },
  "source_release",
  { data_type => "varchar", is_nullable => 1, size => 255 },
  "download",
  {
    data_type => "enum",
    default_value => "Y",
    extra => { list => ["Y", "N"] },
    is_nullable => 1,
  },
  "ordered",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "priority",
  {
    data_type => "integer",
    default_value => 1,
    extra => { unsigned => 1 },
    is_nullable => 1,
  },
  "priority_description",
  { data_type => "varchar", is_nullable => 1, size => 40 },
);

=head1 PRIMARY KEY

=over 4

=item * L</source_id>

=back

=cut

__PACKAGE__->set_primary_key("source_id");

__PACKAGE__->has_many('xrefs', 'Bio::EnsEMBL::Xref::Schema::Result::Xref', 'source_id');
1;
