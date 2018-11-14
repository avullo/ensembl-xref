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
package Bio::EnsEMBL::Xref::Schema::Result::Species;

=head1 NAME

Bio::EnsEMBL::Xref::Schema::Result::Species

=cut

use strict;
use warnings;

use base 'DBIx::Class::Core';

=head1 TABLE: C<species>

=cut

__PACKAGE__->table("species");

=head1 ACCESSORS

=head2 species_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 taxonomy_id

  data_type: 'integer'
  extra: {unsigned => 1}
  is_nullable: 0

=head2 name

  data_type: 'varchar'
  is_nullable: 0
  size: 255

=head2 aliases

  data_type: 'varchar'
  is_nullable: 1
  size: 255

=cut

__PACKAGE__->add_columns(
  "species_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "taxonomy_id",
  { data_type => "integer", extra => { unsigned => 1 }, is_nullable => 0 },
  "name",
  { data_type => "varchar", is_nullable => 0, size => 255 },
  "aliases",
  { data_type => "varchar", is_nullable => 1, size => 255 },
);

=head1 UNIQUE CONSTRAINTS

=head2 C<species_taxonomy_idx>

=over 4

=item * L</species_id>

=item * L</taxonomy_id>

=back

=cut

__PACKAGE__->add_unique_constraint("species_taxonomy_idx", ["species_id", "taxonomy_id"]);

1;
