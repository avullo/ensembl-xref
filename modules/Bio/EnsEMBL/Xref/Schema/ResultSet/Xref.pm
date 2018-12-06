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

package Bio::EnsEMBL::Xref::Schema::ResultSet::Xref;
use strict;
use warnings;

use parent 'DBIx::Class::ResultSet';


=head2 check_direct_xref
  Arg [1]    : HashRef of parameters
                 $params->{accession},
                 $params->{display_label},
                 $params->{description}
  Description:

=cut

sub check_direct_xref {
  my ($self,$params) = @_;

  my $hit = $self->find($params);
  # {
  #   accession => $params->{accession},
  #   label => $params->{display_label},
  #   description => $params->{description}
  # }
  return 1 if defined $hit;
  return;
} ## end sub check_direct_xref

1;
