
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

=head1 NAME

Bio::EnsEMBL::Xref::Parser - Base parser

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Xref::Parser;

use strict;
use warnings;

use Carp;

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

=head1 METHODS

=head2 new

Base constructor

=cut

sub new {
  my ( $caller, %args ) = @_;

  my $class = ref($caller) || $caller;

  my $self = bless { source_id  => $args{source_id},
                     species_id => $args{species_id},
                     species    => $args{species},
                     rel_file   => $args{rel_file},
                     files      => $args{files},
                     xref_dba   => $args{xref_dba},
                     dba        => $args{dba},
                     verbose    => $args{verbose} // 0, }, $class;

  defined $self->{source_id}    and
    defined $self->{species_id} and
    defined $self->{files}      and
    defined $self->{xref_dba} or
    croak "Need to pass (source_id, species_id, files, xref_dba) args";

  # extra necessary param checking
  assert_ref( $self->{files}, 'ARRAY' );
  assert_ref( $self->{xref_dba},
              'Bio::EnsEMBL::Xref::DBSQL::BaseAdaptor' );
  assert_ref( $self->{xref_dba}->dbi, 'DBI::db' );
  assert_ref( $self->{dba},           'Bio::EnsEMBL::DBSQL::DBAdaptor' )
    if defined $self->{dba};

  return $self;
} ## end sub new

=head2 run

Abstract method, implement in parser subclass

=cut

sub run {
  my $self = shift;

  throw( "Cannot call Parser::run abstract method: provide implementation in subclass" );
}

1;
