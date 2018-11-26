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



package Bio::EnsEMBL::Xref::Parser::UniProtParser::Loader;

use strict;
use warnings;

use Carp;


sub new {
  my ( $proto, $arg_ref ) = @_;

  my $self = {
    'batch_size'         => $arg_ref->{'batch_size'}         // 1,
    'checkpoint_seconds' => $arg_ref->{'checkpoint_seconds'} // 0,
    'xref_dba'           => $arg_ref->{'xref_dba'},
    'last_flush'         => time(),
  };
  my $class = ref $proto || $proto;
  bless $self, $class;
  $self->_clear_send_buffer();

  return $self;
}


sub DESTROY {
  my ( $self ) = @_;

  $self->finish();

  return;
}


sub finish {
  my ( $self ) = @_;

  $self->flush();

  return;
}


sub flush {
  my ( $self ) = @_;

  my $xref_dba = $self->{'xref_dba'};

  if ( ! $xref_dba->upload_xref_object_graphs( $self->{'send_buffer'} ) ) {
    confess 'Failed to upload xref object graphs. Check for errors on STDOUT';
  }

  $self->_clear_send_buffer();

  return;
}


sub load {
  my ( $self, $transformed_data ) = @_;

  if ( ! defined $transformed_data ) {
    return;
  }

  $self->_add_to_send_buffer( $transformed_data );

  my $current_time = time();
  if ( ( $self->{'send_backlog'} >= $self->{'batch_size'} )
       || ( $self->{'checkpoint_seconds'} > 0  )
       && ( $current_time - $self->{'last_flush'} > $self->{'checkpoint_seconds'} ) ) {
    $self->flush();
    $self->{'last_flush'} = $current_time;
  }

  return;
}


sub prepare_source_for_dependent_xrefs {
  my ( $self, $source_id ) = @_;

  my $xref_dba = $self->{'xref_dba'};

  $xref_dba->get_dependent_mappings( $source_id );

  return;
}


sub _add_to_send_buffer {
  my ( $self, $entry ) = @_;

  push @{ $self->{'send_buffer'} }, $entry;
  $self->{'send_backlog'}++;

  return;
}


sub _clear_send_buffer {
  my ( $self ) = @_;

  $self->{'send_buffer'}  = [];
  $self->{'send_backlog'} = 0;

  return;
}


1;
