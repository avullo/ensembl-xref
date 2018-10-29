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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

=head1 METHODS

=head2 new

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self =  bless {}, $class;

  return $self;
  
}

# TODO
# sub run {
#   my ($self, %args) = @_;

#   # common stuff to all parsers

#   # dynamic dispatch to either (overridden) run_1 or run_2
#   my $dispatch = {
# 		  opt1 => $self->run_1,
# 		  opt2 => $self->run_2
# 		 }->{(exists $args{db})?'opt2':'opt1'};

#   throw("Unable to get the correct method to call") unless defined $dispatch;
#   $dispatch->(%args);
# }

sub run {
  my ($self, %args) = @_;
  
  my $self->{source_id}    = $args{source_id};
  my $self->{species_id}   = $args{species_id};
  my $self->{species_name} = $args{species};
  my $self->{files}        = $args{files};
  my $self->{verbose}      = $args{verbose} // 0;
  my $self->{dbi}          = $args{dbi} // $self->dbi;

  croak "Need to pass source_id, species_id, files and rel_file as pairs"
    unless defined $self->{source_id} and defined $self->{species_id} and defined $self->{files};

  my $file = shift @{$self->{files}};
  my $self->{io} = $self->get_filehandle($file);
  croak "Could not open $file\n" unless defined $self->{io};

}

sub run_script {
  my ($self, %args) = @_;

  throw("Not yet done");
}

1;
