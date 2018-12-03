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



package Bio::EnsEMBL::Xref::Parser::UniProtParser::Transformer;

use strict;
use warnings;

use Carp;
use Readonly;


Readonly my $PROTEIN_ID_SOURCE_NAME => 'protein_id';
Readonly my $UNIPROT_GN_SOURCE_NAME => 'Uniprot_gn';

Readonly my %whitelisted_crossreference_sources
  => (
      'ChEMBL'                => 1,
      'EMBL'                  => 1,
      'Ensembl'               => 1,
      'MEROPS'                => 1,
      'PDB'                   => 1,
      $PROTEIN_ID_SOURCE_NAME => 1,
      $UNIPROT_GN_SOURCE_NAME => 1,
    );

Readonly my $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD => 2;
Readonly my %source_selection_criteria_for_status
  => (
      'Reviewed'   => [ 'Uniprot/SWISSPROT',
                        sub {
                          return 'sequence_mapped';
                        }, ],
      'Unreviewed' => [ 'Uniprot/SPTREMBL', sub {
                          my ( $level ) = @_;
                          return ( $level <= $MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD ) ?
                            'sequence_mapped' :
                            "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD";
                        }, ],
    );

Readonly my %protein_id_extraction_recipe_for_database
  => (
      'ChEMBL' => \&_get_protein_id_xref_from_embldb_xref,
      'EMBL'   => \&_get_protein_id_xref_from_embldb_xref,
    );
sub _get_protein_id_xref_from_embldb_xref {
  my ( $embldb_extra_info, $linkage_source_id, $source_id ) = @_;

  # For both EMBL and ChEMBL entries protein ID immediately follows
  # their respective accessions i.e will be the first element of extra_info.
  my $protein_id = $embldb_extra_info->[0];

  # Strip the version number, if any, from the protein ID. At the same
  # time, filter out entries with no ID - in which case the ID is a
  # lone hyphen.
  my ( $unversioned_protein_id )
    = ( $protein_id =~ m{
                          \A
                          # Allow hyphens if they are not
                          # the first character
                          ( [^-.] [^.]+ )
                      }msx );

  if ( ! defined $unversioned_protein_id ) {
    return;
  }

  my $xref_link = {
                   'ACCESSION'          => $unversioned_protein_id,
                   'LABEL'              => $protein_id,
                   'LINKAGE_ANNOTATION' => $source_id,
                   'LINKAGE_SOURCE_ID'  => $linkage_source_id,
                   'SOURCE_ID'          => $source_id,
              };
  return $xref_link;
}



sub new {
  my ( $proto, $arg_ref ) = @_;

  my $self = {
              'maps'       => {},
              'species_id' => $arg_ref->{'species_id'},
              'xref_dba'   => $arg_ref->{'xref_dba'},
            };
  my $class = ref $proto || $proto;
  bless $self, $class;

  $self->_load_maps();

  return $self;
}


sub DESTROY {
  my ( $self ) = @_;

  $self->finish();

  return;
}


sub finish {
  my ( $self ) = @_;

  return;
}


# Transforms extracted record into form that can be consumed by
# BaseAdaptor::upload_xref_object_graphs().
sub transform {
  my ( $self, $extracted_record ) = @_;

  $self->{'extracted_record'} = $extracted_record;

  my ( $accession, @synonyms )
    = @{ $extracted_record->{'accession_numbers'} };
  my $source_id = $self->_get_source_id();

  my $xref_graph_node
    = {
       'ACCESSION'     => $accession,
       'DESCRIPTION'   => $extracted_record->{'description'},
       'INFO_TYPE'     => 'SEQUENCE_MATCH',
       'LABEL'         => $accession,
       'SEQUENCE'      => $extracted_record->{'sequence'},
       'SEQUENCE_TYPE' => 'peptide',
       'SOURCE_ID'     => $source_id,
       'SPECIES_ID'    => $self->{'species_id'},
       'STATUS'        => 'experimental',
       'SYNONYMS'      => \@synonyms,
     };

  # UniProt Gene Names links come from the 'gene_names' fields
  my $genename_dependent_xrefs
    = $self->_make_links_from_gene_names( $accession, $source_id );
  # Do not assign an empty array to DEPENDENT_XREFS, current insertion code
  # doesn't like them.
  if ( scalar @{ $genename_dependent_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DEPENDENT_XREFS'} }, @{ $genename_dependent_xrefs };
  }

  # All other xref links come from crossreferences
  my ( $direct_xrefs, $dependent_xrefs )
    = $self->_make_links_from_crossreferences( $source_id );
  # Do not assign empty arrays to FOO_XREFS, current insertion code
  # doesn't like them.
  if ( scalar @{ $direct_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DIRECT_XREFS'} }, @{ $direct_xrefs };
  }
  if ( scalar @{ $dependent_xrefs } > 0 ) {
    push @{ $xref_graph_node->{'DEPENDENT_XREFS'} }, @{ $dependent_xrefs };
  }

  return $xref_graph_node;
}


# Returns a hashref mapping source-name/priority pairs to source
# IDs. Used by the steering code to e.g. handle the setting of release
# numbers on sources.
sub get_source_id_map {
  my ( $self ) = @_;

  # Just in case, even though we presently call _load_maps() in the
  # constructor so it shouldn't be possible to call this method before
  # maps have been loaded
  if ( ! exists $self->{'maps'}->{'named_source_ids'} ) {
    confess 'Source-ID map is missing';
  }

  return $self->{'maps'}->{'named_source_ids'};
}


sub _load_maps {
  my ( $self ) = @_;

  my $xref_dba = $self->{'xref_dba'};

  my $source_id_map
    = {
       'Uniprot/SWISSPROT'
       => {
           'direct'          => undef,
           'sequence_mapped' => undef,
         },
       'Uniprot/SPTREMBL'
       => {
           'direct'            => undef,
           "protein_evidence_gt_$MAX_TREMBL_EVIDENCE_LEVEL_FOR_STANDARD" => undef,
           'sequence_mapped'   => undef,
         },
     };
  while ( my ( $source_name, $pri_ref ) = each %{ $source_id_map } ) {
    foreach my $priority ( keys %{ $pri_ref } ) {
      $pri_ref->{$priority}
        = $xref_dba->get_source_id_for_source_name( $source_name,
                                                    $priority );
    }
  }
  $self->{'maps'}->{'named_source_ids'}
    = $source_id_map;

  my %dependent_source_map
    = $xref_dba->get_xref_sources();
  $self->{'maps'}->{'dependent_sources'}
    = \%dependent_source_map;

  return;
}


# Returns true if the current record describes an entry derived from
# Ensembl, false otherwise.
sub _entry_is_from_ensembl {
  my ( $self ) = @_;

  # The old parser's way of identifying proteins derived from Ensembl
  # was to search for a comment topic "CAUTION" stating "The sequence
  # shown here is derived from an Ensembl". It was confirmed with
  # UniProt in November 2018 that this hadn't been done since late
  # 2017, therefore we no longer check this.

  # As of November 2018, the official way of identifying proteins
  # derived from Ensembl is to look for a special citation identifying
  # the record as such. There are several ways in which such special
  # citations could be identified, it seems the simplest to look at
  # "reference group name" fields to see if any of them says
  # 'Ensembl'.
  my $citation_groups
    = $self->{'extracted_record'}->{'citation_groups'};
  my $from_group_names
    = List::Util::any { $_ eq 'Ensembl' } @{ $citation_groups };

  return $from_group_names;
}


# Translate quality of the extracted entry into the matching Ensembl
# source_id, optionally with an override of priority.
sub _get_source_id {
  my ( $self, $priority_override ) = @_;

  my $source_id_map = $self->{'maps'}->{'named_source_ids'};

  my $entry_quality = $self->{'extracted_record'}->{'quality'};
  my $criteria = $source_selection_criteria_for_status{ $entry_quality->{'status'} };
  my $priority_mapper = $criteria->[1];

  my $source_name = $criteria->[0];
  my $priority = $priority_override
    // $priority_mapper->( $entry_quality->{'evidence_level'} );

  return $source_id_map->{$source_name}->{$priority};
}


# Make xrefs from 'crossreferences' entries in the extracted record,
# in a form suitable to attaching to the main xref's graph node as
# consumed by upload_xref_object_graphs(). Ensembl crossreferences
# become direct xrefs, everything else - dependent ones. If requested
# we additionally generate protein_id dependent xrefs from appropriate
# sources, i.e. EMBL and ChEMBL at present.
sub _make_links_from_crossreferences {
  my ( $self, $xref_source_id ) = @_;

  my $crossreferences = $self->{'extracted_record'}->{'crossreferences'};
  my $dependent_sources = $self->{'maps'}->{'dependent_sources'};

  my @direct_xrefs;
  my @dependent_xrefs;

 REF_SOURCE:
  while ( my ( $source, $entries ) = each %{ $crossreferences } ) {

    if ( ! $whitelisted_crossreference_sources{ $source } ) {
      next REF_SOURCE;
    }

    if ( $source eq 'Ensembl' ) {

    DIRECT_XREF:
      foreach my $direct_ref ( @{ $entries } ) {
        my $xref_link
          = {
             'STABLE_ID'    => $direct_ref->{'id'},
             'ENSEMBL_TYPE' => 'Translation',
             'LINKAGE_TYPE' => 'DIRECT',
             'SOURCE_ID'    => $self->_get_source_id( 'direct' ),
           };
        push @direct_xrefs, $xref_link;
      }

    }
    elsif ( exists $dependent_sources->{$source} ) {
      my $dependent_source_id = $dependent_sources->{$source};

    DEPENDENT_XREF:
      foreach my $dependent_ref ( @{ $entries } ) {
        my $xref_link
          = {
             'ACCESSION'          => $dependent_ref->{'id'},
             'LINKAGE_ANNOTATION' => $dependent_source_id,
             'LINKAGE_SOURCE_ID'  => $xref_source_id,
             'SOURCE_ID'          => $dependent_source_id,
           };
        push @dependent_xrefs, $xref_link;

        my $protein_id_xref_maker
          = $protein_id_extraction_recipe_for_database{ $source };
        if ( $whitelisted_crossreference_sources{ $PROTEIN_ID_SOURCE_NAME }
             && ( defined $protein_id_xref_maker ) ) {

          # Entries for the source 'protein_id' are constructed from
          # crossreferences to other databases
          my $protein_id_xref
            = $protein_id_xref_maker->( $dependent_ref->{'optional_info'},
                                        $xref_source_id,
                                        $dependent_sources->{$PROTEIN_ID_SOURCE_NAME}
                                     );
          if ( defined $protein_id_xref ) {
            push @dependent_xrefs, $protein_id_xref;
          }

        }

      }

    }
  }

  return ( \@direct_xrefs, \@dependent_xrefs );
}


# Make Uniprot_gn dependent xrefs from 'gene_names' entries in the
# extracted record, in a form suitable to attaching to the main xref's
# graph node as consumed by upload_xref_object_graphs().
sub _make_links_from_gene_names {
  my ( $self, $xref_accession, $xref_source_id ) = @_;

  my @genename_xrefs;

  # Are we supposed to process this xref source to begin with?
  if ( ! $whitelisted_crossreference_sources{ $UNIPROT_GN_SOURCE_NAME } ) {
    return [];
  }

  # UniProt Gene Name xrefs are dependent so in order to avoid
  # circular dependencies, do not generate them for proteins derived
  # from Ensembl.
  if ( $self->_entry_is_from_ensembl() ) {
    return [];
  }

  my $gene_names = $self->{'extracted_record'}->{'gene_names'};
  my $dependent_sources = $self->{'maps'}->{'dependent_sources'};
  my $dependent_source_id = $dependent_sources->{$UNIPROT_GN_SOURCE_NAME};

 GN_ENTRY:
  foreach my $gn_entry ( @{ $gene_names } ) {
    if ( ! exists $gn_entry->{'Name'} ) {
      next GN_ENTRY;
    }

    my $name = $gn_entry->{'Name'};
    my $xref = {
                'ACCESSION'          => $xref_accession,
                'LABEL'              => $name,
                'LINKAGE_ANNOTATION' => $dependent_source_id,
                'LINKAGE_SOURCE_ID'  => $xref_source_id,
                'SOURCE_ID'          => $dependent_source_id,
              };

    my $synonyms = $gn_entry->{'Synonyms'};
    if ( defined $synonyms ) {
      push @{ $xref->{'SYNONYMS'} }, @{ $synonyms };
    }

    push @genename_xrefs, $xref;

    # Regardless of how many gene-name entries a protein might have in
    # UniProt, we only want at most *one* such xref - there are bits
    # of Ensembl code which implicitly assume there can only be one
    # xref per unique (accerssion,source_id,species_id) triplet and
    # Uniprot_gn is essentially a last-resort source for when nothing
    # else works. At least we should be able to consistently output
    # the first encountered gene name (and synonyms), the old parser
    # happily ignored the fact there can be multiple gene names and
    # clobbered old names and synonyms every time a new one was
    # encountered - potentially resulting in xrefs with the last
    # encountered name but synonyms from an earlier one.
    last GN_ENTRY;
  }

  return \@genename_xrefs;
}


1;
