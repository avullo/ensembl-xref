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

=head1 DESCRIPTION

Designed to parse the Rat Genome Database download file, historically hosted at
ftp://ftp.rgd.mcw.edu/pub/data_release/GENES_RAT.txt . It comprises 40+ columns in a 
tab-separated format

It contains RGD IDs (which are numeric), and associates them either with Ensembl genes or
RefSeq records (mainly transcripts).

=cut

package Bio::EnsEMBL::Xref::Parser::RGDParser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( Bio::EnsEMBL::Xref::Parser );

=head2 run

Description: Triggers the parsing of the RGD file specified in files parameter
             It uses Text::CSV to consume the source file.

=cut

sub run {

  my ($self) = @_;
  my $source_id    = $self->{source_id};
  my $species_id   = $self->{species_id};
  my $files        = $self->{files};
  my $verbose      = $self->{verbose};
  my $xref_dba     = $self->{xref_dba};

  unless (defined $source_id && defined $species_id && defined $files) {
    confess "Need to pass source_id, species_id and files list";
  }
  $verbose //= 0;

  my $file = @{$files}[0];

  # Used to assign dbIDs for when RGD Xrefs are dependent on RefSeq xrefs
  my (%preloaded_refseq) = %{$xref_dba->get_valid_codes("refseq", $species_id)};
  
  my $rgd_io = $xref_dba->get_filehandle($file);

  my $csv = Text::CSV->new({
    sep => "\t",
    blank_is_undef => 1,
    auto_diag => 1,
    binary => 1,
    allow_loose_quotes => 1
  })
    or confess "Cannot use CSV: ".Text::CSV->error_diag ();
  # WARNING - Text::CSV does not like the GENES-RAT.txt file. It is improperly formatted and contains a non-ASCII character
  # Make sure binary is turned on or it silently fails and you get 1/3rd of the records.
  # strict is turned off to prevent failure on a blank line at the end

  my $line = '#';
  while (substr($line,0,1) eq '#') {
    $line = $rgd_io->getline;
  }
  $csv->parse($line);
  $csv->column_names($csv->fields);
  # Columns we want
  #  GENE_RGD_ID => 0,
  #  SYMBOL => 1,
  #  NAME => 2,
  #  GENBANK_NUCLEOTIDE => 23,
  #  OLD_SYMBOL => 29,
  #  ENSEMBL_ID => 37

  my $count= 0;
  my $ensembl_count = 0;
  my $mismatch = 0;
  my $syn_count = 0;
  my $cols; # Digested columns from CSV
  while ( $cols = $csv->getline_hr($rgd_io) ) {
    next if exists $cols->{GENE_RGD_ID} && ($cols->{GENE_RGD_ID} eq '' || !defined $cols->{GENE_RGD_ID});

    my @nucs;
    @nucs = split(/\;/,$cols->{GENBANK_NUCLEOTIDE}) if defined $cols->{GENBANK_NUCLEOTIDE};
    my $done = 0;
    # @nucs are sorted in the file in alphabetical order. Filter them down
    # to a higher quality subset, then add dependent Xrefs where possible
    foreach my $nuc ($self->sort_refseq_accessions(@nucs)){
      
      if(!$done && exists $preloaded_refseq{$nuc}){
        
        foreach my $xref (@{$preloaded_refseq{$nuc}}){
          my $xref_id = $xref_dba->add_dependent_xref({ 
            master_xref_id => $xref,
            acc            => $cols->{GENE_RGD_ID},
            label          => $cols->{SYMBOL},
            desc           => $cols->{NAME},
            source_id      => $source_id,
            species_id     => $species_id
          });
          $count++;
          $syn_count += $self->process_synonyms($xref_id, $cols->{OLD_SYMBOL});
          $done = 1;
        }
      }
    }

    if (defined $cols->{ENSEMBL_ID}) {
      my @ensembl_ids = split(/\;/, $cols->{ENSEMBL_ID});
      
      foreach my $id (@ensembl_ids) {
        $ensembl_count++;
        $xref_dba->add_to_direct_xrefs({
          stable_id => $id,
          type => 'gene',
          acc => $cols->{GENE_RGD_ID},
          label => $cols->{SYMBOL},
          desc => $cols->{NAME},
          source_id => $source_id,
          species_id => $species_id
        });
        my $xref_id = $xref_dba->get_xref($cols->{GENE_RGD_ID}, $source_id, $species_id);
        $syn_count += $self->process_synonyms($xref_id, $cols->{OLD_SYMBOL});
        $done = 1;
      }
    }
    if(!$done){
      $xref_dba->add_xref({ 
        acc        => $cols->{GENE_RGD_ID},
        label      => $cols->{SYMBOL},
        desc       => $cols->{NAME},
        source_id  => $source_id,
        species_id => $species_id,
        info_type  => "MISC"
      });
      $mismatch++;
    }

  }
  if (! $csv->eof) {
    confess 'Failed to finish parsing RGD file: '.$csv->error_diag();
  }
  $rgd_io->close();

  if($verbose){
    print "$count xrefs succesfully loaded and dependent on refseq\n";
    print "$mismatch xrefs added but with NO dependencies\n";
    print "$ensembl_count direct xrefs successfully loaded\n";
    print "Tried to add $syn_count synonyms, including duplicates\n";
  }
  return;
}

# Predefined importance levels for the most valued RefSeq accession types
my %refseq_priorities = (
  NM => 1,
  NP => 1,
  NR => 1,
  XM => 2,
  XP => 2,
  XR => 2,
);


=head2 sort_refseq_accessions

Arg [1..n]  : Original list of accessions
Description : Filter out any accessions which are not in the "normal" set of
              genomic features. The column in question contains EMBL accessions
              as well as other things, and we don't have the ability to make
              Xrefs to all sources
Returntype  : List of sorted and filtered accessions

=cut
sub sort_refseq_accessions {
  my ($self,@accessions) = @_;
  @accessions = sort {
                  $refseq_priorities{substr $a, 0,2} <=> $refseq_priorities{substr $b, 0,2}
                  || $a cmp $b
                }
                grep { exists $refseq_priorities{substr $_, 0,2} } @accessions;
  return @accessions;
}


=head2 process_synonyms
Arg [1]     : Xref dbID to attach synonyms to
Arg [2]     : Synonym string as read from file
Description : Process the synonym column into potentially many items and add
              them to the synonym table. Synonyms are ';' separated
Returntype  : Int - the count of synonyms added
=cut

sub process_synonyms {
  my ($self, $xref_id, $synonym_string) = @_;
  my $syn_count = 0;
  return $syn_count unless defined $synonym_string && defined $xref_id;

  my @syns = split(/\;/,$synonym_string);
  foreach my $syn (@syns) {
    $self->{xref_dba}->add_synonym($xref_id,$syn);
    $syn_count++;
  }
  return $syn_count;
}

1;
