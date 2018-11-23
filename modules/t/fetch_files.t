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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;
use Test::More;
use Test::Exception;
use Cwd;
use Bio::EnsEMBL::Xref::FetchFiles;

my $client = Bio::EnsEMBL::Xref::FetchFiles->new();

my $uris = [
  'LOCAL:test_file',
  'script:ilovesubshells',
  'ftp://ftp.ensembl.org/pub/current_README', # Testing single FTP file
  'ftp://ftp.ensembl.org/pub/misc-scripts/Variant_effect_predictor_*', # Testing glob-matched FTP files
  'http://www.ensembl.org/', # Test automatic insertion of index.html
  'http://static.ensembl.org/i/species/Mus_musculus.png', # Test http file fetch
  'https://static.ensembl.org/i/species/Homo_sapiens.png', # Test https
];

`echo 'a' > test_file`; # We need a file on disk for the LOCAL: type of input

my $here = getcwd;

my @file_list = $client->fetch_files({
  dest_dir => $here,
  user_uris => $uris,
  deletedownloaded => 0,
  verbose => 1
});

is_deeply(\@file_list, [ 
  'test_file',
  'ilovesubshells',
  map { $here.'/'.$_ } qw/
    current_README
    Variant_effect_predictor_1.0
    Variant_effect_predictor_2.0
    Variant_effect_predictor_2.1
    Variant_effect_predictor_2.2
    Variant_effect_predictor_2.3
    index.html
    Mus_musculus.png
    Homo_sapiens.png
  /],
  'All files reported downloaded'
);

# Remove any downloaded files
foreach my $alleged_file (@file_list) {
  next if $alleged_file eq 'ilovesubshells'; # skip pre-existing files, fetch_files returns the same thing
  ok(-s $alleged_file, $alleged_file.' is present');
  unlink $alleged_file;
}

$uris = ['derp:fake/protocol'];

throws_ok { $client->fetch_files({dest_dir => getcwd(), user_uris => $uris}) }
  qr/Unrecognised fetch method from config file/,
  'Try an invalid protocol';

`echo 'a' > test_file.gz`; # Create a fake gz archive

$uris = ['http://test.url/test_file.gz' ];

throws_ok { $client->fetch_files({dest_dir => getcwd(), user_uris => $uris}) }
  qr/Failed to validate:/,
  'Verify GZIP integrity check operates';

unlink('test_file.gz');

done_testing();