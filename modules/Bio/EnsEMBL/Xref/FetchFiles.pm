
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

package Bio::EnsEMBL::Xref::FetchFiles;

use strict;
use warnings;

# Given one or several FTP or HTTP URIs, download them.  If an URI is
# for a file or MySQL connection, then these will be ignored.  For
# FTP, standard shell file name globbing is allowed (but not regular
# expressions).  HTTP does not allow file name globbing.  The routine
# returns a list of successfully downloaded local files or an empty list
# if there was an error.

use Carp;
use DBI;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;
use POSIX qw(strftime);

use File::Basename;
use File::Spec::Functions;
use IO::File;
use Net::FTP;
use HTTP::Tiny;
use URI;
use URI::file;
use Text::Glob qw( match_glob );
use LWP::UserAgent;

my $base_dir = File::Spec->curdir();

sub new {
  my ($proto) = @_;

  my $class = ref $proto || $proto;
  return bless {}, $class;
}

sub get_ftp {
  my ( $self, $uri, $passive ) = @_;
  my $ftp;

  if ($passive) {
    $ftp = Net::FTP->new( $uri->host(), 'Debug' => 0, Passive => 1 );
  }
  else {
    $ftp = Net::FTP->new( $uri->host(), 'Debug' => 0 );
  }

  if ( !defined($ftp) ) {
    printf( "==> Can not open FTP connection: %s\n", $@ );
    return ();
  }

  if ( !$ftp->login( 'anonymous', '-anonymous@' ) ) {
    printf( "==> Can not log in on FTP host: %s\n", $ftp->message() );
    return ();
  }

  if ( !$ftp->cwd( dirname( $uri->path() ) ) ) {
    printf( "== Can not change directory to '%s': %s\n",
            dirname( $uri->path() ), $ftp->message() );
    return ();
  }

  $ftp->binary();
  return $ftp;
} ## end sub get_ftp

1;
