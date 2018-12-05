
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

=head1 DESCRIPTION

Given one or several FTP or HTTP URIs, download them.  If a URI is
for a file or MySQL connection, then these will be ignored.  For
FTP, standard shell file name globbing is allowed (but not regular
expressions).  HTTP downloads are deferred to Bio::EnsEMB::Utils::Net

=cut

package Bio::EnsEMBL::Xref::FetchFiles;

use strict;
use warnings;

use Carp;
use File::Basename;
use File::Spec::Functions;
use File::Path qw( make_path );
use Net::FTP;
use LWP::UserAgent;
use URI;
use URI::file;
use Text::Glob qw( match_glob );

sub new {
  my ($proto) = @_;

  my $class = ref $proto || $proto;
  return bless {}, $class;
}

=head2 fetch_files

Arg 1      :  Hashref { 
                dest_dir => $path_for_downloads, 
                user_uris => [$uri,...], 
                deletedownloaded => boolean, # Delete any existing files
                verbose => boolean 
              }
Description:  Given a config, this method selects the relevant client and downloads
              the files from the supplies list of URLs. It uses a dispatch table to choose
              which method to run to fetch. It always tries to reuse existing downloaded
              files - set deletedownloaded if you don't want that to happen
Returntype :  List of paths on the local file system

=cut

sub fetch_files {
  my ( $self, $arg_ref ) = @_;
  my $dest_dir         = $arg_ref->{dest_dir};
  my $user_uris        = $arg_ref->{user_uris};
  my $deletedownloaded = $arg_ref->{del_down};
  my $verbose          = $arg_ref->{verbose};

  my @processed_files; # to be returned to caller
  my %dispatch = (
    script => \&script_handler,
    file => \&file_handler,
    ftp => \&ftp_handler,
    http => \&http_handler,
    https => \&http_handler,
    mysql => \&script_handler,
  );

  # Ensure there is a folder to receive downloaded files
  make_path($dest_dir); 

  foreach my $stringy_uri (@$user_uris) {
    # Change old-style 'LOCAL:' URIs into 'file:'.
    $stringy_uri =~ s/^LOCAL:/file:/ix;
    my $uri = URI->new($stringy_uri);

    unless (exists $dispatch{$uri->scheme()}) { 
      croak 'Unrecognised fetch method from config file: '.$uri->scheme().' from '.$stringy_uri;
    }
    my $code = $dispatch{$uri->scheme()};
    my @downloaded_files = $self->$code($stringy_uri, $uri, $dest_dir, $deletedownloaded, $verbose);
    push @processed_files, @downloaded_files;
  }

  if (!@processed_files) { croak "No files were downloaded or existed already"; }
  # Validate any compressed files we downloaded
  my @failures;
  foreach my $file (@processed_files) {
    if ($file =~ /\.(gz|Z)$/x) {
      my $cmd = "gzip -t $file";
      my $return = system($cmd);
      if ($return != 0) {
        push @failures, $file;
      }
    }
  }

  if (@failures) {
    my $error = sprintf "Failed to validate: %s\nThese files have been deleted", join ',',@failures;
    foreach my $doomed (@failures) {
      unlink $doomed;
    }
    confess $error;
  }

  return @processed_files;
}

=head2 script_handler

Description:  A dud function, for the eventuality that we should ever want to
              not do nothing with script-type, mysql-type sources
Caller:       fetch_files

=cut

sub script_handler {
  my ($self, $user_uri) = @_;
  $user_uri =~ s/\Ascript://x;
  return $user_uri;
}

=head2 file_handler

Description:  Handles file system paths, by checking if they're present
Caller:       fetch_files

=cut

sub file_handler {
  my ($self, $user_uri) = @_;

  $user_uri =~ s/\Afile://x;
  if ( -s $user_uri ) {
    return $user_uri;
  } else {
    confess "Claimed file does not exist or is empty: $user_uri";
  }
}

=head2 ftp_handler

Description:  Fetches glob-matchable files from an FTP site. It is able to
              avoid downloading by looking for the expected file.
              Could not be handled in Bio::EnsEMBL::Utils::Net, due to needing
              a manifest of files in order to do glob-matching
Caller:       fetch_files

=cut

sub ftp_handler {
  my ($self, $user_uri, $uri, $dest_dir, $delete, $verbose) = @_;
  
  # In case of single file fetch
  my $path = $self->check_download($dest_dir, basename($uri->path()), $delete, $verbose);
  return $path if defined $path;
  # From here on we are pulling sets of files from the FTP site, either by glob pattern
  # or just getting all of them

  my @ftp_manifest = @{ $self->list_ftp_files($uri) };

  my @download_manifest; # A buffer of file paths
  my $ftp = $self->get_ftp($uri);

  foreach my $remote_file (@ftp_manifest) {
    # Skip if the pattern in the URI does not match individual files?
    # Probably for removing .. and such
    next if ( !match_glob( basename( $uri->path() ), $remote_file ) );

    my $file_path = $self->check_download($dest_dir, basename($remote_file), $delete, $verbose);
    if ($file_path) {
      push @download_manifest, $file_path;
    } else {
      print "Fetching $remote_file to $dest_dir\n" if $verbose;
      my $status = $ftp->get( $remote_file, catfile( $dest_dir, basename($remote_file)) );
      if ($status) {
        push @download_manifest, catfile($dest_dir,basename($remote_file));
      } else {
        confess "Failed to download $remote_file via FTP: ".$ftp->message();
      }
    }
  }
  $ftp->quit();
  
  return @download_manifest;
}

=head2 http_handler

Description:  Fetches a single URL. Will check for existing files by default
Caller:       fetch_files

=cut

sub http_handler {
  my ($self, $uri_string, $uri, $dest_dir, $delete, $verbose) = @_;
  my $remote_file_name;
  if ($uri->path eq '/') { $remote_file_name = "index.html"; } # failover to default webpage
  $remote_file_name //= basename($uri->path());

  my $path = $self->check_download($dest_dir, $remote_file_name, $delete, $verbose);
  return $path if defined $path;

  print "Fetching $uri_string\n" if $verbose;
  my $download_path = catfile($dest_dir,$remote_file_name);
  
  my $ua = LWP::UserAgent->new();
  $ua->ssl_opts( verify_hostname => 0);
  my $response = $ua->get($uri_string, ':content_file' => $download_path);
  if (! $response->is_success) {
    confess "Failed when downloading $uri_string with response ".$response->status_line;

  return $download_path;
}

=head2 get_ftp

Arg 1:        URI::ftp instance
Arg 2:        Passive option for Net::FTP
Arg 3:        FTP user account
Arg 4:        FTP password should it be necessary
Description:  Provides an FTP client configured to download files from the given URI
              It would be a lot more robust if we could use LWP::UserAgent, but this
              code requires an FTP->ls() call, hence a real FTP client is needed
ReturnType:   Net::FTP instance

=cut

sub get_ftp {
  my ($self, $uri, $passive, $user, $pass) = @_;
  my $ftp;
  $passive //= 1; # Default use of passive mode for FTP (better for tight firewalls)
  $user //= 'anonymous';
  $pass //= '-anonymous@';

  $ftp = Net::FTP->new( $uri->host(), Debug => 0, Passive => $passive);

  croak sprintf("Cannot open FTP connection: %s\n", $@) if !defined($ftp);

  my $state = $ftp->login( $user, $pass );
  if (!$state) {
    croak sprintf( "Cannot log in on FTP host: %s\n", $ftp->message() );
  }
  
  my ($filename, $path, $extras) = fileparse( $uri->path() ); # fileparse is a bit safer than dirname
  $state = $ftp->cwd( $path );

  if (!$state) {
    croak sprintf( "== Can not change directory to '%s': %s\n", $path, $ftp->message() );  
  }
  
  $ftp->binary();
  return $ftp;
} ## end sub get_ftp


=head2 check_download

Arg 1:        Base path
Arg 2:        File name
Arg 3:        Delete? (boolean)
Arg 4:        Verbose (boolean)
Description:  Check filesystem for the downloads we expect from this source
              Optionally clean them up so we can replace them
Returntype:   Path to file - A path means the file exists, undef means it does not.
              It may have been deleted to make the result undef, as per Arg 3

=cut

sub check_download {
  my ($self, $base_path, $filename, $delete, $verbose) = @_;

  my $file_path = catfile($base_path,$filename);
  if (-s $file_path) {
    if ($delete) {
      unlink $file_path;
      print "Deleted $file_path\n" if $verbose;
      return;
    }
    return $file_path;
  } else {
    return;
  }
}


=head2 list_ftp_files

Arg 1:        String URI or URI object
Arg 2:        Any file extension you want removed from the list
Description:  Returns a list of files found in the given FTP directory
              Mainly used to see what species data is available in the Xref pipeline
ReturnType :  Listref of file names

=cut

sub list_ftp_files {
  my ($self, $uri, $suffix) = @_;
  
  if (ref($uri) ne 'URI') {
    $uri = URI->new($uri);
  }
  my $ftp_client = $self->get_ftp($uri);
  my @files;
  @files = $ftp_client->ls()
    or confess sprintf "Cannot list content of FTP site %s at %s: %s", $uri->host, $uri->path, $ftp_client->message;
  $ftp_client->quit;

  foreach my $name (@files) {
    $name =~ s/$suffix//x;
  }
  return \@files;
}

1;
