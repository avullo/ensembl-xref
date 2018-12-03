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

=head1 DESCRIPTION

A database class that auto-instantiates a DBIC schema

Convenience methods prevent having to delve into DBIC guts for common activities

=head1 SYNOPSIS

my $db = Bio::EnsEMBL::Xref::DB->new(
  config => {
    host => 'db.com',
    port => 3306,
    user => 'me',
    pass => 'Go on then',
    driver => 'mysql',
    db => 'name_for_db',
    create => 1, # Deploys the schema to the DB for first use
  }
);

$db = Bio::EnsEMBL::Xref::DB->new(
  config_file => 'db.conf' # config options in Config::General format
);

my $dbh = $db->dbh; # $dbh is a DBI database handle borrowed for direct SQL
$dbh->prepare('DROP TABLE dependent_xref');

$db->create_db_row('Xref',{
  xref_id => 1,
  accession => 'YAY',
  description => 'Sample new Xref',
  source_id => 1,
  ...
});

=cut


package Bio::EnsEMBL::Xref::DB;

use strict;
use warnings;

use Moose;
use namespace::autoclean;
use Config::General;
use Carp;
use Bio::EnsEMBL::Xref::Schema;
use DBI;

has schema => (
  isa => 'Bio::EnsEMBL::Xref::Schema',
  is => 'ro',
  builder => '_init_db'
);

has config_file => (
  isa => 'Str',
  is => 'rw',
  builder => '_guess_config'
);

has config => (
  isa => 'HashRef',
  is => 'rw'
);


=head2 new
  Arg [1]    : HashRef of configuation parameters (driver, db, host, port, user, pass)
  Description: Initialise the core database.
  Return type: schema
  Caller     : internal

=cut

sub _init_db {
  my $self = shift;
  print STDERR "Setting up schema\n";
  $self->_init_config if ! defined $self->config;
  $self->_validate_config($self->config);
  my %conf = %{ $self->config };
  my %opts;
  $opts{mysql_enable_utf8} = 1 if ($conf{driver} eq 'mysql');
  $opts{mysql_auto_reconnect} = 1 if ($conf{driver} eq 'mysql');
  $opts{sqlite_unicode} = 1 if($conf{driver} eq 'SQLite');
  my $dsn;
  if ($conf{driver} eq 'SQLite') {
    $dsn = sprintf("dbi:%s:database=%s",$conf{driver},$conf{file});
  } else {
    $dsn = sprintf("dbi:%s:database=%s;host=%s;port=%s",$conf{driver},$conf{db},$conf{host},$conf{port});
  }

  my %deploy_opts = ();
  # Example deploy option $deploy_opts{add_drop_table} = 1;
  print STDERR 'Connecting: '.$dsn."\n";
  my $schema = Bio::EnsEMBL::Xref::Schema->connect($dsn, $conf{user},$conf{pass}, \%opts);

  if ($conf{create} == 1 && $conf{driver} eq 'mysql') {
    my $dbh = DBI->connect(sprintf("DBI:%s:database=;host=%s;port=%s",$conf{driver},$conf{host},$conf{port}),$conf{user},$conf{pass}, \%opts);
    $dbh->do('CREATE DATABASE '.$conf{db}.';');
    $dbh->disconnect;
  }
  $schema->deploy(\%deploy_opts) if $conf{create} == 1;

  return $schema;
}


=head2 new
  Description: Don't want production use to guess at least at the moment.
               This mainly exists so TestDB can override and replace with a
               useful default
  Return type: undef
  Caller     : internal

=cut

sub _guess_config {
  return;
}


=head2 new
  Arg [1]    : HashRef of configuation parameters (driver, db, host, port, user, pass)
  Description: Initialisae the loading of the configuration file.
  Return type: HashRef - $self->config
  Caller     : internal

=cut

sub _init_config {
  my $self = shift;

  if (defined $self->config_file) {
    my $conf = Config::General->new($self->config_file);
    my %opts = $conf->getall();
    $self->config(\%opts);
  } else {
    confess 'No config or config_file provided to new(). Cannot execute';
  }

  return $self->config;
}


=head2 new
  Arg [1]    : HashRef of configuation parameters (driver, db, host, port, user, pass)
  Description: Configuration file parameter validation
  Return type: DBI database handle
  Caller     : internal

=cut

sub _validate_config {
  my ($self,$config) = @_;
  my @required_keys = qw/driver/;
  if ($config->{driver} eq 'mysql') {
    push @required_keys, qw/db host port user pass/;
  } elsif ($config->{driver} eq 'SQLite') {
    push @required_keys, qw/file/;
  } else {
    confess "TestDB config requires parameter 'driver' with value mysql or SQLite";
  }
  my @errors;
  foreach my $constraint (@required_keys) {
    if (! exists $config->{$constraint}) {
      push @errors, "Missing argument '$constraint'";
    }
  }
  if (scalar @errors > 0) {
    confess sprintf "%s \n%s",
      ($self->config_file) ? 'Missing options in '.$self->config_file. ': ' : 'Missing options in supplied config: ',
      join ';',@errors;
  }
}


=head2 new
  Description: Shortcut for accessing a database handle directly. I get the
               impression we might be doing this a lot.
  Return type: DBI database handle
  Caller     : internal

=cut

sub dbh {
  my $self = shift;
  return $self->schema->storage->dbh;
}


=head2 new
  Arg [1]    : model
  Arg [2]    : arguments : These should be key-value pairs matching the rows in
                           the table
  Description: Shortcut for creating things on the fly
  Return type:
  Caller     : internal

=cut

sub create_db_row {
  my ($self,$model, $params) = @_;
  my $source = $self->schema->resultset($model)->create(
    $params
  );
  return $source;
}

__PACKAGE__->meta->make_immutable;

1;
