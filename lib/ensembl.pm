package ensembl;

our $VERSION = '0.4.1';

=head1 NAME

   ensembl.pm - module to bundle commonly used ensembl functions

=head1 SYNOPSIS

  use ensembl;
  my $ens = ensembl->new();                   # create new ensembl object
  $ens->species("Saccharomyces_cerevisiae");  # specify species
  my $registry = $ens->connect();             # connect to Ensembl database

=head1 DESCRIPTION

Some aspects of the Ensembl API are a bit fiddly and it's best to bundle them in one place, rather than keep re-implement them. The main (only?) one being which database to connect to given your species of interest.

=cut

use Moose;
use Carp;
use DBI;
use Bio::EnsEMBL::Registry;

## internal attribute defining whether given species is from
## main Ensembl database or Ensembl Genomes
has '_isEnsemblMain' => (
   isa => 'Int',
   is => 'rw',
   default => 1
);

=item B<species()>

Setter/Getter method for 'species' attribute. Needs to be a valid Ensembl/Ensembl Genomes species.

Default: 'human'
Returns: string

=cut

has 'species' => (
   isa => 'Str',
   is => 'rw',
   default => 'human',
   required => 1,
   trigger => \&_checkSpecies
);

=item B<genomeBuild()>

Setter/Getter method for the 'genomeBuild' attribute. Define which genome build to use. Only of relevance to human.

Default: 'GRCh38' 
Returns: string

=cut

has 'genomeBuild' => (
   isa => 'Str',
   is => 'rw',
   default => 'GRCh38',
   trigger => \&_checkBuild
);

=item B<VERBOSE()>

Setter/Getter method for 'VERBOSE' attribute. 

Default: 1
Returns: int

=cut

has 'VERBOSE' => (
   isa => 'Int',
   is => 'rw',
   default => 1
);


=item B<connect()>

Method to connect to Ensembl. Requires species to be set.

Returns: Ensembl registry object

=cut

sub connect {
   my $self = shift;
   
   my $registry = 'Bio::EnsEMBL::Registry';
   
   my $port = 3306;
   $port = 3337 if ($self->genomeBuild eq 'GRCh37' && $self->species eq 'human');
   
   if ($self->_isEnsemblMain) {  # this is for the main API species
      print "Connect to main Ensembl API...\n" if $self->VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'ensembldb.ensembl.org',
          -user => 'anonymous',
          -port => $port
      );
      
   } else {  # this is for the Ensemble Genomes API species
      print "Connecting to Ensembl Genomes API...\n" if $self->VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'mysql-eg-publicsql.ebi.ac.uk',
          -user => 'anonymous',
          -port => 4157,
      );
      
   }
   return($registry);
}

# method to check whether a species in the Ensembl main table or not
# Can't really check whether it is a valid species in Ensembl as there
# are too many species to check, would be really slow.
sub _checkSpecies {
   my $self = shift;
   my $species = shift;
   
   my %mainSpecies;
   $mainSpecies{$_}++ foreach (qw/chicken Gallus_gallus human Homo_sapiens mouse Mus_musculus Saccharomyces_cerevisiae/);
   
   if (exists($mainSpecies{$species})) {
      $self->_isEnsemblMain(1);
   } else {
      $self->_isEnsemblMain(0);
   }
}

## internal method checking the defined build is valid or appropriate
sub _checkBuild {
   my $self = shift;
   my $build = shift;
   
   if ($self->species eq 'human') {
      unless ($build eq 'GRCh38' or $build eq 'GRCh37') {
         carp "Unrecognised human genome build '$build'. Defaulting to 'GRCh38'\n";
         $self->genomeBuild('GRCh38');
      }
   } else {
      carp "Specifying a genome build is only supported for human. Can only use the default build for this species.\n"
   }
}

1;

=head1 AUTHOR

Chris Cole <c.cole@dundee.ac.uk>

=head1 COPYRIGHT

Copyright 2016, Chris Cole. All rights reserved.

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut