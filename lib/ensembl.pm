package ensembl;

our $VERSION = '0.1';

=head1 NAME

   ensembl.pm - module to bundle commonly used ensembl functions

=head1 SYNOPSIS

  use ensembl;
  my $ens = ensembl->new();

=head1 DESCRIPTION

Some aspects of the Ensembl API are a bit fiddly and it's best to bundle them in one place, rather than keep re-implement them. The main (only?) one being which database to connect to given your species of interest.

=cut

use Moose;
use Carp;
use DBI;
use Bio::EnsEMBL::Registry;


=item B<VERBOSE()>

Accessor method for 'VERBOSE' attribute. 

=cut

has '_isEnsemblMain' => (
   isa => 'Int',
   is => 'rw',
   default => 1
);

has 'species' => (
   isa => 'Str',
   is => 'rw',
   default => 'human',
   required => 1,
   trigger => \&_checkSpecies
);

has 'VERBOSE' => (
   isa => 'Int',
   is => 'rw',
   default => 1
);


## this is required in order to pick the correct
## connection parameters to the Ensembl API as 
## species from the Ensembl Genomes projects differ from the main API
sub connect {
   my $self = shift;
   my $species = shift;
   
   my $registry = 'Bio::EnsEMBL::Registry';
   
   if ($self->_isEnsemblMain) {  # this is for the main API species
      print "Connect to main Ensembl API...\n" if $self->VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'ensembldb.ensembl.org',
          -user => 'anonymous'
      );
      
   } else {  # this is for the Ensemble Genomes API species
      print "Connecting to Ensembl Genomes API...\n" if $self->VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'mysql.ebi.ac.uk',
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
   $mainSpecies{$_}++ foreach (qw/chicken Gallus_gallus human Homo_sapiens mouse Mus_musculus yeast Saccharomyces_cerevisae/);
   
   if (exists($mainSpecies{$species})) {
      $self->_isEnsemblMain(1);
   } else {
      $self->_isEnsemblMain(0);
   }
}

1;

=head1 AUTHOR

Chris Cole <c.cole@dundee.ac.uk>

=head1 COPYRIGHT

Copyright 2016, Chris Cole. All rights reserved.

This library is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut