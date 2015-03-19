#!/usr/bin/perl

=head1 NAME

annotate_with_gene_name.pl - annotate an Ensembl fasta with a gene name for each entry

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;

use lib '/opt/perl/bioperl-live';
## correct path to the below module required via -I on commandline
use Bio::EnsEMBL::Registry;

$| = 1;
my $fasta;
my $species;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;

GetOptions (
   'in=s'      => \$fasta,
   'species=s' => \$species,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') if (!$fasta or !-e $fasta);
pod2usage(-msg => 'Please supply a species.') if (!$species);

my $registry = connectEnsemblRegistry($species);
my $gene_adaptor = $registry->get_adaptor($species, 'Core', 'Gene');
die "ERROR - failed to get gene adaptor for '$species'. Check spelling and that it's a valid Ensembl species. Or check that you're using the correct API.\n" unless (defined($gene_adaptor));
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($registry->version_check($registry->get_DBAdaptor($species, 'core')));

print "Annotating entries...\n" if $VERBOSE;
my $out = "$fasta.tmp";
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
open(my $IN, "<", $fasta) or die "ERROR - unable to open '$fasta': ${!}\nDied";
my $c = 0;
while(<$IN>) {
   if (/^>.* gene:(\S+)/) {  # extract gene ID from description line
      ++$c;
      print "Fetching gene '$1' from Ensembl...\n" if $DEBUG;
      printf "\r%5s", $c unless $DEBUG;
      my $g = $gene_adaptor->fetch_by_stable_id($1); # fetch it from ensembl
      my $name;
      if (!defined($g)) {  # check gene exists
         warn "Warning - gene '$1' not found\n";
         $name = 'unknown';
      } else {
         $name = $g->external_name(); # get the name
      }
      if (defined($name)) { 
         s/ / name:$name /;   # add name to first space in desc line
      } else {
         s/ / name:unknown /; # or if no name put 'unknown'
      }
   }
   
   print $OUT $_;
}
close($IN);
if ($c == 0) {
   warn "Warning - no gene names found. Nothing changed\n";
   unlink($out);
   exit;
}
rename($out, $fasta) or die "Failed to rename '$out' to '$fasta'\n";
print "\nDone!\n";

## this is required in order to pick the correct
## connection parameters to the Ensembl API as 
## species from the Ensembl Genomes projects differ from the main API
sub connectEnsemblRegistry {
   my $species = shift;
   
   my $registry = 'Bio::EnsEMBL::Registry';
   my %main;
   $main{$_}++ foreach (qw/chicken human mouse/);
   
   if ($main{$species}) {  # this is for the main API species
      print "Connect to main Ensembl API...\n" if $VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'ensembldb.ensembl.org',
          -user => 'anonymous'
      );
      
   } else {  # this is for the Ensemble Genomes API species
      print "Connecting to Ensembl Genomes API...\n" if $VERBOSE;
      
      $registry->load_registry_from_db(
          -host => 'mysql.ebi.ac.uk',
          -user => 'anonymous',
          -port => 4157,
      );
      
   }
   return($registry);
}

=head1 SYNOPSIS

perl -I<path API modules> annotate_with_gene_name.pl --in <file> --species <name> [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Fasta files downloaded from Ensembl have a lot of information on the defline for each entry, but none of it human friendly: it's all ID numbers. Generally what people want is to know what the transcript/protein is.

This script rectifies this by adding the gene name as the second field in the defline.

B<NB> owing to the current situation with Ensembl requiring different APIs for different species you need to supply the path to the correct modules for the species your searching for. This will probably change in future.

=head1 OPTIONS

=over 5

=item B<--in>

Input fasta file as retrieved from Ensembl.

=item B<--species>

Species name as understood by Ensembl.

=item B<--verbose|--no-verbose>

Toggle verbosity. [default:none]

=item B<--debug|--no-debug>

Toggle debugging output. [default:none]

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=head1 AUTHOR

Chris Cole <christian@cole.name>

=cut