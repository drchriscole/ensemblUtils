#!/usr/bin/perl

=head1 NAME

add_gene_name_column.pl - given a tab delimited file with Ensembl IDs in the first col, add a gene name.

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;

our $VERSION = '1.2';

my $file;
my $species = '';
my $feature = 'Gene';
my $desc = 0;
my $coords = 0;
my $delim = "\t";
my $out = 'ensembl_annotated.csv';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;

GetOptions (
   'in=s'      => \$file,
   'species=s' => \$species,
   'feature-type=s' => \$feature,
   'desc!'     => \$desc,
   'coords!'   => \$coords,
   'delim=s'   => \$delim,
   'out=s'     => \$out,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') unless ($file && -e $file);
pod2usage(-msg => 'Please supply a species name.') unless ($species);

$desc = 1 if ($coords); # enable descriptions if co-ordinates are set.

die "ERROR - invalid feature-type '$feature'. Use 'Gene' or 'Transcript'." unless ($feature =~ /^(gene|transcript)$/i); 

printf "NOTE: using Ensembl API version %s\n", software_version();
my $registry = connectEnsemblRegistry($species);
my $adaptor = $registry->get_adaptor($species, 'Core', 'Gene');
die "ERROR - failed to get adaptor for '$species'. Check spelling and that it's a valid Ensembl species. Or check that you're using the correct API.\n" unless (defined($adaptor));
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($registry->version_check($registry->get_DBAdaptor($species, 'core')));

print "Annotating entries...\n" if $VERBOSE;
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
open(my $IN, "<", $file) or die "ERROR - unable to open '$file': ${!}\nDied";
my $c = 0;
while(<$IN>) {
   chomp;
   s/\"//g; # remove quotes;
   my @F = split(/$delim/, $_);
   my $id = shift @F;
   unless ($id =~ /^ENS/) { # print header line
      if ($desc) {
         print $OUT "$id${delim}GeneName${delim}Description$delim";
         print $OUT "Chromosome${delim}Start${delim}End${delim}Strand$delim" if ($coords);
         print $OUT join("$delim",@F),"\n";
      } else {
         print $OUT "$id${delim}GeneName${delim}",join("$delim",@F),"\n";
      }
      next;
   }
   ++$c;
   print "Fetching gene '$id' from Ensembl...\n" if $DEBUG;
   my $g;
   if ($feature =~ /transcript/i) {
      $g = $adaptor->fetch_by_transcript_stable_id($id)
   } else {
      $g = $adaptor->fetch_by_stable_id($id); # fetch it from ensembl
   }
   my $name;
   my $gDesc;
   my $chr = '-';
   my $start = '-';
   my $end = '-';
   my $strand = '-';
   
   if (!defined($g)) {  # check gene exists
      warn "Warning - gene '$id' not found\n";
      $name = 'unknown';
      $gDesc = 'none';
   } else {
      $name = $g->external_name(); # get the name
      $name = 'unknown' unless ($name);
      if ($desc) {
         $gDesc = $g->description();
         $gDesc = 'none' unless ($gDesc);
         $gDesc =~ s/\[Source.*\]//;
      }
      if ($coords) {
         $chr = $g->slice->seq_region_name();
         $start = $g->start();
         $end = $g->end();
         $strand = $g->strand();
      }
      
   }
   if ($desc) {
      print $OUT "$id$delim$name$delim$gDesc$delim";
      print $OUT "$chr$delim$start$delim$end$delim$strand$delim" if ($coords);
      print $OUT join("$delim",@F),"\n";
   } else {
      print $OUT "$id$delim$name$delim",join("$delim",@F),"\n";
   }
   
}
close($IN);
if ($c == 0) {
   warn "Warning - no gene names found. Nothing changed\n";
   unlink($out);
   exit;
}
#rename($out, $fasta) or die "Failed to rename '$out' to '$fasta'\n";
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

add_gene_name_column.pl --in <file> --species <name> [--feature-type <string>] [--desc|--no-desc] [--coords|--no-coords] [--delim <string>] [--out <file>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Annotate delimited datafiles with addtional gene information from ensembl.

This script is dependent on the ensembl perl API and requires the path to be set in PERL5LIB e.g.

   export PERL5LIB=/sw/opt/ensembl-api/72/ensembl/modules:/homes/ccole/lib:/sw/opt/ensembl-api/bioperl-live


=head1 OPTIONS

=over 5

=item B<--in>

Input tab-delimited file. 1st column must be an ensembl ID.

=item B<--species>

Species name.

=item B<--feature-type>

Specify whether the Ensembl IDs correspond to genes or transcripts. [default: Genes]

=item B<--desc|--no-desc>

Toggle whether to add gene description as well as gene name. [default: off]

=item B<--coords|--no-coords>

Toggle whether to include genomic coordinates (assumes --desc is set). [default: off]

=item B<--delim>

Specify the field delimiter in the input and output files. [default: tab]

=item B<--out>

Output filename. [default: STDOUT]

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