#!/usr/bin/perl

=head1 NAME

overlapping_genes.pl - count number of overlapping genes

=cut

use strict;
use warnings;

use Getopt::Long qw(:config auto_version);
use Pod::Usage;
use FindBin qw($Bin);
use lib "$Bin/lib";
use ensembl;
use Bio::EnsEMBL::ApiVersion;

my $species = 'human';
my $chromosome = '21';
my $out;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '0.2';

GetOptions (
   'species=s' => \$species,
   'chromosome=s' => \$chromosome,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);


# load ensembl object
printf "NOTE: using Ensembl API version %s\n", software_version() if $VERBOSE;
my $ens = ensembl->new(species => $species, VERBOSE => $VERBOSE);
print "Species: ", $ens->species, "\n" if $VERBOSE;

# connect to ensembl and do some checks
my $registry = $ens->connect();
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($registry->version_check($registry->get_DBAdaptor($species, 'core')));


## clean chromosome definition
$chromosome =~ s/chr//;

## load ensembl adaptors
my $gene_adaptor  = $registry->get_adaptor( $species, 'Core', 'Gene' );
die "ERROR - unable to connect to Ensembl for genes\n" unless ($gene_adaptor);
my $exon_adaptor  = $registry->get_adaptor( $species, 'Core', 'Exon' );
die "ERROR - unable to connect to Ensembl for exons\n" unless ($exon_adaptor);
my $slice_adaptor  = $registry->get_adaptor( $species, 'Core', 'slice' );
die "ERROR - unable to connect to Ensembl for slice\n" unless ($slice_adaptor);

## fetch chromosome slice
my $slice = $slice_adaptor->fetch_by_region('chromosome',$chromosome) or die "ERROR - failed to fetch slice for chromosome: $chromosome";
   
# get genes on slice
my $genes = $gene_adaptor->fetch_all_by_Slice($slice);
my $nGenes = scalar @$genes;
print "Found $nGenes genes in Ensembl for chromosome $chromosome\n" if $VERBOSE;

# now for each gene...
my %geneOverlaps;
my $sameStrand = 0;
for (my $i = 0; $i < $nGenes; ++$i) {
   for (my $j = $i + 1; $j < $nGenes; ++$j) {
      if ($genes->[$i]->overlaps($genes->[$j])) { # see if it overlaps...
         printf "%s overlaps with %s\n", $genes->[$i]->stable_id, $genes->[$j]->stable_id if $DEBUG;
         next unless (doExonsOverlap($genes->[$i],$genes->[$j])); # check that at least one exon overlaps...
         
         # store counts for genes which overlap an the exon level
         $geneOverlaps{$genes->[$i]->stable_id}++;
         $geneOverlaps{$genes->[$j]->stable_id}++;
         ++$sameStrand if ($genes->[$i]->strand eq $genes->[$j]->strand);
      }
   }
}

printf "Found %d genes with overlaps. $sameStrand are on same strand\n", scalar keys %geneOverlaps;

# take two ensembl gene objects and compare their exons to see if
# any pair overlaps.
#
# return True at first overlapping pair
# return False if no exon pairs overlap
#
sub doExonsOverlap {
   my $g1 = shift;
   my $g2 = shift;
   
   # get exons for each gene and see if they overlap
   my $e1s = $g1->get_all_Exons();
   my $e2s = $g2->get_all_Exons();
   
   foreach my $e1 (@$e1s) {
      foreach my $e2 (@$e2s) {
         if ($e1->overlaps($e2)) {
            printf "%s overlaps %s\n", $e1->stable_id(), $e2->stable_id() if $DEBUG;
            return(1);
         }
      }
   }
   print "No exon overlaps\n" if $DEBUG;
   return(0);
}

=head1 SYNOPSIS

overlapping_genes.pl [--species <name>] [--chromosome <str>] [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Quick script to count overlapping genes for a given chromosome and species. 

Overlaps are defined by having at least one exon in each gene which overlaps by whatever criterion is used by Ensembl's API Bio::EnsEMBL::Feature::overlaps() method.   

=head1 OPTIONS

=over 5

=item B<--species>

Species name as understood by Ensembl core database. [default: human]

=item B<--chromosome>

Chromosome name. [default: 21]

=item B<--out>

Output filename. [default: STDOUT]

=item B<--version>

Report version information and exit

=item B<--verbose|--no-verbose>

Toggle verbosity. [default:none]

=item B<--debug|--no-debug>

Toggle debugging output. [default:none]

=item B<--help>

Brief help.

=item B<--man>

Full manpage of program.

=back

=head1 TODO

Short list of immediate improvements to make.

   - Return the list of genes (and exons) which are overlapping. Maybe as BED file?
   - Provide option to do whole genome rather than single chromosomes only
   - Expand to use Ensembl Genomes species

=back

=head1 AUTHOR

Chris Cole <c.cole@dundee.ac.uk>

=head1 COPYRIGHT

Copyright 2016, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut