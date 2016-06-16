#!/usr/bin/perl

=head1 NAME

annotate_promoters.pl - annotate promoters of given genes with putative binding sites

=cut

use strict;
use warnings;

use Getopt::Long qw(:config auto_version);
use Pod::Usage;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/lib";
use ensembl;
use Bio::EnsEMBL::ApiVersion;

my $file;
my $qMotif = '[AG]CGTG'; # HRE motif
my $promoterLength = 20;
my $genome = 'GRCh37';
my $out = 'motif_matches.bed';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '1.0';

GetOptions (
   'in=s'      => \$file,
   'motif=s'   => \$qMotif,
   'length=i'  => \$promoterLength,
   'genome=s'  => \$genome,
   'out=s'     => \$out,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') unless ($file && -s $file);

# get gene list - HGNC symbols, one gene per line
my @genes = getGeneList($file);
die "ERROR - no genes found\n" unless (scalar @genes);
printf "Read %d genes in file '$file'\n", scalar @genes if $VERBOSE;

# load ensembl object
printf "NOTE: using Ensembl API version %s\n", software_version() if $VERBOSE;
my $ens = ensembl->new(genomeBuild => $genome, VERBOSE => $VERBOSE);
print "Species: ", $ens->species, "\n" if $VERBOSE;

# connect to ensembl and do some checks
my $registry = $ens->connect();
print "Genome version: ".$registry->get_adaptor('human', 'core', 'genomecontainer')->get_version()."\n" if $DEBUG;
my $sliceAdaptor = $registry->get_adaptor('human', 'core', 'slice');
my $geneAdaptor = $registry->get_adaptor('human', 'core', 'gene');


# for each gene get upstream sequence via ensembl
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
print $OUT "track name=\"Motif_BS\"\n";
foreach my $g (@genes) {
   my ($slice, $strand) = getPromoterSeq($geneAdaptor->fetch_by_display_label($g), $promoterLength);   
   
   # check strand and revComp if on reverse strand
   # only look at promoter part of seq - ensembl returns whole gene + up/down stream region
   my $seq = '';
   if ($strand eq '-1' ) {
      $seq = substr(reverseComplement($slice->seq()),0,$promoterLength);
   } else {
      $seq = substr($slice->seq(),0,$promoterLength);
   }
   
   # capture all motif sequence regions
   while ($seq =~ /($qMotif)/g) {
      my $motif = $1;
      my $len = length($motif);
      print "Match $1 found for gene '$g'\n" if $VERBOSE;
      my $idx = index($seq,$motif); # find where the match is
      
      # report matches is BED format: chr start end name score strand
      if ($strand eq '-1') {
         printf $OUT "%s\t%d\t%d\t$g:$motif\t100\t-\n", $slice->seq_region_name(), $slice->end()-$idx-$len, $slice->end()-$idx
      } else {
         printf $OUT "%s\t%d\t%d\t$g:$motif\t100\t+\n", $slice->seq_region_name(), $slice->start()+$idx, $slice->start()+$idx+$len
      }
   } 
}
close($OUT);

## for a given ensembl gene object and +/- padding region,
## return slice object and strand
sub getPromoterSeq {
   my $gene = shift;
   my $length = shift;
   
   die "ERROR - no ensembl gene found for '$gene'\n" unless (defined $gene);
   printf "Found %s stableID for '$gene'\n", $gene->stable_id() if $DEBUG;
   
   my $slice = $sliceAdaptor->fetch_by_gene_stable_id($gene->stable_id(), $length);
   printf "Locus for '$gene' is %s:%d:%d\n", $slice->seq_region_name(), $slice->start(), $slice->end() if $DEBUG;
   return($slice, $gene->strand());
}

## parse gene list 
sub getGeneList {
   my $file = shift;
   
   my @data;
   open(my $fh, "<", $file) or die "ERROR - unable to open '$file': ${!}\nDied";
   while(<$fh>) {
      chomp();
      next unless(length($_));
      push @data, $_;
   }
   close($fh);
   
   return(@data);
}

## generate reverse complement sequence 
sub reverseComplement {
   my $seq = shift;
   
   my $rev = reverse $seq;
   $rev =~ tr/ACGT/TGCA/;
   return($rev);
}

=head1 SYNOPSIS

annotate_promoters.pl --in <file> [--length <int>] [--genome <str>] [--out <file>] [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

This script identifies putative motifs & binding sites in the upstream/promoter region for a given list of genes.

All identified sites are reported in BED format.

=head1 KNOWN ISSUES

When using a version of the ensembl perl API which is more than two versions older than the most recent, all output will be GRCh38 specific regardless of which assembly specified via I<--version>.

=head1 OPTIONS

=over 5

=item B<--in>

Input file with list of genes.

=item B<--length>

Length of upstream region to search.

=item B<--genome>

Specify the genome build to use (GRCh37|GRCh38). [default: GRCh37]

=item B<--out>

Output filename. [default: motif_matches.bed]

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

=head1 AUTHOR

Chris Cole <c.cole@dundee.ac.uk>

=head1 COPYRIGHT

Copyright 2015, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut