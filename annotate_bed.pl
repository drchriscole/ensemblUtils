#!/usr/bin/perl

=head1 NAME

annotate_bed.pl - annotate a BED file with gene name

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
my $species = 'human';
my $build = 'GRCh38';
my $out = 'annotated.bed';
my $unique = 1;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '0.4';

GetOptions (
   'in=s'      => \$file,
   'species=s' => \$species,
   'genome-build=s' => \$build,
   'unique-names!' => \$unique,
   'out=s'     => \$out,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') unless ($file && -s $file);

# load ensembl object
my $ens = ensembl->new(species => $species, genomeBuild => $build, VERBOSE => $VERBOSE);
print "Species: ", $ens->species, "\n" if $VERBOSE;

# connect to ensembl and do some checks
my $reg = $ens->connect();
printf "NOTE: using Ensembl API version %s\n", software_version() if $VERBOSE;
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($reg->version_check($reg->get_DBAdaptor($species, 'core')));

# load adaptors
my $slice_adaptor = $reg->get_adaptor($species, "core", "slice");
my $gene_adaptor = $reg->get_adaptor($species, "core", "gene");

# read in BED file
my %geneNames;
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
open(my $BED, "<", $file) or die "ERROR - unable to open '$file': ${!}\nDied";
while(<$BED>) {
   ## skip track and browser lines
   if (/^track/) {
      print $OUT $_;
      next;
   }
   if (/^browser/) {
      print $OUT $_;
      next;
   }
   chomp();
   my ($chrom, $chromStart, $chromEnd, @rest) = split (/\t/);
   die "ERROR - missing required column 'chromEnd'. Check '$file' is a valid BED file\n" unless ($chromEnd);
   
   $chrom =~ s/chr//;  # remove any chromosome prefix
   
   # assume fwd strand unless otherwise defined
   my $strand = 1;
   if (defined($rest[2]) && $rest[2] eq '-') {
      $strand = -1
   } 
   # now define ensembl slice and retrieve genes
   my $slice = $slice_adaptor->fetch_by_region('chromosome', $chrom, $chromStart, $chromEnd, $strand);
   warn "Warning - region '' at line $. is unrecognised\n" unless (defined($slice));
   my $genes = $gene_adaptor->fetch_all_by_Slice($slice);
   my $geneStr = 'Intergenic';
   if (scalar @$genes) {
      my @gids;
      foreach my $g (@$genes) {
         push @gids, $g->external_name();
      }
      $geneStr = join(":",@gids);
   }
   $geneNames{$geneStr}++;
   
   # Write out
   print $OUT "$chrom\t$chromStart\t$chromEnd";
   
   # uniquify the gene names (or not)
   if ($unique) {
      printf $OUT "\t${geneStr}_%d", $geneNames{$geneStr};
   } else {
      print $OUT "\t$geneStr";
   }

   # Overwrite existing name column (if there is one), but keep other columns
   if (scalar @rest > 1) {
      my $null = pop @rest;
      printf $OUT "\t%s", $geneNames{$geneStr}, join("\t",@rest);
   }
   print $OUT "\n";
}
close($BED);
close($OUT);


=head1 SYNOPSIS

annotate_bed.pl --in <file> [--species <speices>] [--genome-build <build>] [--out <file>] [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Script to take a BED file and annotate the name (4th) column with the gene(s) that each line in the bed file overlaps with in the genome.

NB: any existing name column information will be overwritten

=head1 OPTIONS

=over 5

=item B<--in>

Input file in BED format.

=item B<--species>

Species name as recognisable by ensembl [default: human]

=item B<--genome-build>

Genome build - only relevant for human [ default: GRCh38]

=item B<--out>

Output filename. [default: annotate.bed]

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

Copyright 2016, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut