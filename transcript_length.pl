#!/usr/bin/perl

=head1 NAME

transcript_length.pl - output the transcript lengths for the longest transcript in each gene

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;

use FindBin qw($Bin);
use lib "$Bin/lib";
use ensembl;
use Bio::EnsEMBL::ApiVersion;

my $out = 'length.csv';
my $species;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '0.9';

GetOptions (
   'out=s'     => \$out,
   'species=s' => \$species,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a species name.') unless ($species);

# load ensembl object
my $ens = ensembl->new(species => $species, VERBOSE => $VERBOSE);
print "Species: ", $ens->species, "\n" if $VERBOSE;

# connect to ensembl and do some checks
my $reg = $ens->connect();
printf "NOTE: using Ensembl API version %s\n", software_version() if $VERBOSE;
my $gene_adaptor = $reg->get_adaptor($species, 'Core', 'Gene');
my $trans_adaptor = $reg->get_adaptor($species, 'Core', 'Transcript');
die "ERROR - failed to get adaptor for '$species'. Check spelling and that it's a valid Ensembl species. Or check that you're using the correct API.\n" unless (defined($gene_adaptor));
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($reg->version_check($reg->get_DBAdaptor($species, 'core')));

my $genes = $gene_adaptor->fetch_all();
printf "Found %d genes from '$species'\n", scalar @$genes if $VERBOSE;

open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
foreach my $g (@$genes) {
   printf "Gene: %s\n", $g->stable_id() if $DEBUG;
   my $transcripts = $trans_adaptor->fetch_all_by_Gene($g);
   printf "Found %d transcripts\n", scalar @$transcripts if $DEBUG;
   my $max = 0;
   foreach my $t (@$transcripts) {
      printf "Transcript: %s\n", $t->stable_id() if $DEBUG;
      $max = $t->length() if ($t->length() > $max)
   }
   printf $OUT "%s\t$max\n", $g->stable_id();
}



=head1 SYNOPSIS

transcript_length.pl --species <file> [--out <file>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION


=head1 OPTIONS

=over 5

=item B<--species>

Species to query.

=item B<--out>

Output filename. [default: 'length.csv']

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