#!/usr/bin/perl

=head1 NAME

get_full_transcripts.pl - retrieve complete set of full-length transcript sequences from ensembl

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


my $species = 'human';
my $out = 'transcripts.fasta';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '0.3';

GetOptions (
   'species=s' => \$species,
   'out=s'     => \$out,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);

my $ens = ensembl->new(species => $species);

print "Species: ", $ens->species, "\n";

my $reg = $ens->connect();
printf "NOTE: using Ensembl API version %s\n", software_version();
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($reg->version_check($reg->get_DBAdaptor($species, 'core')));

my $trans_adaptor = $reg->get_adaptor($species, 'Core', 'Transcript');
die "ERROR - failed to get adaptor for '$species'. Check spelling and that it's a valid Ensembl species. Or check that you're using the correct API.\n" unless (defined($trans_adaptor));

my $trans = $trans_adaptor->fetch_all();
printf "Found %d transcripts\n", scalar @$trans;

open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
my $n = 0;
foreach my $t (@$trans) {
   last if ($n == 10);
   my $bseq = $t->seq();
   printf $OUT ">%s %s:%s:%s:%s:%s:%s\n%s\n", $t->stable_id, $t->external_name(), $t->biotype, $t->seqname, $t->start, $t->end, $t->strand,$bseq->seq;
   ++$n;
}
close($OUT);

=head1 SYNOPSIS

get_full_transcripts.pl [--species <name>] [--out <file>] [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Ensembl helpfully provide pre-generated files containing sequence datasets, but they don't provide one with full-length transcripts. The nearest they do is cDNAs which only include protein-coding genes, but doesn't include the UTR sequences.

When performing a transcript quantification experiment (e.g. with RNA-seq) a dataset complete with UTRs and non-protein coding genes is required. This scripts fills this gap.

=head1 OPTIONS

=over 5

=item B<--species>

Species name as understood by Ensembl. [Default: 'human']

=item B<--out>

Output filename. [default: 'transcripts.fasta']

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

Chris Cole <christian@cole.name>

=head1 COPYRIGHT

Copyright 2012, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut