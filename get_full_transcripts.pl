#!/usr/bin/perl

=head1 NAME

get_full_transcripts.pl - retrieve complete set of full-length transcript sequences from ensembl

=cut

use strict;
use warnings;

use Time::HiRes qw(gettimeofday tv_interval);
use Getopt::Long qw(:config auto_version);
use Pod::Usage;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/lib";
use ensembl;
use Bio::EnsEMBL::ApiVersion;

$| = 1;
my $species = 'human';
my $out = 'transcripts.fasta';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;
our $VERSION = '0.5';

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

# load ensembl object
my $ens = ensembl->new(species => $species, VERBOSE => $VERBOSE);
print "Species: ", $ens->species, "\n" if $VERBOSE;

# connect to ensembl and do some checks
my $reg = $ens->connect();
printf "NOTE: using Ensembl API version %s\n", software_version() if $VERBOSE;
warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($reg->version_check($reg->get_DBAdaptor($species, 'core')));

# get transcript set
my $trans_adaptor = $reg->get_adaptor($species, 'Core', 'Transcript');
die "ERROR - failed to get adaptor for '$species'. Check spelling and that it's a valid Ensembl species. Or check that you're using the correct API.\n" unless (defined($trans_adaptor));
my $trans = $trans_adaptor->fetch_all();
my $tot = scalar @$trans;
die "ERROR - no transcripts found\n", unless ($tot);
printf "Found %d transcripts\n", $tot if $VERBOSE;

# iterate through transcripts and print them out
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
my $n = 0;
my $start_time = [gettimeofday];
foreach my $t (@$trans) {
   last if ($DEBUG && $n == 10);
   
   TimeRemaining($start_time, $n / $tot) if $VERBOSE;
   
   my $bseq = $t->seq();
   printf $OUT ">%s %s:%s:%s:%s:%s:%s\n%s\n", $t->stable_id, $t->external_name(), $t->biotype, $t->seqname, $t->start, $t->end, $t->strand,$bseq->seq;
   ++$n;
}
close($OUT);

# code from Marek to report time elapsed and remaining.
sub TimeRemaining {
  my ($start, $frac) = @_;

  my $elapsed = tv_interval($start);
  my $remaining = ($frac > 0) ? $elapsed * (1 - $frac) / $frac : 0;

  my $ela = hms($elapsed);
  my $rem = hms($remaining);

  my $s = sprintf "%5.1f%%  Elapsed %8s  Remaining %8s", 100*$frac, $ela, $rem;
  my $back = "\b" x length($s);
  print $back . $s;
}

# convert seconds into hours, minutes, seconds
sub hms {
  my $time = shift;
  my $h = int($time / 3600);
  my $m = int(($time - 3600*$h) / 60);
  my $s = int($time - 3600*$h - 60*$m);
  sprintf "%02d:%02d:%02d", $h, $m, $s;
}

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

Copyright 2016, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut