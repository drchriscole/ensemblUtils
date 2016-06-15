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
our $VERSION = '0.1';

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

my $n = 0;
foreach my $t (@$trans) {
   last if ($n == 10);
   my $bseq = $t->seq();
   printf ">%s %s:%s:%s:%s\n%s\n", $t->stable_id, $t->biotype, $t->seqname, $t->start, $t->end,$bseq->seq;
   ++$n;
}


=head1 SYNOPSIS

get_full_transcripts.pl --in <file> [--out <file>] [--verbose|--no-verbose] [--version] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION


=head1 OPTIONS

=over 5

=item B<--in>

Input file.

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

=head1 AUTHOR

Chris Cole <christian@cole.name>

=head1 COPYRIGHT

Copyright 2012, Chris Cole. All rights reserved.

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut