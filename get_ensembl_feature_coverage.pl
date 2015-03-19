#!/usr/bin/perl

=head1 NAME

get_ensembl_feature_coverage.pl - for a given species and chromosome generate a file of bases covered by Exons and Genes as defined by Ensembl

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $species = 'Human';
my $chromosome;
my $out;
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;

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
pod2usage(-msg => 'Please supply a chromosome.') unless ($chromosome);

print "Searching $species for chromosome $chromosome data in Ensembl....\n" if $VERBOSE;

my $coverage = ensemblFeatureCoverage($chromosome, $species);
die "ERROR - no data found for $species $chromosome. Are these valid?\n" unless (scalar keys %$coverage);

$out = "${species}_${chromosome}_ensembl_coverage.csv";
print "Printing out data to $out...\n" if $VERBOSE;
open(my $OUT, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
print $OUT "Strand\t\tPosition\tExon\n";
foreach my $id (keys %$coverage) {
   if ($id =~ /^(\w{3})(\d+)$/) {
      print $OUT "$1\t$2";
   } else {
      die "ERROR - position id '$id' is not valid!\nDied "
   }
   
   if ($coverage->{$id}) {
      print $OUT "\ty"
   } else {
      print $OUT "\tn"
   }
   print $OUT "\n";
}
printf "%d bases covered by genes or exon definitions.\n", scalar keys %$coverage;
close($OUT);


## create a hash with all positions covered by genes
## or exons. Then this can be used to determine whether
## a peak is within a feature or not.
sub ensemblFeatureCoverage {
   #### preload the EnsEMBL database ####
   # this is where the API code is installed
   # on nanna and the cluster.
   # Using the registry simplifies access 
   # to the API greatly.
   use lib '/opt/perl/bioperl-live';
   use lib '/opt/perl/ensembl/modules';
   use Bio::EnsEMBL::Registry;
   
   my $registry = 'Bio::EnsEMBL::Registry';
   
   $registry->load_registry_from_db(
       -host => 'ensembldb.ensembl.org',
       -user => 'anonymous'
   );
   ######################################

   my $chr = shift;
   my $species = shift;
   $species = 'Human' unless defined($species);
   
   $chr =~ s/chr//;
   
   my $gene_adaptor  = $registry->get_adaptor( $species, 'Core', 'Gene' );
   die "ERROR - unable to connect to Ensembl for genes\n" unless ($gene_adaptor);
   my $exon_adaptor  = $registry->get_adaptor( $species, 'Core', 'Exon' );
   die "ERROR - unable to connect to Ensembl for exons\n" unless ($exon_adaptor);
   my $slice_adaptor  = $registry->get_adaptor( $species, 'Core', 'slice' );
   die "ERROR - unable to connect to Ensembl for slice\n" unless ($slice_adaptor);
   
   my $slice = $slice_adaptor->fetch_by_region('chromosome',$chr) or die "ERROR - failed to fetch slice for chromosome: $chr";
   
   
   
   # First get gene data...
   my $genes = $gene_adaptor->fetch_all_by_Slice($slice);
   printf "Found %d genes in Ensembl for chromosome $chr\n", scalar @$genes if $DEBUG;
   
   my $fwd = 0;
   my $rev = 0;
   my %data;
   foreach my $g (@$genes) {
      my $id = $g->stable_id();
      die "ERROR - no ID found for gene\n", unless ($id);
      my $strand = $g->strand();
      die "ERROR - no strand found for gene '$id'\n", unless ($strand);
      my $start = $g->start();
      die "ERROR - no start found for gene '$id'\n", unless ($start);
      my $end = $g->end();
      die "ERROR - no end found for gene '$id'\n", unless ($end);
      
      # convert strand definition to local form
      if ($strand == 1) {
         $strand = 'fwd';
         ++$fwd;
      } else {
         $strand = 'rev';
         ++$rev;
      }
      
      my $pos = $start;
      while($pos <= $end) {
         my $uid = $strand.$pos;
         $data{$uid} = 0;
         ++$pos;
      } 
   }
   print "Found $fwd forward and $rev reverse genes\n" if $DEBUG;
   
   # Then, get exon data...
   my $exons = $exon_adaptor->fetch_all_by_Slice($slice);
   printf "Found %d exons in Ensembl for chromosome $chr\n", scalar @$exons if $DEBUG;
   
   $fwd = 0;
   $rev = 0;
   foreach my $e (@$exons) {
      my $id = $e->stable_id();
      die "ERROR - no ID found for gene\n", unless ($id);
      my $strand = $e->strand();
      die "ERROR - no strand found for gene '$id'\n", unless ($strand);
      my $start = $e->start();
      die "ERROR - no start found for gene '$id'\n", unless ($start);
      my $end = $e->end();
      die "ERROR - no end found for gene '$id'\n", unless ($end);
      
      # convert strand definition to local form
      if ($strand == 1) {
         $strand = 'fwd';
         ++$fwd;
      } else {
         $strand = 'rev';
         ++$rev;
      }
      
      my $pos = $start;
      while($pos <= $end) {
         my $uid = $strand.$pos;
         $data{$uid}= 1;
         ++$pos;
      } 
   }
   print "Found $fwd forward and $rev reverse exons\n" if $DEBUG;
   return(\%data);
   
}


=head1 SYNOPSIS

get_ensembl_feature_coverage.pl --chromosome <string> [--species <string>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION


=head1 OPTIONS

=over 5

=item B<--chromosome>

Specify which chromosome to search.

=item B<--species>

Specify a species. [default: Human]

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