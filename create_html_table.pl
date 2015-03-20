#!/usr/bin/perl

=head1 NAME

create_html_table.pl - given a tab-delimited file with ensembl identifiers create an HTML table with links to Ensembl

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Basename;

# requires that ensembl perl API is in $PATH
use Bio::EnsEMBL::Registry;

our $VERSION = '1.1_json.2';

my $file;
my $title = 'Ensembl Data';
my $out = 'out.html';
my $check = 0;
my $ensURL;
my $desc = 0;
my $species = 'human';
my $cols = '';
my $VERBOSE = 1;
my $DEBUG = 0;
my $help;
my $man;

GetOptions (
   'in=s'      => \$file,
   'out=s'     => \$out,
   'title=s'   => \$title,
   'check-ensembl!' => \$check,
   'species=s' => \$species,
   'pval-cols=s' => \$cols,
   'verbose!'  => \$VERBOSE,
   'debug!'    => \$DEBUG,
   'man'       => \$man,
   'help|?'    => \$help,
) or pod2usage();

pod2usage(-verbose => 2) if ($man);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please supply a valid filename.') if (!$file or !-e $file);

my $gene_adaptor;
if ($check) {
   my $registry = connectEnsemblRegistry($species);
   $gene_adaptor = $registry->get_adaptor($species, 'Core', 'Gene');
   die "ERROR - failed to get gene adaptor for '$species'. Check spelling and that it's a valid Ensembl species. Or check that you're using the correct API.\n" unless (defined($gene_adaptor));
   warn "Warning - API version check has failed. You probably need to update your local install.\n" unless ($registry->version_check($registry->get_DBAdaptor($species, 'core')));
}

# generate a URL for links to ensembl
$ensURL = constructURL($species);

print "Parsing input file...\n" if $VERBOSE;
my ($headers, $data);
if ($check) {
	($headers, $data) = parseTsv($file, $gene_adaptor)
} else {
	($headers, $data) = parseTsv($file);
}
printf "Found %d columns and %d data rows\n", scalar @$headers, scalar @$data if $VERBOSE;

print "Generating HTML table...\n" if $VERBOSE;
## open HTML file and print headers
## Includes funky JSON, jQuery and dataTables javascript goodness...
## The javascript allows better formatting of table data via pagination with
## user-selectable number of rows to display at a time. Also, the table is 
## searchable and sortable. By default the table is shown in its original 
## state; no additional sorting is done.
open(my $html, ">", $out) or die "ERROR - unable to open '$out' for write: ${!}\nDied";
print $html "<?xml version='1.0' encoding='utf-8' ?>
<!DOCTYPE html PUBLIC '-//W3C//DTD XHTML 1.0 Strict//EN' 'http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd'>
<html xmlns='http://www.w3.org/1999/xhtml'>
<head>
  <title>Ensembl data</title>

  <script type = 'text/javascript' src='http://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js'> </script>
  <script type = 'text/javascript' src='http://cdn.datatables.net/1.10.5/js/jquery.dataTables.min.js'> </script>
  <script type = 'text/javascript'>
     \$(document).ready(function() {
        
           /* two functions for sorting p-value data sensibly */
         jQuery.fn.dataTableExt.oSort['allnumeric-asc']  = function(a,b) {
            var x = parseFloat(a);
            var y = parseFloat(b);
            return ((x < y) ? -1 : ((x > y) ?  1 : 0));
         };
       
         jQuery.fn.dataTableExt.oSort['allnumeric-desc']  = function(a,b) {
            var x = parseFloat(a);
            var y = parseFloat(b);
            return ((x < y) ? 1 : ((x > y) ?  -1 : 0));
         };
         
         /* json data for tabulating */
         var myData = [
";

## contruct json data structure here.
foreach my $row (@$data) {
	printf $html "            ['%s'],\n", join("','", @$row);
}
          
print $html "         ];

         \$('#myTable').dataTable(
	       {
	          'data': myData,
	          'columns': [
";

## add column headers here
foreach my $head (@$headers) {
   print $html "	             { 'title': '$head' },\n";
}

## print rest of javascript and HTML
print $html "	          ],
	       	  'aoColumnDefs': [
	       	      { 'sType': 'allnumeric', 'aTargets': [ $cols ] }
	       	  ],
	       	  'aaSorting': [],
	           'aLengthMenu': [[25, 50, 100, 200, -1], [25, 50, 100, 200, 'All']],
	           'iDisplayLength': 25,
	           'sPaginationType': 'full_numbers'
	       }
	
	    );
	 } );
  </script>
  <link rel='stylesheet' type='text/css' href='http://www.compbio.dundee.ac.uk/user/ccole/css/datatable.css'/>
  <style type='text/css'>
     h1 {
        text-align: center;
     }
     #main {
        float: left;
        padding:  0px 30px 0px 30px;
     }
     caption {
        text-align: left;
     }
  </style>
</head>
<body>
<h1>$title</h1>
<div id='main'>
<table class = 'display' id = 'myTable' cellpadding = '5' width = '100%'>
<caption>Click column headers to sort rows</caption>
</table>\n</div>\n</body>\n</html>
";

close($html);
print "Done.\n" if $VERBOSE;
exit;

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

## the ensembl URL is dependent on species, so generate it here
sub constructURL {
   my $species = shift;
   
   $species =~ tr/A-Z/a-z/;
   my %known = (
      'human' => 'Homo_sapiens',
      'mouse' => 'Mus_musculus',
      'chicken' => 'Gallus_gallus'
   );
   
   if ($known{$species}) {
      return("http://www.ensembl.org/$known{$species}/Gene/Summary?g=")
   } else {
      warn "Warning - species '$species' is not recognised. Links to Ensembl won't work.\n          This is a minor problem and can be solved easily in the contructURL() function.\n";
      return('');
   }
}

## parse data
sub parseTsv {
	my $file = shift;
	my $ens_adaptor= shift;
	
	## read source data file and append HTML table data with it
	my @colHeaders;
	my @data;
	my $i = 0;
	open(my $IN, "<", $file) or die "ERROR - unable to open '$file': ${!}\nDied";
	while(<$IN>) {
		chomp;
		s/\"//g; # remove quotes;
		print join(":",split(/\t/, $_))."\n" if $DEBUG;
		if ($. == 1) { # header line - assume it exists
			@colHeaders = split(/\t/, $_);
		} else {
			my @F = split(/\t/, $_);
			my $id = $F[0];
			die "ERROR - wrong number of data columns at line $.\n" unless (scalar @F == scalar @colHeaders);
			
			my $validID = 1;
			if ($ens_adaptor) { # if adapter provided, check geneID is valid
				print "Fetching gene '$id' from Ensembl...\n" if $DEBUG;
				my $g = $ens_adaptor->fetch_by_stable_id($id);
				warn "Warning - gene ID '$id' is not known by ensembl\n";
				$validID = 0;
			}
			
			# if the geneID is valid turn it into an HTML <a> tag
			$F[0] = "<a href=\\'$ensURL$id\\'>$id</a>" if ($validID);
			
			# store data
			push @{$data[$i]},  @F;
			++$i;
		}
	}
	close($IN);
	return(\@colHeaders, \@data);
}

=head1 SYNOPSIS

create_html_table.pl --in <file> [--check-ensembl|--no-check-ensembl] [--title <string>] [--species <name>] [--pval-cols <string>] [--out <file>] [--verbose|--no-verbose] [--debug|--no-debug] [--man] [--help]

=head1 DESCRIPTION

Sometimes when sharing data with others you want more than just a simple tab-delimited file. You want the ability to sort, search and have links - especially to ensembl geneIDs.

This script solves this problem.

Take any tab-delimited file with ensembl geneIDs as the first column and creates a jQuery dataTable in a HTML page.

=head1 REQUIREMENTS

Requires the ensembl Perl API.

=head1 OPTIONS

=over 5

=item B<--in>

Input file. Tab-delimitted with Ensembl identifiers in first column

=item B<--check-ensembl|--no-check-ensembl>

Check identifiers against Ensembl. [default: off]

=item B<--title>

Title for web page.

=item B<--species>

Specify species. [default: human]

=item B<--pval-cols>

Provide a comma-separated list of column indexes (starting at 0) which have scientific notation data. [default: none]

=item B<--out>

Output file. [default: out.html]

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
