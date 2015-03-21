# ensemblUtils
Utility scripts which use the Ensembl perl-API

Herein are a bunch of scripts that I've used over the years to perform 
various tasks with the Ensembl API.

I'm putting them here as they may be useful to others.

All require the Ensembl API is already installed and available via 
PERL5LIB. More information regarding each script is available with
the --help and --man arguments e.g.

    perl add_gene_name_column.pl --man

## Annotate Data with Ensembl Information

Scenario: you a have a data file with ensembl genes, but you want more
information on the genes rather than just their IDs.

The add\_gene\_name_column.pl script will do this. By default, it will
add gene names as an extra column immediately to the right of the gene
IDs (in the first column). All other columns are shifted to the right.

    perl add_gene_name_column.pl --in data.tsv --species human --desc --out annotated.tsv

This example adds an additional column with the text description of the
gene.

## Make a More Dynamic Version of Data Files

Tab-delimited files are a perfect way to work or share data. However,
they are not the easiest to explore the data unless you are handy in R.
You could view them in Excel, but it has well-documented problems when 
working with gene IDs.

The create\_html_table.pl script creates a dynamic, searchable HTML 
table enriched with jQuery and dataTables goodness, including direct 
links to ensembl for known gene IDs (first column).

    perl create_html_table.pl --in data.tsv --title "My dynamic data" --check-ensembl --out "data.html"

This example enforces a check of the validity of ensembl gene IDs. Any
that are not valid are reported and no link to ensembl provided in the
final output.
