#!/usr/bin/env perl
# EPN, Mon May 11 10:32:44 2015
#
# dnaorg_summarize_parsed_ftable.pl: companion script for parsing output
# of dnaorg_parse_ftable.pl.
# 
use strict;
use warnings;
use Getopt::Long;

## Example 'parsed ftable', created by dnaorg_parse_ftable.pl:
## NC_010435.CDS.tbl: 
##
##################################################
##full-accession	accession	coords	strand	min-coord	product	transl_table	protein_id	db_xref	gene
#gb|EU822322.1|	EU822322.1	118..465	+	118	AV2 protein	1	gb|ACF42108.1|	-	-
#gb|EU822322.1|	EU822322.1	complement(269..586)	-	269	AC5 protein	1	gb|ACF42109.1|	-	-
#gb|EU822322.1|	EU822322.1	278..1048	+	278	coat protein	1	gb|ACF42110.1|	-	-
#gb|EU822322.1|	EU822322.1	complement(1051..1455)	-	1051	replication enhancer protein	1	gb|ACF42111.1|	-	-
#gb|EU822322.1|	EU822322.1	complement(1148..1600)	-	1148	transcription activator protein	1	gb|ACF42112.1|	-	-
#gb|EU822322.1|	EU822322.1	complement(1497..2588)	-	1497	replication initiator protein	1	gb|ACF42113.1|	-	-
#gb|EU822322.1|	EU822322.1	complement(2138..2437)	-	2138	AC4 protein	1	gb|ACF42114.1|	-	-
#gb|EU822321.1|	EU822321.1	118..465	+	118	AV2 protein	1	gb|ACF42101.1|	-	-
#gb|EU822321.1|	EU822321.1	complement(269..586)	-	269	AC5 protein	1	gb|ACF42102.1|	-	-
#gb|EU822321.1|	EU822321.1	278..1048	+	278	coat protein	1	gb|ACF42103.1|	-	-
##################################################
#
# This script will 'summarize' a parsed ftable by giving a report on the number of 
# non-empty rows for each column.
#
# For the above parsed ftable:
#
##################################################
##qualifier      count  fract
#full-accession     10  1.000
#accession          10  1.000
#coords             10  1.000
#strand              4  0.400
#min-coord          10  1.000
#product            10  1.000
#transl_table       10  1.000
#protein_id         10  1.000
#db_xref             0  0.000
#gene                0  0.000
##################################################
##
# NOTE: you'd never see a '0' count in a real full parsed ftable output, but here I've truncated
# one so their are qualifiers that have 0 counts.
my $usage = "\ndnaorg_summarize_parsed_ftable.pl <.tbl file output from dnaorg_parse_ftable.pl>\n";
$usage .= "\n"; 

my $line = <>;
if($line !~ m/^\#/) { 
  die "ERROR, first line does not start with #";
}
$line =~ s/^\#//;
chomp $line;
my @qualifiers_A = split(/\s+/, $line);
my $nqualifiers = scalar(@qualifiers_A);
my @ct_A = ();
my $w_qualifier = length("qualifier") + 1;
for(my $i = 0; $i < $nqualifiers; $i++) { 
  $ct_A[$i] = 0; 
  if(length($qualifiers_A[$i]) > $w_qualifier) { $w_qualifier = length($qualifiers_A[$i]); }
}

my $nrows = 0;
while($line = <>) { 
  if($line !~ m/^\#/) { 
    chomp $line;
    my @values_A = split(/\t/, $line);
    if(scalar(@values_A) != $nqualifiers) { die "ERROR, didn't find $nqualifiers tab-delimited tokens in line: $line\n"; }
    for(my $i = 0; $i < $nqualifiers; $i++) { 
      if($values_A[$i] ne "-") { $ct_A[$i]++; }
    }
    $nrows++;
  }
}

printf("#%-*s  %5s  %5s\n", $w_qualifier-1, "qualifier", "count", "fract");
for(my $i = 0; $i < $nqualifiers; $i++) { 
  printf("%-*s  %5d  %5.3f\n", $w_qualifier, $qualifiers_A[$i], $ct_A[$i], ($ct_A[$i] / $nrows));
}
  
