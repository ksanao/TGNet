#!/usr/bin/env perl -w
# 
# input: blat output (query ESTs against assembly)
# discard header: awk 'FNR>5' file.psl
# output: filtered blat output by %of match length/EST length & max number of mismatched bases
#
# OR 100318
# OR 100319 added opt, added filtering by EST size
# OR 1004.. add filter by exon size and number
# OR 100409 make a possibility to filter by match bp length or match fraction
# OR 100510 added filter by average gap in genome (column 8/column 7)
# OR 100511 changed deafault limit value (now set to 0, limit deactivated)
# TESTED 100511 on a file test_blat_filterResults

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max);

#default option values
my $min_match_length=0;
my $min_match_fraction=0;
my $max_mismatch=0.05;
my $max_gap=20000;
my $min_EST_length=300;
my $min_exon_size=30;
my $min_exon_nb=1;
# a limit for duplicate matches, must be within this distance from the end/beginning of a scaffold, e.g 20000
my $limit=0;



# optionName|alternativeName
my $result=GetOptions(
	"matchlen|l=i"      => \$min_match_length,  
	"matchfraction|f=f" => \$min_match_fraction,  
	"mismatch|m=f"      => \$max_mismatch,
	"maxgap|g=i"        => \$max_gap,
	"estlen|e=i"        => \$min_EST_length,
	"minexon|x=i"       => \$min_exon_size,
	"minexonNb|b=i"     => \$min_exon_nb,
	"limit|t=i"         => \$limit,
	"help"              => sub { HelpMessage() },
	);

while(<>) 
{ 
	chomp $_;         
	my @values = split(/\t/,$_); 

	my $avg_gap=0;
	$avg_gap=$values[7]/$values[6], if $values[6]>0;
			
	if($values[0]>=$min_match_fraction*$values[10]
	&& $values[0]>=$min_match_length
	&& $values[1]/$values[0]<=$max_mismatch 
	&& $values[10]>=$min_EST_length
	&& $avg_gap<=$max_gap)
	{
		#keep only matches with at least 1 exon bigger than threshold and with a min nb of exons
		if(max(split(/,/, $values[18]))>=$min_exon_size 
		&& $values[17]>=$min_exon_nb)
		{
			#keep only matches on extremities, defined by $limit
			if(!$limit || ($values[15]<=$limit || ($values[14]-$values[16])<=$limit))
			{
				print $_."\n";
			}
		}
	}

}

#col	Value (-1 for an array position)
#1 match
#2 mismatch
#3 matches in repeat regions
#4 N bases
#5 EST gaps (insert number)
#6 EST bases in gaps
#7 Scaf/contig gaps (insert number)
#8 Scaf/contig bases in gaps
#9 strand
#10 EST name
#11 EST size
#12 EST start
#13 EST end
#14 Scaf/contig name
#15 Scaf/contig size
#16 Scaf/contig start
#17 Scaf/contig end
#18 block count
#19 Comma-separated list of sizes of each block 
#20 EST block start positions
#21 Scaf/contig block start positions

#   1. matches - Number of bases that match that aren't repeats
#   2. misMatches - Number of bases that don't match
#   3. repMatches - Number of bases that match but are part of repeats
#   4. nCount - Number of 'N' bases
#   5. qNumInsert - Number of inserts in query
#   6. qBaseInsert - Number of bases inserted in query
#   7. tNumInsert - Number of inserts in target
#   8. tBaseInsert - Number of bases inserted in target
#   9. strand - '+' or '-' for query strand. In mouse, second '+'or '-' is for genomic strand
#  10. qName - Query sequence name
#  11. qSize - Query sequence size
#  12. qStart - Alignment start position in query
#  13. qEnd - Alignment end position in query
#  14. tName - Target sequence name
#  15. tSize - Target sequence size
#  16. tStart - Alignment start position in target
#  17. tEnd - Alignment end position in target
#  18. blockCount - Number of blocks in the alignment
#  19. blockSizes - Comma-separated list of sizes of each block
#  20. qStarts - Comma-separated list of starting positions of each block in query
#  21. tStarts - Comma-separated list of starting positions of each block in target 
sub HelpMessage
{
	print "\nUsage: \n\tawk 'FNR>5' file.psl | blat_filterResults.pl [options] \n\t--matchlen      | -l : minimum match length (bp, default 0)\n\t\t\t\trefers to a number of bases that aren't repeats\n\t--matchfraction | -f : minimum match length as a fraction of EST size (floating, default 0)\n\t\t\t\trefers to a number of bases that aren't repeats\n\t--mismatch      | -m : maximum allowed mismatch residues as a fraction of match length (floating, default 0.05)\n\t--maxgap        | -g : maximum average gap in genomic sequence (bp, default 20kb) \n\t--estlen        | -e : minimum length of macthed EST (integer, default 300)\n\t--minexon       | -x : at least 1 exon should be bigger than this threshold (bp, default 30)\n\t--minexonNb     | -b : set the minimum number of exons per EST match (integer, default 1)\n\t--limit         | -t : duplicate matches must be within this distance from the end/beginning of a scaffold \n\t\t\t\t(bp, default deactivated, limit set to 0)\n\n";
	exit 0;
}

exit 0;
