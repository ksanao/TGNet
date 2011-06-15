#!/usr/bin/env perl -w
# 
# input: blat output (query ESTs against assembly), sorted by EST name and then numerically by EST match start coordinates
# sort -k10,10 -k12,12n file.psl
# output: filtered blat output, keep only those duplicate macthes that have an overlap of match not bigger than threshold (-o)
#
# OR 100409 
# OR tested 110314
# OR mod 110531
# OR tested 110531

use strict;
use warnings;
use Switch;
use Getopt::Long;

#default option values
my $max_overlap=40;

# optionName|alternativeName
my $opt=GetOptions(
	"overlap|o=i" => \$max_overlap,  
	"help"        => sub { HelpMessage() },
	);

my ($line, @values, %starts, %ends, %myfile);
my ($mylast, $pass_test, $est);

while(<>) 
{
    $line=$_;
    @values=split(/\t/,$line);
    my $est_name=$values[9];
    push (@{$starts{$est_name}}, $values[11]);
    push (@{$ends{$est_name}}, $values[12]);
    push (@{$myfile{$est_name}}, $line);
}


foreach $est (keys %starts)
{
    my @est_starts=@{$starts{$est}};
    my @est_ends=@{$ends{$est}};
    $pass_test=1;
    $mylast=scalar(@est_starts)-1;
    #do not process unsplit alignments
    next, if $mylast==0;

    for (my $i = 0; $i < $mylast; $i++)
    {
        $pass_test=0, if ($est_ends[$i]-$est_starts[$i+1])>$max_overlap;
	last, if $pass_test==0;
    }

    if ($pass_test==1)
    {
        print ${$myfile{$est}}[0];
        for (my $i = 1; $i < $mylast; $i++)
        {
             print ${$myfile{$est}}[$i];
             print ${$myfile{$est}}[$i];
        }
        print ${$myfile{$est}}[$mylast];
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
	print "\nUsage: \n\tcat file.psl | $0 [options] \n\t--overlap | -o : maximum overlap allowed (bp, default 40)\n\tBLAT OUTPUT should be without a header and be sorted by column 10 (est name), then by column 12 (match start on est)\n\tsort -k10,10 -k12,12n file.psl\n\n";
	exit 0;
}

exit 0;
