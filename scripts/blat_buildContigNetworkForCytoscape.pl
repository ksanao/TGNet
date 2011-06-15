#!/usr/bin/env perl
# 
# input: blat output (query ESTs against assembly) only 2xmatches
# output: filtered blat output, keep only those duplicate macthes that have an overlap of match not bigger than 40bp
#
# OR 100413 
# OR 100709 added comments

use strict;
use warnings;
use Switch;
use Getopt::Long;
use List::Util qw(max min);

#default option names
my $testit=1;

# optionName|alternativeName
my $opt=GetOptions(
    "infile|i=s"  => \my $infile,
    "test|t"      => \$testit,
    "help"        => sub { HelpMessage() },
    );

my ($i, %scaffolding, @values1, @values2);

while(<>) 
{ 
    chomp ($_);
    my @values=split(/\t/,$_);
    #1 line contains 4 elements: scaffold start end contig
    #contig is a key, other elements are values array assigned as a reference
    #values array [0]=scaffold [1]=start [2]=end
    my $key=pop(@values);
    $scaffolding{$key}=\@values;
}

open (IN, "<", $infile)  || die "fail open $! $infile\n";
$infile=~s/.psl$//;
my $file=$infile.".est_network";
open (OUT, ">", $file)  || die "fail open $! $file\n";
print OUT "EST\tESTsize\tstrand1\tESTmatch1_start\tESTmatch1_end\tcontig_name1\tscaf1\tscaf1start\tscaf1end\tstrand2\tESTmatch2_start\tESTmatch2_end\tcontig_name2\tscaf2\tscaf2start\tscaf2end\tdistanceOnScaffold\ttag\n";

while (!eof(IN))
{
    my $line=<IN>;
    chomp ($line);
    #get every two lines
    $i=$.%2;
    
    switch ($i)
    {
        case (1) 
        {
            @values1 = split(/\t/,$line);
        }
        case (0) 
        {
            @values2 = split(/\t/,$line);
            CompareScaffolding(\@values1, \@values2, \%scaffolding, 0);
        }
    }
}

close IN;
close OUT;



###FUNCTIONS

sub CompareScaffolding
{
    my ($values1_ref, $values2_ref, $scaffolding_ref, $do_test)=(shift, shift, shift, shift);
#   @{$values1_ref}
#   ${$hash_ref}{KEY}
    my @print_array;
    
    #check if this is a duplicate match 
    #column 10/array position 9 should contain same EST name
    die "Fail to get duplicate match lines (EST name is not the same)!\n", unless ${$values1_ref}[9] eq ${$values2_ref}[9];
    
    #check if contigs are used in scaffolding
    #check if duplicate match occurs on different contigs
    if (exists ${$scaffolding_ref}{${$values1_ref}[13]} && exists ${$scaffolding_ref}{${$values2_ref}[13]} && ${$values1_ref}[13] ne ${$values2_ref}[13])
    {
        my @scaffolding1=@{${$scaffolding_ref}{${$values1_ref}[13]}};
        my @scaffolding2=@{${$scaffolding_ref}{${$values2_ref}[13]}};
        
        #                   EST(same in 1&2)    ESTsize(same in 1&2) strand1(ESTmatch)    ESTmatch1_start     ESTmatch1_end        contig_name1        scaff1/start/end   strand2(ESTmatch)  ESTmatch2_start     ESTmatch2_end        contig_name2  scaff2/start/end      scaf2_size
        push (@print_array, ${$values1_ref}[9], ${$values1_ref}[10], ${$values1_ref}[8], ${$values1_ref}[11], ${$values1_ref}[12], ${$values1_ref}[13], @scaffolding1,     ${$values2_ref}[8], ${$values2_ref}[11], ${$values2_ref}[12], ${$values2_ref}[13], @scaffolding2);
        
        #are contigs assigned ot the same scaffold?
        my $same_scaffold=$scaffolding1[0] eq $scaffolding2[0];
        #do matches have the same orientation?
        my $same_strand=${$values1_ref}[8] eq ${$values2_ref}[8];
        #is the order of contigs correct?
        #if + strand, then max end coordinate of match should correspond to a maximum coordinate of scaffold
        #if - strand then max end coordinate of match should correspond to a minimum coordinate of scaffold
        my $correct_order=0;
        my $distanceOnScaf=0;
        if ($same_scaffold && $same_strand)
        {
            if (${$values1_ref}[8] eq "+" 
             && ((${$values1_ref}[12]>${$values2_ref}[12] && $scaffolding1[2]>$scaffolding2[2])
             ||  (${$values1_ref}[12]<${$values2_ref}[12] && $scaffolding1[2]<$scaffolding2[2])))
            {
                $correct_order=1;
            }
            elsif (${$values1_ref}[8] eq "-" 
               && ((${$values1_ref}[12]>${$values2_ref}[12] && $scaffolding1[2]<$scaffolding2[2])
               ||  (${$values1_ref}[12]<${$values2_ref}[12] && $scaffolding1[2]>$scaffolding2[2])))
            {
                $correct_order=1;
            }

            #get the distance between last match on Contig1 and first match on Contig2
            if ($correct_order)
            {
                if ($scaffolding1[2]>$scaffolding2[2])
                {
                    #!! blat start coordinates are -1 adjsuted (start on 1 is noted as 0)
                    my $first=$scaffolding2[1]+${$values2_ref}[16];
                    my $last=$scaffolding1[1]+${$values1_ref}[15];
                    $distanceOnScaf=$last-$first;
                }
                else
                {
                    #!! blat start coordinates are -1 adjsuted (start on 1 is noted as 0)
                    my $first=$scaffolding1[1]+${$values1_ref}[16];
                    my $last=$scaffolding2[1]+${$values2_ref}[15];
                    $distanceOnScaf=$last-$first;
                }   
            }
        }
        
        
        #scaffolding is confirmed if 
        #1:contigs are assigned to same scaffold,
        #2:match has same orientation 
        #3:contigs are positioned in a coorect order
        my $error_type="NA";
        
        if ($same_scaffold && $same_strand && $correct_order)
        {           
            print OUT join("\t", @print_array)."\t$distanceOnScaf\tOK\n", if !$do_test;
        }
        else
        {   
            if ($same_scaffold)
            {
                if ($same_strand)
                {
                    $error_type="ORDER";
                }
                else
                {
                    $error_type="ORIENTATION";
                }
            }
            else
            {
                $error_type="SCAFFOLD";
            }
    
            print OUT join("\t", @print_array)."\t$distanceOnScaf\t$error_type\n", if !$do_test;                
        }
        
        if ($do_test)
        {
            push(@print_array, $error_type, $distanceOnScaf);
            return @print_array;
        }
        
    }
    
}


#col    Value (-1 for an array position)
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
    print "\nUsage: \n\tcut -f1,2,3,6 -s file.agp | $0 -i file.psl \n\n";
    print "Based on duplicate blat matches to contigs checks if scaffolding is consistent with cDNA alignment.\n";
    print "BLAT OUTPUT should contain only duplicate (exactly 2x) matches without a header sorted by cDNA name.\n";
    print "Requires scaffold-begin-end-contig columns from file.agp\n\n";
    print "\t-i | --infile : name of the input file.psl\n\t-t | --test : flag to turn testing on/off (default testing is done)\n\n";
    print "Output: \n\tinfile.est_network with EST match links between scaffolded contigs\n"; 
    print "\tagreement condition: exons on same scaffold, same strand and in correct order \n\tif 1 or more conditions for agreement are not met, error tags are assigned accordingly (SCAFFOLD, ORDER or ORIENTATION)\n\tdisagreement cases could be gap filling and scaffold joining cases (do scaffold network to assess)\n\tcolumns: EST ESTsize strand1(ESTmatch) ESTmatch1_start ESTmatch1_end contig_name1 scaff1 scaff1start scaff1end strand2(ESTmatch) ESTmatch2_start ESTmatch2_end contig_name2 scaff2 scaff2start scaff2end distanceOnScaffold tag\n\n";
    exit;
}

exit 0;

