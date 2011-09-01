#!/usr/bin/env perl -w
#requires AGP file, network file (nodes, edges, edge attributes), length attribute file for each sequence
#network and attribute file as produced by makeCytoscapeNetwork.sh
#OR 101020
#OR 101215 corrected and tested for partial merge cases
#OR 110304 added help message

# It assumes that scaffolds are as W:N:W blocks, W stands for WGS contig, 
#    could be other sequence element, except gap as defined in AGP format:
#        smallest_possible_scaffold=W [1 contig element, no gaps], 
#        bigger always having W atr extremities and Ns in between each 2 W

# In network file, the first scaffold listed has the start of EST and the second listed has the end of EST

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max maxstr sum min);

#default valuies
my $intron_limit=20000;
my $maxSizeDifference=4000;

my $setOptions=GetOptions(
    "network|n=s" => \my $nw_file,
    "length|l=s"  => \my $length_file,
    "agp|a=s"     => \my $agp_file,
    "intron|i=i"  => \$intron_limit,
    "diff|d=i"      => \$maxSizeDifference,
    "help|h"      => sub { HelpMessage() },
    );

my (@scaf1, @scaf2, @coord_scaf1, @coord_scaf2);

#exit if in files are not prvided
if (!$agp_file || !$nw_file || !$length_file)
{
	HelpMessage();
}

#print header
print  "source\ttarget\tinteraction\tmatchOnSource\tmatchOnTarget\tmatchStartOnSource\tmatchEndOnSource\tmatchStartOnTarget\tmatchEndOnTarget\tstatus\ttype\testLength\tDownStream\tUpStream\n";

#length file has no header as is temp file
my(@len)=split(/\t/, `awk '{ORS="\t"; print}' $length_file`);
my(%SQ_len)=(@len);

my(@NW) = split(/;/, `awk 'FNR>1{ORS=";"; print}' $nw_file`);
	
foreach my $line (@NW)
   {
      my @values=split(/\t/,$line);
      if ($values[10] eq "join")
      {
         print "$line\tNA\tNA\n";
      }
      else
      {
         @scaf1=split(/\t/,`awk -v myScaf=$values[0] '\$1==myScaf{ORS="\t";if(\$5=="N"){print \$6} else {print \$8}}' $agp_file`);
         @scaf2=split(/\t/,`awk -v myScaf=$values[1] '\$1==myScaf{ORS="\t";if(\$5=="N"){print \$6} else {print \$8}}' $agp_file`);
			
         #if sequence is not in AGP file, it is considered as unscaffolded contig
         $scaf1[0]=$SQ_len{$values[0]}, if (scalar(@scaf1)==0);
         $scaf2[0]=$SQ_len{$values[1]}, if (scalar(@scaf2)==0);

         #blat PSL poutput uses 0-based start, which means that start coordinates are offset by -1
         #see http://genomewiki.ucsc.edu/index.php/Coordinate_Transforms
         #therefore adjust the starts by +1
         @coord_scaf1=@values[5,6];
         $coord_scaf1[0]+=1;
         @coord_scaf2=@values[7,8];
         $coord_scaf2[0]+=1;

         #if match on scaffold is on - strand, then reverse the fragments and update the coordinates (new_coord=total-coord+1)
         if ($values[3]=~/-/)
         {
            @scaf1=reverse(@scaf1);
            @coord_scaf1=reverse(@coord_scaf1);
            $coord_scaf1[0]=$SQ_len{$values[0]}-$coord_scaf1[0]+1;
            $coord_scaf1[1]=$SQ_len{$values[0]}-$coord_scaf1[1]+1;
         }			

         if ($values[4]=~/-/)
         {
            @scaf2=reverse(@scaf2);
            @coord_scaf2=reverse(@coord_scaf2);
            $coord_scaf2[0]=$SQ_len{$values[1]}-$coord_scaf2[0]+1;
            $coord_scaf2[1]=$SQ_len{$values[1]}-$coord_scaf2[1]+1;
         }


         my $lastMatchElement_scaf1=findArrayElement(\@scaf1, max(@coord_scaf1));
         my $firstMatchElement_scaf2=findArrayElement(\@scaf2, min(@coord_scaf2));
#print "scaf1 el no $lastMatchElement_scaf1 scaf2 el no $firstMatchElement_scaf2 of length $scaf2[$firstMatchElement_scaf2]\n";

			#if scaffolds are of one element then depending if they contain end or start of EST 
			#then the up is compirsed in down stream calculation or vice versa
			#same if scaffold 1 match element is the last one or scaffold 2 match element is the first one
			my ($statusDown, $statusUp);
			if (scalar(@scaf2)==1 || $firstMatchElement_scaf2==0)
			{
				$statusDown=1;
			}
			else
			{ 
				$statusDown=compareSizesDownStream(\@scaf1, \@scaf2, $lastMatchElement_scaf1, ($firstMatchElement_scaf2-1), $maxSizeDifference);
                        }

			if (scalar(@scaf1)==1 || (scalar(@scaf1)-$lastMatchElement_scaf1)==1)
			{
				$statusUp=1;
			}
			else
			{ 
				$statusUp=compareSizesUpStream(\@scaf1, \@scaf2, ($lastMatchElement_scaf1+1), $firstMatchElement_scaf2, $maxSizeDifference);
			}
			
			#print $values[0]."\t".$values[1]."\t$statusDown\t$statusUp\n";
			if ($statusDown && $statusUp)
			{$values[9]="OK";}
			else
                        {$values[9]="Pb";}

			if ($statusDown)
			{
				$values[12]="fits";
			}
			else
			{
				$values[12]="conflict";
			}

                        if ($statusUp)
         {
            $values[13]="fits";
         }
         else
         {
            $values[13]="conflict";
         }

         my $newline=join("\t", @values);
         print "$newline\n";

      }
   }

sub findArrayElement
{
	my ($arrayRef, $position)=@_;
	my @array=@{$arrayRef};
	
	my $i;
	my $sum=0;

#print "position is $position\n";
        for ($i=0; $i<scalar(@array); $i++)
        {
		$sum+=$array[$i];
#print "$i element $array[$i] sum  $sum\n";
		last if ($sum>=$position);
        }
	return ($i);
}

sub compareSizesDownStream
{
	my ($arrayRef1, $arrayRef2, $nb1, $nb2, $cutoff)=@_;
        my @array1=@{$arrayRef1};
        my @array2=@{$arrayRef2};

        my $i;
	my($tag)=1;
        for ($i=0; $i<=min($nb1, $nb2); $i++)
        {  
		my $elementNo1=$nb1-$i;
		my $elementNo2=$nb2-$i;
		my $case1=$elementNo1%2;
		my $case2=$elementNo2%2;
		my $element1=$array1[$elementNo1];
		my $element2=$array2[$elementNo2];
		#we expect to get scaffolds as W-N-W-N-W
		#threrefore all divisible by 2 are gaps
		#others are contiguous sequences (W or other as in AGP format)
		if (!$case1)
		{
#print "scaf1 is contig $element1 & scaf2 is gap $element2\n";
			$tag=0, if (($element1-$element2)>$cutoff);
		}
		if (!$case2)
		{
#print "scaf1 is gap $element1 & scaf2 is contig $element2\n";
			$tag=0, if (($element2-$element1)>$cutoff);
		}		
 
                last if (!$tag);
        }

	return ($tag);
}


sub compareSizesUpStream
{
        my ($arrayRef1, $arrayRef2, $nb1, $nb2, $cutoff)=@_;
        my @array1=@{$arrayRef1};
        my @array2=@{$arrayRef2};

        my $i;
        my $delta;
        my($tag)=1;
        for ($i=0; $i<min((scalar(@array1)-$nb1), (scalar(@array2)-$nb2)); $i++)
        {
                my $elementNo1=$nb1+$i;
                my $elementNo2=$nb2+$i;
                my $case1=$elementNo1%2;
                my $case2=$elementNo2%2;
                my $element1=$array1[$elementNo1];
                my $element2=$array2[$elementNo2];
                #we expect to get scaffolds as W-N-W-N-W
                #threrefore all divisible by 2 are N's
                if (!$case1) 
                {
#print "scaf1 is contig $element1 & scaf2 is gap $element2\n";
                	$tag=0, if (($element1-$element2)>$cutoff);
		}
                if (!$case2)
                {
#print "scaf1 is gap $element1 & scaf2 is contig $element2\n";
                        $tag=0, if (($element2-$element1)>$cutoff);
                }

                last if (!$tag);
        }

        return ($tag);
}

sub HelpMessage
{
    print "\nusage $0 [options]\n\n";
    print "--network|-n file with (pre-)network generated by makeCytoscapeNetwork.sh\n";
    print "--length|-l file with contigs and scaffold lengths (sequenceID <tab> length)\n";
    print "--agp|-a scaffolding file in agp format\n";
    print "--intron|-i maximum intron size (bp)\n";
    print "--diff|-d maximum difference between a gap and a contig to fill in the gap (bp)\n";
    print "--help|-h help message\n\n";
    exit;
}

exit 0;

