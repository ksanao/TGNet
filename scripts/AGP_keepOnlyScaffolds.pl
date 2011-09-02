#!/usr/bin/env perl -w
#OR 110812
#IN: AGP
#OUT: AGP with features containing more than one element (SCAFFOLDS)
 
use strict;
use warnings;
use Switch;
use Getopt::Long;

# optionName|alternativeName
my $opt=GetOptions(
        "help"        => sub { HelpMessage() },
        );

my (%counts, %agp);

while (<>)
{
    my $line=$_;
    #scaffold name is the first field in tab delimited line
    my @values=split(/\t/, $line);
    $counts{$values[0]}+=1;
    $agp{$values[0]}.=$line;
}

foreach my $feature ( sort keys %agp )
{
    if ($counts{$feature}>1)
    {
        print $agp{$feature};
    }
}

sub HelpMessage
{
    print "\n$0 file.agp\nOutputs AGP features that have more then 1 element (scaffolds)\n\n";
    exit 0;
}

exit 0;

