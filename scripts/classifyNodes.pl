#!/usr/bin/env perl -w
#Defines whether a node is a scaffold or contig using pre-processed agp file
#if a node is not in agp file then it is considered a contg

#OR 110812

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max maxstr sum min);

my $setOptions=GetOptions(
    "node|n=s" => \my $node_attr,
    "agp|a=s"     => \my $agp_file,
    "help|h"      => sub { HelpMessage() },
    );


my (%node_in_AGP);

open (AGP, "<", $agp_file)  || die "fail open $agp_file $!\n";
while (<AGP>)
{
    my $line=$_;
    my @values=split(/\t/, $line);
    $node_in_AGP{$values[0]}=1, if !exists($node_in_AGP{$values[0]});
}
close AGP;

open (NOD, "<", $node_attr)  || die "fail open $node_attr $!\n";

while (<NOD>)
{
    my $line=$_;
    chomp $line;
    my @values=split(/\t/, $line);
    if (exists($node_in_AGP{$values[0]}))
    {
        print "$line\tscaffold\n";
    }
    else
    {
        print "$line\tcontig\n";
    }
}

close NOD;

sub HelpMessage
{
    print "\n$0 -a file.agp -n file.node.attributes\nOutputs node attributes along with node state in the assembly (contig or scaffold)\n\n";
    exit 0;
}

exit 0;
