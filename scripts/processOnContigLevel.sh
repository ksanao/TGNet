#!/usr/bin/env bash
#OR 110202

outPrefix=$1
myAGP=$2
scriptDIR=$3

myPSL="$outPrefix"_filtered.psl	

#build network of contigs linked by transcript alignments
myTag=contig.nw
echo -e "Building transcript alignment network"
cut -f1,2,3,6 -s $myAGP | perl $scriptDIR/blat_buildContigNetworkForCytoscape.pl -i $myPSL

#make simplified format for cytoscape (matches pre-defined visual properties VizMap file)
printf "contig_name1\tinteraction\tcontig_name2\t" > "$outPrefix"."$myTag"
printf "interactionType\tstrand1strand2\tESTsize\t" >> "$outPrefix"."$myTag"
printf "distanceOnScaffold\ttag\n" >> "$outPrefix"."$myTag"
awk 'FNR>1{if ($18=="ORDER" || $18=="ORIENTATION"){myType="alignment_Pb"}else{myType="alignment"}; print $6"\t"$1"\t"$13"\t"myType"\t"$3$10"\t"$2"\t"$17"\t"$18}' "$outPrefix"_filtered.est_network >> "$outPrefix"."$myTag"

#build scaffolding network for contigs present in transcript alignments network
echo Building scaffolding network
awk 'FNR>1{print $6"\t"$7"\n"$13"\t"$14}' "$outPrefix"_filtered.est_network | sort -k1,1 | uniq > temp.list
#put all scaffolds in a nonredundant array
myS=(`awk 'FNR>1{print $7"\n"$14}' "$outPrefix"_filtered.est_network | sort | uniq`)
#for each scaffold define one to one connections for its contigs with transcript alignments
for oneS in ${myS[@]}
do
    #printf "    ..processing "$oneS"\n"
    grep $oneS temp.list | sort | uniq > temp
    myL=`cat temp | wc -l`
    if [ "$myL" -gt 1 ]
    then 
        awk -v myVar="$myL" 'FNR<myVar{print $1}' temp > temp1
        awk 'FNR>1{print $1}' temp > temp2
        awk 'FNR>1{print $2}' temp > temp3
        paste temp1 temp2 temp3 | awk '{print $1"\t"$3"\t"$2"\tscaffolding\t\t\t\t"}' >> "$outPrefix"."$myTag"
    fi
done

#generate node attributes: length of contigs, scaffold assignment
myTag=contig.attr
echo Generating node attributes
#I want to discard split EST matches to same contig and not output those to final node attributes
printf "contig\tlength\tscaffold\n" > "$outPrefix"."$myTag"
awk '$6!=$13' "$outPrefix"_filtered.est_network | awk 'FNR>1{myL1=$9-$8+1; myL2=$16-$15+1; print $6"\t"myL1"\t"$7"\n"$13"\t"myL2"\t"$14}' | sort | uniq >> "$outPrefix"."$myTag"

