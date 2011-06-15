#!/usr/bin/env bash
#OR 110202

outPrefix=$1
scriptDIR=$2
myMaxIntron=$3
myMaxDifference=$4
myAGP=$5

myPSL="$outPrefix"_filtered.psl

#generate node attributes: length of contigs
echo generating node attributes
printf "node\tlength\n" > "$outPrefix".scaffold.attr
awk '{print $14"\t"$15}' $myPSL | sort | uniq >> "$outPrefix".scaffold.attr

#generate EST alignment network
echo building transcript alignment network
myEST=(`cut -f10 $myPSL | sort | uniq`)

printf "" > "$outPrefix"_GnLocation
printf "" > "$outPrefix"_EstLocation

# we will compute max intron distance from each scaffold end, thus actual distance can be $myMaxIntron*2

for oneE in ${myEST[@]}
do
    awk -v myEst="$oneE" '$10==myEst' $myPSL |
    awk -v myDist=$myMaxIntron '{Fdist=$16-1; Tdist=$15-$17; \
            if ( ( NR % 2 ) == 0 ) {ORS="\n"} else {ORS="\t"}
            if (Fdist<=myDist && Tdist<=myDist) {print $14"\tboth"$9"\t"$16"\t"$17}
       else if (Fdist<=myDist && Tdist>myDist)  {print $14"\tfive"$9"\t"$16"\t"$17} 
       else if (Fdist>myDist && Tdist<=myDist)  {print $14"\tthree"$9"\t"$16"\t"$17}
                                          else  {print $14"\tmid"$9"\t"$16"\t"$17}}' >> "$outPrefix"_GnLocation

# printf "\n" >> "$outPrefix"_GnLocation

    awk -v myEst="$oneE" '$10==myEst' $myPSL | 
    awk '{if ( ( NR % 2 ) == 0 ) {ORS="\n"} else {ORS="\t"}
         print $14"\t"$10"\t"$11"\t"$12}' | 
    awk '{if ($4>$8) {print $1"\t"$5"\t"$2"\te\ts\t"$3}
                else {print $1"\t"$5"\t"$2"\ts\te\t"$3}}'>> "$outPrefix"_EstLocation
done

#source has starting EST exons
#target has ending EST exons
printf "source\ttarget\tinteraction\t\
matchOnSource\tmatchOnTarget\t\
matchStartOnSource\tmatchEndOnSource\t\
matchStartOnTarget\tmatchEndOnTarget\t\
status\ttype\testLength\n" > "$outPrefix".network.temp

paste "$outPrefix"_GnLocation "$outPrefix"_EstLocation | 
awk '{if($12=="s"){print $1"\t"$5"\t"$11"\t"$2"\t"$6"\t"$3"\t"$4"\t"$7"\t"$8"\t"$14}
              else{print $5"\t"$1"\t"$11"\t"$6"\t"$2"\t"$7"\t"$8"\t"$3"\t"$4"\t"$14}}' | 
awk '$1!=$2' | 
awk '{if(($4=="three+" || $4=="five-" || $4~/both/) && ($5=="three-" || $5=="five+" || $5~/both/)){mySTATUS="OK"; myCASE="join"}
      else{mySTATUS="Pb"; myCASE="merge"}
      print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"mySTATUS"\t"myCASE"\t"$10}' >> "$outPrefix".network.temp
# do not print if empty
paste "$outPrefix"_GnLocation "$outPrefix"_EstLocation | 
awk '{if($12=="s"){print $1"\t"$5"\t"$11"\t"$2"\t"$6"\t"$3"\t"$4"\t"$7"\t"$8"\t"$14}
    else{print $5"\t"$1"\t"$11"\t"$6"\t"$2"\t"$7"\t"$8"\t"$3"\t"$4"\t"$14}}' | 
awk '$1==$2' > "$outPrefix".withinOneScaffold

myNW="$outPrefix".network.temp

if [[ -e $myAGP ]]
then
   perl $scriptDIR/generateScaffoldingAttributes.pl -a "$myAGP" -n $myNW -l "$outPrefix".scaffold.attr -i $myMaxIntron -d $myMaxDifference > "$outPrefix".network.processed.temp
   myNW="$outPrefix".network.processed.temp
fi

#make network simpler and match mizual mapping properties in VizMap.props file
printf "contig_name1\tinteraction\tcontig_name2\t\
interactionType\tstrand1strand2\tESTsize\t\
type\ttag\n" > "$outPrefix".scaffold.nw

awk 'FNR>1{if($10=="Pb"){myType="alignment_Pb"}
   else{myType="alignment"};
   if($4~/\-/){strand1="-"};
   if($4~/\+/){strand1="+"};
   if($5~/\-/){strand2="-"};
   if($5~/\+/){strand2="+"};
   print $1"\t"$3"\t"$2"\t"myType"\t"strand1 strand2"\t"$12"\t"$11"\t"$10}' $myNW >> "$outPrefix".scaffold.nw

