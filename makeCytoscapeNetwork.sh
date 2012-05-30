#!/usr/bin/env bash
#OR 110202

#IN: assembly fasta, transcript fatsa, AGP scaffolding file
#OUT: network files with interactions and properties, (and optionally files generated in process like blat output)
#DO: build network definition

ERROR_ARGUMENTS=65
ERROR_NOFILE=66

USAGE="\nUsage: \n\
`basename $0` -g genome.fasta -t transcript.fasta -a file.agp [-chklinoqsw]\n\
`basename $0` -p blat.alignment.psl -a file.agp [-chklinoqsw]\n\n\
-a file specifying scaffolding way in agp format\n\
-c flag to build network on contig level (required if contig.fasta is submitted)\n\
-g fasta file with scaffolded genome contigs (requires -c flag) or\n\
\tfasta file with genomic scaffolds (can include unscaffolded contigs)\n\
-h print help\n\
-k flag to keep blat alignment files\n\
\tdefault: do not keep\n\
-l match length (bp, minimum transcript match length for blat output filterting)\n\
\tdefault=200 \n\
-n network.name (string to be used for out files)\n\
\tdefault=out\n\
-o overlap (bp, maximum overlap of split transcript alignments)\n\
\tdefault=50 \n\
-p transcript to genome alignment file in psl format (can be submitted instead fasta files)\n\
-q flag to keep folder with temporary files\n\
-i intron size threshold (bp)\n\
\tdefault=20000 \n\
-s maximum allowed size difference between a gap and a sequence to fill in the gap (bp, required only for scaffold level) \n\
\tdefault=4000 \n\
-t fasta file with transcript sequences\n\
-w mismatch (floating, percentage of allowed mismatch bases)\n\
\tdefault=0.05\n"


# Print usage: scriptname -options, if script invoked with no command-line args
if [ $# -eq 0 ]  
then
	echo -e $USAGE
	exit $ERROR_ARGUMENTS
fi  

SCNversion=1.0.0

# Options followed by : expect an argument 
while getopts "a:g:i:l:n:o:p:s:t:w:chkqv" Option
do
  case $Option in
    a ) inAGP=$OPTARG;;
    c ) cflag=1;;
    g ) inAssembly=$OPTARG;;
    h ) echo -e $USAGE; exit 1;;
    i ) myMaxIntron=$OPTARG;;
    k ) kflag=1;;
    l ) matchLength=$OPTARG;;
    n ) outPrefix=$OPTARG;;
    o ) maxOverlap=$OPTARG;;
    p ) inUserPSL=$OPTARG;;
    q ) qflag=1;;
    s ) myMaxDifference=$OPTARG;;
    t ) inTranscripts=$OPTARG;;
    w ) misMatch=$OPTARG;;
    v ) echo version: $SCNversion; exit;;
    #if no argument is supplied to the option requiring it, it will fall into default
    * ) echo "Unimplemented option chosen.";  exit $ERROR_ARGUMENTS;;   # DEFAULT
  esac
done

# Decrements the argument pointer so it points to next argument.
shift $(($OPTIND - 1))

#store start time
startTime=`echo $(date) |awk '{split($4,myS,":"); print myS[1]*360+myS[2]*60+myS[3]}'`

#### CHECKING SCRIPTS in RELATIVE SCRIPT DIRECTORY
scriptDIR="$(cd ${0%/*}; echo "$PWD")"
scriptDIR=$scriptDIR/scripts
requiredScripts=(AGP_keepOnlyScaffolds.pl
                 blat_filterResults.pl
                 blat_filterDuplicatesWithOverlap.pl
                 blat_buildContigNetworkForCytoscape.pl
                 generateScaffoldingAttributes.pl)

for myscript in ${requiredScripts[@]}
do
        if [[ ! -r $scriptDIR/$myscript ]]
        then
                echo Cannot open $scriptDIR/$myscript
                exit $ERROR_NOFILE
        fi
done

#### DEFAULT VALUES
maxOverlap_DEFAULT=50
outPrefix_DEFAULT="out"
matchLength_DEFAULT=200
misMatch_DEFAULT=0.05
myMaxIntron_DEFAULT=20000
myMaxDifference_DEFAULT=4000

#### SETTING UNDEFINED OPTIONS TO DEFAULTS
: ${maxOverlap=$maxOverlap_DEFAULT}
: ${outPrefix=$outPrefix_DEFAULT}
: ${matchLength=$matchLength_DEFAULT}
: ${misMatch=$misMatch_DEFAULT}
: ${myMaxIntron=$myMaxIntron_DEFAULT}
: ${myMaxDifference=$myMaxDifference_DEFAULT}

#### SETUP A TEMPORARY WORKING DIRECTORY
currentDir=$PWD
tempDIR=temp"$(date +%Y%m%d_%H%M%S)"
mkdir $tempDIR
if [[ -z $inUserPSL ]]
then
    #use symlinks in temp dir
    inFILES=($inAssembly $inTranscripts)
    myAssembly=${inAssembly##*/}
    myTranscripts=${inTranscripts##*/}
else
    inFILES=($inUserPSL)
    myUserPSL=${inUserPSL##*/}
fi

#make symlinks for input files (with absoulute paths)
for onefile in ${inFILES[@]}
do
    ln -s "$(cd `dirname $onefile`; echo "$PWD")"/"${onefile##*/}" $tempDIR/
done

#if AGP is input, process it to retain only scaffolds (more than 1 element features)
if [[ ! -z $inAGP ]]
then
    echo Retrieving scaffold features from AGP
    perl $scriptDIR/AGP_keepOnlyScaffolds.pl $inAGP > $tempDIR/scaffold_only.agp
    myAGP=scaffold_only.agp
fi

#go to RUN directory
cd $tempDIR

##### STEP1: RUN BLAT
if [[ -z $inUserPSL ]]
then
    if [[ ! -x `which blat` ]]; then
        echo Cannot find blat or is not an executable
        exit $ERROR_NOFILE 
    else
        myPSL=$outPrefix.psl
        BLATcall="blat $myAssembly $myTranscripts $myPSL"
        echo Aligning: $BLATcall
        echo `$BLATcall`
    fi
else
    myPSL=$myUserPSL
fi

##### STEP2: FILTER BLAT OUTPUT (PSL FORMAT)
echo Filtering psl
mytag=l"$matchLength"f0m"$misMatch"g"$myMaxIntron"e0x30b1t0
awk 'FNR>5' "$myPSL" | perl $scriptDIR/blat_filterResults.pl -l $matchLength -e 0 -m $misMatch -g "$myMaxIntron" > "$outPrefix"_"$mytag".psl
sort -k10,10 -k12,12n "$outPrefix"_"$mytag".psl | perl $scriptDIR/blat_filterDuplicatesWithOverlap.pl -o $maxOverlap > "$outPrefix"_filtered.psl

##### STEP 3 generate alignment statistics
if [[ -z "$cflag" ]]
then
    mytag=scaffold
else
    mytag=contig
fi
myOutFile=$outPrefix.$mytag.blat.stat
printf "RunName\t"$outPrefix"\n" > $myOutFile
printf "unfilteredMatches\t" >> $myOutFile
awk 'FNR>5{print $10}' "$myPSL" | sort | uniq | wc -l >> $myOutFile
printf "Matches_max_mismatch_$misMatch\t" >> $myOutFile
mytag=l"$matchLength"f0m"$misMatch"g"$myMaxIntron"e0x30b1t0
cut -f10 "$outPrefix"_"$mytag".psl | sort | uniq | wc -l >> $myOutFile
myFraction=(0.50 0.90 0.95 0.99)
for oneFraction in ${myFraction[@]}
do
   printf $oneFraction"LengthMatches\t" >> $myOutFile
   mytag=l"$matchLength"f"$oneFraction"m"$misMatch"g"$myMaxIntron"e0x30b1t0
   awk 'FNR>5' "$myPSL" | perl $scriptDIR/blat_filterResults.pl -l $matchLength -f $oneFraction -e 0 -m $misMatch -g "$myMaxIntron" > "$outPrefix"_"$mytag".psl
   cut -f10 "$outPrefix"_"$mytag".psl | sort | uniq | wc -l >> $myOutFile
done

###### STEP4 generate cytoscape network
#Swith between contig and scaffold level
unset -v mytag
if [[ -z "$cflag" ]]
then
	##### SCAFFOLD-LEVEL
        # can process (reduced level) without agp
	echo Processing on scaffold-level
	bash $scriptDIR/processOnScaffoldLevel.sh $outPrefix $scriptDIR $myMaxIntron $myMaxDifference $myAGP
	myOutFiles=("$outPrefix".scaffold.nw "$outPrefix".scaffold.attr "$outPrefix".withinOneScaffold $outPrefix.scaffold.blat.stat)
        mytag=scaffold
else	
	#### CONTIG LEVEL
        # agp file is required
	echo Processing on contig-level
        bash $scriptDIR/processOnContigLevel.sh $outPrefix $myAGP $scriptDIR
	myOutFiles=("$outPrefix".contig.nw "$outPrefix".contig.attr $outPrefix.contig.blat.stat)
        mytag=contig
fi

mv ${myOutFiles[@]} ..

#print options and in/out files used/generated
echo -e "\nRUN SUMMARY\n\
Input_files: ${inFILES[@]}\n\
Output_directory: $currentDir"
printf "Output_files: "
printf "%s " "${myOutFiles[@]}"

#remove alignment files, unless -k flag was set
if [[ ! -z "$kflag" ]]
then
    outname="$outPrefix"_filtered."$mytag".psl
    mv "$outPrefix"_filtered.psl ../"$outname"
    printf "%s " "$outname"

    if [[ -z "$inUserPSL" ]]
    then 
        outname="$outPrefix"."$mytag".psl
        mv "$outPrefix".psl ../"$outname"
        printf "%s " "$outname"
    fi
fi

cd ..

if [[ -z "$qflag" ]]
then
    rm -r $tempDIR
else
    outname="$outPrefix"_"$mytag"_"$tempDIR"
    mv $tempDIR "$outname"
    printf "\nTemporary_directory: $outname"
fi

#continue printing
echo -e "\n\
Options: maximum overlap between aligned transcript parts = $maxOverlap bp\n\
         minimum match length = $matchLength bp\n\
         maximum mismatch fraction of alignment = $misMatch\n\
         intron threshold = $myMaxIntron bp"
if [[ -z "$cflag" ]]
then
echo -e "         maximum difference between gap and scaffold to fit in = $myMaxDifference"
fi

endTime=`echo $(date) |awk '{split($4,myS,":"); print myS[1]*360+myS[2]*60+myS[3]}'`
printf $endTime"\t"$startTime | awk '{myTime=($1-$2)/60; print "Processing time: "myTime" minutes."}'

exit 0
