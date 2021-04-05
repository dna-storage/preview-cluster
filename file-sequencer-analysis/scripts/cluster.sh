#!/bin/bash
#this script launches a set of jobs to use starcode clustering analysis on fastq file
outputDir=""
strippedFastQDir=""
biasDir=""
biasconfig=""
lsf=false

SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin                                                                                                                                                                        
SCRIPTPATH=$(dirname "$SCRIPT")

while getopts ":o:s:b:c:l" opt; do
    case ${opt} in 
	o )
	    outputDir=$OPTARG
	    ;;
	s )
	    strippedFastQDir=$OPTARG
	    ;;
	b )
	    biasDir=$OPTARG
	    ;;
	c )
	    biasconfig=$OPTARG
	    ;;
	l )
	    lsf=true
	    ;;
	\? )
	    echo "Usage [-o] <outputdir> [-s] <strippedfastqdir> [-b] <biasdir> [-c] <biasconfigfile>" 1>&2
	    exit 1
	    ;;
	: )
	    echo "Invalid option $OPTARG requires an argument" 1>&2
	    exit 1
	    ;;
    esac
done


ls $strippedFastQDir | cat > fastQList.txt

if [ ! -e $outputDir ]; then
    mkdir $outputDir
fi

while read line; do
    echo $line

    sample="$(cut -d'_' -f1 <<<$line)"
    
    echo "Sample ID: " $sample
    
    if [ ! -e $outputDir/$sample ]; then
	mkdir $outputDir/$sample
    fi
    
    
    if [ "$lsf" = true ] ; then
	bsub -W 480 -N -C 0 -n 4 -o stdout_cluster_$line.txt -e sterr_cluster_$line.txt -P DNA_CLUSTER -R "rusage[mem=8000]" python3 $SCRIPTPATH/../real_analysis/main.py --input $strippedFastQDir/$line --biasconfig $biasconfig --output $outputDir/$sample/$sample.csv --editdistance 8 --biascount 20 --biasdirectory $biasDir/ --scdir $SCRIPTPATH/../starcode/
    fi
    
    if [ "$lsf" = false ] ; then
	python3 $SCRIPTPATH/../real_analysis/main.py --input $strippedFastQDir/$line --biasconfig $biasconfig --output $outputDir/$sample/$sample.csv --editdistance 8 --biascount 20 --biasdirectory $biasDir/ --scdir $SCRIPTPATH/../starcode/
	
    fi
done < fastQList.txt

rm fastQList.txt

