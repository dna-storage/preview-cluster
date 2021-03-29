#!/bin/sh
#this script launches a set of jobs to use starcode clustering analysis on fastq file

#strippedFastQDir=StrippedFastQ_12_16_2020_Template
#outputDir=mapped_strands_12_16_2020_template
#biasDir=/tuck_data/kvolkel/file-sequencer-analysis/bias_lib/base_preview_lib
#biasconfig=/tuck_data/kvolkel/file-sequencer-analysis/config_lib/config_preview_2020/bias_preview_12_16_2020_template.json                                                                           

#strippedFastQDir=StrippedFastQ_12_16_2020_AccessSets
#outputDir=mapped_strands_12_16_2020_access_sets
#biasDir=/tuck_data/kvolkel/file-sequencer-analysis/bias_lib/preview_0HD
#biasconfig=/tuck_data/kvolkel/file-sequencer-analysis/config_lib/config_preview_2020/bias_preview_12_16_2020_access_sets_bg.json

#strippedFastQDir=StrippedFastQ_11_23_2020
#outputDir=mapped_strands_11_23_2020
#biasDir=./bias_lib/preview_error_primer
#biasconfig=./config_lib/config_preview_2020/bias_preview_11_23_2020_error.json


while getopts ":o:s:b:c:" opt; do
    case ${opt} in 
	o )
	    outputDir=$OPTARG
	    ;;
	s )
	    strippedfastQDir=$OPTARG
	    ;;
	b )
	    biasDir=$OPTARG
	    ;;
	c )
	    biasconfig=$OPTARG
	    ;;
	\? )
	    echo "Usage [-o] <outputdir> [-s] <strippedfastqdir> [-b] <biasdir> [-c] <biasconfigfile>"
	    exit 1
	    ;;
	: )
	    echo "Invalid option $OPTARG requires an argument" 1>&2
	    exit 1
	    ;;

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

 bsub -W 480 -N -C 0 -n 4 -o stdout_cluster_$line.txt -e sterr_cluster_$line.txt -P DNA_CLUSTER -R "rusage[mem=8000]" python real_analysis/main.py --input $strippedFastQDir/$line --biasconfig $biasconfig --output $outputDir/$sample/$sample.csv --editdistance 8 --biascount 20 --biasdirectory $biasDir/ --scdir starcode/

done < fastQList.txt

rm fastQList.txt

