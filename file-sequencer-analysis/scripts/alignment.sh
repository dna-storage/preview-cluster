#!/bin/sh
#this script launches a set of jobs to perform alignment in parallel on raw fastq files
fastQDir=""
outputDir=""
alignmentFile=""
lsf=false

while getopts ":o:f:a:l" opt; do
    case ${opt} in 
	o )
	    outputDir=$OPTARG
	    ;;
	f )
	    fastQDir=$OPTARG
	    ;;
	a )
	    alignmentFile=$OPTARG
	    ;;
	l )
	    lsf=true
	    ;;
	\? )
	    echo "Usage [-o] <outputdir> [-f] <fastqdir> [-a] <alignmentfile>" 1>&2
	    exit 1
	    ;;
	: )
	    echo "Invalid option $OPTARG requires an argument" 1>&2
	    exit 1
	    ;;
    esac
done


ls $fastQDir | cat > fastQList.txt

if [ ! -e $outputDir ]; then
    mkdir $outputDir
fi

COUNTER=1
while read line; do
    if [ ${line: -6} == ".fastq" ] ; then
       PREFIX="s${COUNTER}_"
       echo $line
       echo $PREFIX
       if [ "$lsf" = true ] ; then
	   bsub -W 480  -N -C 0 -n 4 -o stdout_align_$line.txt -e sterr_align_$line.txt -P DNA_ALIGN -R "rusage[mem=4000]" python real_analysis/Sequence_Alignment.py --input $fastQDir/$line --config $alignmentFile --output $outputDir --skip_alignment --prefix $PREFIX
       fi
       if [ "$lsf" = false ] ; then
	   python real_analysis/Sequence_Alignment.py --input $fastQDir/$line --config $alignmentFile --output $outputDir --skip_alignment --prefix $PREFIX
       fi
       COUNTER=$((COUNTER+1))
    fi
done < fastQList.txt
rm fastQList.txt
