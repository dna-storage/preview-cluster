#!/bin/sh
#this script launches a set of jobs to perform alignment in parallel on raw fastq files

#fastQDir=/tuck_data/kvolkel/preview_12_26_2020_template
#outputDir=StrippedFastq_12_16_2020_Template
#alignmentFile=/tuck_data/kvolkel/file-sequencer-analysis/config_lib/config_preview_2020/alignment_preview_12_16_2020_template.config


#fastQDir=/tuck_data/kvolkel/preview_12_26_2020_access_sets
#outputDir=StrippedFastQ_12_16_2020_AccessSets
#alignmentFile=/tuck_data/kvolkel/file-sequencer-analysis/config_lib/config_preview_2020/alignment_preview_12_16_2020_access_sets.config

fastQDir=../FastQFiles_Preview_2020
outputDir=StrippedFastQ_11_23_2020
alignmentFile=./config_lib/config_preview_2020/alignment_preview_11_23_2020.config






ls $fastQDir | cat > fastQList.txt

if [ ! -e $outputDir ]; then
    mkdir $outputDir
fi

while read line; do
echo $line

 bsub -W 480  -N -C 0 -n 4 -o stdout_align_$line.txt -e sterr_align_$line.txt -P DNA_ALIGN -R "rusage[mem=4000]" python real_analysis/Sequence_Alignment.py --input $fastQDir/$line --config $alignmentFile --output $outputDir --skip_alignment
done < fastQList.txt
rm fastQList.txt
