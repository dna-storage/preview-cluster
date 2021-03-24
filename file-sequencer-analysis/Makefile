.PHONY: clean help init align cluster preview_1_set_12_16 preview_2_set_12_16 preview_3_set_12_16 preview_4_12_16

SCRIPT_DIR=scripts



cluster:
	$(SCRIPT_DIR)/cluster.sh

align:
	$(SCRIPT_DIR)/alignment.sh


init: #run init steps
	git clone https://github.com/gui11aume/starcode
	make -C starcode

test:
        #! /usr/bin/python2.7.15
	nosetests tests

clean:
	rm -rf build
	rm -rf dist
	rm -rf dnastorage.egg-info
	rm -rf generate.egg-info


#commands used to process results for mapped strands for preview 12_16_2020 sequencing data
preview_1_set_12_16:
	bsub -W 480 -N -C 0 -n 4 -o stdout_analysis_1.txt -e sterr_analysis_1.txt -P Analysis_1 -R "rusage[mem=8000]" python real_analysis/mapped_strands_analysis_starcode.py --original_encoding_directory bias_lib/preview_0HD/ --stripped_fastq_directory preview_12_16_2020_access_sets_1_Stripped/ --primer_partition_top_directory mapped_strands_12_16_2020_access_sets/ --sample_distributions --output_book preview_sets_distributions_1.xlsx --sample_range 1-5

preview_2_set_12_16:
	bsub -W 480 -N -C 0 -n 4 -o stdout_analysis_2.txt -e sterr_analysis_2.txt -P Analysis_1 -R "rusage[mem=8000]" python real_analysis/mapped_strands_analysis_starcode.py --original_encoding_directory bias_lib/preview_0HD/ --stripped_fastq_directory preview_12_16_2020_access_sets_2_Stripped/ --primer_partition_top_directory mapped_strands_12_16_2020_access_sets/ --sample_distributions --output_book preview_sets_distributions_2.xlsx --sample_range 6-10

preview_3_set_12_16:
	bsub -W 480 -N -C 0 -n 4 -o stdout_analysis_3.txt -e sterr_analysis_3.txt -P Analysis_1 -R "rusage[mem=8000]" python real_analysis/mapped_strands_analysis_starcode.py --original_encoding_directory bias_lib/preview_0HD/ --stripped_fastq_directory preview_12_16_2020_access_sets_3_Stripped/ --primer_partition_top_directory mapped_strands_12_16_2020_access_sets/ --sample_distributions --output_book preview_sets_distributions_3.xlsx --sample_range 11-15

preview_4_set_12_16:
	bsub -W 480 -N -C 0 -n 4 -o stdout_analysis_4.txt -e sterr_analysis_4.txt -P Analysis_1 -R "rusage[mem=8000]" python real_analysis/mapped_strands_analysis_starcode.py --original_encoding_directory bias_lib/preview_0HD/ --stripped_fastq_directory preview_12_16_2020_access_sets_4_Stripped/ --primer_partition_top_directory mapped_strands_12_16_2020_access_sets/ --sample_distributions --output_book preview_sets_distributions_4.xlsx --sample_range 16-20

install:
	python setup.py install --user

