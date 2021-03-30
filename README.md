
# preview-cluster

# System Requirements

The clustering analysis has been verified to work with python 3. And starcode was compiled with gcc version 4.8.5. If the user has an LSF compatible parallel computing environment they can toggle its use in the commands explained below.

# Clustering Analysis Overview

This repository can be used as a tool to aid in the analysis of sequencing data, specifically sequencing data that represents DNA strands that encode data for an actual file. The main mechanism for the analysis is based on the starcode clustering algorithm (https://github.com/gui11aume/starcode). At a high level the clustering algorithm groups together strands that are within some small edit distance of each other. This is based on the idea that each file's encoded logical strands will have many copies, but the copies may be slightly corrupted due to insertions, deletions, and substitutions. However,it is assumed that the non-corrupted strands will be sequenced at a much higher rate than those that have errors, and the clustering finds centroids that have high amounts of reads compared to strands that have a close edit distance match to it. The close edit distance mathces are then included in the cluster around that centroid. The analysis in this repository leverages the clustering mechanism by seeding the clustering pass with strands we expect to see. For example, with an encoded file, we take the encoded strands without error, replicate them several times, and then include the replicated encoded strands with the fastq sequencing file. This establishes centroids automatically at the strands that we are interested in, and the clustering algorithm will group the sequencing reads to the centroids with the closest edit distance. After placing each sequencing strand in a cluster, a mapping file is created to store which file strand a sequencing strand belongs to. This mapping information can then be used for analysis in a subsequent script to develop read distributions, error rates, etc.

# Clustering Analysis Set up and Tools

To set up the analysis, first starcode needs to be downloaded and compiled. This can be done with the following command:

`make init`

After compiling this, the files that indicate the canonical strands that should form the centroids need to be set up, and configuration files need to be created. There are examples of configuration files and canonical strand files found in the directories `./config_lib` and `./bias_lib`, respectively. Canonical strand files should just be those files that represent the encoded strands for the library, and the configuration files are used to indicate to the clustering analysis which files from the library will be used to set up the centroids, e.g. you can chose a configuration that creates centroids only using strands from 1 file in the library. There are two specific configuration files that need to be created. One is an alignment configuration that is used to pair a *sequencing index DNA sequence* with each *experiment* in your data set. Whether this actually gets used on not depends on how the sequencing machine/service prepares the data. If it *does* remove the *sequencing index*, this configuration file is ignored and the scripts/alignment.sh script needs to be configured to reflect that. If the index does not get removed the alignment scripts will pick out the best matching parts of each sequencing strand that match to the index, and the remaining bases between the index regions are taken and output to a stripped fastq directory. Once an alignment configuration script is made that reflects the experiments being stripped, the following command can be run to strip fastq files. Currently, analysis assumes these indexes are removed.

`./scripts/alignment.sh [-o] <outputdir> [-f] <fastqdir> [-a] <alignmentfile> [-l]`

Where `<outputdir>` is the name of the directory the stripped fastq files will be placed, `<alignmentfile>` provides information on the sequencing index used, `<fastqdir>` is the path to raw fastq files to be processed, and `-l` is used to toggle whether or not lsf should be used. If present, the python command will be submitted to a processor through the lsf interface. 

There should now be a directory populated with stripped fastq files. The next step is to run the clustering process to get a mapping of the sequencing strands. This is where the second config file (known as biasing files) is needed, which indicates which file should be used for seeding. An example of such a config can be found under `./bias_lib/config_preview_2020`. With a config file made, the following command can be run to pair sequencing reads with the expected encoded strands.

`./scripts/cluster.sh [-o] <outputdir> [-s] <strippedfastqdir> [-b] <biasdir> [-c] <biasconfigfile> [-l]`

Where `<outputdir>` is the output directory name where clustering results will be stored, `<strippedfastqdir>` is an input directory name of the path of stripped fastq files to analyze, `<biasdir>` is the path name for the directory that has files for representitive strands that we are searching for and are used to bias the clustering analysis, `<biasconfigfile>` is a configuration file path name for the json file that will be used to specify which files from `<biasdir>` should be included in the biasing process.

Now the results should be in the directory that was set up with the `<outputdir>`, and can be processed by the script at the path `./real_analysis/mapped_strands_analysis_starcode.py` to generate erorr rate analysis and read distribution analyses. The list of options can be generated with:

`python real_analysis/mapped_strands_analysis_starcode.py --help`

The options that set up the analysis with the proper information are explained in more detail as follows:

`--sample_range`: The range of experiments to include in one spread sheet. This is based on the assumption that all sequencing fastq files start with *sX* where *X* is an arbitrary ID number for the fastq file. This can be achieved witha simple renaming before performing all of the analysis. One should check the directory in which mapped strands have been dumped to see which range of sample ID's have been generated since one could generated disjoint ranges 1 at a time, e.g samples 1-5 and 10-15 without 6-9.

`--original_encoding_directory`: path to the directory that represents the original encoded strands. Should be the same as the `<biasdir>` path previously used.

`--stripped_fastq_directory`: path to the stripped fastq data generated by the `alignment.sh` script.

`--primer_partition_top_directory`: path to the top level directory of mapped strands. This should be the same as the `<outputdir>` variable in the `cluster.sh` script.

`--output_book`: name of the excel spread sheet that data will be dumped to. In the book, there will be a separate sheet for each sample identified by `--sample_range`.


