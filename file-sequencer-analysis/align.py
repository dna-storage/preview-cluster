if __name__ == "__main__":
    import argparse
    import os
    import subprocess
    import shutil
    
    parser = argparse.ArgumentParser(description="Parse file paths")

    # stripped file to analyze
    parser.add_argument("--fastq_set", nargs = '+',help="Array of fastq directories to be stripped",required=True)
    parser.add_argument("--align_config",nargs="+",help="Corresponding alignment config paths for each fastq set", required=True) # config file
    parser.add_argument("--skip_indexing", action="store_true", default=False, help="Skip the indexing strip phase") 
    args = parser.parse_args()

    for dir_index, fastq_dir in enumerate(args.fastq_set):
        base_dir=os.path.basename(fastq_dir) #base name
        dump_dir=base_dir+'_Stripped'
        if os.path.exists(dump_dir):
            shutil.rmtree(dump_dir)
        os.mkdir(dump_dir)
        if not os.path.exists(fastq_dir):
            raise ValueError("Directory not found {}".format(fastq_dir))
        if not os.path.isdir(fastq_dir):
            raise ValueError("Path {} is not a directory".format(fastq_dir))
        #get a list of raw fastq file paths
        fastq_list=os.listdir(fastq_dir)
        fastq_list=[_ for _ in fastq_list if ".fastq" in _]
        fastq_path_list = [os.path.join(fastq_dir,_) for _ in fastq_list]
        #at this point have the list of fastq files under the top level fastq_dir, launch bsub command
        for fastq_file_index, fastq_file in enumerate(fastq_list):
            print("Submitting File {} for stripping".format(fastq_file))
            command="bsub -W 480 -C 0 -n 4"
            stdout_path="_".join(["stdout","align",fastq_file])+".txt"
            stdout_opt="-o "+stdout_path
            stderr_path="_".join(["sterr","align",fastq_file])+".txt"
            stderr_opt="-e "+stderr_path
            project_name_opt="-P DNA_ALIGN"
            mem_usage_opt="-R \"rusage[mem=4000]\""
            python_command_opt="python real_analysis/Sequence_Alignment.py"
            input_fastq_opt = "--input "+fastq_path_list[fastq_file_index]
            config_opt="--config "+args.align_config[dir_index]
            output_dir_opt = "--output " + dump_dir
            skip_indexing_opt = ""
            if args.skip_indexing:
                skip_indexing_opt="--skip_alignment"
            command=" ".join([command,stdout_opt,stderr_opt,project_name_opt,mem_usage_opt,python_command_opt,input_fastq_opt,config_opt,output_dir_opt,skip_indexing_opt])
            #launch the command
            print("Launching command {} \n".format(command))
            bsub_process=os.system(command)
    
            
