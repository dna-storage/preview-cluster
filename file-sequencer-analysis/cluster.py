if __name__ == "__main__":
    import argparse
    import os
    import subprocess
    import shutil
    
    parser = argparse.ArgumentParser(description="Parse file paths")

    # stripped file to analyze
    parser.add_argument("--stripped_fastq_set", nargs = '+',help="Array of fastq (stripped) directories to be stripped",required=True)
    parser.add_argument("--bias_config",nargs="+",help="Corresponding bias config paths for each fastq set", required=True) # config file
    parser.add_argument("--bias_dir",type=str,help="Corresponding bias directory that houses the strands we are expecting",required=True)
    parser.add_argument("--out_dir",help="Top level output directory that the mapping for each of the samples will be placed")
    args = parser.parse_args()

    if len(args.bias_dir)==1:
        bias_dir=args.bias_dir #global bias directory to be used for each set of fastq files


    if os.path.exists(args.out_dir):
        shutil.rmtree(args.out_dir)
    os.mkdir(args.out_dir)

    for fastq_dir_index,fastq_dir in enumerate(args.stripped_fastq_set):
        if not os.path.exists(fastq_dir):
            raise ValueError("Directory not found {}".format(fastq_dir))
        if not os.path.isdir(fastq_dir):
            raise ValueError("Path {} is not a directory".format(fastq_dir))
        stripped_fastq_list=os.listdir(fastq_dir)
        stripped_fastq_list=[_ for _ in stripped_fastq_list if ".txt" in _]
        stripped_path_list=[os.path.join(fastq_dir,_) for _ in stripped_fastq_list]
        for stripped_file_index, fastq_file in enumerate(stripped_fastq_list):
            print("Submitting File {} for clustering".format(fastq_file))
            sample_ID=fastq_file.split("_")[0] #get the sample ID: sXX
            map_out_dir=os.path.join(args.out_dir,sample_ID)
            if os.path.exists(map_out_dir):
                shutil.rmtree(map_out_dir)
            os.mkdir(map_out_dir)
            output_file_path=os.path.join(map_out_dir,sample_ID+".csv")
            command="bsub -W 480 -C 0 -n 4"
            stdout_path="_".join(["stdout","align",fastq_file])+".txt"
            stdout_opt="-o "+stdout_path
            stderr_path="_".join(["sterr","align",fastq_file])+".txt"
            stderr_opt="-e "+stderr_path
            project_name_opt="-P DNA_Cluster"
            mem_usage_opt="-R \"rusage[mem=4000]\""
            python_command_opt="python real_analysis/main.py"
            input_opt="--input "+stripped_path_list[stripped_file_index]
            biasconfig_opt = "--biasconfig "+args.bias_config[fastq_dir_index]
            output_opt="--output " + output_file_path
            ed_opt="--editdistance 8"
            biascount_opt="--biascount 20"
            bias_dir_opt="--biasdirectory " + args.bias_dir
            scdir_opt="--scdir starcode/"
            bsub_command=" ".join([command,stdout_opt,stderr_opt,project_name_opt,mem_usage_opt,python_command_opt,input_opt,biasconfig_opt,output_opt,ed_opt,biascount_opt,bias_dir_opt,scdir_opt])
            print("Launching command {} \n".format(bsub_command))
            os.system(bsub_command)
