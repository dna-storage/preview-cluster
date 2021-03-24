if __name__ == "__main__":
    import argparse
    import Analysis
    import Strand_Pool_Biasing
    import subprocess

    parser = argparse.ArgumentParser(description="Parse file paths")

    # stripped file to analyze
    parser.add_argument("--input", help="Stripped Input")
    parser.add_argument("--biasconfig", help="Config File")  # config file
    parser.add_argument("--biasdirectory", default="",
                        help="Bias Strand Directory")  # dir of original strand
    parser.add_argument("--output", help="Spreadsheet Output")  # .csv output
    parser.add_argument("--editdistance", type=int,
                        help="Starcode Edit Distance")
    parser.add_argument("--biascount", type=int, help="Bias Count")
    parser.add_argument("--scdir", help="Starcode Directory")

    args = parser.parse_args()

    strand_dict = {}

    print "Biasing...\n"
    bias_array = Strand_Pool_Biasing.biasing(args.biasdirectory, args.input,
                                             args.biasconfig, strand_dict,
                                             args.biascount)
    print "Complete!\n"
    print "Now running starcode...\n"

    process = subprocess.Popen([args.scdir + "./starcode", "-d" +
                               str(args.editdistance), "-q", "--seq-id"],
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)

    starout = process.communicate('\n'.join(bias_array))
    clustered = starout[0]

    with open("Clusters.DEBUG", 'w') as w:
        w.write(clustered)

    print "Starcode Complete!\n"
    print "Starting File Analysis...\n"

    Analysis.File_Analysis(clustered, strand_dict, bias_array, args.output)

    print "File Analysis Complete!\n"
    print "Program finished!\n"
