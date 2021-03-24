'''
This program takes the output from the starcode on a biased txt file and the
.p strand dictionary file to output a .csv with information on each strand
for further analysis.
'''


def Edit_Distance(sampleindex, originalstrand, canonicalstrand):
    from Levenshtein import distance
    return distance(originalstrand[sampleindex].rstrip(),
                    canonicalstrand.rstrip())


def File_Analysis(clusters, strand_dict, bias_array, fileoutput):
    import csv
    with open(fileoutput, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        clusters = clusters.splitlines()
        for line in clusters:
            strand = line.split('\t')
            # Starcode separates each group via tab & each index by ','
            sampleindices = strand[2].strip('\n').split(',')
            # print int(strand[1]) < 5 or len(strand[0]) < 10
            if int(strand[1]) < 5 or len(strand[0]) < 10:
                # Enter null values for each index if < 5
                for index in sampleindices:
                    writer.writerow([index, "-1", "N/A", "N/A", "NC",
                                    "N/A", "N/A", "N/A"])
            else:
                for index in sampleindices:
                    if int(index) < strand_dict["startindex"]:
                        try:  # If there is a dict error, null the entry.
                            output = strand_dict[strand[0]]
                        except:
                            writer.writerow([int(index) - 1, "-1", "N/A",
                                             "N/A", "NC", "N/A", "N/A",
                                             "N/A"])
                        else:
                            distance = Edit_Distance(int(index) - 1,
                                                     bias_array,
                                                     strand[0])
                            writer.writerow([int(index) - 1,
                                            output["index"],
                                            output["filename"], distance,
                                            output["coding"], 0,
                                            len(strand[0]),
                                            output["hierarchy"]])
