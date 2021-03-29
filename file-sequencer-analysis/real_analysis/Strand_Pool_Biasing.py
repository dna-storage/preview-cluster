'''
This function biases several encoded files according to the configuration file
and dumps the biased strands into a txt file.
'''


def file_length(filename):
    index = 1  # This ensures the index is recorded where the biasing starts.
    with open(filename) as f:
        for i in f:
            index += 1
    return index


def remove_strand(strand, primerlength):
    # Removes the first x bases from the strand
    return strand[primerlength:]


def compliment(base):
    # Converts all bases to their compliment
    # (Should have used a dictionary in hindsight)
    if base == 'A':
        return 'T'
    if base == 'T':
        return 'A'
    if base == 'G':
        return 'C'
    if base == 'C':
        return 'G'
    if base == '\n':
        return ''


def reverse_compliment(strand):  # Takes reverse compliment of the strand
    strandreverse = strand[::-1]  # Reverses the strand
    strandnc = list(strandreverse)  # Converts the strand into a list
    for base in range(len(strandreverse)):
        # Iterates through and returns the base compliment
        strandnc[base] = compliment(strandnc[base])
    return "".join(strandnc)  # Reconverts the list back to a string


'''
This procedure takes all of the configurations for the specific file and adds
the resulting biased strands (regular and reverse compliment) to the array,
performing a hierarchy operation if needed as well.
'''


def bias_strand(fileinput, bias_array, hierarchy, primerlength, biascount,
                filename, strand_dict):
    if hierarchy:  # If true, run without hierarchy first.
        bias_strand(fileinput, bias_array, False, primerlength,
                    biascount, filename, strand_dict)
    with open(fileinput, 'r') as f:
        index = 0  # For file index tracking (0-based)
        for strand in f:
            if hierarchy:  # If true, remove the primer.
                strand = remove_strand(strand, primerlength)
            dictionary_addition(strand.rstrip(), strand_dict, index, hierarchy,
                                filename, "C")  # Adds first strand to dict
            for i in range(biascount):
                # Append the array with the strand by the number of times
                # specified by bias count
                bias_array.append(strand.rstrip('\n'))
            strandnc = reverse_compliment(strand)  # Convert reverse compliment
            for i in range(biascount):  # Same idea here
                bias_array.append(strandnc)
            dictionary_addition(strandnc.rstrip(), strand_dict, index,
                                hierarchy, filename, "NC")
            # Adds nc strand to dict
            index += 1  # Increment index for next strand


def dictionary_addition(strand, strand_dict, index, hierarchy, filename,
                        coding):
    strand_dict[strand] = {"index": index, "filename": filename,
                           "coding": coding, "hierarchy": hierarchy}


def biasing(biasdirectory, biasfile, config, strand_dict, biascount):
    import argparse
    import json

    print(json.load(open(config)))
    print("biasdirectory: {}".format(biasdirectory))
    bias_array = []  # Sets up empty bias array
    strand_dict["startindex"] = file_length(biasfile)

    with open(biasfile, 'r') as f:
        bias_array = [strand.rstrip('\n') for strand in f]

    [bias_strand(biasdirectory + _file["input_file"], bias_array,
                 _file["hierarchy"], _file["primer_length"], biascount,
                 _file["input_file"], strand_dict)
        for _file in json.load(open(config))["configs"]]

    return bias_array
