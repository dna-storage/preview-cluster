'''
Filename: preview_strip.py
Description: Take in a file that represents a graduated preview architecture and strip it down to the canonical amplification primer strands for use with clustering algorithms.

File Assumption: Strands are ordered with respect to the hamming distance graduation. E.g. The first strand is the 0HD.

Author: Kevin Volkel
'''
import Levenshtein as ld

#calculate how much end match there is between two strands, foward is used to toggle where the 3 primer end refernece is
def end_match(forward,string1,string2):
    _string2=string2
    _string1=string1
    if forward==False:
        _string2=string2[::-1]
        _string1=string1[::-1]
    #padding out strings so their 3' ends match up
    if len(_string1)>len(_string2):
        _string2=["x"]*(len(_string1)-len(_string2))
        _string2="".join(_string2)
        _string2=_string2+string2
        assert(len(_string1)==len(_string2))
    elif len(_string2)>len(_string1):
        _string1=["x"]*len(_string2)-len(_string1)
        _string1="".join(_string1)
        _string1=_string1+string1
    for moves in range(0,len(string1)):
        if _string1[moves::]==_string2[moves::]:
            return len(string1)-moves
    return 0

def calculate_HD(string1,string2):
    assert len(string1)==len(string2)
    count=0
    for i in range(len(string1)):
        if string1[i]!=string2[i]:
            count+=1
    return count

#slide the primer to find candidates
def slide_primer(strand,primer,forward,HD_cutoff,prime_cutoff):
    mod_list=[]
    assert len(strand)>len(primer)
    for j in range((len(strand)-len(primer))+1):
        window=strand[j:j+len(primer)]
        prime_match=end_match(forward,window,primer)
        HD=calculate_HD(window,primer)
        if prime_match>=prime_cutoff or HD<=HD_cutoff:
            candidate_tuple=(j,j+len(primer),prime_match,HD) #(start,end (non-inclusive),prime_match score, HD score)
            mod_list.append(candidate_tuple)
    return mod_list

#return True if strand in list with ed<=cutoff
def check_ed(strand_list,strand,cutoff):
    for s in strand_list:
        if ld.distance(s,strand)<=cutoff: return True
    return False


import os
import pandas as pd
if __name__=="__main__":
    import argparse
    
    parser=argparse.ArgumentParser(description="File that modifies preview encoded files to form that of what we expect after hamming distance amplification")
    parser.add_argument('--strip_path',action="store",type=str,default=None,help="File to strip")
    parser.add_argument('--out_dir',action="store",type=str,default=None,help="Directory to place the stripped file")
    parser.add_argument('--end_prime_cutoff',action="store",type=int,default=6,help="If at least these many three prime matches, mark as possible amplification spot")
    parser.add_argument('--HD_cutoff',action="store",type=int,default=8,help="If less than or equal to this HD mark as possible amplification spot")
    parser.add_argument('--ed_cutoff',action="store",type=int,default=8,help="If there is a strand with this ED or less, don't add it")
    parser.add_argument('--keep_na',action="store_true",default=False,help="If True, keep n/a strands")
    args=parser.parse_args()

    if args.strip_path==None:
        assert 0
    if args.out_dir==None:
        assert 0

    if not os.path.exists(args.out_dir):
        os.mkdir(args.out_dir)
    HD=args.HD_cutoff
    end_match_value=args.end_prime_cutoff
    assert os.path.exists(args.out_dir) and os.path.exists(args.strip_path)
    preview_file=open(args.strip_path,'r')
    strand_lines=preview_file.readlines()
    stripped_list=[]
    problematic_list=[]
    canonical_forward=None
    canonical_backward=None
    region_counter=0
    current_region=None
    
    seed_strand_stats={"canonical_ID":[],"is_canonical":[],"HD-F":[],"HD-R":[],"3-P-Match-F":[],"3-P-Match-R":[]}
    for strand_index, strand_string in enumerate(strand_lines):
        strand=strand_string.strip()
        error_strand=None
        if strand_index==0:
            #get the canonical 0HD primers
            canonical_forward=strand[20:40]
            canonical_backward=strand[-40:-20]

        strand_backward=strand[-40:-20]
        strand_forward=strand[20:40]
        assert canonical_backward!=None and canonical_forward!=None
        new_strand=canonical_forward+strand[40:-40]+canonical_backward
        stripped_list.append(new_strand)
        #get some stats on the canonical strand
        seed_strand_stats["canonical_ID"].append(strand_index)
        seed_strand_stats["is_canonical"].append(True)
        seed_strand_stats["HD-F"].append(calculate_HD(strand_forward,canonical_forward))
        seed_strand_stats["HD-R"].append(calculate_HD(strand_backward,canonical_backward))
        seed_strand_stats["3-P-Match-F"].append(end_match(True,strand_forward,canonical_forward))
        seed_strand_stats["3-P-Match-R"].append(end_match(False,strand_backward,canonical_backward))

    #go back over the strand list to find non-canonical outcomes
    for strand_index, strand_string in enumerate(strand_lines):
        strand=strand_string.strip("\n")
        forward_outcome=slide_primer(strand,canonical_forward,True,HD,end_match_value)
        reversed_outcome=slide_primer(strand,canonical_backward,False,HD,end_match_value)
        #First apply just the forward and reverse outcomes
        if args.keep_na:
            for (forward_tuple,reversed_tuple) in zip(forward_outcome,reversed_outcome):
                f_start,f_end,f_prime,f_HD=forward_tuple
                r_start,r_end,r_prime,r_HD=reversed_tuple
                f_new_strand=canonical_forward+strand[f_end::]
                if f_new_strand not in stripped_list and not check_ed(stripped_list,f_new_strand,args.ed_cutoff): 
                    stripped_list.append(f_new_strand)
                    seed_strand_stats["canonical_ID"].append(strand_index)
                    seed_strand_stats["is_canonical"].append(False)
                    seed_strand_stats["HD-F"].append(f_HD)
                    seed_strand_stats["HD-R"].append("n/a")
                    seed_strand_stats["3-P-Match-F"].append(f_prime)
                    seed_strand_stats["3-P-Match-R"].append("n/a")
                r_new_strand=strand[0:r_start]+canonical_backward
                if r_new_strand in stripped_list or check_ed(stripped_list,r_new_strand,args.ed_cutoff): continue 
                stripped_list.append(r_new_strand)
                seed_strand_stats["canonical_ID"].append(strand_index)
                seed_strand_stats["is_canonical"].append(False)
                seed_strand_stats["HD-F"].append("n/a")
                seed_strand_stats["HD-R"].append(r_HD)
                seed_strand_stats["3-P-Match-F"].append("n/a")
                seed_strand_stats["3-P-Match-R"].append(r_prime)

        #now apply all combos together to get the effect of two amplifications in both off target directions
        for forward_tuple in forward_outcome:
            f_start,f_end,f_prime,f_HD=forward_tuple
            for reversed_tuple in reversed_outcome:
                r_start,r_end,r_prime,r_HD=reversed_tuple
                if r_start<f_end: continue #regions overlap
                new_strand=canonical_forward+strand[f_end:r_start]+canonical_backward
                if new_strand in stripped_list or check_ed(stripped_list,new_strand,args.ed_cutoff): continue
                stripped_list.append(new_strand)
                seed_strand_stats["canonical_ID"].append(strand_index)
                seed_strand_stats["is_canonical"].append(False)
                seed_strand_stats["HD-F"].append(f_HD)
                seed_strand_stats["HD-R"].append(r_HD)
                seed_strand_stats["3-P-Match-F"].append(f_prime)
                seed_strand_stats["3-P-Match-R"].append(r_prime)
    #dump stripped list out to a file
    preview_file.close()
    stripped_list=stripped_list
    base_name=os.path.basename(args.strip_path)
    base_name_no_extension=base_name.split(".")[0]
    dump_path=os.path.join(args.out_dir,base_name)
    dump_file=open(dump_path,'w+')
    for new_strand in stripped_list:
        dump_file.write(new_strand+"\n")
    dump_file.close()
    #open pandas file for dumping stats
    stats_path=os.path.join(args.out_dir,"stats_"+base_name_no_extension+".csv")
    stats_dataframe=pd.DataFrame(seed_strand_stats)
    stats_dataframe.to_csv(stats_path,index=False)
    
