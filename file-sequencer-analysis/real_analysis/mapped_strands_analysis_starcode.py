#this script provides the ability to go through mapped samples and calculate statistics
import os
import xlsxwriter
import csv
import Levenshtein as lv
import random

#compute the reverse compliment for strands
reverse_table={'A':'T','T':'A','G':'C','C':'G','N':'N'}
def calculate_reverse_compliment(strand):
    reverse=strand[::-1]
    return ''.join([reverse_table[nuc] for nuc in reverse])


def build_sample_pool(sample_path):
    strands=[line.rstrip('\n') for line in open(sample_path)]
    return strands


def parse_range(range_):
    r_lower,r_upper=range_.split('-')
    return int(r_lower),int(r_upper)


#this function goes through all of the mapping files for each of the samples specified
#will create a workbook of spreadsheets, 1 sheet per sample which breaks down the distribution of reads
#for each strand in a file in a sample
def analyze_sample_distributions(path_dictionary,output_book,encoding_directory,fastq_path_dict):
    workbook=xlsxwriter.Workbook(output_book)
    file_length_dict={}
    for sample in sorted(path_dict,key=lambda x: int(x.split('s')[1])):
        analysis_dict={}
        separationFiles={}
        file_pool_dict={}
        sample_pool=build_sample_pool(fastq_path_dict[sample]) #build an array of strands from the original sample
        for mapped_data_path in path_dictionary[sample]:
            mapped_file=open(mapped_data_path,'r')
            mapped_data=csv.reader(mapped_file,delimiter=',')
            for row in mapped_data:
                #if this strand actually maps, proceed
                if int(row[1]) >=0:
                    file_name=row[2]
                    ed=row[3]
                    coding=row[4]
                    file_index=int(row[1])
                    sample_index=int(row[0])
                    #write out the strands associated with each file 
                    if os.path.exists(sample+'_'+file_name) == False:
                        separationFiles[sample+'_'+file_name]=open(sample+'_'+file_name,'w')
                    if file_name not in analysis_dict:
                        file_path=encoding_directory+'/'+file_name
                        analysis_dict[file_name]={}
                        #get the length (number of strands) of the original file
                        file_length_dict[file_name]=len([line.rstrip('\n') for line in open(file_path)])
                        analysis_dict[file_name]["total"]=[0]*file_length_dict[file_name]
                        analysis_dict[file_name]["indexes"]=range(0,file_length_dict[file_name])
                        file_pool_dict[file_name]=[line.rstrip('\n') for line in open(file_path)]
                    if coding=="C":
                        separationFiles[sample+'_'+file_name].write(sample_pool[sample_index]+','+file_pool_dict[file_name][file_index]+'\n')
                    else:
                        separationFiles[sample+'_'+file_name].write(sample_pool[sample_index]+','+calculate_reverse_compliment(file_pool_dict[file_name][file_index])+'\n')
                    if ed not in analysis_dict[file_name]:
                        #initialize an array for the files strands
                        analysis_dict[file_name][ed]=[0]*file_length_dict[file_name]
                    #increment and count for the file index
                    analysis_dict[file_name][ed][file_index]+=1
                    analysis_dict[file_name]["total"][file_index]+=1
        #at this point have a bunch of data in the analysis dictionary, need to dump to a sheet for the sample
        #open a new worksheet
        sheet=workbook.add_worksheet(sample)
        #row and column coordinates
        row=0
        col=0

        
        #write out files to the worksheet
        for found_file in sorted(analysis_dict, key=lambda x: x[0]):
            assert(sheet.write(row,col,found_file) == 0)
            row+=1
            #need to sort out the headers 
            sorted_headers=[0]*len(analysis_dict[found_file])
            ed_index=1
            for ed in analysis_dict[found_file]:
                if ed == "indexes":
                    sorted_headers[0]="indexes"
                elif ed == "total":
                    sorted_headers[-1]="total"
                else:
                    sorted_headers[ed_index]=ed
                    ed_index+=1
            sorted_headers=[sorted_headers[0]]+sorted(sorted_headers[1:-1], key=lambda x: int(x))+[sorted_headers[-1]]
            for ed in sorted_headers:
                sheet.write(row,col,ed)
                row+=1
                for data in analysis_dict[found_file][ed]:
                    assert(sheet.write(row,col,data)==0)
                    row+=1
                col+=1
                row=1
            row=0
            col+=2
        for sFile in separationFiles:
            separationFiles[sFile].close() # close files that hold the separation data
    workbook.close()


#analyze the error rates for each file in a sample and across the whole sample
def analyze_sample_error_rates(path_dictionary,output_book,encoding_directory,fastq_path_dict):
    workbook=xlsxwriter.Workbook(output_book)
    for sample in sorted(path_dict,key=lambda x: int(x.split('s')[1])):
        analysis_dict={}
        file_pool_dict={}
        sample_pool=build_sample_pool(fastq_path_dict[sample])
        for mapped_data_path in path_dictionary[sample]:
            mapped_file=open(mapped_data_path,'r')
            mapped_data=csv.reader(mapped_file,delimiter=',')
            for row in mapped_data:
                #if this strand actually maps, proceed
                if int(row[1]) >=0:
                    sample_index=int(row[0])
                    file_index=int(row[1])
                    file_name=row[2]
                    #print file_name
                    ed=row[3]
                    file_index=int(row[1])
                    coding=row[4]
                    sample_start=int(row[5])
                    sample_end=int(row[6])
                    if file_name not in analysis_dict:
                        file_path=encoding_directory+'/'+file_name
                        analysis_dict[file_name]={}
                        analysis_dict[file_name]["replace"]=[0]*200
                        analysis_dict[file_name]["delete"]=[0]*200
                        analysis_dict[file_name]["insert"]=[0]*200
                        #total counts the total number of mapped reads for a strand
                        analysis_dict[file_name]["total_long"]=0
                        analysis_dict[file_name]["total_short"]=0
                        file_pool_dict[file_name]=[line.rstrip('\n') for line in open(file_path)]
                    short=0
                    #if result is not clean do a editops calculation
                    if ed != '0':
                        #make sure we get the right length file strand for hierarchy
                        if row[7]=="True":
                            file_strand=file_pool_dict[file_name][file_index][20:]
                            analysis_dict[file_name]["total_short"]+=1
                            short=1
                        else:
                            file_strand=file_pool_dict[file_name][file_index][0:]
                            analysis_dict[file_name]["total_long"]+=1
                        if coding == "C":
                            sample_strand=sample_pool[sample_index][sample_start:sample_end]
                        else:
                            sample_strand=calculate_reverse_compliment(sample_pool[sample_index])[sample_start:sample_end]
                        edit_operations=lv.editops(file_strand,sample_strand)
                        #print "file strand: {} sample_strand: {} edit_ops: {}\n".format(file_strand,sample_strand,edit_operations)
                        #print "len file strand: {} len sample_strand: {}".format(len(file_strand),len(sample_strand))
                        for edit in edit_operations:
                            #if the strand is short, we ought to offset it by 20
                            if short==0:
                                analysis_dict[file_name][edit[0]][edit[1]-1]+=1
                            else:
                                analysis_dict[file_name][edit[0]][edit[1]+20-1]+=1
                    else:
                        #increment total counters for the edit distance = 0 case
                        if row[7]=="True":
                            analysis_dict[file_name]["total_short"]+=1
                        else:
                            analysis_dict[file_name]["total_long"]+=1
        #at this point have a bunch of data in the analysis dictionary, need to dump to a sheet for the sample
        #open a new worksheet
        sheet=workbook.add_worksheet(sample)
        #row and column coordinates
        row=0
        col=0
        sorted_headers=["insert","delete","replace","overall"]
        #calculate the overall error rate for the files
        #calculate the total rate across all files in the sample
        analysis_dict["6_total"]={}
        analysis_dict["6_total"]["replace"]=[0]*200
        analysis_dict["6_total"]["delete"]=[0]*200
        analysis_dict["6_total"]["insert"]=[0]*200
        analysis_dict["6_total"]["overall"]=[0]*200
        analysis_dict["6_total"]["total_short"]=0
        analysis_dict["6_total"]["total_long"]=0
        for file_ in analysis_dict:
            if file_!="6_total":
                analysis_dict[file_]["overall"]=[0]*200
                for error_type in analysis_dict[file_]:
                    if error_type!="overall" and error_type!="total_short" and error_type!="total_long":
                        for index, data in enumerate(analysis_dict[file_]["overall"]):
                            analysis_dict[file_]["overall"][index]+=analysis_dict[file_][error_type][index]
                            analysis_dict["6_total"][error_type][index]+=analysis_dict[file_][error_type][index]
                analysis_dict["6_total"]["total_short"]+=analysis_dict[file_]["total_short"]
                analysis_dict["6_total"]["total_long"]+=analysis_dict[file_]["total_long"]
        #calculate the overall error rate across the sample
        for error_type in analysis_dict["6_total"]:
            if error_type!="overall" and error_type!="total_short" and error_type!="total_long":
                for index, data in enumerate(analysis_dict["6_total"]["overall"]):
                    analysis_dict["6_total"]["overall"][index]+=analysis_dict["6_total"][error_type][index]

        #write out files to the worksheet
        for found_file in sorted(analysis_dict, key=lambda x: x[0]):
            assert(sheet.write(row,col,found_file) == 0)
            row+=1
            for error_type in sorted_headers:
                sheet.write(row,col,error_type)
                row+=1
                for index,data in enumerate(analysis_dict[found_file][error_type]):
                    if index<20:
                        divisor=float(analysis_dict[found_file]["total_long"])
                    else:
                        divisor=float(analysis_dict[found_file]["total_long"]+analysis_dict[found_file]["total_short"])
                    if divisor!=0:
                        normalized_data=float(data)/divisor
                    else:
                        normalized_data=float(0)
                    assert(sheet.write(row,col,normalized_data)==0)
                    row+=1
                col+=1
                row=1
            row=0
            col+=2
    workbook.close()



#analyze the error rates for each file in a sample and across the whole sample
def analyze_strand_rate(path_dictionary,output_book,encoding_directory,fastq_path_dict):
    workbook=xlsxwriter.Workbook(output_book)
    for sample in sorted(path_dict,key=lambda x: int(x.split('s')[1])):
        analysis_dict={}
        file_pool_dict={}
        sample_pool=build_sample_pool(fastq_path_dict[sample])
        for mapped_data_path in path_dictionary[sample]:
            mapped_file=open(mapped_data_path,'r')
            mapped_data=csv.reader(mapped_file,delimiter=',')
            for row in mapped_data:
                #if this strand actually maps, proceed
                if int(row[1]) >=0:
                    sample_index=int(row[0])
                    file_index=int(row[1])
                    file_name=row[2]
                    ed=row[3]
                    file_index=int(row[1])
                    coding=row[4]
                    sample_start=int(row[5])
                    sample_end=int(row[6])
                    if file_name not in analysis_dict:
                        file_path=encoding_directory+'/'+file_name
                        file_pool_dict[file_name]=[line.rstrip('\n') for line in open(file_path)]
                        analysis_dict[file_name]={}
                        analysis_dict[file_name]["replace"]=[0]*len(file_pool_dict[file_name])
                        analysis_dict[file_name]["delete"]=[0]*len(file_pool_dict[file_name])
                        analysis_dict[file_name]["insert"]=[0]*len(file_pool_dict[file_name])
                        analysis_dict[file_name]["reads"]=[0]*len(file_pool_dict[file_name])
                    
                    short=0
                    analysis_dict[file_name]["reads"][file_index]+=1
                    #if result is not clean do a editops calculation
                    if ed != '0':
                        #make sure we get the right length file strand for hierarchy
                        if row[7]=="True":
                            file_strand=file_pool_dict[file_name][file_index][20:]
                            short=1
                        else:
                            file_strand=file_pool_dict[file_name][file_index][0:]
                        if coding == "C":
                            sample_strand=sample_pool[sample_index][sample_start:sample_end]
                        else:
                            sample_strand=calculate_reverse_compliment(sample_pool[sample_index])[sample_start:sample_end]
                        edit_operations=lv.editops(file_strand,sample_strand)
                        for op in edit_operations:
                            #accumalate counts for the file index
                            analysis_dict[file_name][op[0]][file_index]+=1
                            
        #at this point have a bunch of data in the analysis dictionary, need to dump to a sheet for the sample
        #open a new worksheet
        sheet=workbook.add_worksheet(sample)
        #row and column coordinates
        row=0
        col=0
        sorted_headers=["insert","delete","replace"]
        #write out files to the worksheet
        for found_file in sorted(analysis_dict, key=lambda x: x[0]):
            assert(sheet.write(row,col,found_file) == 0)
            row+=1
            for error_type in sorted_headers:
                sheet.write(row,col,error_type)
                row+=1
                for index,data in enumerate(analysis_dict[found_file][error_type]):
                    divisor=float(analysis_dict[found_file]["reads"][index])
                    if divisor!=0:
                        normalized_data=float(data)/divisor
                    else:
                        #mark down a missing strand
                        normalized_data=float(-1)
                    assert(sheet.write(row,col,normalized_data)==0)
                    row+=1
                col+=1
                row=1
            row=0
            col+=2
    workbook.close()

#convert the read space file to a dictionary indexed by 'sX' where 'X' is the sample number
def convert_read_space(read_space):
    read_space_dict={}
    read_space_file=open(read_space,'r')
    read_space_csv=csv.reader(read_space,delimiter=',')
    for sample_number, read_count in enumerate(read_space_csv[0]):
        read_space_dict['s'+str(sample_number+1)]=int(read_count)
    return read_space_dict
#function to down sample the read space to down_count
def random_sample(read_space,down_count):
    #if the actual read space is less than the down_count, just return the read space
    if read_space<=down_count:
        return range(0,read_space)
    else:
        return random.sample(range(0,read_space),down_count)
    
#function that down samples the samples to some read count by randomly "selecting" a lower amount of strands that were read 
def analyze_down_sample_distribution(path_dictionary,output_book,encoding_directory,down_count,read_space):
    workbook=xlsxwriter.Workbook(output_book)
    file_length_dict={}
    read_space_dict=convert_read_space(read_space)
    for sample in sorted(path_dict,key=lambda x: int(x.split('s')[1])):
        analysis_dict={}
        for mapped_data_path in path_dictionary[sample]:
            mapped_file=open(mapped_data_path,'r')
            mapped_data=csv.reader(mapped_file,delimiter=',')
            sampled_array=random_sample(read_space_dict[sample],down_count)
            for row in mapped_data:
                #if this strand actually maps and is in the sampled_array, proceed
                if int(row[1])>=0 and int(row[0]) in sampled_array:
                    file_name=row[2]
                    ed=row[3]
                    file_index=int(row[1])
                    if file_name not in analysis_dict:
                        file_path=encoding_directory+'/'+file_name
                        analysis_dict[file_name]={}
                        #get the length (number of strands) of the original file
                        file_length_dict[file_name]=len([line.rstrip('\n') for line in open(file_path)])
                        analysis_dict[file_name]["total"]=[0]*file_length_dict[file_name]
                        analysis_dict[file_name]["indexes"]=range(0,file_length_dict[file_name])
                    if ed not in analysis_dict[file_name]:
                        #initialize an array for the files strands
                        analysis_dict[file_name][ed]=[0]*file_length_dict[file_name]
                    #increment and count for the file index
                    analysis_dict[file_name][ed][file_index]+=1
                    analysis_dict[file_name]["total"][file_index]+=1
        #at this point have a bunch of data in the analysis dictionary, need to dump to a sheet for the sample
        #open a new worksheet
        sheet=workbook.add_worksheet(sample)
        #row and column coordinates
        row=0
        col=0


def analyze_total_reads(path_dictionary,output_book,encoding_directory):
    workbook=xlsxwriter.Workbook(output_book)
    file_length_dict={}
    row_start=0
    sheet=workbook.add_worksheet('total reads found')
    for sample in sorted(path_dict,key=lambda x: int(x.split('s')[1])):
        analysis_dict={}
        for mapped_data_path in path_dictionary[sample]:
            mapped_file=open(mapped_data_path,'r')
            mapped_data=csv.reader(mapped_file,delimiter=',')
            for row in mapped_data:
                #if this strand actually maps, proceed
                if int(row[1]) >=0:
                    file_name=row[2]
                    ed=row[3]
                    file_index=int(row[1])
                    if file_name not in analysis_dict:
                        analysis_dict[file_name]=0;
                    #increment and count for the file index
                    analysis_dict[file_name]+=1
        #at this point have a bunch of data in the analysis dictionary, need to dump to a sheet for the sample
        #row and column coordinates
        col=0
        sheet.write(row_start,col,sample)
        col+=1

        row=row_start
        #write out data to the worksheet
        for found_file in sorted(analysis_dict, key=lambda x: x[0]):
            if(row_start==0):
                assert(sheet.write(row,col,found_file) == 0)
                row+=1
            #need to sort out the headers 
            sheet.write(row,col,analysis_dict[found_file])
            row=row_start
            col+=2
        if(row_start==0):
            row_start+=2
        else:
            row_start+=1
    workbook.close()





   
if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Perform analysis on mapping results")
    parser.add_argument('--sample_range',dest="sample_range",action="store",default="1-10", help="Range for samples.")
    parser.add_argument('--down_count',dest='down_count',action='store',default=87900,type=int,help='Number to down sample the read space to')
    
    #what action the script should make 
    parser.add_argument('--sample_distributions',dest='sample_dist',action='store_true',default=False,help='Calculates the distribution of strand reads for each sample, puts them into a excel spreadsheet')
    parser.add_argument('--error_rates',dest='error_rate',action='store_true',default=False,help='Calculates the error rates for each base and each error type for a sample')
    parser.add_argument('--strand_rate',dest='strand_rate',action='store_true',default=False,help='Calculates the rate for strands in each file in each sample')
    parser.add_argument('--down_sample_distribution',dest="down_sample",action='store_true',default=False,help='Calculate distributions that are down sampled')
    parser.add_argument('--total_reads',dest='total_reads',action='store_true',default=False,help="Calculate the total reads for each file in a sample")
    #Directory arguments
    parser.add_argument('--primer_partition_top_directory',dest='part_dir',action="store", default="mapped_strands",help="path to top partition directory")
    parser.add_argument('--original_encoding_directory',dest='enc_dir',action="store",default="original_encodings", help="path to original file encodings")
    parser.add_argument('--stripped_fastq_directory',dest='stripped_fastq',action="store",default="StrippedFastQ", help="path to stripped FASTQ files")
    parser.add_argument('--output_book',dest='book_name',action="store",default="distribution.xlsx",help="name of the excel workbook to be generated")
    parser.add_argument('--read_space_file',dest='read_space',action="store",default="NGS_20180716_reads.csv",help="name of csv that holds the read space counts for each sample")
    args = parser.parse_args()

    sample_lower,sample_upper=parse_range(args.sample_range)
        

    
    sample_names=[]
    path_dict={}
    
    for i in range(sample_lower,sample_upper+1):
        sample_names.append('s'+str(i))

    #create paths of relevant partition files 
    for sample in sample_names:
        path_dict[sample]=[]
        path_dict[sample].append('/'.join([args.part_dir,sample,sample+'.csv']))
    
    #create paths to relevant sample files
    sample_nums = ['s'+str(x)+'_' for x in range(sample_lower,sample_upper+1)]
    sample_files=[]
    #extract only the sample files that are desired
    for file_ in os.listdir(args.stripped_fastq):
        for sample in sample_nums:
            if sample in file_:
                sample_files.append(file_)
    fastq_path_dict={}
    for sample_file_name in sample_files:
        fastq_path_dict[sample_file_name.split('_')[0]]='/'.join([args.stripped_fastq,sample_file_name])

        
    if(args.sample_dist):
        analyze_sample_distributions(path_dict,args.book_name,args.enc_dir,fastq_path_dict)
    if(args.error_rate):
        analyze_sample_error_rates(path_dict,args.book_name,args.enc_dir,fastq_path_dict)
    if(args.strand_rate):
        analyze_strand_rate(path_dict,args.book_name,args.enc_dir,fastq_path_dict)
    if(args.down_sample):
        analyze_down_sample_distribution(path_dict,args.book_name,args.enc_dir,args.down_count,args.read_space)
    if(args.total_reads):
        analyze_total_reads(path_dict,args.book_name,args.enc_dir)
