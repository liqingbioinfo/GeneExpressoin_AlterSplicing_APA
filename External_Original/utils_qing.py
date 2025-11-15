import pandas as pd
import numpy as np
import os
import sys
import re

def getGenesNameFromGtf(gftFileName):
    GenesNames = []
    with open(gftFileName,"r") as f:
        line=f.readline().strip()
        while line:
            gene_info = line.split("\t")[9]
            gene_info_arr = gene_info.split(";")
            for item in gene_info_arr:
                if item.startswith(" gene_name "):
                    gene_name = item.split(" ")[-1]
                    if gene_name not in GenesNames:
                        GenesNames.append(gene_name)
                    break
            line=f.readline().strip()
    return(set(GenesNames))

def getENSGNameFromGtf(gftFileName):
    GenesNames = []
    with open(gftFileName,"r") as f:
        line=f.readline().strip()
        while line:
            gene_info = line.split("\t")[9]
            gene_info_arr = gene_info.split(";")
            for item in gene_info_arr:
                if item.startswith(" gene_id "):
                    gene_name = item.split(" ")[-1]
                    if gene_name not in GenesNames:
                        GenesNames.append(gene_name)
                    break
            line=f.readline().strip()
    return(set(GenesNames))


def getSampleIDsTfam(tfam_file):
    ##tfam is not seperated by \t, it is seperated by " "
    SampleIDs=[]
    with open(tfam_file,"r") as f:
        line=f.readline().strip()
        while line:
            line_arr=line.split(" ")
            if len(line_arr)==6:
                SampleIDs.append(line_arr[0])
            else:
                print("Error in spliting lines, please check")
                exit()
            line=f.readline().strip()
    return SampleIDs
    
def tped2tsv(tped_file, tfam_file, output_file,header=False, ref_alt="NA"):
    ##tped is not seperated by \t, it is seperated by " "
    SampleIDs=getSampleIDsTfam(tfam_file)
    output_file_bw=open(output_file,"w")
    output_header="\"varID\""
    for SampleID in SampleIDs:
        output_header+="\t\""+SampleID+"\""
    output_file_bw.write(output_header+"\n")
    
    with open(tped_file,"r") as f:
        if header:
            f.readline()
        line=f.readline().strip()
        ###So line_arr[1] contains reference and alternative allele info
        if ref_alt=="NA":
            while line:
                line_arr = line.split(" ")
                chr_,pos,ref,alt,build=line_arr[1].split("_")
                chr_num=re.search(r'\d{1,2}',chr_).group()
                s=str(chr_num)+"_"+pos+"_"+ref+"_"+alt+"_"+build
                for i in range(4,len(line_arr),2):
                    if(line_arr[i]==ref and line_arr[i+1]==ref):
                        s+="\t0"
                    elif(line_arr[i]==alt and line_arr[i+1]==alt):
                        s+="\t2"
                    else:
                        s+="\t1"
                output_file_bw.write(s+"\n")
                line=f.readline().strip()                    
        else:
            print("TBF")
    output_file_bw.close()

            
# def two2one(snp):
    # result="N"
    # if(snp=="A A" or snp=="C C" or snp=="G G" or snp=="T T"):
        # result=snp.split(" ")[0]
    # else:
        # print("HET!")
        # if(snp=="A C" or snp=="C A"): result= "M"
        # elif(snp=="A G" or snp=="G A"): result="R"
        # elif(snp=="A T" or snp=="T A"): result="W"
        # elif(snp=="C G" or snp=="G C"): result="S"
        # elif(snp=="C T" or snp=="T C"): result="Y"
        # elif(snp=="G T" or snp=="T G"): result="K"
        # else: print(snp+": Not a correct SNP!")
    # return result

##Stop and return the index of a value that nearest larger than target value
def binary_search(a, x, lo=0, hi=None):
    if hi is None:
        hi = len(a)
    while lo < hi:
        mid = (lo+hi)//2
        midval = a[mid]
        if midval < x:
            lo = mid+1
        elif midval > x: 
            hi = mid
        else:
            return mid
    return mid
    
def get_startList_and_startEndDict(bedfile,header=True):
    chr_pos_start_lists=[]
    chr_pos_start_lists_sorted=[]
    chr_pos_start_end_dict={}
    for i in range(22):
        chr_pos_start_lists.append(list())

    f=open(bedfile,"r")
    if header:
        f.readline()
    line=f.readline().strip()
    while line:
        if "chrX" not in line and "chr23" not in line and "23" not in line:
            line_arr = line.split("\t")
            if "chr" in line_arr[0] or "CHR" in line_arr[0]:
                chr_=int(re.search(r'\d{1,2}',line_arr[0]).group())
            else:
                chr_=int(line_arr[0])
            chr_pos_start_lists[chr_-1].append(int(line_arr[1]))
            key=str(chr_)+"_"+line_arr[1]
            chr_pos_start_end_dict[key]=key+"_"+line_arr[2]
        line=f.readline().strip()
    f.close()
    
    for i in range(22):
        chr_pos_start_lists_sorted.append(sorted(chr_pos_start_lists[i]))
        
    return chr_pos_start_lists_sorted,chr_pos_start_end_dict


def number_SNPs_in_out_of_func_regions(target_SNPs_set,bedfile):

    flank=0
    SNP_in_func_regions=[]
    SNP_outof_func_regions=[]
    
    chr_pos_start_lists_sorted,chr_pos_start_end_dict = get_startList_and_startEndDict(bedfile)
    for snp in target_SNPs_set:
        snp_arr = snp.split("_")
        if "chr" in snp_arr[0] or "CHR" in snp_arr[0]:
            re_res = re.search(r'\d{1,2}',snp_arr[0])
            if re_res:
                chr_=str(re_res.group())
            else: ##This snp does not from autosomes, skip. 
                continue
        else:
            chr_=str(snp_arr[0])
        
        pos = snp_arr[1]
        if not chr_=="23":
            chr_tmp=chr_pos_start_lists_sorted[int(chr_)-1]
            pos_index = binary_search(chr_tmp,int(pos))
            if pos_index >0:
                front_start = chr_tmp[pos_index-1]
            else:
                front_start = chr_tmp[pos_index]
            behind_start = chr_tmp[pos_index]

            front_start_end = chr_pos_start_end_dict.get(chr_+"_"+str(front_start))
            behind_start_end = chr_pos_start_end_dict.get(chr_+"_"+str(behind_start))

            ##check if the pos is located in the front or behind gene
            _,fs,fe=front_start_end.split("_")
            _,bs,be=behind_start_end.split("_")

            if (int(pos)<=int(fe) and int(pos)>=int(fs)) or (int(pos)<=int(be) and int(pos)>=int(bs)):
                SNP_in_func_regions.append(snp)
#                 print(snp,fs,fe,bs,be)
            else:
                SNP_outof_func_regions.append(snp)  
#         else:
#             print(snp)
    return SNP_in_func_regions,SNP_outof_func_regions    

                
                
                
                
        