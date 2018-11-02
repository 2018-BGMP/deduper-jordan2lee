#!/usr/bin/env python

#SBATCH --partition=short  
#SBATCH --job-name=deduper_data2
#SBATCH --time=1-00:00:00   
#SBATCH --nodes=1         
#SBATCH --ntasks-per-node=14    
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jordanl.7890@gmail.com

file="/projects/bgmp/jlee33/bi624_gen_lab/3_deduper/data/Dataset1.sam"
umi_file="/projects/bgmp/jlee33/bi624_gen_lab/3_deduper/data/STL96.txt"


##########################################
############ Problem: #######################
##########################################


#PCR duplicates need to be removed from 
#SAM files before downstream analysis can continue

# input: SAM file of uniquely mapped reads
# full files on Talapas /projects/bgmp/shared/deduper/

#Tips on getting right most position POS with cigar string (minus strand)
		#https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
		
##########################################
########## Arg Parse #########################
##########################################

# import argparse

# def get_arguments():
    # parser = argparse.ArgumentParser(description='Purpose: remove PCR duplicates from sam file. Input must be sorted by chromosome then by POS, only single end reads, cigar string can only have M, N, D, I, or S')
    # parser.add_argument("-f", "--file", help ="absolute path to sorted sam file", required=True, type=str)
    # parser.add_argument("-p", "--paired", help ="file contains paired end reads", required=False)
    # parser.add_argument("-u", "--umi_file", help ="designates file containing the list of UMIs (unset if randomers instead of UMIs)", required=True, type=str)
    # return parser.parse_args()
    
# args = get_arguments()
# file = args.file
# paired = args.paired
# umi_file = args.umi_file


# This commented out section was used only for debugging purposes. Ignore.
# file = "/mnt/c/Users/Jordan/UO_Docs/fall_bi624_genomics_lab/deduper-jordan2lee/data/test_files/pairedend1_test.sam"
#file = "/mnt/c/Users/Jordan/UO_Docs/fall_bi624_genomics_lab/deduper-jordan2lee/data/test_files/singleend1_testSORTED.sam"
#umi_file = "/mnt/c/Users/Jordan/UO_Docs/fall_bi624_genomics_lab/deduper-jordan2lee/data/STL96.txt"


#########################################
##############  Functions ####################
#########################################

import re


def check_paired(flag):
    '''input SAM file FLAG and returns True (paired end data) or False (single end data)'''
    strand = ((flag & 1) == 1)    #means 1 in that position
    return strand

def get_UMI(QNAME):
	'''input QNAME and will return string of UMI sequence'''
	QNAME = QNAME.split(":")
	UMI = QNAME[7]
	return UMI

def max_cigar(CIGAR):
    '''Input cigar of minus strand and output max value of sum(M+N+I+D+S2). For determining sliding window size of minus strand'''
    # Grab softclipping if occur last on cigar --> output 10S
    S2 = re.findall('\d*S$',CIGAR)
    # Grab M, N, I, D --> output 10M10N...
    cig = re.findall('\d*[MNID]',CIGAR)
    # Convert to str
    S2 = "".join(S2)
    cig = "".join(cig)
    # combine into one max cigar str
    max = cig + S2
    # grab just numeric --> output 10
    max = re.findall(r'\d+', max)
    # sum of max
    max = sum([int(i) for i in max])
    return max
	
def expect_UMI(UMI):
	'''input UMI seq and output True (if seq matches with any expected UMI seqs) or False (if seq doesn't match with any expected UMI seqs)'''
	for key in umi_dict:
		#now will enter this "if" statement for first key in dict, 2nd iteration will check the next key in dict, etc
		if UMI == key:
			#UMI is in dict --> return true + exit function
			return UMI == key
		else:
			#if UMI isn't in dict (doesn't match with any expected UMI seqs)
			#and continue through fucntion by repeating above "if" for the next key in dict
			pass
	#returns false if UMI seq doesn't match with any seq in dict
	return False
	
def check_strand(flag):
    '''input SAM file FLAG and returns plus (plus strand) or minus (minus strand)'''
    strand = ((flag & 16) != 16)    #means 1 in that position
    if strand == True:
        return "plus"
    return "minus"

def adjust_POS_plus(CIGAR, POS):
    '''Input cigar of plus strand and output adjusted POS. Based on adjPOS=POS-S1. Output as str'''
    # Grab softclipping if occur FIRST on cigar --> output 10S
    softclip = re.search(r"^(\d+)S",CIGAR)
    # if no softclippping at start then adjPOS is original POS
    if softclip is None:
        adjPOS = POS
    # if there is softclipping at start then...
    else:
	    # grab the int (10S --> 10)
        S1 = int(softclip.group(1))
        adjPOS = POS - S1
    return adjPOS

def adjust_POS_minus(CIGAR, POS):
	'''Input cigar of minus strand and output adjusted POS. based on adjPOS = (POS -1) + M + D + S2 + N'''
	#grab softclipping if occur last on cigar 
	softclip2 = re.findall(r'\d*S$',CIGAR)
	#if there is softclipping
	if softclip2 is not None:
		#S2 = softclipping at end of cigar
		S2 = softclip2
		cig = re.findall('\d*[MND]',CIGAR)
		# Combine all into one str
		S2 = "".join(S2)
		cig = "".join(cig)
		total = S2 + cig
		# Grab just numeric values (10M -->10)
		total = re.findall(r'\d+', total)
		# Sum total
		total = sum([int(i) for i in total])
		# Calc adjPOS = (POS -1) + M + D + S2 + N
		adjPOS = POS - 1 + total
	return adjPOS

	
##########################################
############# Dedup   Algorithm ################
##########################################

# Create dict of plus and minus strand of entries already seen. Removing entries as they are no longer a possible dup
plus_dict = {}  # for strand +     key=concat of "umi_rname_adjPOS" value= full read line
minus_dict = {}   # for strand -     key=concat of "umi_rname_adjPOS" value= full read line

			
# Create umi_dict with expected UMI sequences from UMI file inputed into argparse
		#key = UMI seq   value= 0
umi_dict = {}
with open(umi_file, "r") as umi_fh:
	for line in umi_fh:
		line = line.strip()
		umi_dict[line] = 0

############ Read SAM First Time ############################
############ Purpose: Define Size of Sliding Window ################

# input must be SAM tool sorted file by chromosome then by POS		
with open(file, "r") as sam, open(file[:-4] + "_deduped.sam", "w") as out_fh:
	for line in sam:
		if line.startswith("@"):
			#do nothing and continue for header lines (aka start with @)
			continue
		
		# Now only have alignment lines.....so
		#split line --- so can index specific columns
		line=line.split("\t")
		# Important info to determine sliding window: cigar, minus strand, rname(chrom)
				#rename obj for readability
		flag = int(line[1])    	#tells if minus/plus strand ...plus more info if needed
		rname = line[2]  		#tells what chromosome on 
		cigar = line[5]           #tells how many matches, mismatch, deletions, splicing, etc
		
		# Double check user hasn't entered paired end data
		if check_paired(flag) == True:   #found paired end data
			print("ERROR: Paired end data found. Input only single ended data")
			break


# print("Plus strand sliding window size: ", plus_wind)
# print("Minus strand sliding window size: ", minus_wind_dict)

############ Read SAM Second Time ############################
############ Purpose: Remove PCR Duplicates ######################


plus_ct = 1  #tells first line has value 1
minus_ct = 1 #tell first line has value of 1
past_plus="" #initialize previous
current_plus=""
past_minus=""
current_minus="" 

# input must be SAM tool sorted file by chromosome then by POS		
with open(file, "r") as sam, open(file[:-4] + "_deduped.sam", "w") as out_fh:
	for line in sam:
		#write header lines directly to output file
		if line.startswith("@"):
			out_fh.write(line)

		else:
			#split line --- so can index specific columns
			line=line.split("\t")
			#rename obj for readability
			qname = line[0]
			flag = int(line[1])
			rname = line[2]
			pos = int(line[3])
			cigar = str(line[5])
			num = line[6]  #just for debugging purposes
			umi = get_UMI(qname)


			# Check that UMI seq is one of the expected UMIs
			if expect_UMI(umi) == False:
				# print("Unexpected UMI found on read", line)
				continue #skip this entry and restart with next entry in SAM file
			
			# Now: Separate potential pcr duplicates by strand+ and strand- ##
			#by incorporating sliding windows (plus_wind or minus_wind_dict)
			elif expect_UMI(umi) == True:
				# First create a str to compare against keys in plus_dict/minus_dict
				info = [get_UMI(qname), str(rname), str(adjust_POS_plus(cigar, pos))] 
				#concat info into one str
				info = "_".join(info)   # info = UMI_chrom_adjPOS
				
				###### 1. PLUS stranded reads ######
				if check_strand(flag) == "plus":							
					#First check chromosomes
					
					# if first line for that chrom
					if plus_ct == 1:
						#set previous = rname
						past_plus = rname	

					# if same chrom     as seen in past
					if rname == past_plus:
						# Check for duplicates
						
						#if duplicate found
						if info in plus_dict:
							#increment counter
							plus_ct+=1
						
						#if no duplicate found  (info not in plus_dict)
						else:
							#add info to plus dict
							plus_dict[info] = "none"
							#write to outputfile
							line ="\t".join(line)
							out_fh.write(line) 
							#increment counter
							plus_ct+=1
							
							
					# if different chrom    as seen in past
					elif rname != past_plus:
						#clear plus_dict
						plus_dict = {}
						#add info to plus_dict
						plus_dict[info] = "none"
						#update past_plus to current chrom
						past_plus = rname
						#write to outputfile
						line ="\t".join(line)
						out_fh.write(line) 
						#increment counter
						plus_ct+=1

						
				###### 2. MINUS stranded reads	########			
				elif check_strand(flag) == "minus":
					#First check chromosomes
					
					# if first line for that chrom
					if minus_ct == 1:
						#set previous = rname
						past_minus = rname	

					# if same chrom     as seen in past
					if rname == past_minus:
						# Check for duplicates
						
						#if duplicate found
						if info in minus_dict:
							#increment counter
							minus_ct+=1
						
						#if no duplicate found  (info not in minus_dict)
						else:
							#add info to plus dict
							minus_dict[info] = "none"
							#write to outputfile
							line ="\t".join(line)
							out_fh.write(line) 
							#increment counter
							minus_ct+=1
							
							
					# if different chrom    as seen in past
					elif rname != past_minus:
						#clear plus_dict
						minus_dict = {}
						#add info to plus_dict
						minus_dict[info] = "none"
						#update past_minus to current chrom
						past_minus = rname
						#write to outputfile
						line ="\t".join(line)
						out_fh.write(line) 
						#increment counter
						minus_ct+=1



