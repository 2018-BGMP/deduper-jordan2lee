#!/usr/bin/env python3

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


# Include flags: 
#-f, --file: required arg, absolute file path
#-p, --paired: optional arg, designates file is paired end (not single-end)
#-u, --umi: optional arg, designates file containing the list of UMIs (unset if randomers instead of UMIs)
#-h, --help: optional arg, prints a USEFUL help message (see argparse docs)

# Plus flag for if input paired end ---> return an error message and quit

# --help message along the lines of:
#input data must not have any quality trimming done
#data can only have insertions to the ref (I), deletions from the ref (D), alignment match/mismatch (M), splicing (N)
#input data must have already had Samtools to sort by chromosome
#no randomers and only single end reads (no paired end reads)
#cigar string can only have M, N, D, I, or S


#########################################
##############  Functions ####################
#########################################


# Create check_strand function
		#input: FLAG
		#body: takes the FLAG (col 2) and checks if entry is strand+ (0) or strand- (1)
		#strand-_strd = ((flag & 16) == 16)    #means 1 in that position
		#return: strand-_strd == True
				#meaning: False (strand+), True (strand-)

# Create check_UMI function 
		#input: UMI seq (as str)
		#body: grab QNAME by using regex
		#return: True (UMI matches with expected UMIs), False (UMI doesn't match with expected UMIs)
		
#function: adjust_strand+_POS()
	#input cigar string for strand+ only, and adjust POS value to leftmost position
	#body: adjusted = POS - S1  
			#(where S1 is integer of softclipping occuring before any M) 
			#note: ignore M, D, N, and I 
	#return: integer of new adjusted POS value 

#function: adjust_strand-_POS()
	#input cigar string for strand- strand only, and adjust POS value to rightmost position
	#body: adjusted = POS +S2 + M +N +I + D 
			#(where S2 is integer of softclipping occuring after any M)
			#note: S1
	#return: integer of new adjusted POS value

# Create check_paired()
	#input: FLAG 
	#body: paired = ((flag & 1?) == 1?)
	#return: paired == True
		#meaning: true (paired end data), false (single stranded data)

		
##########################################
############# Dedup   Algorithm ################
##########################################

## input must be SAM tool sorted file by chromosome


## Check that data is paired end (not single end) (FLAG in SAM col 2)
if check_paired == True     #contains pairded end data:
	write error message "paired end data found, only single end data allowed"
	break # exit algorithm

	
## Determine size of sliding window
# create empty dictionaries
dict+ = {}  # for strand +     key=cigarstring letter (S1,S2,M,N,D,I)   value=max corresponding number
dict- = {}   # for strand -     key=cigarstring letter (S1,S2,M,N,D,I)   value=max corresonding number
if strand+:
	check if corresponding number of the cigarStrLetter is greater than value in dict:
		update dict with that value
if strand-:
	check if corresponding number of the cigarStrLetter is greater than value in dict:
		update dict with that value 
#after iterate thru entire file: save the max corresponding number as an object
S1_+ = value of dict+
...etc
S1_- = value of dict_
...etc
# create size of sliding window as per objects printed out
		#thus this will be different for each input sam file

		
# Set up output SAM file
elif see line start with @:
	write header lines to deduper.sam
	
## Check for N in UMI seq (different than N in CIGAR)
# if N found in UMI seq
elif check_unknown_seq == True:
	continue  # stop looking at that entry
	
## Check that UMI seq is one of the expected UMIs
#if UMI seq not one we expected
elif check_UMI == FALSE:
	continue # stop looking at that entry
	

## Now: Separate potential pcr duplicates by strand+ and strand- ##

#incorporate size of sliding window to read through file

# strand +: potential PCR duplicates
elif same UMI {QNAME col 1} and chromosome (RNAME col 3) and strand+ (via check_strand ==False):
	# calcu adjusted POS value -- aka take into account softclippping, insertions, deletions, splicing
	if adjust_strand+_POS() of one entry == adjust_strand+_POS() of any other entry:
		#Found duplicate
		write first entry to deduper.sam (and do nothing with duplicates)
		
		
# strand-: potential PCR duplicates
elif same UMI (QNAME col 1) and chromosome (RNAME col 3) and strand- (via check_strand ==True):
	# calcu adjusted POS value  -- aka take into account softclippping, insertions, deletions, splicing
	if adjust_strand-_POS() of one entry == adjust_strand-_POS() of any other entry:
		#Found duplicate
		write first entry to deduper.sam (and do nothing with duplicates)

# Entries with no potential PCR duplicates
else:
	# not a duplicate so write to good file
	write entry to deduper.sam