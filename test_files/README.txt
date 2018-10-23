################################################
## Message:

test1.sam is created testfile and 
correct_testoutput_deduper.sam is the expected file output if deduper algorithm performs correctly

#Tips on getting right most position POS with cigar string (minus strand)
		#https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/
		
#################################################
## TEST file need to contain:
-unexpected UMI
-UMI with N
-chromosomes
-strand+ and strand-
-strand+  = soft clipping before M, M
-strand-  = soft clippping after M, M, insert, delete, splicing
-leftmost start POS


unsorted. max transcript length is 101 (max length of Illumina read) but file contains trimmed  (potentially quality, adapter, and index trimmed)
#################################################


## Test file contains entry of 
1. paired end read (to test error) 
2. UMI has N in seq
		CTGTNCAC, rest can be anything e.g. strand+, chrom=3, POS=10, 97M
3. UMI unexpected

[dups are 4,5, 6,7,8]
4. strand+ , UMI matches with #5, chrom matches with #5, same POS as #5, cigarstr
		 strand+, CTGTTCAC, chrom=3, POS = 10, 100M
		[perfect M]
		
5. strand+, UMI matches with #4, chrom matches with #4, same POS as #4, citgarstr
		strand+, CTGTTCAC, chrom=3, POS = 10, 50M 40S
		[all same but S1]
		
6. strand+, UMI matches with #4, chrom matches with #4, differ POS as #4, cigarstr
		strand+, CTGTTCAC, chrom=3, POS = 60, 50S 50M
		[differ POS but S1 allows for same adjPOS]
		
7. strand+, UMI matches with #4, chrom matches with #4, differ POS as #4, cigarstr
		strand+, CTGTTCAC, chrom=3, POS = 60, 50S 10M 10N 30M
		[differ POS but S1 allowd ror same adjPOS and has N  in cigar]
		
8. strand+, UMI matches with #4, chrom matches with #4, differ POS as #4, cigarstr
		strand+, CTGTTCAC, chrom=3, POS = 60, 50S 10M 20I 10D 10M
		[differ POS but S1 allows for same adjPOS and has I/D in cigar]
		
9. strand+ that has no potential duplicates
		strand+, ATCCATGG, chrom=3, POS=10, 60M 10D 10M
		[good strand+ read with no potential duplicates]
		
10. strand+, UMI matches with #5, chrom doesnt match with #5
		forward, CTGTTCAC, chrom=2 [e.g. POS=10, 50M 40S]

		
[dups are 11, 12,13, 14]
11. strand-, UMI matches with 12, chrom matches with 12
		strand-, GTTCACCT, chrom=2, POS=10, 100M
		[perfect M] adjPOS=110
		
12. strand-, UMI matches with 11, chrom matches with 11, same POS as #11, cigarstr
		strand-, GTTCACCT, chrom=2, POS=20, 10S 90M
		[S1 present but doesn't do anything]   adjPOS=110bp

13. strand-, UMI matches with 11, chrom matches with 11, same POS as #11, cigarstr
		strand-, GTTCACCT, chrom=2, POS=10, 90M 10S
		[differ POS but S2 allows for same adjPOS]  adjPOS=110bp

14. strand-, UMI matches with 11, chrom matches with 11, same POS as #11, cigarstr
		strand-, GTTCACCT, chrom=2, POS=10, 10S 10M 10N 20D 10I  50M
		[contains S1, M, N, D and I]  adjPOS=110

15. strand-, UMI unqiue, chrom matches with 11, same POS as #11, cigarstr
		strand-, TTCGCCTA, chrom=2, POS=10, 100M
		adjuPOS=110
		
16. strand- , same as #4 but minus strand
		 strand-, CTGTTCAC, chrom=3, POS = 10, 100M
	
	
	