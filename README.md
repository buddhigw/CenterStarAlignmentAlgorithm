# CenterStarAlignmentAlgorithm
Center Star Alignment Method
Idea: Find a sequence which is the most similar to all the rest and then it use as the center of the star to align all the other sequences to it.

S1= ATTCGGATT
S2= ATCCGGATT
S3= ATGGAATTTT
S4= ATGTTGTT
S5= AGTCAGG
1 Step: Compute all the pairwise scores. 

In this part, 2 sequences taken from the array.
From the for loops call one by one character in these two sequences. And the sequences must not equal one to another
These characters store in tow char arrays as sequences as align1 and align 2
Then these two called in the method scoreTable
scoreTable calculate the sum of pair score of each pair.

The TcoreTable function, create the score table of each pair of sequences. 
After that, get the traceback result of each pair and finally get the finalize result of each pair and store them in 2D matrix as scorematrix. 

2 Step: Calculate the row wise summation and find the maximum. The maximum summation sequence is the center sequence of the alignment.



3 Step: Calculate the optimal pairwise alignments between S1 and all other sequences.
	S1: A T T C G G A T T
	S2: A T C C G G A T T

	S1: A T T C G G A – T T - -
	S3: A T - - G G A A T T T T 

	S1: A T – T C G G A T T
	S4: A T G T T G - - T T

	S1: A T T C – G G A T T
	S5: A G T C A G G - - -

3 Step: Find the maximum lengthened among the aligned sequences and merge all the alignment using “once a gap always a gap” principle. Then add gaps in the end of the sequences until the length of the sequence equal to maximum length.


Finding max lengthened sequence

Start with S1 and S2, then add S3, then two spaces need to be added to the ends of S1 and S2.
S2: A T C C G G A T T - -
S1: A T T C G G A T T - -
S3: A T G – G A A T T T T

Then add S4 and S5
S2: A T C C G G A T T - -
S1: A T T C G G A T T - -
S3: A T G – G A A T T T T
S4: A T G T T G – T T - -
S5: A G T C A G G - - - -

Running Time Analysis:
Suppose all sequences in S are of length O(n).
Step 1: Computing the optimal distance D(Si,Sj) for all i,j: O(k2n2)
Step 2. Computing the center string Sc: O(k2)
Step 3: Generating the k pairwise alignments with Sc: O(kn2)
Step 4: Insert spaces into Sc in order to satisfy all multiple alignments of step 3 simultaneously: O(k2n)
Total running time: O(k2n2)

Examples with different sequence sets.

1. "AGTG","ATCC","ATCG","TCCT","TCGA","TGCG"
Result:-
Maximum score: 90
Center Sequence: ATCG
 Max len seq: AGT-G
AGT-G
ATCG-
ATCC-
-TCCT
-TCGA
-TGCG

2. "ATCGTGGTACTG"
"CCGGAGAACTAG"
"AACGTGCTACTG"
"ATGGTGAAGTG"
"CCGGAAAACTTG"
"TGGCCCTGTATC"
Result:-
Maximum score: 265
Center Sequence: ATCGTGGTACTG
 Max len seq: -TGGCCCT-GTA-TC
-TGGCCCT-GTA-TC
ATCGTGGTACTG---
-CCG-GAGAACTAG-
AACGTGCTACTG---
ATGGT-GAAGTG---
-CCG-GAAAACTTG-

Here I attached the algrithm which was built using Python and C programming languages
