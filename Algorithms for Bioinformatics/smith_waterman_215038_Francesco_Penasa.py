#!/usr/bin/env python
"""
Implementation of Smith Watermann algorithm in python3
@author: Francesco Penasa 215038

Input: sequence 1 and sequence 2 
Output: input sequences and best local alignment, with score

Input: you can use any solution 
    (e.g. the sys.argv function to deal with the 2 input sequences and the score values, 
    the argparse module to create a user-friendly interface)
Output: the minimum output printed to screen is the best sub-alignment, 
    with the corresponding score. 
    More information (e.g. length of the subalignment, number of matches, number of mismatches, number of gaps) is facultative 
Evaluation: the script must return the correct solution 
    (the best local alignment, given 2 sequences and given a score system). 
    Additional points are given if the script is user-friendly and well documented
"""


# === imports === #
import numpy as np
import sys
# =============== #


# === functions === #
def init_scoring_matrix(seq1, seq2):
	"""
	Initialize the scoring matrix all to zero, 
	the scoring matrix will have the shape of 
	len(seq1) + 1 as columns and 
	len(seq2) + 1 as rows.
	The original algorithm ask to initialize only the first row and column,
	but it results more easy to just initialize all the matrix to 0.
	"""

	scoring_matrix = np.zeros(((len(seq1) + 1,len(seq2) + 1,)))

	# this is the step that why should have been done to init the first row and column to zero 
	# if the matrix was not initialized to zero.
	#for i in range(len(seq1)):
	#	scoring_matrix[i,0] = 0;
	#for j in range(len(seq2)):
	#	scoring_matrix[0,j] = 0;

	return scoring_matrix

def scoring(scoring_matrix, seq1, seq2, match, mismatch, gap):
	"""
	Update the scoring_matrix row after row with the new values calculated 
	through match score, mismatch score and gap penality.
	"""

	for i in range(1, len(seq1) + 1):		# i = row
		for j in range(1, len(seq2) + 1):	# j = column
			tmp_m = 0; # where to store match or mismatch result		
			a = seq1[i-1] # i-1 because the matrix rows are == seq1 + 1
			b = seq2[j-1] # j-1 because the matrix columns are == seq2 + 1

			if(a == b): # match
				tmp_m = match + scoring_matrix[i-1, j-1]	# match case
			else: # mismatch  
				tmp_m = mismatch + scoring_matrix[i-1, j-1]	# mismatch case
		
			# gaps
			gap_r = gap + scoring_matrix[i-1, j] # gap from the left element
			gap_c = gap + scoring_matrix[i, j-1] # gap from the upper element

			# assign the max value between match (or mismatch) gap (on row or on column) and zero.
			scoring_matrix[i, j] = max(tmp_m, gap_r, gap_c, 0)

	return scoring_matrix

def traceback(scoring_matrix, r,c, match,mismatch,gap, seq1,seq2):
    return traceback_rec(scoring_matrix, r,c, match,mismatch,gap, seq1,seq2, [""],[""])

def traceback_rec(scoring_matrix, i,j, match,mismatch,gap, seq1,seq2, align1,align2):
    """
	traceback function done in a recursive way to discover all the possible alignments (dynamic programming technique).
	starting at the element with the score equal to the one sent in input, traceback until a 0 is encountered.
	"""
    if scoring_matrix[i,j] == 0:
        return align1,align2
    else:
        all_alignments1, all_alignments2 = [],[]
        tmp_list1, tmp_list2 = [],[]

        if ((scoring_matrix[i-1,j-1]+match == scoring_matrix[i,j] and seq1[i-1] == seq2[j-1]) # match
        or (scoring_matrix[i-1,j-1]+mismatch == scoring_matrix[i,j] and seq1[i-1] != seq2[j-1])): # mismatch
            a1,a2 = align1.copy(),align2.copy()
            a1[-1] = seq1[i-1]+a1[-1]
            a2[-1] = seq2[j-1]+a2[-1]
            tmp_list1, tmp_list2 = traceback_rec(scoring_matrix, i-1,j-1, match,mismatch,gap, seq1,seq2, a1,a2)
            all_alignments1, all_alignments2 = all_alignments1+tmp_list1, all_alignments2+tmp_list2

        if scoring_matrix[i, j-1]+gap == scoring_matrix[i, j]: # gap left       
            a1,a2 = align1.copy(),align2.copy()
            a1[-1] = "-"+a1[-1]
            a2[-1] = seq2[j-1]+a2[-1]
            tmp_list1, tmp_list2 = traceback_rec(scoring_matrix, i,j-1, match,mismatch,gap, seq1,seq2, a1,a2)
            all_alignments1, all_alignments2 = all_alignments1+tmp_list1, all_alignments2+tmp_list2

        if scoring_matrix[i-1, j]+gap == scoring_matrix[i, j]: # gap top
            a1,a2 = align1.copy(),align2.copy()
            a1[-1] = seq1[i-1] + a1[-1]
            a2[-1] = "-" + a2[-1]
            tmp_list1, tmp_list2 = traceback_rec(scoring_matrix, i-1,j, match,mismatch,gap, seq1,seq2, a1,a2)
            all_alignments1, all_alignments2 = all_alignments1+tmp_list1, all_alignments2+tmp_list2
        
        return all_alignments1, all_alignments2 


def find_top_score(scoring_matrix):
    """
    find the max score in the matrix and 
    return such score and the indexes of all the cells with such score
    """
    score = np.amax(scoring_matrix)
    
    if score == 0:
        print("There are no possible alignment")
        sys.exit(3) # code error 3 = there are no alignments

    r,c = np.where(scoring_matrix == score)
    if (not r[0]):
        r = [r,]
    if (not c[0]):
        c = [c,]
    return score, r, c

def help():
    """
    Show the usage manual
    """
    print("Minimum Usage:")
    print("python smith_waterman.py seq1 seq2")
    print("There are some additional options that you can use.")
    print("-i --ifile > To specify an input file use: python smith_waterman.py -i <filename> ")
    print("-o --ofile > To specify an output file use: python smith_waterman.py -o <filename> ")
    print("-int --interactive > To specify the sequence during runtime: python smith_waterman.py -int ")
    sys.exit(1)

def parse_input_file(file):
    """
    Parse the content of a file assuming the FASTA format
    """
    fasta = {}
    with open(file) as f:
        fasta = {line.strip(">\n"):next(f).rstrip() for line in f}
    items = []
    for elem in fasta.values():
        items.append(elem)
    
    return items[0], items[1]

def check_input(argv):
    """
    Check input arguments
    """
    inputfile = ""
    outputfile = ""
    seq1 = ""
    seq2 = ""

    if len(sys.argv) < 2:
        print("Invalid number of input arguments:")
        help() # show the manual
        

    for i in range(1, len(argv)):
        if argv[i] == '-h' or argv[i] == '--help':
            help()
            
        elif argv[i] == '-i' or argv[i] == '--ifile':
            inputfile = argv[i+1]
            seq1, seq2 = parse_input_file(argv[i+1])

        elif argv[i] == '-o' or argv[i] == '--ofile':
            outputfile = argv[i+1]
            sys.stdout = open(outputfile, "w") # redirect std output of prints to a file
            
        elif argv[i] == '-int' or argv[i] == '--interactive':
            seq1 = input("Insert the first sequence: ")
            seq2 = input("Insert the second sequence: ")
            
        elif i == 1: # if no option have been used
            seq1 = argv[i]
            seq2 = argv[i+1]

    if len(seq1) < 1 or len(seq2) < 1:
        print("Error: parse of sequence failed!")
        
    print ('Input file is "', inputfile)
    print ('Output file is "', outputfile)
    print(sys.argv)
    return seq1, seq2

        
        
if __name__ == '__main__':

    # INPUT
    seq1, seq2 = check_input(sys.argv)

    print("Sequence_1 = \"" + seq1 + "\"")  # on the row of the matrix
    print("Sequence_2 = \"" + seq2 + "\"")  # on the column of the matrix

    # MATCH SCORE, MISMATCH SCORE AND GAP PENALTY
    match = 3
    mismatch = -3
    gap = -2

    # INITIALIZE THE SCORING MATRIX 
    scoring_matrix = init_scoring_matrix(seq1, seq2)

    # SCORING
    scoring_matrix = scoring(scoring_matrix, seq1, seq2, match, mismatch, gap)
    #print("Scoring matrix:")
    #print(scoring_matrix)

    # find top score
    score, rows, columns = find_top_score(scoring_matrix)
    print("Score:")
    print(score)

    # TRACEBACK
    seq1_alignments = []
    seq2_alignments = []

    for i in range(len(rows)):
        alignments1, alignments2 = traceback(scoring_matrix, rows[i],columns[i], match,mismatch,gap, seq1,seq2)
        seq1_alignments = seq1_alignments + alignments1
        seq2_alignments = seq2_alignments + alignments2

    print("Results:")
    print("There are ",len(seq1_alignments), " different alignments")
    print("")

    for i in range(len(seq1_alignments)):
        print(seq1_alignments[i])
        print(seq2_alignments[i])
        print("")