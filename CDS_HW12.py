# Problems 
# 1. Create a random sequence and copy it. In the copy remove a couple of letters at different places. Use NW to align these two sequences. 
import dynprog as dyn
import blosum as b
seq1='ASDFASDFASDFASDF'
seq2='ASFASDFASDASDF'
subvals=dyn.FastSubValues(b.BLOSUM50,b.PBET,seq1,seq2)
scormat, arrow = dyn.FastNW(subvals,seq1,seq2)
print(scormat, arrow)

# 2. Repeat the above problem but change the gap penalty. Does the alignment change if the gap penalty is -16? Does it change if it is -2? 
import dynprog as dyn
import blosum as b
seq1='ASDFASDFASDFASDF'
seq2='ASFASDFASDASDF'
subvals=dyn.FastSubValues(b.BLOSUM50,b.PBET,seq1,seq2)
scormat, arrow = dyn.FastNW(subvals,seq1,seq2)
print(scormat, arrow)
scormat, arrow = dyn.FastNW(subvals,seq1,seq2, -16)
print(scormat, arrow)
scormat, arrow = dyn.FastNW(subvals,seq1,seq2, -2)
print(scormat, arrow)

# 3. Create a scoring matrix which is, M[i, j] = (5 i = j −1 i 6= j . 
# Align two amino acid sequences (of at least 100 characters) using the BLOSUM50 matrix and the above M matrix.
# Are the alignments significantly different? 


# 4. Modify the BlosumScore algorithm 
# to align DNA strings such that the 3-rd element in each codon 
# is weighted half as much as the other two.
 
 
# 5. Create a string with the form XAXBX where 
# X is a set of random letters and A and B are specific strings designed by the user. 
# Each X can be a different length. 
# Create a second string with the form YAYBY where Y is a different set of random letters and each Y can have a different length. 
# Align the sequences using Smith-Waterman. 
# The scoring matrix will have two major maximum for the alignments of the A and B regions. 
# Modify the program to extract both alignments.



# import numpy as np
# import random as rd
# import string

# N=12
# r1 = "abc-ABC"
# r2 = "abd-ABd"
# r3 = ''.join(rd.choices(string.ascii_uppercase + string.ascii_lowercase, k=N))
# r4 = ''.join(rd.choices(string.ascii_uppercase + string.ascii_lowercase, k=N))
# r5 = ''.join(rd.choices(string.ascii_uppercase + string.ascii_lowercase, k=N))
# r6 = ''.join(rd.choices(string.ascii_uppercase + string.ascii_lowercase, k=N))

# def SimpleScore ( s1, s2):
#     a1 = np.array(list(s1))
#     a2 = np.array(list(s2))
#     score = (a1==a2).astype(int).sum()
#     score = score -(a1!=a2 ).astype(int).sum()
#     ngaps = s1.count( '-' ) + s2.count('-')
#     score = score - ngaps
#     print(score)
#     return score

# SimpleScore (r1,r2)
# SimpleScore (r3,r4)
# SimpleScore (r5,r6)

# # # CDS_HW12.py
# import simplealign as sal
# sal.SimpleScore( 'AGTCGATCGATT', 'AGTCGATCGATT')
# sal.SimpleScore( 'AGTCGATCGATT', 'AGTCGATCGAAT')
# sal.SimpleScore( 'AGTCGATCGATT', 'AGTCGATCGA-T')


#simplealign .py 2 
# def SimpleScore (s1, s2): 
#   3 a1 = np. array (list (s1)) 4 a2 = np. array (list (s2)) 5 score = (a1 == a2). astype (int). sum () 6 score = score -( a1 != a2). astype (int). sum () 7 ngaps = s1. count (✬-✬) + s2. count (✬-✬) 8 score = score - ngaps 9 return score 10 11 
# >>> import simplealign as sal 
# 12 >>> sal. SimpleScore (✬ AGTCGATCGATT ✬, ✬ AGTCGATCGATT ✬) 
# 13 12 14 >>> sal. SimpleScore (✬ AGTCGATCGATT ✬, ✬ AGTCGATCGAAT ✬) 
# 15 10 16 >>> sal. SimpleScore (✬ AGTCGATCGATT ✬, ✬ AGTCGATCGA-T✬)

# Kinser, Jason. Computational Methods for Bioinformatics: Python 3.4 (Page 374).  . Kindle Edition. 


# 3	Problem
# Load and run the code presented in slides 
#   16 and 18 of the video (Code 25.2 and 25.4 in the text).

# import blosum
# blosum.BLOSUM50[:4 ,:4]
# # array([[5, -2, -1, -2], 4 [-2, 7, -1, -2], 5 [-1, -1, 7, 2], 6 [-2, -2, 2, 8]])
# blosum.PBET
# # Kinser, Jason. Computational Methods for Bioinformatics: Python 3.4 (Page 37).  . Kindle Edition. 
# from Bio import pairwise2
# from Bio.SubsMat import MatrixInfo as matlist
# from Bio.SubsMat.MatrixInfo import blosum62 as blosum
# matrix = matlist.blosum50
# print(blosum.BLOSUM50[:4 ,:4])
# print(blosum.PBET)
# print(matrix)
# print(type(matrix))
