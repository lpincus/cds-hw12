# Problems 

def prob1():
    """ 1. Create a random sequence and copy it. In the copy remove a couple of letters at different places. Use NW to align these two sequences. 
    """
    import dynprog as dyn
    import blosum as b
    seq1='ASDFASDFASDFASDF'
    seq2='ASFASDFASDASDF'
    subvals=dyn.FastSubValues(b.BLOSUM50,b.PBET,seq1,seq2)
    scormat, arrow = dyn.FastNW(subvals,seq1,seq2)
    print(scormat, arrow)
# prob1()

def prob2():
    """2. Repeat the above problem but change the gap penalty. Does the alignment change if the gap penalty is -16? Does it change if it is -2? 
    """
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
# prob2()

# 3. Create a scoring matrix which is, 
#     M[i, j] = (5  i = j 
#               (−1 i != j 
#   Align two amino acid sequences (of at least 100 characters) 
#   using the BLOSUM50 matrix and the above M matrix.
#   Are the alignments significantly different? 
def prob3():
    from blosum import BLOSUM50, PBET
    import dynprog as dpg
    import random as rd
    N=101
    rd.seed(23)
    seq1 = ''.join(rd.choices(PBET, k=N))
    seq2 = ''.join(rd.choices(PBET, k=N))
    subvals = dpg.FastSubValues(BLOSUM50, PBET, seq1, seq2)
    schmat, arrow = dpg.FastNW(subvals, seq1, seq2)
    t1, t2 = dpg.SWBacktrace(schmat, arrow, seq1, seq2)
    print(t1, t2)


# 4. Modify the BlosumScore algorithm 
# to align DNA strings such that the 3-rd element in each codon 
# is weighted half as much as the other two.
def BlosumScore2( mat, abet, s1, s2, gap=-8 ):
    sc = 0
    n = min( [len(s1), len(s2)] )
    for i in range( n ):
        if s1[i] == '-' or s2[i] == '-' and s1[i] != s2[i]:
            if i // 3 == 0:
                sc += gap * 0.5
            sc += gap
        elif s1[i] == '.' or s2[i] == '.':
            pass
        else:
            n1 = abet.index( s1[i] )
            n2 = abet.index( s2[i] )
            if i // 3 == 0:
                sc += mat[n1,n2] * 0.5
            sc += mat[n1,n2]
    return sc
def prob4():
    from blosum import BLOSUM50, PBET
    from simplealign import BlosumScore
    import dynprog as dpg
    import random as rd
    N=101
    rd.seed(23)
    seq1 = ''.join(rd.choices(PBET, k=N))
    seq2 = ''.join(rd.choices(PBET, k=N))
    score = BlosumScore(BLOSUM50, PBET, seq1, seq2)
    score2 = BlosumScore2(BLOSUM50, PBET, seq1, seq2)
    print(score, score2)
# prob4()   

# 5. Create a string with the form XAXBX where 
# - X is a set of random letters and 
# - A and B are specific strings designed by the user. 
# Each X can be a different length. 
def stringXAXBX(a, b):
    import string 
    import random as rd
    from blosum import PBET
    rd.seed(123)
    N = max(len(a), len(b))
    Nx = rd.randint(1, N)
    X = ''.join(rd.choices(PBET, k=Nx))
    Ny = rd.randint(1, N)
    Y = ''.join(rd.choices(PBET, k=Ny))
    Nz = rd.randint(1, N)
    Z = ''.join(rd.choices(PBET, k=Nz))
    return '{}{}{}{}{}'.format(X,a,Y,b,Z)

# Create a second string with the form YAYBY where 
# - Y is a different set of random letters 
#   - each Y can have a different length. 
# Align the sequences using Smith-Waterman. 
# The scoring matrix will have two major maximum for the alignments of the A and B regions. 
# Modify the program to extract both alignments.
import dynprog as dpg
from blosum import BLOSUM50, PBET
import random as rd
a = ''.join(rd.choices(PBET, k=12))
b = ''.join(rd.choices(PBET, k=12))
Xstr = stringXAXBX(a, b)
Ystr = stringXAXBX(a, b)

subvals = dpg.FastSubValues(BLOSUM50, PBET, Xstr, Ystr)
scmat, arrow = dpg.FastSW(subvals, Xstr, Ystr)
t1, t2, = dpg.SWBacktrace(scmat, arrow, Xstr, Ystr)
print(t1, t2)


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
