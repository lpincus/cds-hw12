# Problems 
# 1. Create a random sequence and copy it. In the copy remove a couple of letters at different places. Use NW to align these two sequences. 
# 2. Repeat the above problem but change the gap penalty. Does the alignment change if the gap penalty is -16? Does it change if it is -2? 
# 3. Create a scoring matrix which is, M[i, j] = (5 i = j −1 i 6= j . Align two amino acid sequences (of at least 100 characters) using the BLOSUM50 matrix and the above M matrix. Are the alignments significantly different? 
# 4. Modify the BlosumScore algorithm to align DNA strings such that the 3-rd element in each codon is weighted half as much as the other two.
# 5. Create a string with the form XAXBX where X is a set of random letters and A and B are specific strings designed by the user. Each X can be a different length. Create a second string with the form YAYBY where Y is a different set of random letters and each Y can have a different length. Align the sequences using Smith-Waterman. The scoring matrix will have two major maximum for the alignments of the A and B regions. Modify the program to extract both alignments.
import numpy as np
import random as rd
import string

N=12
r1 = "abc-ABC"
r2 = r1
r3 = ''.join(rd.choices(string.ascii_uppercase + string.ascii_lowercase, k=N))
r4 = r3
r5 = ''.join(rd.choices(string.ascii_uppercase + string.ascii_lowercase, k=N))
r6 = r5

def SimpleScore ( s1, s2):
    a1 = np.array(list(s1))
    a2 = np.array(list(s2))
    score = (a1==a2).astype(int).sum()
    score = score -(a1!=a2 ).astype(int).sum()
    ngaps = s1.count( '-' ) + s2.count('-')
    score = score - ngaps
    print(score)
    return score

SimpleScore (r1,r2)
SimpleScore (r3,r4)
SimpleScore (r5,r6)