# simplealign.py
# Python 3.4
# Python scripts are provided as an educational tool. They are offered without guarantee of effectiveness or accuracy. Python scripts composed by the author may not be used for commercial uses without the author's explicit written permission.
# (c) 2016 Jason M. Kinser

import numpy as np

def SimpleScore( s1, s2 ):
    a1 = np.array(list(s1))
    a2 = np.array(list(s2))
    # count matches
    score = (a1==a2).astype(int).sum()
    # count mismatches
    score = score -(a1!=a2 ).astype(int).sum()
    # gaps
    ngaps = s1.count( '-' ) + s2.count('-')
    score = score - ngaps
    return score

def BlosumScore( mat, abet, s1, s2, gap=-8 ):
    sc = 0
    n = min( [len(s1), len(s2)] )
    for i in range( n ):
        if s1[i] == '-' or s2[i] == '-' and s1[i] != s2[i]:
            sc += gap
        elif s1[i] == '.' or s2[i] == '.':
            pass
        else:
            n1 = abet.index( s1[i] )
            n2 = abet.index( s2[i] )
            sc += mat[n1,n2]
    return sc

def BruteForceSlide( mat, abet, seq1, seq2 ):
    # length of strings
    l1, l2  = len( seq1 ), len( seq2 )
    # make new string with leader
    t1 = len( seq2) * '.' + seq1
    lt = len( t1 )
    answ = np.zeros( lt, int )
    for i in range( lt ):
        answ[i] = BlosumScore( mat, abet, t1[i:], seq2 )
    return answ
