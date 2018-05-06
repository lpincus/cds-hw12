import numpy as np 

#easyalign .py 2 
def BruteForceSlide (mat, abet, seq1, seq2): 
    l1, l2 = len (seq1), len (seq2) 
    t1 = len (seq2) * '.' + seq1 
    lt = len (t1) 
    answ = np. zeros (lt, int) 
    for i in range (lt): 
        answ [i] = BlosumScore (mat, abet, t1 [i:], seq2) 
    return answ
