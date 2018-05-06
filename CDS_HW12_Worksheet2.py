import simplealign as sal
import blosum 


# s1 = 'RNDKPKFSTARN'
# s2 = 'AAAAARNQKPKWWTATN'
# v = sal.BruteForceSlide(blosum.BLOSUM50, blosum.PBET, s1, s2)
# print(len(s2) - v.argmax())
# s1 = 'AAAAAAARNDKPKFSTARN'
# s2 = 'RNQKPKWWTATN'
# v = sal.BruteForceSlide(blosum.BLOSUM50, blosum.PBET, s1, s2)
# print(len(s2) - v.argmax())
# print(s1)
# print(7* '.' +s2)

import dynprog as dpg 
s1 = 'IQIFSFIFRQEWNDA'
s2 = 'QIFFFFRMSVEWND'
scormat, arrow = dpg. ScoringMatrix (blosum. BLOSUM50, blosum. PBET, s1, s2)

def Backtrace (arrow, seq1, seq2): 
    st1, st2 = '','' 
    v, h = arrow. shape 
    ok = 1 
    v-=1 
    h-=1 
    while ok: 
        if arrow [v, h] == 0: 
            st1 += seq1 [v-1] 
            st2 += '-'
            v -= 1 
        elif arrow [v, h] == 1: 
            st1 += '-'
            st2 += seq2 [h-1] 
            h -= 1 
        elif arrow [v, h] == 2: 
            st1 += seq1 [v-1] 
            st2 += seq2 [h-1] 
            v -= 1 
            h -= 1 
        if v ==0 and h ==0: 
            ok = 0 
#reverse the strings 
    st1 = st1 [::-1] 
    st2 = st2 [::-1] 
    return st1, st2 
st1, st2 = dpg. Backtrace (arrow, s1, s2) 
print(st1) 
print(st2)

import blosum
sq1 = 'KMTIFFMILK'
sq2 = 'NQTIFF'
subvals = dpg. FastSubValues (blosum.BLOSUM50, blosum.PBET, sq1, sq2)
scmat, arrow = dpg. FastSW (subvals, sq1, sq2)
t1, t2 = dpg. SWBacktrace (scmat, arrow, sq1, sq2)

print(subvals, scmat, arrow, t1, t2)