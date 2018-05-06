# dynprog.py
# Python 3.4
# Python scripts are provided as an educational tool. They are offered without guarantee of effectiveness or accuracy. Python scripts composed by the author may not be used for commercial uses without the author's explicit written permission.
# (c) 2016 Jason M. Kinser

import numpy as np

def ScoringMatrix( mat, abet, seq1, seq2, gap=-8 ):
    l1, l2 = len( seq1), len(seq2)
    scormat = np.zeros( (l1+1,l2+1), int )
    arrow = np.zeros( (l1+1,l2+1), int )
    # create first row and first column
    scormat[0] = np.arange(l2+1)* gap
    scormat[:,0] = np.arange( l1+1)* gap
    arrow[0] = np.ones(l2+1)
    for i in range( 1, l1+1 ):
        for j in range( 1, l2+1 ):
            f = np.zeros( 3 )
            f[0] = scormat[i-1,j] + gap
            f[1] = scormat[i,j-1] + gap
            n1 = abet.index( seq1[i-1] )
            n2 = abet.index( seq2[j-1] )
            f[2] = scormat[i-1,j-1] + mat[n1,n2]
            scormat[i,j] = f.max()
            arrow[i,j] = f.argmax()
    return scormat, arrow

def Backtrace( arrow, seq1, seq2 ):
    st1, st2 = '',''
    v,h = arrow.shape
    ok = 1
    v-=1
    h-=1
    while ok:
        if arrow[v,h] == 0:
            st1 += seq1[v-1]
            st2 += '-'
            v -= 1
        elif arrow[v,h] == 1:
            st1 += '-'
            st2 += seq2[h-1]
            h -= 1
        elif arrow[v,h] == 2:
            st1 += seq1[v-1]
            st2 += seq2[h-1]
            v -= 1
            h -= 1
        if v==0 and h==0:
            ok = 0
    # reverse the strings
    st1 = st1[::-1]
    st2 = st2[::-1]
    return st1, st2

# create a fast substitution matrix
def FastSubMatrix( mat, abet, seq1, seq2 ):
    l1, l2 = len( seq1), len(seq2)
    LM = len( mat )
    # convert the sequences to numbers
    si1 = np.zeros( l1, int )
    si2 = np.zeros( l2, int )
    for i in range( l1 ):
        si1[i] = abet.index( seq1[i] )
    for i in range( l2 ):
        si2[i] = abet.index( seq2[i] )
    # create a matrix that contains the blosum values
    br = mat.ravel()
    bvndx = np.add.outer( LM*si1, si2 )
    bv = np.take( br, ravel( bvndx ))
    bv = bv.reshape( (l1,l2) )
    bvt = np.zeros( (l1+1,l2+1), int )
    bvt[1:,1:] = bv + 0
    bv = bvt.ravel()
    return bv


def FastSubValues( mat, abet, seq1, seq2 ):
    l1, l2 = len( seq1 ), len(seq2)
    subvals = np.zeros( (l1+1,l2+1), int )
    # convert the sequences to numbers
    si1 = np.zeros( l1, int )
    si2 = np.zeros( l2, int )
    for i in range( l1 ):
        si1[i] = abet.index( seq1[i] )
    for i in range( l2 ):
        si2[i] = abet.index( seq2[i] )
    for i in range( 1, l1+1 ):
        subvals[i,1:] = mat[ [si1[i-1]]*l2,si2]
    return subvals
    
def CreateIlist( l1, l2 ):
    ilist = []
    for i in range( l1 + l2 -1 ):
        st1 = min( i+1, l1 )
        sp1 = max( 1, i-l2+2 )
        st2 = max( 1, i-l1+2 )
        sp2 = min( i+1, l2 )
        ilist.append( (np.arange(st1,sp1-1,-1),np.arange(st2,sp2+1)))
    return ilist

def FastNW( subvals, seq1, seq2, gap=-8 ) :
    l1, l2 = len( seq1), len(seq2)
    scormat = np.zeros( (l1+1,l2+1), int )
    arrow = np.zeros( (l1+1,l2+1), int )
    # create first row and first column
    scormat[0] = np.arange(l2+1)* gap
    scormat[:,0] = np.arange( l1+1)* gap
    arrow[0] = np.ones( l2+1 )
    # compute the ilist
    ilist = CreateIlist( l1, l2 )
    # fill in the matrix
    for i in ilist:
        LI = len( i[0] )
        f = np.zeros( (3,LI), float )
        x,y = i[0]-1, i[1]+0
        f[0] = scormat[x,y] + gap
        x,y = i[0]+0,i[1]-1
        f[1] = scormat[x,y] + gap
        x,y = i[0]-1,i[1]-1
        f[2] = scormat[x,y] + subvals[i]
        f += 0.1 * np.sign(f) * np.random.ranf( f.shape ) 
        mx = (f.max(0)).astype(int)   # best values
        maxpos = f.argmax( 0 )
        scormat[i] = mx + 0
        arrow[i] = maxpos + 0
    return scormat, arrow
        
def FastSW( subvals, seq1, seq2, gap=-8 ) :
    l1, l2 = len( seq1), len(seq2)
    scormat = np.zeros( (l1+1,l2+1), int )
    arrow = np.zeros( (l1+1,l2+1), int )
    # create first row and first column
    arrow[0] = np.ones( l2+1 )
    # compute the ilist
    ilist = CreateIlist( l1, l2 )
    # fill in the matrix
    for i in ilist:
        LI = len( i[0] )
        f = np.zeros( (4,LI), float )
        x,y = i[0]-1, i[1]+0
        f[0] = scormat[x,y] + gap
        x,y = i[0]+0,i[1]-1
        f[1] = scormat[x,y] + gap
        x,y = i[0]-1,i[1]-1
        f[2] = scormat[x,y] + subvals[i]
        f += 0.1 * np.sign(f) * np.random.ranf( f.shape ) 
        mx = (f.max(0)).astype(int)   # best values
        maxpos = f.argmax( 0 )
        scormat[i] = mx + 0
        arrow[i] = maxpos + 0
    return scormat, arrow


def SWBacktrace( scormat, arrow, seq1, seq2 ):
    st1, st2 = '',''
    v,h = arrow.shape
    ok = 1
    #arrow[0] = ones( h )
    #mx = amax(amax(scormat ))
    v,h = divmod( scormat.argmax(), len(seq2)+1 )
    while ok:
        #print v,h,arrow[v,h], scormat[v,h],'    ',
        if arrow[v,h] == 0:
            st1 += seq1[v-1]
            st2 += '-'
            v -= 1
        elif arrow[v,h] == 1:
            st1 += '-'
            st2 += seq2[h-1]
            h -= 1
        elif arrow[v,h] == 2:
            st1 += seq1[v-1]
            st2 += seq2[h-1]
            v -= 1
            h -= 1
        elif arrow[v,h] == 3:
            ok = 0
        if (v==0 and h==0) or scormat[v,h]==0:
            ok = 0
    # reverse the strings
    st1 = st1[::-1]
    st2 = st2[::-1]
    return st1, st2
