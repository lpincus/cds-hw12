# -*- coding: utf-8 -*-
# brutecompare.py
"""
Created on Thu Jun 29 16:14:40 2017

@author: jkinser
"""

import genbank as gb
import simplealign
import blosum

def DumpSequences( gbname ):
    data = gb.ReadFile( gbname )
    klocs = gb.FindKeywordLocs(data)
    N = len( klocs )
    genes = []
    for i in range( N ):
        g = gb.Translation( data, klocs[i])
        genes.append( g )
    return genes
    
def CompareAll( seq, genes, bmat, abet ):
    scores = []
    NG = len( genes )
    for i in range( NG ):
        print(i)
        a = simplealign.BruteForceSlide(bmat, abet, seq, genes[i])
        scores.append(a.max()/len(seq))
    return scores
    
# load all sequences from genome 1
# load all sequences from genome 2
# using brute force slide. compare one sequence to all in a genome. 
# return the best scoring match.
# repeat for all in genome 1.