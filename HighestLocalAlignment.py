#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 19:36:36 2019

@author: Michael linderman
"""


import sys 
from enum import IntEnum
import numpy as np
from Bio.SubsMat.MatrixInfo import pam250




def readFile (fileName1, fileName2):
    """
    Opening file, returning two sequences.
    """
    file1 = open(fileName1, 'r')
    file1 = file1.readlines()
    word1 = ''
    for i in range(1, len(file1)):
        word1 += file1[i].strip()

    file2 = open(fileName2, 'r')
    file2 = file2.readlines()
    word2 = ''
    for j in range(1, len(file2)):
        word2 += file2[j].strip()

    return word1, word2



class Back(IntEnum):
    """Enum of edges"""
    MAT=1
    VRT=2
    HRZ=3
    END=4  # Edge from beginning

def local_alignment(seqV, seqW, scoring, indel):
    """Compute local alignment with score
    Args:
        seq1, seq2: Strings of sequence to be aligned
        scoring: SubstitutionMatrix
        indel: Indel score (negative is penalty)
        
    Returns:
        Tuple of (score, aligned seqV, aligned seqW)"""
        
    n, m = len(seqV)+1, len(seqW)+1
    scores = np.zeros((n,m), dtype=np.int)
    back = np.full((n,m), Back.END, dtype=np.int)
    # Treating seq1 as the "reference", so vertical move (i axis) is "insertion"
    # and horizontal move (j-axis) is a "deletion"
    if indel > 0:
    # If indel is penalty then we take edge from 0,0
        for i in range(1, n):
            scores[i,0] = scores[i-1,0] + indel
            back[i,0] = Back.VRT
        for j in range(1, m):
            scores[0,j] = scores[0,j-1] + indel
            back[0,j] = Back.HRZ
    
    for i in range(1, n):
        for j in range(1, m):
            # Needs to be in same order as Back enum above
            try: 
                cost = scoring[seqV[i-1], seqW[j-1]]
            except:
                cost = scoring[seqW[j-1],seqV[i-1]]
            incoming = (scores[i-1,j-1] + cost,
                        scores[i-1,j] + indel,
                        scores[i,j-1] + indel,
                        0)
            index = np.argmax(incoming)
            back[i, j] = index + 1  # To account for enums starting at 1
            scores[i, j] = incoming[index]
    
    np.set_printoptions(threshold=np.inf)# Find maximum value in matrix (taxi edge to the end)
    i, j = np.unravel_index(np.argmax(scores), scores.shape)
    
    align_score, alignV, alignW = scores[i, j], "", ""
    
    while i > 0 or j > 0:
        if back[i, j] == Back.MAT:
            alignV = seqV[i-1] + alignV
            alignW = seqW[j-1] + alignW
            i -= 1
            j -= 1
        elif back[i, j] == Back.VRT:
            alignV = seqV[i-1] + alignV
            alignW = "-" + alignW
            i -= 1
        
        elif back[i, j] == Back.HRZ:
            alignV = "-" + alignV
            alignW = seqW[j-1] + alignW
            j -= 1
        elif back[i, j] == Back.END:
            break
        
    return align_score, alignV, alignW




if __name__ == "__main__":
    
    
    
    seqV, seqW = readFile(sys.argv[1], sys.argv[2])
    
    
    print(seqV)
    
    
    
    score, alignV, alignW = local_alignment(seqV, seqW, pam250, -5)
    print(score)
    print()
    print()
    print(alignV)
    print(len(seqV))
    print(len(alignV))
    print()
    print()
    print(alignW)
    print(len(seqW))
    print(len(alignW))
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            