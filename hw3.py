# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 23:06:29 2021
@author: dillonchewwx

Algorithms for DNA Sequencing - Homework 3

Adapt the editDistance function we saw in practical to answer questions 1 and 2 below. Your function should take arguments p (pattern), t (text) and should return the edit distance of the match between P and T with the fewest edits.
"""
def main():
    # Q1. What is the edit distance of the best match between pattern GCTGATCGATCGTACG and the excerpt of human chromosome 1?  (Don't consider reverse complements.)
    p = "GCTGATCGATCGTACG"
    t = readGenome("Data/chr1.GRCh38.excerpt.fasta")
    print(editDistance(p,t))
    # Q2. What is the edit distance of the best match between pattern GATTTACCAGATTGAG and the excerpt of human chromosome 1?  (Don't consider reverse complements.)
    p2 = "GATTTACCAGATTGAG"
    print(editDistance(p2, t))
    # Q3. Download and parse the read sequences from the provided Phi-X FASTQ file. Next, find all pairs of reads with an exact suffix/prefix match of length at least 30. Don't overlap a read with itself; if a read has a suffix/prefix match to itself, ignore that match. Picture the overlap graph corresponding to the overlaps just calculated.  How many edges are in the graph?  In other words, how many distinct pairs of reads overlap?
    phi_seq, phi_qual = readFastq("Data/ERR266411_1.for_asm.fastq")
    map = overlapMap(phi_seq, 30)
    print(len(map))
    # Q4. Picture the overlap graph corresponding to the overlaps computed for the previous question. How many nodes in this graph have at least one outgoing edge?  (In other words, how many reads have a suffix involved in an overlap?)
    outgoing_nodes = set()
    for pair in map:
        outgoing_nodes.add(pair[0])
    print(len(outgoing_nodes))
    
def overlap(a, b, min_length=3):
    """ 
    Return length of longest suffix of 'a' matching
    a prefix of 'b' that is at least 'min_length'
    characters long.  If no such overlap exists,
    return 0. 
    """
    start = 0  # start all the way at the left
    while True:
        start = a.find(prefix(b, min_length), start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def overlapMap(reads, k):
    # Let every k-mer in the dataset have an associated Python set object.
    # Use dictionary to associate each k-mer with its corresponding set
    kmers_set = {}
    for read in reads:
        kmers = getKmer(read, k)
        for kmer in kmers:
            if not kmer in kmers_set.keys():
                kmers_set[kmer] = set()
            kmers_set[kmer].add(read) # For every k-mer in a read, we add the read to the set object corresponding to that k-mer.
            
    # For each read a, find all overlaps involving a suffix of a by taking a's length-k suffix and find all reads containing that k-mer and call overlap for each.
    overlaps = []
    for a in reads:
        for b in kmers_set[suffix(a, k)]:
            if a != b and overlap(a, b, k): # ignore same reads and ensure at least 30 overlaps
                overlaps.append((a, b))
    return overlaps

def getKmer(read, k):
    kmer = set()
    for i in range(len(read)-k+1):
        kmer.add(read[i:i+k])
    return kmer
        
def prefix(sequence, length):
    return sequence[:length]

def suffix(sequence, length):
    return sequence[-length:]
    
def editDistance(p, t):
    """
    Rows are labeled with bases from p and columns with bases from t
    """
    # Create distance matrix
    D = []
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))
    # Initialize first row and column of matrix
    for i in range(len(p)+1):
        D[i][0] = i # Elements in the first column are set to 0, 1, 2...
    for i in range(len(t)+1):
        D[0][i] = 0 # Elements in the first row are set to 0
    # Fill in the rest of the matrix
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance of the closest match between p and t is the minimal value in the bottom row.
    return min(D[-1])

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

if __name__ == "__main__":
    main()