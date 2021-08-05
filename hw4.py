# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 02:28:28 2021
@author: dillonchewwx
"""
def main():
    # Q1. What is the length of the shortest common superstring of the following strings? "CCT", "CTT", "TGC", "TGG", "GAT", "ATT"
    reads = ["CCT", "CTT", "TGC", "TGG", "GAT", "ATT"]
    shortest_sub, shortest_all = scs_all(reads)
    print(len(shortest_sub))
    # Q2. How many different shortest common superstrings are there for the input strings given in the previous question?
    print(len(shortest_all))
    # Q3. How many As are there in the full, assembled genome?
    virus_seq, virus_qual = readFastq("Data/ads1_week4_reads.fq")    
    assembled = greedy_scs(virus_seq, 30)
    print(len(assembled))
    print(assembled.count("A"))
    # Q4. How many Ts are there in the full, assembled genome?
    print(assembled.count("T"))

def pick_maximal_overlap(reads, k):
    """
    Returns a pair of reads from the list with a maximal suffix/prefix overlap >= k. Returns overlap length 0 if there are no such overlaps 
    """
    read_a, read_b = None, None
    best_overlap_length = 0
    
    # make index 
    kmers_set = {}
    for read in reads:
        kmers = getKmer(read, k)
        for kmer in kmers:
            if not kmer in kmers_set.keys():
                kmers_set[kmer] = set()
            kmers_set[kmer].add(read)
    
    for a in reads:
        for b in kmers_set[suffix(a, k)]:
            if a != b:
                overlap_length = overlap(a, b, k)
                if overlap_length > best_overlap_length:
                    read_a, read_b, best_overlap_length = a, b, overlap_length
    return read_a, read_b, best_overlap_length

def greedy_scs(reads, k):
    read_a, read_b, overlap_length = pick_maximal_overlap(reads, k)
    while overlap_length > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[overlap_length:])
        read_a, read_b, overlap_length = pick_maximal_overlap(reads, k)
    return "".join(reads)

def getKmer(read, k):
    kmer = set()
    for i in range(len(read)-k+1):
        kmer.add(read[i:i+k])
    return kmer
        
def prefix(sequence, length):
    return sequence[:length]

def suffix(sequence, length):
    return sequence[-length:]

def overlap(a, b, min_length=3):
    """ 
    Return length of longest suffix of 'a' matching
    a prefix of 'b' that is at least 'min_length'
    characters long.  If no such overlap exists,
    return 0. 
    """
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's suffx in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def scs(ss):
    """ Returns shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    import itertools
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
    return shortest_sup  # return shortest

def scs_all(ss):
    """ Returns all shortest common superstring of given
        strings, which must be the same length """
    shortest_sup = None
    shortest_all = set()
    import itertools
    for ssperm in itertools.permutations(ss):
        sup = ssperm[0]  # superstring starts as first string
        for i in range(len(ss)-1):
            # overlap adjacent strings A and B in the permutation
            olen = overlap(ssperm[i], ssperm[i+1], min_length=1)
            # add non-overlapping portion of B to superstring
            sup += ssperm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup  # found shorter superstring
        elif len(sup) == len(shortest_sup):
            shortest_all.add(sup)
    return shortest_sup, shortest_all

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