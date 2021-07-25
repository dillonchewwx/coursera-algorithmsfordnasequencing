# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 11:45:58 2021
@author: dillonchewwx

Algorithms for DNA Sequencing - Homework 2

Implement versions of the naive exact matching and Boyer-Moore algorithms that additionally count and return (a) the number of character comparisons performed and (b) the number of alignments tried. Roughly speaking, these measure how much work the two different algorithms are doing.
"""

def main():
    # Q1. How many alignments does the naive exact matching algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1?  (Don't consider reverse complements.)
    p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"
    t = readGenome("Data/chr1.GRCh38.excerpt.fasta")
    occurrences, num_alignments, num_character_comparisons = naive_with_counts(p, t)
    print(num_alignments)
    # Q2. How many character comparisons does the naive exact matching algorithm try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1?  (Don't consider reverse complements.)
    print(num_character_comparisons)
    # Q3. How many alignments does Boyer-Moore try when matching the string GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG (derived from human Alu sequences) to the excerpt of human chromosome 1?  (Don't consider reverse complements.)
    from bm_preproc import BoyerMoore
    p_bm = BoyerMoore(p, "ACGT")
    bm_occurrences, bm_num_alignments, bm_num_character_comparisons = boyer_moore_with_counts(p, p_bm, t)
    print(bm_num_alignments)
    # Q4. Implement the pigeonhole principle using Index to find exact matches for the partitions. Assume P always has length 24, and that we are looking for approximate matches with up to 2 mismatches (substitutions). We will use an 8-mer index. Write a function that, given a length-24 pattern P and given an Index object built on 8-mers, finds all approximate occurrences of P within T with up to 2 mismatches. Insertions and deletions are not allowed. Don't consider any reverse complements. How many times does the string GGCGCGGTGGCTCACGCCTGTAAT, which is derived from a human Alu sequence, occur with up to 2 substitutions in the excerpt of human chromosome 1?  (Don't consider reverse complements here.)
    p_4 = "GGCGCGGTGGCTCACGCCTGTAAT"
    occurrences, index_hits = approximate_match_Index_v2(p_4, t, 2, 8)
    print(occurrences)
    # Q5. Using the instructions given in Question 4, how many total index hits are there when searching for occurrences of GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1?
    print(index_hits)
    # Q6. Write a function that, given a length-24 pattern P and given a SubseqIndex object built with k = 8 and ival = 3, finds all approximate occurrences of P within T with up to 2 mismatches. When using this function, how many total index hits are there when searching for GGCGCGGTGGCTCACGCCTGTAAT with up to 2 substitutions in the excerpt of human chromosome 1?  (Again, don't consider reverse complements.)
    occurrences_subseq, index_hits_subseq = approximate_match_SubseqIndex(p_4, t, 2, 8, 3)
    print(index_hits_subseq)
    
def approximate_match_SubseqIndex(p, t, n, k, ival):
    all_matches = set()
    from kmer_subseq_index import SubseqIndex
    index = SubseqIndex(t, k, ival)
    index_hits = 0
    for i in range(n+1): # query first three kmers with respect to p.
        matches = index.query(p[i:])
        for match in matches:
            index_hits += 1
            offset = match-i
            if offset < 0  or offset + len(p) > len(t):
                continue
            if HammingDistance(p, t[offset:offset+len(p)]) <= n:
                all_matches.add(offset)
    occurrences = list(all_matches)
    return occurrences, index_hits
    
def approximate_match_Index_v2(p, t, n, k):
    segment_length = int(round(len(p)/(n+1))) # divide p into n+1 segments where at least one of those segments must match perfectly against t.
    all_matches = set()
    from kmer_index import Index
    index = Index(t, k)
    index_hits = 0
    for i in range(n+1): # iterate across segments
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = index.query(p[start:end])    
        for match in matches: # iterate across matches
            index_hits += 1    
            offset = match-start
            if offset < 0  or offset + len(p) > len(t):
                continue
            if HammingDistance(p[:start], t[offset:offset+start]) + HammingDistance(p[end:], t[offset+end:offset+len(p)]) <= n:
                all_matches.add(offset)
    occurrences = list(all_matches)
    return occurrences, index_hits
            
def approximate_match_Index(p, t, n, k):
    segment_length = int(round(len(p)/(n+1))) # divide p into n+1 segments where at least one of those segments must match perfectly against t.
    all_matches = set()
    from kmer_index import Index
    index = Index(t, k)
    for i in range(n+1): # iterate across segments
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        matches = index.query(p[start:end])
        for match in matches:
            offset = match-start
            if offset < 0 or offset + len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if p[j] != t[offset+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end, len(p)):
                if p[j] != t[offset+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(offset)
    return list(all_matches)

def boyer_moore_with_counts(p, p_bm, t):
    """ 
    Do Boyer-Moore matching. p=pattern, t=text,
    p_bm = BoyerMoore object for p 
    """
    i = 0
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    while i < len(t) - len(p) + 1: # loop over alignments
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1): # loop over characters
            num_character_comparisons += 1
            if p[j] != t[i+j]: # compare characters
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
        num_alignments += 1
    return occurrences, num_alignments, num_character_comparisons

def naive_with_counts(p, t):
    """
    Does naive exact matching given p=pattern and t=text.
    Does not consider reverse complements.
    Returns occurrences, number of character comparisons, and number of alignments tried.
    """
    occurrences = []
    num_alignments = 0
    num_character_comparisons = 0
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            num_character_comparisons += 1
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
        num_alignments += 1
    return occurrences, num_alignments, num_character_comparisons

def HammingDistance(p, q):
    return sum(a != b for a, b in zip(p,q))

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome
    
if __name__ == "__main__":
    main()