""" 
Algorithms for DNA Sequencing - Homework 1

First, implement a version of the naive exact matching algorithm that is strand-aware. That is, instead of looking only for occurrences of P in T, additionally look for occurrences of the reverse complement of P in T. If P is "ACT", your function should find occurrences of both "ACT" and its reverse complement "AGT" in T.

If P and its reverse complement are identical (e.g. "AACGTT"), then a given match offset should be reported only once. So if your new function is called naive_with_rc, then the old naive function and your new naive_with_rc function should return the same results when P equals its reverse complement.
"""
def main():
    # Q1. How many times does AGGT or its reverse complement ACCT occur in the lambda virus genome?
    lambda_genome = readGenome("Data/lambda_virus.fa")
    print(len(naive_with_rc("AGGT", lambda_genome)))
    # Q2. How many times does TTAA or its reverse complement occur in the lambda virus genome?
    print(len(naive_with_rc("TTAA", lambda_genome)))
    # Q3. What is the offset of the leftmost occurrence of ACTAAGT or its reverse complement in the Lambda virus genome? 
    print(min(naive_with_rc("ACTAAGT", lambda_genome)))
    # Q4. What is the offset of the leftmost occurrence of AGTCGA or its reverse complement in the Lambda virus genome?
    print(min(naive_with_rc("AGTCGA", lambda_genome)))
    # Q5. How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 mismatches? 
    print(len(naive_2mm("TTCAAGCC", lambda_genome)))
    # Q6. What is the offset of the leftmost occurrence of AGGAGGTT in the Lambda virus genome when allowing up to 2 mismatches?
    print(min(naive_2mm("AGGAGGTT", lambda_genome)))
    # # Q7. Report which sequencing cycle has the problem.
    sequences, qualities = readFastq("Data/ERR037900_1.first1000.fastq")
    import numpy as np
    meanScores = np.array(qualityScores(qualities)).mean(axis=0)
    print(list(meanScores).index(min(meanScores)))

def qualityScores(qualities):
    output = []
    for qual in qualities:
        phredscore = []
        for phred in qual:
            q = phred33ToQ(phred)
            phredscore.append(q)
        output.append(phredscore)
    return output

def phred33ToQ(qual):
    return ord(qual) - 33    

def HammingDistance(p, q):
    return sum(a != b for a, b in zip(p,q))

def naive_2mm(p,t):
    occurrences = []
    # don't consider reverse complement here
    for i in range(len(t) - len(p) + 1): # loop over alignments
        if HammingDistance(p, t[i:i+len(p)]) <= 2: # compare forward strand
            occurrences.append(i)
    return occurrences

def naive_with_rc(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1): # loop over alignments
        if t[i:i+len(p)] == p: # compare forward strand
            occurrences.append(i)
        elif t[i:i+len(p)] == reverseComplement(p): # compare reverse strand
            occurrences.append(i)
    return occurrences

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

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
