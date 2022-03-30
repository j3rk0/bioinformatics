import numpy as np


def skew(dna_string):
    """return an array with the skew for each position
    in dna string given as input"""
    curr_skew = 0
    skew = []
    for base in dna_string:
        if base == 'G':  # if the curr base is G increase skew
            curr_skew += 1
        elif base == 'C':  # if curr base is C decrease skew
            curr_skew -= 1
        skew.append(curr_skew)  # save the skew at current position
    return np.array(skew)


def reverse_complement(dna_string):
    """return reverse complement of dna_string in input"""
    rev_comp = []
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}  # map each base to it's complement
    for base in dna_string[::-1]:  # iterate trough reversed dna string
        rev_comp.append(comp[base])
    return "".join(rev_comp)


def frequent_words(dna_string, k):
    """return a tuple with most frequent kmers and their frequency"""
    frequency_table = {}  # dictionary of kmer frequency
    max_frequency = 0  # maximum frequency of a kmer
    frequent_patterns = []  # kmers with maximum frequency

    for i in range(len(dna_string) - k):
        pattern = dna_string[i:i + k]  # slide window over dna string

        # update frequency
        if pattern in frequency_table.keys():
            frequency_table[pattern] += 1
        else:
            frequency_table[pattern] = 1

        # update max
        if frequency_table[pattern] > max_frequency:
            max_frequency = frequency_table[pattern]
            frequent_patterns = [pattern]
        elif frequency_table[pattern] == max_frequency:
            frequent_patterns.append(pattern)

    return frequent_patterns, max_frequency


def compare_with_mismatch(s1, s2, d):
    for i in range(len(s1)):
        if not s1[i] == s2[i]:
            if d > 0:  # reduce the number of mismatches allowed
                d -= 1
            else:  # too many mismatch!!
                return False
    return True
