import numpy as np


# Input: list of dna strings of the same length
# Output: score of the motifs
def calculate_score(motifs):
    string_size = len(motifs[0])
    n_strings = len(motifs)
    ret = 0
    for i in range(string_size):  # for each position
        # initialize scores and frequency
        curr_score = {'A': 0, 'C': 0, 'T': 0, 'G': 0}
        max_freq = 0
        for j in range(n_strings):  # for each string
            curr_score[motifs[j][i]] += 1  # update frequency
            if curr_score[motifs[j][i]] > max_freq:  # update max frequency
                max_freq = curr_score[motifs[j][i]]
        ret += n_strings - max_freq  # update score
    return ret


def calculate_motifs(dna, k):
    motifs_number = len(dna)

    for i in range(motifs_number):
        pass


# Input: a dna string and an integer k
# Output: a list of k-mers
def get_k_mers(dna_string, k) -> []:
    k_mers = list()
    for i in range(len(dna_string) - k):
        k_mers.append(dna_string[i:i+k])

    return k_mers


# Input: A collection of strings Dna and an integer k.
# Output: A collection Motifs of k-mers, one from each string in Dna, minimizing
# SCORE(Motifs) among all possible choices of k-mers.
def brute_force_motifs_search(dna, k) -> []:
    k_motifs = np.array([[]])

    for dna_string in dna:
        k_motifs = np.append(k_motifs, get_k_mers(dna_string, k), 0)
