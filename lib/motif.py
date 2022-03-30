import numpy as np
from random import randrange


def kmer_prob(kmer, prof):
    """given a kmer and a profile matrix dict calculate
     the probability of kmer given profile"""
    pr = 1
    for i in range(len(kmer)):
        pr *= prof[kmer[i]][i]
    return pr


def get_motifs(prof, dna_strings):
    """given a collection of dna strings of the same size and a
     profile matrix dict this function returns the collection of
      kmers formed by the profile-most probable kmers in each
      dna string"""
    k = prof['A'].shape[0]
    string_size = len(dna_strings[0])
    motifs = []
    for dna_string in dna_strings:  # for each string
        max_prob = -1
        best_kmer = None
        for i in range(string_size - k):  # for each kmer
            curr_kmer = dna_string[i:i + k]
            curr_prob = kmer_prob(curr_kmer, prof)
            if max_prob < curr_prob:  # update most probable
                best_kmer = curr_kmer
                max_prob = curr_prob
        motifs.append(best_kmer)
    return motifs


def profile(motifs):
    """ given a list of dna-strings of the same length calculate the
     profile matrix as a dictionary"""
    n_strings = len(motifs)
    string_size = len(motifs[0])
    # initialize occurence count to one (laplace rule)
    ret = {'A': np.ones(len(motifs[0])),
           'C': np.ones(len(motifs[0])),
           'G': np.ones(len(motifs[0])),
           'T': np.ones(len(motifs[0]))}

    # count occurences
    for i in range(string_size):  # for each position of motifs
        for j in range(n_strings):  # for each char in motif string
            ret[motifs[j][i]][i] += 1  # get count of the base

    # convert occurence counts to frequency
    ret['A'] /= len(motifs)+4
    ret['C'] /= len(motifs)+4
    ret['T'] /= len(motifs)+4
    ret['G'] /= len(motifs)+4
    return ret


def profile_most_probable(dna_string, prof, k):
    """ get the most probable k-mer in a dna string given a profile """
    most_prob = None
    prob = -1
    for i in range(len(dna_string) - k):  # for each kmer
        curr_kmer = dna_string[i:i + k]
        curr_prob = kmer_prob(curr_kmer, prof)  # calculate kmer probabilities
        if curr_prob > prob:  # update most probable
            most_prob = curr_kmer
            prob = curr_prob
    return most_prob


def score(motifs):
    """given a list of dna string of the same length calculate the
    score"""
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


def randomized_motif_search(dna_strings, k, t):
    string_size = len(dna_strings)
    motifs = []
    for i in range(t):  # choose randomly first motifs
        a = randrange(0, string_size - k)
        motifs.append(dna_strings[i][a:a + k])
    best_motifs = motifs
    while True:
        prof = profile(motifs)  # update profile
        motifs = get_motifs(prof, dna_strings)  # update motifs
        if score(motifs) < score(best_motifs): # we can improve
            best_motifs = motifs
        else:  # best motif found
            return best_motifs


def greedy_motif_search(dna_strings, k, t):
    """greedy algorithm that, given a list of string of the same length
    and two integer k and t returns a list of motif"""
    string_size = len(dna_strings[0])
    # initialize best motif with first kmer of each string
    best_motif = [s[:k] for s in dna_strings]
    for i in range(string_size - k):
        # initialize motifs with the i-th one of the first string
        motifs = [dna_strings[0][i:i + k]]
        for j in range(1, t):
            prof = profile(motifs)  # calculate profile for motifs
            most_prob = profile_most_probable(dna_strings[j], prof, k)  # get the most probable kmer in j-th string
            motifs.append(most_prob)  # add new motif
        if score(motifs) < score(best_motif):  # if score is improved update best motifs
            best_motif = motifs
    return best_motif
