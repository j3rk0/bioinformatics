import random
import numpy as np
from MotifFindingProblem import calculate_score, get_k_mers


# Input: a dictionary of k mers prob.
# Output: a biased random k mer
def biased_random_k_mer_selector(probabilities):
    probabilities_sum = sum(probabilities.values())
    k_mers = dict()
    for k_mer, prob in probabilities.items():
        k_mers[k_mer] = (prob / probabilities_sum)

    return max(k_mers)


# Input: a list of k mers and a profile matrix
# Output: conditional prob. of k mers
def get_conditional_probabilities(k_mers, profile) -> dict:
    probabilities = dict()
    for k_mer in k_mers:
        probability = 1
        for nucleotides in k_mer:
            probability *= profile.get(nucleotides)
        probabilities[k_mer] = probability
    return probabilities


# Input: a string, a profile matrix, an integer k
# Output: a biased random k mer from profile matrix
def profile_random_generator_k_mer(deleted_string, profile, k) -> str:
    k_mers = get_k_mers(deleted_string, k)
    probabilities = get_conditional_probabilities(k_mers, profile)
    selected_k_mer = biased_random_k_mer_selector(probabilities)
    return selected_k_mer


# Input: a list of motifs, a string k mer, an integer i
# Output: a list of motifs
def new_motifs_from(motifs, k_mer, index) -> list:
    new_motifs = list()
    for i in range(len(motifs)):
        if i != index:
            new_motifs.append(motifs[i])
        else:
            new_motifs.append(k_mer)
    return new_motifs


# Input: A collection of strings Dna, an integer k, an integer t and an integer n.
# Output: the best motifs and the score
def gibbs_sampler(dna, k, t, n):
    best_motifs = random_motifs_selector(dna, k)
    for j in range(n):
        i = random.randint(0, t - 1)
        profile = get_profile_except_i_th_string(best_motifs, i)
        k_mer = profile_random_generator_k_mer(dna[i], profile, k)
        new_motifs = new_motifs_from(best_motifs, k_mer, i)
        if calculate_score(best_motifs) < calculate_score(new_motifs):
            best_motifs = new_motifs

    return best_motifs, calculate_score(best_motifs)


# Input: list of dna-strings of the same length and an integer index
# Output: profile matrix as a dictionary
def get_profile_except_i_th_string(motifs, index) -> dict:
    n_strings = len(motifs)
    string_size = len(motifs[0])
    # initialize occurrence count to one (laplace rule)
    # change to np.zeros if you don't wont laplace rule
    ret = {'A': np.ones(len(motifs[0])),
           'C': np.ones(len(motifs[0])),
           'G': np.ones(len(motifs[0])),
           'T': np.ones(len(motifs[0]))}

    # count occurences
    for i in range(string_size):  # for each position of motifs
        if i != (index - 1):
            for j in range(n_strings):  # for each char in motif string
                ret[motifs[j][i]][i] += 1  # get count of the base

    # convert occurence counts to frequency
    # remove +4 if you don't wont laplace
    ret['A'] /= (len(motifs) - 1) + 4
    ret['C'] /= (len(motifs) - 1) + 4
    ret['T'] /= (len(motifs) - 1) + 4
    ret['G'] /= (len(motifs) - 1) + 4
    return ret


# Input: a list of strings dna, an integer k
# Output: a list of motifs
def random_motifs_selector(dna, k) -> list:
    motifs = list()
    for dna_string in dna:
        dna_string_length = len(dna_string)
        start_index = random.randint(0, dna_string_length - k)
        motifs.append(dna_string[start_index:start_index + k])
    return motifs


def test_gibbs_sampler():
    dna = ["TTACCTTAAC", "GATGTCTGTC", "CCGGCGTTAG", "CACTAACGAG", "CGTCAGAGGT"]
    k = 4
    iteration_number = 10
    print(gibbs_sampler(dna, k, len(dna), iteration_number))


if __name__ == '__main__':
    test_gibbs_sampler()
