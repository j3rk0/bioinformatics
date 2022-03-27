import sys

def get_hamming_distance(pattern, pattern2) -> int:
    distance = 0
    for i in range(len(pattern)):
        if pattern[i] != pattern2[i]:
            distance += 1

    return distance


def minimum_hamming_distance(pattern, text) -> int:
    k = len(pattern)
    min_distance = sys.maxsize
    for i in range(len(text) - k):
        distance = get_hamming_distance(pattern, text[i:i+k])
        if min_distance >= distance:
            min_distance = distance

    return min_distance


def calculate_distance(pattern, dna) -> int:
    distance = 0
    for dna_string in dna:
        distance += minimum_hamming_distance(pattern, dna_string)

    return distance


def generate_all_possible_k_mers(alphabet, prefix, k, alphabet_length, result_set):
    if k == 0:
        return result_set.add(prefix)

    for i in range(alphabet_length):
        new_prefix = prefix + alphabet[i]
        generate_all_possible_k_mers(alphabet, new_prefix, k - 1, alphabet_length, result_set)


# Input: A collection of strings Dna and an integer k.
# Output: A k-mer Pattern minimizing d(Pattern, Dna) among all k-mers Pattern.
def median_string(dna, k):
    min_distance = sys.maxsize
    median = ""

    all_possible_k_mers = set()
    generate_all_possible_k_mers(["G", "A", "C", "T"], "", k, 4, all_possible_k_mers)

    for k_pattern in all_possible_k_mers:
        distance = calculate_distance(k_pattern, dna)
        if min_distance >= distance:
            min_distance = distance
            median = k_pattern

    return median

if __name__ == '__main__':
    dataset = ["AAATTGACGCAT", "GACGACCACGTT", "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTACGGGACAG"]
    print(median_string(dataset, 3))

