def skew(genome_string):
    g_count = 0
    c_count = 0
    skew = []
    for base in genome_string:
        if base == 'G':
            g_count += 1
        elif base == 'C':
            c_count += 1
        skew.append(g_count - c_count)
    return skew


def reverse_complement(genome_string):
    rev_comp = []
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_genome = genome_string[::-1]
    for base in reverse_genome:
        rev_comp.append(comp[base])
    return "".join(rev_comp)


def frequency_table(Text, k):
    frequent_patterns = {}
    for i in range(len(Text) - k):
        pattern = Text[i:i + k]
        if pattern in frequent_patterns.keys():
            frequent_patterns[pattern] += 1
        else:
            frequent_patterns[pattern] = 1
    return frequent_patterns


def frequent_words(genome_string, k):
    frequent_patterns = []
    freq_map = frequency_table(genome_string, k)
    max_frequency = freq_map[max(freq_map, key=freq_map.get)]

    for pattern in freq_map.keys():
        if freq_map[pattern] == max_frequency:
            frequent_patterns.append(pattern)

    return frequent_patterns, max_frequency
