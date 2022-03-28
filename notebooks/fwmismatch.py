# d number of mishmatches
def frequent_words_with_mismatch(genome_string, k, d):
    frequentPatterns = {}
    close = []
    frequencyArray = []
    for i in range(0,(pow(4,k) - 1)):
        close[i] = 0
        frequencyArray = 0
    for i in range(0,(len(genome_string) - k)):
        neighborhood = neighbors(genome_string(i,k), d)
        for pattern in neighborhood:
            index = pattern_to_number(pattern)
            close[index] = 1
    for i in range(0,(pow(4,k) - 1)):
        if close[i] == 1:
            pattern = number_to_pattern(i,k)
            frequencyArray[i] = approximate_pattern_count(genome_string, pattern, d)
        maxCount = max(frequencyArray)
        for i in range(0,(pow(4,k) - 1)):
            if frequencyArray[i] == maxCount:
                pattern = number_to_pattern(i, k)
                frequentPatterns[i] = pattern #TODO: controllare indice
        return frequentPatterns

def neighbors(pattern, d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}
    neighborhood = []
    suffixNeighbors = neighbors(suffix(pattern), d)
    for text in suffixNeighbors:
        if hamming_distance(suffix(pattern), text) < d:
            for x_nucleotide in pattern:
                concatenation = text.join(x_nucleotide)
                neighborhood.append(concatenation)
        else:
            firstSymbol = first_symbol(pattern)
            neighborhood.append(text.join(firstSymbol))
    return neighborhood

#returns hamming distance
def hamming_distance(string1, string2):
    # Return the Hamming distance between equal-length sequences
    if len(string1) != len(string2):
        raise ValueError("Undefined for sequences of unequal length")
    else:
        return sum(ch1 != ch2 for ch1, ch2 in zip(string1, string2))

#removing the first symbol of pattern we will obtain a (k-1)-mer
def suffix(pattern):
    return pattern[1:]

#removing the last symbol of pattern we will obtain a (k-1)-mer
def prefix(pattern):
    return pattern[:-1]

#returns the first symbol of the pattern
def first_symbol(pattern):
    return pattern[0]

#returns the last symbol of the pattern
def last_symbol(pattern):
    string_len = len(pattern)
    return pattern[string_len - 1]

#TODO: capire che fa
def pattern_to_number(pattern):
    if not pattern:
        return 0
    symbol = last_symbol(pattern)
    pref = prefix(pattern)
    return 4*(pattern_to_number(pref)+symbol_to_number(symbol))

#transforms the symbols A,C,G,T to the respective integers 0,1,2,3
def symbol_to_number(string):
    string = string.replace("A", "0")
    string = string.replace("C", "1")
    string = string.replace("G", "2")
    string = string.replace("T", "3")
    return int(string)

#transforms the integers 0,1,2,3 to the respective symbols A,C,G,T
def number_to_symbol(index):
    string = str(index)
    string = string.replace("0", "A")
    string = string.replace("1", "C")
    string = string.replace("2", "G")
    string = string.replace("3", "T")
    return string

def number_to_pattern(index, k):
    if k == 1:
        return number_to_symbol(index)
    prefixIndex = quotient(index, 4)
    r = remainder(index, 4)
    symbol = number_to_symbol(r)
    prefixPattern = number_to_pattern(prefixIndex, k-1)
    return prefixPattern + symbol

def quotient(n, m):
    return n // m

def remainder(n, m):
    return n % m

def approximate_pattern_count(text, pattern,d):
    count = 0
    for i in range(0, (len(text) - len(pattern))):
        pattern1 = text(i, len(pattern))
        if hamming_distance(pattern, pattern1) <= d:
            count += 1
    return count