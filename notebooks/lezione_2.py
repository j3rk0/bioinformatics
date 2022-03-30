from lib.motif import *

k = 3
t = 5
dna_strings = ["GGCGTTCAGGCA",
               "AAGAATCAGTCA",
               "CAAGGAGTTCGC",
               "CACGTCAATCAC",
               "CAATAATATTCG"]

morifs = greedy_motif_search(dna_strings, k, t)

# %%
from lib.motif import *

k = 3
t = 5
s = ["GGCGTTCAGGCA",
     "AAGAATCAGTCA",
     "CAAGGAGTTCGC",
     "CACGTCAATCAC",
     "CAATAATATTCG"]

motifs = greedy_motif_search(s, k=k, t=t)

print(motifs)

# %%
s = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
     "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
     "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
     "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
     "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

k = 8
t  =5

motifs = randomized_motif_search(s, k=k, t = 5)
