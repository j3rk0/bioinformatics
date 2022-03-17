from lib.motif import *

f = open("data/genoma_colera.txt", "r")
colera = f.read()
f = open("data/reverse_genoma_colera.txt", "r")
rc_colera = f.read()
# %%
# let's try reverse complement
x = reverse_complement(colera)
print(x == rc_colera)
# %%
# let's plot skew
import numpy as np
import matplotlib.pyplot as plt

sk = np.array(skew(colera))
plt.scatter(x=np.array(range(sk.shape[0])), y=sk)
plt.show()

#%%

from lib.motif import *

string = "atcaatgatcaacgtaagcttctaagcatgatcaaggtgctcacacagtttatccacaac\
ctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgacca\
cggaaagatgatcaagagaggatgatttcttggccatatcgcaatgaatacttgtgactt\
gtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggatt\
acgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttagga\
tagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaat\
tgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaag\
atcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtt\
tccttaaccctctattttttacggaagaatgatcaagctgctgctcttgatcatcgtttc"

max_frequency, max_kmers = frequent_words(string,5)

print(f" max_frequency:{max_frequency} kmers:{max_kmers}")

