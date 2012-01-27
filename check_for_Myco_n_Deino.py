import os,sys,random
from miscBowTie import BowTieReader, BowTieWriter

def check_seq(r, seq):
	matched = 0
	l = len(seq)
	for i in xrange(10):
		s = random.randint(0, l-20)
		matched += seq[s:s+20] in r['seq']
	return matched >= 5

myco  = 'AGGGAATTTTTCACAATGAGCGAAAGCTTGATGGAGCAATGCCGCGTGAACGATGAAGGTCTTTAAGATTGTAAAGTTCTTTTATTTGGGAAGAATGACTTTAGCAGGTAATGGCTAGAGTTTGACTGTACCATTTTGAATAAGTGACGACTAACTAT'
deino = 'TAGGAATCTTCCACAATGGGCGCAAGCCTGATGGAGCGACGCCGCGTGAGGGATGAAGGTTTTCGGATCGTAAACCTCTGAATCTGGGACGAAAGAGCCTTAGGGCAGATGACGGTACCAGAGTAATAGCACCGGCTAACTCC' 

input = sys.argv[1] # XXXX.aligned.composite.gz
h1 = BowTieWriter(input+'.alien_MP') #matches to mycoplasma pneumoniae
h2 = BowTieWriter(input+'.alien_DR') #matches to deinococcus radiodurans

for r in BowTieReader(input, False):
	if check_seq(r, myco):
		h1.write(r)
	elif check_seq(r, deino):
		h2.write(r)
h1.close()
h2.close()
