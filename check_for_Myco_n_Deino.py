import os,sys,random
from miscBowTie import BowTieReader, BowTieWriter
from Bio.pairwise2 import align

def check_seq(r, seq):
	matched = 0
	l = len(seq)
	for i in xrange(10):
		s = random.randint(0, l-20)
		matched += seq[s:s+20] in r['seq']
	return matched >= 5

def id_to_template(template, read):
	"""
	Align a possible MP or DR to the real sequence locally.
	If the alignment is bad, it must have a score (or id%) very low,
	which we will consider as a non-MP/DR.
	Returns: (local) id% to the template, local id% to the template IGNORING N,
	          comma-separated list of <0-based pos>:X>Y
	ex: 1:A>N means at pos 1 in template it was A but in read was N
	"""
	# NOTE: the alignment is 0-based start, 1-based end!!
	# this means aligned_temp[start:end] is the aligned region
	aligned_temp, aligned_read, score, start, end = \
			align.localms(template, read, 1, 0, -10**5, -10**5)[0]
	#print aligned_temp, len(aligned_temp)
	#print aligned_read
	#print score, start, end
	mismatch = []
	Ncount = 0

	for i in xrange(start, end):
		if aligned_temp[i]!=aligned_read[i]:
			if aligned_read[i] == 'N': Ncount += 1
			mismatch.append("{0}:{1}>{2}".format(i, aligned_temp[i], aligned_read[i]))
	#print mismatch
	return score*1./(end-start), (score-Ncount)*1./(end-start), ",".join(mismatch)

# this is MP and DR, just V3, primer removed
myco  = 'AGGGAATTTTTCACAATGAGCGAAAGCTTGATGGAGCAATGCCGCGTGAACGATGAAGGTCTTTAAGATTGTAAAGTTCTTTTATTTGGGAAGAATGACTTTAGCAGGTAATGGCTAGAGTTTGACTGTACCATTTTGAATAAGTGACGACTAACTAT'
deino = 'TAGGAATCTTCCACAATGGGCGCAAGCCTGATGGAGCGACGCCGCGTGAGGGATGAAGGTTTTCGGATCGTAAACCTCTGAATCTGGGACGAAAGAGCCTTAGGGCAGATGACGGTACCAGAGTAATAGCACCGGCTAACTCC' 

input = sys.argv[1] # DSXXXXX.aligned.composite.gz.primer_good.gz
h1 = BowTieWriter(input+'.alien_MP') #matches to mycoplasma pneumoniae
h2 = BowTieWriter(input+'.alien_DR') #matches to deinococcus radiodurans
h3 = open(input+'.alien.summary', 'w')
h3.write("SEQID\tLEN\tMATCH\tLOWEST_QUAL\tSIM\tSIM_ignoreN\tMISMATCHES\n")

for r in BowTieReader(input, False):
	if check_seq(r, myco):
		id, idN, mism = id_to_template(myco, r['seq'])
		lq = min(ord(q)-33 for q in r['qual'])
		#print id, lq, r['ID']
		if id >= 0.9:
			h1.write(r)
			h3.write("{seq}\t{len}\tMP\t{lq}\t{sim}\t{sim2}\t{mism}\n".format(\
					seq=r['ID'],len=len(r['seq']),lq=lq,sim=id, sim2=idN,mism=mism)) 
	elif check_seq(r, deino):
		id, idN, mism = id_to_template(deino, r['seq'])
		lq = min(ord(q)-33 for q in r['qual'])
		#print id, lq, r['ID']
		if id >= 0.9:
			h2.write(r)
			h3.write("{seq}\t{len}\tDR\t{lq}\t{sim}\t{sim2}\t{mism}\n".format(\
					seq=r['ID'],len=len(r['seq']),lq=lq,sim=id, sim2=idN,mism=mism)) 

h1.close()
h2.close()
h3.close()
