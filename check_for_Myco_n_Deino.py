import sys,random
import itertools
from miscBowTie import BowTieReader, BowTieWriter
from Bio.Seq import Seq
#from Bio.pairwise2 import align

def diff_seq(s, offset, refseq):
	"""
	Returns: list of {globalpos}:{real}>{observed}
	"""
	mm = [] 
	for i,x in enumerate(s):
		if refseq[i]!=x: mm.append("{0}:{1}>{2}".format(i+offset, refseq[i], x))
	return mm

def check_seq(r, seq):
	"""
	Check for matching to MP or DR by randomly matching perfect k-mers
	Returns: True (is match) if 5 out of 10 k-mers match!
	"""
	k = 10
	matched = 0
	l = len(r['seq'])
	for i in xrange(10):
		s = random.randint(0, l-k)
		matched += seq[s:s+k] in r['seq']
	return matched >= 5

#def id_to_template(template, read):
#	"""
#	Align a possible MP or DR to the real sequence locally.
#	If the alignment is bad, it must have a score (or id%) very low,
#	which we will consider as a non-MP/DR.
#	Returns: (local) id% to the template, local id% to the template IGNORING N,
#	          comma-separated list of <0-based pos>:X>Y
#	ex: 1:A>N means at pos 1 in template it was A but in read was N
#	"""
#	# NOTE: the alignment is 0-based start, 1-based end!!
#	# this means aligned_temp[start:end] is the aligned region
#	aligned_temp, aligned_read, score, start, end = \
#			align.localms(template, read, 1, 0, -10**5, -10**5)[0]
#	#print aligned_temp, len(aligned_temp)
#	#print aligned_read
#	#print score, start, end
#	mismatch = []
#	Ncount = 0
#	for i in xrange(start, end):
#		if aligned_temp[i]!=aligned_read[i]:
#			if aligned_read[i] == 'N': Ncount += 1
#			mismatch.append("{0}:{1}>{2}".format(i, aligned_temp[i], aligned_read[i]))
#	#print mismatch
#	return score*1./(end-start), (score-Ncount)*1./(end-start), ",".join(mismatch)

# this is MP and DR, just V3, primer removed
myco  = 'AGGGAATTTTTCACAATGAGCGAAAGCTTGATGGAGCAATGCCGCGTGAACGATGAAGGTCTTTAAGATTGTAAAGTTCTTTTATTTGGGAAGAATGACTTTAGCAGGTAATGGCTAGAGTTTGACTGTACCATTTTGAATAAGTGACGACTAACTAT'
deino = 'TAGGAATCTTCCACAATGGGCGCAAGCCTGATGGAGCGACGCCGCGTGAGGGATGAAGGTTTTCGGATCGTAAACCTCTGAATCTGGGACGAAAGAGCCTTAGGGCAGATGACGGTACCAGAGTAATAGCACCGGCTAACTCC' 

myco_r = Seq(myco).reverse_complement().tostring()
deino_r = Seq(deino).reverse_complement().tostring()

input1 = sys.argv[1] # ex: DSXXXXX.aligned.composite.gz.primer_good.gz, BowTie format, gzipped
input2 = sys.argv[2]
sample = sys.argv[3] # ex: DSXXXXX

h3 = open(sample+'.alien.summary', 'w')
h3.write('Sample\tID\tForR\tPRIMERlen\tMATCH\tSEQ_sansPrimer\tQUAL_sansPrimer\tMISMATCHES_localPos\n')

for r1,r2 in itertools.izip(BowTieReader(input1, False), BowTieReader(input2, False)):
	match = None
	if check_seq(r1, myco) and check_seq(r2, myco_r):
		mm1 = diff_seq(r1['seq'], 0, myco)
		mm2 = diff_seq(r2['seq'], 0, myco_r)
		if len(mm1)*1./len(r1['seq']) < .1 and len(mm2)*1./len(r2['seq']) < .1:
			match = 'MP'
	elif check_seq(r1, deino) and check_seq(r2, deino_r):
		mm1 = diff_seq(r1['seq'], 0, deino)
		mm2 = diff_seq(r2['seq'], 0, deino_r)
		if len(mm1)*1./len(r1['seq']) < .1 and len(mm2)*1./len(r2['seq']) < .1:
			match = 'DR'
	if match is not None:			
		h3.write(sample + '\t')
		h3.write(r1['ID'] + '\t')
		h3.write('F' + '\t')
		h3.write(str(r1['offset']) + '\t')
		h3.write(match + '\t')
		h3.write(r1['seq'] + '\t')
		h3.write(r1['qual'] + '\t')
		h3.write(",".join(mm1) + '\n')
		
		h3.write(sample + '\t')
		h3.write(r2['ID'] + '\t')
		h3.write('R' + '\t')
		h3.write(str(r2['offset']) + '\t')
		h3.write(match + '\t')
		h3.write(r2['seq'] + '\t')
		h3.write(r2['qual'] + '\t')
		h3.write(",".join(mm2) + '\n')

h3.close()
