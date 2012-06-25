import os,sys
from csv import DictReader
from Bio import Seq
from miscBowTie import BowTieReader
sys.path.append('/shared/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/EBB_Illumina/data.htseq.org/FHCRC_Hullar/FCD05FJ/IlluminaPE')
import remove_primers
find_primer=remove_primers.match_primer_len

"""
Because I did MP/DR checking on the composite-primer-good reads first, 
now I'm going back and picking just these MP/DR reads from .aligned.gz
where I can look at individual sequence errors in the forward/reverse read.

In this way, I'm ignoring potential MP/DR that have either 
(a) failed alignment OR (b) failed primer matching.

USAGE: python xxx.py <sample_name>

ex: python xxx.py DS21057 results in output file: DS21057.aligned.gz.alien_MPnDR.summary
with fields SEQID, Forward or Reverse, primer len, MP or DR, seq w/ primer removed
(and reverse-complemented if is reverse read), qual w/ primer removed,
and list of errors with {global position w/ primer}:{real base}>{observed base}
"""



f_primer = 'ACTCCTACGGGAGGCAGCAGT'
r_primer = 'GTGCCAGCAGCCGCGGTAATAC'
myco  = 'AGGGAATTTTTCACAATGAGCGAAAGCTTGATGGAGCAATGCCGCGTGAACGATGAAGGTCTTTAAGATTGTAAAGTTCTTTTATTTGGGAAGAATGACTTTAGCAGGTAATGGCTAGAGTTTGACTGTACCATTTTGAATAAGTGACGACTAACTAT'
deino = 'TAGGAATCTTCCACAATGGGCGCAAGCCTGATGGAGCGACGCCGCGTGAGGGATGAAGGTTTTCGGATCGTAAACCTCTGAATCTGGGACGAAAGAGCCTTAGGGCAGATGACGGTACCAGAGTAATAGCACCGGCTAACTCC'

myco_r = Seq.Seq(myco).reverse_complement().tostring()
deino_r = Seq.Seq(deino).reverse_complement().tostring()

sample = sys.argv[1]

SUMMARY    = sample + '.aligned.composite.gz.primer_good.gz.alien.summary2'
ALIGNED_GZ = sample + '.aligned.gz'
OUTPUT     = sample + '.aligned.gz.alien_MPnDR.summary'

def list_mm(real, obs, offset):
	mm = ''
	for i in xrange(len(obs)):
		if obs[i]!=real[i]: 
			mm += "{0}:{1}>{2},".format(i+offset,real[i],obs[i])
	return mm[:-1]

d = {} # SEQID --> {MATCH, MISMATCHES}

for x in DictReader(open(SUMMARY),delimiter='\t'):
	id = x['SEQID']
	id = id[:id.find('/')]
	d[id] = x

f = open(OUTPUT,'w')
f.write("ID\tForR\tPRIMERlen\tMATCH\tSEQ_sansPrimer\tQUAL_sansPrimer\tMISMATCHES_globalPos\n")
    
for x1,x2 in BowTieReader(ALIGNED_GZ, True):
	id = x1['ID'][:-2]
	if id in d:
		if d[id]['MATCH'] == 'MP':
			real,real_r = myco,myco_r
		else:
			real,real_r = deino,deino_r
		i = find_primer(x1['seq'], f_primer, 2, 10, False)
		f.write("{id}\tF\t{primer}\t{match}\t{seq}\t{qual}\t{mm}\n".format(\
				id=id, primer=i, match=d[id]['MATCH'], seq=x1['seq'][i:],\
				qual=x1['qual'][i:], mm=list_mm(real, x1['seq'][i:], i)))

		i = find_primer(x2['seq'], r_primer, 2, 10, True)
		x2_r = Seq.Seq(x2['seq'][:-i]).reverse_complement().tostring()
		q2_r = x2['qual'][:-i]
		q2_r = q2_r[::-1]
		f.write("{id}\tR\t{primer}\t{match}\t{seq}\t{qual}\t{mm}\n".format(\
				id=id, primer=i, match=d[id]['MATCH'], seq=x2_r,\
				qual=q2_r, mm=list_mm(real_r, x2_r, i)))

