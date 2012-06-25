import os,sys
from miscBowTie import FastqReader, FastqWriter, BowTieReader

# this is E.coli V6 forward SANS primer
ecoli1 = 'TGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGA'
ecoli2 = 'TGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTT'

def diff(s, offset, ecoli):
	mm = [] 
	for i,x in enumerate(s):
		if ecoli[i]!=x: mm.append("{0}:{1}>{2}".format(i+offset, ecoli[i], x))
	return mm

def main(input_filename):
	f = open(input_filename + '.summary', 'w')
	f.write('ID\tForR\tPRIMERlen\tMATCH\tSEQ_sansPrimer\tQUAL_sansPrimer\tMISMATCHES_globalPos\n')
	for r in BowTieReader(input_filename, False):
		mm1 = diff(r['seq'], int(r['offset']), ecoli1)
		mm2 = diff(r['seq'], int(r['offset']), ecoli2)
		mm = mm1 if len(mm1) < len(mm2) else mm2
		if len(mm) > 5:
			print >> sys.stderr, "MORE than 5 errors. Discard!!!!"
			continue
		f.write(r['ID'] + '\t')
		f.write('F\t' + str(r['offset']) + '\tECOLI\t')
		f.write(r['seq'] + '\t')
		f.write(r['qual'] + '\t')
		f.write(",".join(mm) + '\n')
	f.close()

if __name__ == "__main__":
	main('primered/run1_AGC.Fprimer_good.gz.experror_good.gz')
