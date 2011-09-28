import os,sys,gzip
sys.path.append('/home/etseng/Dropbox/SchoolWork/GitCode/gbpML/data/evaluations/simerr/')
from miscBowTie import BowTieReader
from collections import defaultdict

def uniquify_bowtie_output_to_fastq(filename):
	"""
	<filename> is in BowTie format (either pre- or post-composite), gzipped
	read through it, ignore the qual scores, and simply output:
	(1) <filename>.unique.fasta.gz
	(2) <filename>.unique.count.gz
	"""
	seen_seq = {} # sequence --> list of ids
	for r in BowTieReader(filename, False):
		if r['seq'] in seen_seq:
			seen_seq[r['seq']].append(r['ID'])
		else:
			seen_seq[r['seq']] = [r['ID']]
	
	f1 = gzip.open(filename+'.unique.fasta.gz', 'w')
	f2 = gzip.open(filename+'.unique.count.gz', 'w')
	items = seen_seq.items()
	items.sort(key=lambda x: len(x[1]), reverse=True)
	for seq, ids in items:
		f1.write(">{0}\n{1}\n".format(ids[0], seq))
		f2.write("{0}\t{1}\t{2}\n".format(ids[0], len(ids), ",".join(ids)))
	f1.close()
	f2.close()

if __name__ == "__main__":
	uniquify_bowtie_output_to_fastq(sys.argv[1])
