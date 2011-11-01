import os,sys,gzip
sys.path.append('/shared/silo_researcher/Lampe_J/Gut_Bugs/UW_Larry/gbpML/data/evaluations/simerr/')
from miscBowTie import BowTieReader, gziplines

def filter_low_qual_seqs(gz_filename, phred_cutoff=20):
	"""
	Takes a BowTie-style gzipped file (ex: .aligned.composite.gz)
	and retain only seqs that have every base phred >= <cutoff>
	also uniquifies the sequences

	Outputs: .phred<cutoff>_passed.unique.fasta
	         .phred<cutoff>_passed.unique.count
	"""
	seen = {} # seq --> {'ids':list of IDs, 'index':associated cluster index}
	index = 0 # for tracking unique clusters
	bad = 0
	good = 0
	f = open(gz_filename + ".phred{0}_passed.unique.fasta".format(phred_cutoff), 'w')
	for r in BowTieReader(gz_filename, False):
		if all(ord(x)-33 >= phred_cutoff for x in r['qual']):
			good += 1
			if r['seq'] in seen:
				seen[r['seq']]['ids'].append(r['ID'])
			else: 
				seen[r['seq']] = {'index':index, 'ids':[r['ID']]}
				f.write('>' + str(index) + '\n')
				f.write(r['seq'] + '\n')
				index += 1
		else:
			bad += 1
	f.close()
	with open(gz_filename + ".phred{0}_passed.unique.count".format(phred_cutoff), 'w') as f:
		for d in seen.itervalues():
			f.write("{0}\t{1}\n".format(d['index'], "\t".join(d['ids'])))
	with open(gz_filename + ".phred{0}_passed.unique.log".format(phred_cutoff), 'w') as f:
		f.write("Running filter_low_qual_seqs\n")
		f.write("Input: " + gz_filename + '\n')
		f.write("PhredCutoff: " + str(phred_cutoff) + '\n')
		f.write("RemovedDueToLowQual: " + str(bad) + '\n')
		f.write("RemainingTotal: " + str(good) + '\n')
		f.write("RemainingUnique: " + str(len(seen)) + '\n')

def filter_low_count_low_qual_seqs(gz_filename):
	"""
	Reads in a bowtie gzip file (ex: DS19187.aligned.composite.gz)
	which must have a corressponding .unique.count.gz file 

	Outputs a text file which denotes seqs that should be REMOVED
	(ex: via Qiime's filter_seqs.py) because it has
	(1) only 1 count
	and
	(2) has 1 or more phred-2 bases ('#')

	Output is written to .unique.count1phred2.filter.txt
	"""
	failed = {}
	for r in BowTieReader(gz_filename, False):
	    if r['qual'].count('#') >= 1: failed[r['ID']] = 1
	print >> sys.stderr, "finished reading gz"

	with open(gz_filename + '.unique.count1phred2.filter.txt', 'w') as f:
		for line in gziplines(gz_filename + '.unique.count.gz'):
			a = line.strip().split('\t')
			if int(a[1]) == 1 and a[0] in failed:
				del failed[a[0]]
				f.write(a[0] + '\n')
#				print >> sys.stderr, "removing", a[0]

if __name__ == "__main__":
#	filter_low_count_low_qual_seqs(sys.argv[1])
	filter_low_qual_seqs(sys.argv[1])