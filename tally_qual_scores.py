import os, sys
sys.path.append('/home/etseng/Dropbox/SchoolWork/GitCode/gbpML/data/evaluations/simerr/')
from collections import defaultdict
from miscBowTie import BowTieReader

def tally_qual_scores(filename, output):
	quals_at = defaultdict(lambda: defaultdict(lambda: 0))
	with open(filename) as f:
			f.readline()
			f.readline()
			f.readline()
			line = f.readline().strip()
			#assert len(line) == 100
			for pos,q in enumerate(line):
				assert 41 >= ord(q) - 33 >= 0
				quals_at[pos][ord(q) - 33] += 1

	poses = quals_at.keys()
	poses.sort()

	with open(output, 'w') as f:
		f.write("POS," + ",".join([str(x) for x in xrange(42)]) + '\n')
		for pos in poses:
			f.write(str(pos) + ',' + ",".join([str(quals_at[pos][x]) for x in xrange(42)]) + '\n')

def tally_qual_scores_gz_bowtie(gz_filename, output, strand, reverse_pos):
	quals_at = defaultdict(lambda: defaultdict(lambda: 0))
	f = BowTieReader(gz_filename, False)
	max_phred_seen = 42
	count = 0
	for r in f:
		if strand is not None and r['strand']!=strand:
			continue
		count += 1
		r_qual = r['qual']
		if reverse_pos:
			r_qual = r_qual[::-1]
		for pos,q in enumerate(r_qual):
			assert ord(q) - 33 >= 0 # for combined reads it is possible to go above 41
			quals_at[pos][ord(q) - 33] += 1
			max_phred_seen = max(max_phred_seen, ord(q) - 33)

	# sanity check
	for pos in quals_at:
		sum(quals_at[pos]) == count

	poses = quals_at.keys()
	poses.sort()
	print >> sys.stderr, "{0} reads used".format(count)

	with open(output, 'w') as f:
		f.write("POS," + ",".join([str(x) for x in xrange(max_phred_seen)]) + '\n')
		for pos in poses:
			f.write(str(pos) + ',' + ",".join([str(quals_at[pos][x]) for x in xrange(max_phred_seen)]) + '\n')

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", dest="input", required=True, type=str, help="input aligned/composite filename (.gz)")
	parser.add_argument("-s", "--strand", dest="strand", default=None, choices=['+','-'], help="minus/plus strand, don't use for composite which should be None")
	parser.add_argument("-r", "--reverse", dest="reverse", default=False, action="store_true", help="Reverse strand (default False)")
	stuff = parser.parse_args()
	stuff.output = stuff.input
	if stuff.strand is not None:
		stuff.output += ".plus" if stuff.strand == '+' else ".minus"
	if stuff.reverse:
		stuff.output += ".reverse"
	stuff.output += ".phred_scores.txt"
	#print stuff
	print >> sys.stderr, "Output will be written to", stuff.output
	tally_qual_scores_gz_bowtie(stuff.input, stuff.output, stuff.strand, stuff.reverse)
