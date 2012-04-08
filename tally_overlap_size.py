import os, sys
sys.path.append('/home/etseng/Dropbox/SchoolWork/GitCode/gbpML/data/evaluations/simerr/')
from collections import defaultdict
from miscBowTie import BowTieReader

def tally_overlap_size_from_count_file(count_filename):
	"""
	count files have format:
	<cluster id> \t <tab-delimited seq IDs>
	seq IDs are format: HWI-ST700693:182:D0MGFACXX:1:1102:20792:55431/1 COMPOSED/23
	"""
	tally = defaultdict(lambda: 0)
	with open(count_filename) as f:
		for line in f:
			for b in line.strip().split('\t')[1:]:
				overlap = int(b.split()[-1][len('COMPOSED/'):])
				tally[overlap] += 1
	_min = min(tally)
	_max = max(tally)
	print("OVERLAP," + ",".join(map(str,xrange(_min, _max+1))))
	print("COUNT," + ",".join(str(tally[i]) for i in xrange(_min, _max+1)))
	return tally
		
def tally_overlap_size(composite_gz_filename):
	tally = defaultdict(lambda: 0)
	for r in BowTieReader(composite_gz_filename, False):
		overlap = int(r['ID'].split()[-1][len('COMPOSED/'):])
		tally[overlap] += 1
	_min = min(tally)
	_max = max(tally)
	print("OVERLAP," + ",".join(map(str,xrange(_min, _max+1))))
	print("COUNT," + ",".join(str(tally[i]) for i in xrange(_min, _max+1)))
	return tally

if __name__ == "__main__":
	tally_overlap_size(sys.argv[1])
	#tally_overlap_size_from_count_file(sys.argv[1])
