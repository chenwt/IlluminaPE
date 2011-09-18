import os, sys
sys.path.append('/home/etseng/Dropbox/SchoolWork/GitCode/gbpML/data/evaluations/simerr/')
from collections import defaultdict
from miscBowTie import BowTieReader

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
