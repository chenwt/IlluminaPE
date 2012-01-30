import os,sys,time
from miscBowTie import BowTieReader, BowTieWriter

def mismatch_exceeded(s1, s2, max_mm_allowed):
	mm = 0
	for i in xrange(min(len(s1),len(s2))):
		if s1[i]!=s2[i]:
			mm += 1
			if mm > max_mm_allowed:
				return True
	return False

def match_primer_len(seq, primer, max_mm_allowed, max_offset_allowed, is_reverse):
	l = len(primer)
	if not is_reverse: # forward primer
		for i in xrange(max_offset_allowed):
			if not mismatch_exceeded(primer[i:], seq, max_mm_allowed): # found a match!
				return l-i
	else:
		for i in xrange(max_offset_allowed):
			if not mismatch_exceeded(primer[:l-i], seq[-(l-i):], max_mm_allowed):
				return l-i
	return 0

def remove_primers(input, f_primer, r_primer, max_mm_allowed=2, max_offset_allowed=10):
	os.system("rm {0}.primer_*".format(input))
	hgood = BowTieWriter(input + '.primer_good')
	hbad  = BowTieWriter(input + '.primer_bad')
	hlog = open(input + '.primer.log', 'w')
	start_t = time.time()
	good, bad = 0,0
	for r in BowTieReader(input, False):
		match_f_len = match_primer_len(r['seq'], f_primer, max_mm_allowed, max_offset_allowed, False)
		match_r_len = match_primer_len(r['seq'], r_primer, max_mm_allowed, max_offset_allowed, True)
		if match_f_len > 0 and match_r_len > 0:
			r['seq'] = r['seq'][match_f_len+1:-match_r_len]
			hgood.write(r)
			good += 1
		else:
			hbad.write(r)
			bad += 1
	hlog.write("Primer detection and removal took {0} sec.\n".format(time.time()-start_t))
	hlog.write("# of original reads: {0}\n".format(good+bad))
	hlog.write("# of reads removed: {0} ({1:.2f})\n".format(bad,bad*1./(good+bad)))
	hlog.write("# of reads remaining: {0} ({1:.2f})\n".format(good,good*1./(good+bad)))

	hgood.close()
	hbad.close()
	hlog.close()
	os.system("gzip " + hgood.f.name)
	os.system("gzip " + hbad.f.name)

if __name__ == "__main__":
	import argparse
	
	f_primer = 'ACTCCTACGGGAGGCAGCAGT'
	r_primer = 'GTGCCAGCAGCCGCGGTAATAC'

	parser = argparse.ArgumentParser(description='Remove complete or partial forward & reverse primers')
	parser.add_argument("-i", dest="input", required=True, help="Input file (BowTie format, gzipped)")

	args = parser.parse_args()

	remove_primers(args.input, f_primer, r_primer)


