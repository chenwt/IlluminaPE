import os,sys,time
from miscBowTie import BowTieReader, BowTieWriter

def remove_high_expected_error(input, max_expected_error):
	"""
	Remove all reads where the expected error (sum of err probs from phred scores)
	exceeds <max_expected_error>
	"""
	os.system("rm {0}.experror_*".format(input))
	hgood = BowTieWriter(input + '.experror_good')
	hbad  = BowTieWriter(input + '.experror_bad')
	hlog = open(input + '.experror.log', 'w')
	start_t = time.time()
	good, bad = 0,0
	for r in BowTieReader(input, False):
		if sum(10**-((ord(x)-33)/10.) for x in r['qual']) <= max_expected_error:
			hgood.write(r)
			good += 1
		else:
			hbad.write(r)
			bad += 1
	hlog.write("Expected error filtering took {0} sec.\n".format(time.time()-start_t))
	hlog.write("Max allowed expected error: {0}\n".format(max_expected_error))
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

	parser = argparse.ArgumentParser(description='Remove reads with high expected errors')
	parser.add_argument("-i", dest="input", required=True, help="Input file (BowTie format, gzipped)")
	parser.add_argument("-e", dest="maxerr", required=True, type=int, help="Max expected error")
	args = parser.parse_args()

	remove_high_expected_error(args.input, args.maxerr)


