import os, time
import itertools
from PEassembly.miscBowTie import BowTieReader, BowTieWriter

#def remove_high_expected_error(input_prefix, max_expected_error):
#	"""
#	Remove all reads where the expected error (sum of err probs from phred scores)
#	exceeds <max_expected_error>
#	"""
#	os.system("rm {0}.experror_*".format(input_prefix))
#	hgood = BowTieWriter(input_prefix + '.experror_good')
#	hbad  = BowTieWriter(input_prefix + '.experror_bad')
#	hlog = open(input_prefix + '.experror.log', 'w')
#	start_t = time.time()
#	good, bad = 0,0
#	for r in BowTieReader(input_prefix, False):
#		if sum(10**-((ord(x)-33)/10.) for x in r['qual']) <= max_expected_error:
#			hgood.write(r)
#			good += 1
#		else:
#			hbad.write(r)
#			bad += 1
#	hlog.write("Expected error filtering took {0} sec.\n".format(time.time()-start_t))
#	hlog.write("Max allowed expected error: {0}\n".format(max_expected_error))
#	hlog.write("# of original reads: {0}\n".format(good+bad))
#	hlog.write("# of reads removed: {0} ({1:.2f})\n".format(bad,bad*1./(good+bad)))
#	hlog.write("# of reads remaining: {0} ({1:.2f})\n".format(good,good*1./(good+bad)))
#
#	hgood.close()
#	hbad.close()
#	hlog.close()
#	os.system("gzip " + hgood.f.name)
#	os.system("gzip " + hbad.f.name)
	
def remove_high_expected_error_PE(file1, file2, max_expected_error):
	"""
	Remove all reads where the expected error (sum of err probs from phred scores)
	exceeds <max_expected_error>
	"""
	assert os.path.exists(file1) and os.path.exists(file2)
	os.system("rm {0}.experror_*".format(file1))
	os.system("rm {0}.experror_*".format(file2))
	hgood1 = BowTieWriter(file1 + '.experror_good')
	hgood2 = BowTieWriter(file2 + '.experror_good')
	hbad1  = BowTieWriter(file1 + '.experror_bad')
	hbad2  = BowTieWriter(file2 + '.experror_bad')
	hlog = open(file1 + '.experror.log', 'w')
	start_t = time.time()
	good, bad = 0,0
	for r1, r2 in itertools.izip(BowTieReader(file1, False), BowTieReader(file2, False)):
		if sum(10**-((ord(x)-33)/10.) for x in r1['qual']) <= max_expected_error and \
		sum(10**-((ord(x)-33)/10.) for x in r2['qual']) <= max_expected_error:
			hgood1.write(r1)
			hgood2.write(r2)
			good += 1
		else:
			hbad1.write(r1)
			hbad2.write(r2)
			bad += 1
	hlog.write("Expected error filtering took {0} sec.\n".format(time.time()-start_t))
	hlog.write("Max allowed expected error: {0}\n".format(max_expected_error))
	hlog.write("# of original reads: {0}\n".format(good+bad))
	hlog.write("# of reads removed: {0} ({1:.2f})\n".format(bad,bad*1./(good+bad)))
	hlog.write("# of reads remaining: {0} ({1:.2f})\n".format(good,good*1./(good+bad)))

	hgood1.close()
	hgood2.close()
	hbad1.close()
	hbad2.close()
	hlog.close()
	os.system("gzip " + hgood1.f.name)
	os.system("gzip " + hgood2.f.name)
	os.system("gzip " + hbad1.f.name)
	os.system("gzip " + hbad2.f.name)

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(description='Remove reads with high expected errors')
	parser.add_argument("-1", dest="file1", required=True, help="Input file 1 (BowTie format, gzipped)")
	parser.add_argument("-2", dest="file2", required=True, help="Input file 2 (BowTie format, gzipped)")
	parser.add_argument("-e", dest="maxerr", required=True, type=int, help="Max expected error")
	args = parser.parse_args()

	remove_high_expected_error_PE(args.file1, args.file2, args.maxerr)


