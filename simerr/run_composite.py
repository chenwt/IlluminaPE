if __name__ == "__main__":
	import os, sys
	sys.path.append(os.path.join(os.environ['EBB'],'IlluminaPE/simerr'))
	from miscBowTie import BowTieReader, BowTieWriter
	from cPickle import *
	from optparse import OptionParser
	from c_composite import compose2 # the cython version of composite
	# from composite import compose2 # this is the python version, comment out if using cython version

	parser = OptionParser()
	parser.add_option("--input", dest="input", help="Input bowtie aligned file (gzipped)")
	parser.add_option("--output", dest="output", help="Output composite read filename")
	
	options, args = parser.parse_args()

	reader = BowTieReader(options.input, is_paired=True)
	print >> sys.stderr, "calculating base frequencies"
	base_freq_pickle = options.input + ".base_freq.pickle"
	if os.path.exists(base_freq_pickle):
		base_freq = load(open(base_freq_pickle))
	else:
		base_freq = reader.get_base_frequency()
		with open(options.input + ".base_freq.pickle", 'w') as f:
			dump(base_freq, f)
	print >> sys.stderr, "reading bowtie aligned file..."
	writer = BowTieWriter(options.output)
	for r1, r2 in reader:
		seq, qual, overlap = compose2(r1, r2, base_freq)
		writer.write_composite(r1, r2, seq, qual, overlap)
