if __name__ == "__main__":
	import os, sys
	from collections import defaultdict
	sys.path.append(os.path.join(os.environ['EBB'],'IlluminaPE/simerr'))
	from miscBowTie import BowTieReader, BowTieWriter
	from cPickle import *
	from optparse import OptionParser
	from c_composite import compose2 # the cython version of composite
	#from composite import compose2 # this is the python version, comment out if using cython version

	parser = OptionParser()
	parser.add_option("--input", dest="input", help="Input bowtie aligned file (gzipped)")
	parser.add_option("--output", dest="output",  help="Output composite read filename")
	parser.add_option("--mutate", dest="mutate_freq", help="Mutate base frequencies pickle")
	
	options, args = parser.parse_args()
	
	if options.mutate_freq is None:
		mutate_freq = {'A':defaultdict(lambda:1/3.),\
					'T':defaultdict(lambda:1/3.),\
					'C':defaultdict(lambda:1/3.),\
					'G':defaultdict(lambda:1/3.)}
	else:
		mutate_freq = load(open(options.mutate_freq))
		for k, v in mutate_freq.iteritems():
			assert sum(v.itervalues()) == 1.

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
		# remove /1 from r1 if there is one, same with /2 from r2
		while r1['ID'].endswith('/1') and r2['ID'].endswith('/2'):
			r1['ID'] = r1['ID'][:-2]
			r2['ID'] = r2['ID'][:-2]
		seq, qual, overlap = compose2(r1, r2, base_freq, mutate_freq)
		writer.write_composite(r1, r2, seq, qual, overlap)
