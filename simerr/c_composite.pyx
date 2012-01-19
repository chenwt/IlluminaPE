import cython
import math
from libc.stdlib cimport free, malloc

cdef extern from "stdlib.h":
	int atoi(char* x)

cdef extern from "math.h":
	double log10(double x)
	double pow(double x, double y)
	double floor(double x)

cdef object compose(bytes b1, bytes b2, double e1, double e2, object base_freq):
	"""
	b1, b2 -- base calls
	e1, e2 -- associated error probabilities
	base_freq --- base frequencies

	Returns: the composite base (as Python bytes) 
	         and the error probability (not yet in phred scores!)
	"""
	cdef double err2 = 10**(-2/10.) # the error probability for phred score 2
	cdef double p, _sum, fy
	cdef bytes n, y

	assert e1 < 1.
	assert e2 < 1.
	assert sum(base_freq.itervalues()) == 1.

	if b1 == 'N' or b2 == 'N':
		if b1 == b2:
			n = b1
			p = err2
		elif b1 != 'N':
			n = b1
			p = e1
		else: # b2 != 'N'
			n = b2
			p = e2
		return n, p
	elif b1 == b2:
		n = b1
		x = base_freq[b1]
		p = (1-e1)*(1-e2)*x
		_sum = 0
		for y,fy in base_freq.iteritems():
			if y==n: _sum += fy*(1-e1)*(1-e2)
			else: _sum += fy*e1*e2/9.
		p /= _sum
	elif e1 < e2:
		n = b1
		x = base_freq[b1]
		p = (1-e1)*(e2/3.)*x
		_sum = 0
		for y,fy in base_freq.iteritems():
			if y==b1: _sum += fy*(1-e1)*e2/3.
			elif y==b2: _sum += fy*e1/3.*(1-e2)
			else: _sum += fy*e1*e2/9.
		p /= _sum
	else:
		n = b2
		p = (e1/3.)*(1-e2)*base_freq[b2]
		_sum = 0
		for y,fy in base_freq.iteritems():
			if y==b1: _sum += fy*(1-e1)*e2/3.
			elif y==b2: _sum += fy*e1/3.*(1-e2)
			else: _sum += fy*e1*e2/9.
		p /= _sum
	return n, 1-p

cdef object compose3(char* seq1, char* seq2, char* qual1, char* qual2, int offset, int N, object base_freq):
	"""
	Should be called by the Python callable compose2

	seq1, seq2 --- the paired read seqs
	qual1, qual2 --- the paired qual strings  ex: '##BAC' in ASCII-33
	offset --- the overlap is at seq1[offset:offset+N] and seq2[:N]
	N --- overlap length
	base_freq --- dict for base freq of 'A','T','C','G'

	Returns: composite read string, composite qual string
	"""
	cdef double q1, q2
	cdef int i
	cdef bytes n
	cdef double e
	out_qual = ''
	out_seq  = ''
	for i in xrange(N):
		q1 = pow(10, -(qual1[offset+i] - 33)/10.) # don't need atoi cuz a char is kept is int in cython
		q2 = pow(10, -(qual2[i] - 33)/10.) 
		n, e = compose(<bytes>seq1[offset+i], <bytes>seq2[i], q1, q2, base_freq)
		out_seq += n
		out_qual += <bytes><int>floor((-10*log10(e) + 33)) # NOTE: occassionally off-by-1 with the Python version because of very small numerical differences....
	return out_seq, out_qual

def compose2(r1, r2, base_freq):
	"""
	Takes the paired end records (see miscBowTie.BowTieReader)
	and returns the composite read

	Returns: composite read string, composite qual string, overlap_size
	"""
	cdef int N, insert, delta
	seq = ''
	qual = ''
	#print compose3(r1['seq'], r2['seq'], r1['qual'], r2['qual'], len(r1['seq']), base_freq)
	N = len(r1['seq'])
	delta = r2['offset'] - r1['offset']
	insert = N - delta
	if insert <= 0: # there is no overlap! put a bunch of Ns in the middle
		seq = r1['seq'] + 'N'*(-insert) + r2['seq']
		qual = r1['qual'] + '#'*(-insert) + r2['qual']
	else:
		seq = r1['seq'][:delta]
		qual = r1['qual'][:delta]
		i_seq, i_qual = compose3(r1['seq'], r2['seq'], r1['qual'], r2['qual'], \
				delta, insert, base_freq)
		seq += i_seq
		qual += i_qual
		seq += r2['seq'][insert:]
		qual += r2['qual'][insert:]
	return seq, qual, insert

#if __name__ == "__main__":
#	from miscBowTie import BowTieReader, BowTieWriter
#	from cPickle import *
#	from optparse import OptionParser
#
#	parser = OptionParser()
#	parser.add_option("--input", dest="input", help="Input bowtie aligned file (gzipped)")
#	parser.add_option("--output", dest="output", help="Output composite read filename")
#	
#	options, args = parser.parse_args()
#
#	reader = BowTieReader(options.input, is_paired=True)
#	print >> sys.stderr, "calculating base frequencies"
#	base_freq_pickle = options.input + ".base_freq.pickle"
#	if os.path.exists(base_freq_pickle):
#		base_freq = load(open(base_freq_pickle))
#	else:
#		base_freq = reader.get_base_frequency()
#		with open(options.input + ".base_freq.pickle", 'w') as f:
#			dump(base_freq, f)
#	print >> sys.stderr, "reading bowtie aligned file..."
#	writer = BowTieWriter(options.output)
#	for r1, r2 in reader:
#		seq, qual, overlap = compose2(r1, r2, base_freq)
#		writer.write_composite(r1, r2, seq, qual, overlap)
