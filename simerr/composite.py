import os, sys
import random
import math
sys.path.append('/home/etseng/Dropbox/SchoolWork/GitCode/gbpML/data/evaluations/simerr/')

def simread(length):
	true = ''
	seq1 = ''
	seq2 = ''
	qual1 = ''
	qual2 = ''
	for i in xrange(length):
		phred1 = random.randint(1,41)
		err1 = 10**(-phred1/10.)
		qual1 += chr(phred1 + 33)
		phred2 = random.randint(1,41)
		err2 = 10**(-phred2/10.)
		qual2 += chr(phred2 + 33)
		p = random.choice(['A','T','C','G'])
		true += p
		if random.random() <= err1:
			# give an errorneous base
			seq1 += random.choice([x for x in ['A','T','C','G'] if x!=p])
		else:
			seq1 += p
		if random.random() <= err2:
			seq2 += random.choice([x for x in ['A','T','C','G'] if x!=p])
		else:
			seq2 += p

	return true,seq1,seq2,qual1,qual2

def compose(b1, b2, e1, e2, base_freq):
	"""
	b1, b2 -- base calls
	e1, e2 -- associated error probabilities
	base_freq --- base frequencies

	Returns: the composite base and the error probability (not yet in phred scores!)
	"""
	err2 = 10**(-2/10.) # the error probability for phred score 2
	if e1 >= 1.: # this is probably phred score, convert!
		#print >> sys.stderr, "smartly converting phred scores to probability"
		e1 = 10**(-e1/10.)
	if e2 >= 1.:
		#print >> sys.stderr, "smartly converting phred scores to probability"
		e2 = 10**(-e2/10.)
	assert sum(base_freq.itervalues()) == 1.
	if b1 == 'N' or b2 == 'N':
		if b1 == b2 == 'N':
			return 'N', err2
		elif b1 != 'N':
			return b1, e1
		else: # b2 != 'N'
			return b2, e2
	elif b1==b2:
		n = b1
		x = base_freq[n]
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

def compose2(r1, r2, base_freq):
	"""
	Takes the paired end records (see miscBowTie.BowTieReader)
	and returns the composite read

	Returns: composite read, composite qual (in phred), overlap_size
	"""
	log10 = lambda x: math.log(x, 10)
	N = len(r1['seq'])
	delta = r2['offset'] - r1['offset'] # the overlap is N-delta
	if N <= delta: # there is no overlap! put a bunch of ?s in the middle
		insert = delta-N
		seq = r1['seq'] + '?'*insert + r2['seq']
		qual = r1['qual'] + '#'*insert + r2['qual']
	else:
		seq = r1['seq'][:delta]
		qual = r1['qual'][:delta]
		for i in xrange(N-delta): 
			c, e = compose(r1['seq'][delta+i], r2['seq'][i], \
					ord(r1['qual'][delta+i])-33, ord(r2['qual'][i])-33, base_freq)
			seq += c
			qual += chr(int(-10*log10(e)+33))
		seq += r2['seq'][N-delta:]
		qual += r2['qual'][N-delta:]
	return seq, qual, N-delta

if __name__ == "__main__":
	from miscBowTie import BowTieReader, BowTieWriter
	from cPickle import *
	from optparse import OptionParser

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
