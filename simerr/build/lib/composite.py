import random
from math import log

def simread(length):
	true = []
	seq1 = []
	seq2 = []
	qual1 = []
	qual2 = []
	for i in xrange(length):
		phred1 = random.randint(1,41)
		err1 = 10**(-phred1/10.)
		qual1.append(err1)
		phred2 = random.randint(1,41)
		err2 = 10**(-phred2/10.)
		qual2.append(err2)
		p = random.choice(['A','T','C','G'])
		true.append(p)
		if random.random() <= err1:
			# give an errorneous base
			seq1.append(random.choice([x for x in ['A','T','C','G'] if x!=p]))
		else:
			seq1.append(p)
		if random.random() <= err2:
			seq2.append(random.choice([x for x in ['A','T','C','G'] if x!=p]))
		else:
			seq2.append(p)

	return true,seq1,seq2,qual1,qual2

def compose(b1, b2, e1, e2, base_freq):
	"""
	b1, b2 -- base calls
	e1, e2 -- associated error probabilities
	base_freq --- base frequencies

	Returns: the composite base and the error probability (not yet in phred scores!)
	"""
	assert sum(base_freq.itervalues()) == 1.
	if b1==b2:
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

if __name__ == "__main__":
	true,seq1,seq2,qual1,qual2 = simread(10**6)

	c = [compose(seq1[i],seq2[i],qual1[i],qual2[i],{'A':.25,'T':.25,'C':.25,'G':.25}) for i in xrange(len(seq1))]
