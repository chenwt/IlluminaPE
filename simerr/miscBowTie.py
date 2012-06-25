import os, sys
from subprocess import Popen, PIPE
import cStringIO
io_method = cStringIO.StringIO


def gziplines(fname):
	f = Popen(['zcat', fname], stdout=PIPE)
	for line in f.stdout:
		yield line
		
class FastqWriter:
	def __init__(self, filename):
		assert not os.path.exists(filename)
		self.f = open(filename, 'w')
	
	def close(self):
		self.f.close()
		
	def write(self, r, id=None):
		self.f.write("@{id}\n{seq}\n+\n{qual}\n".format(\
			id=r['ID'] if id is None else id, seq=r['seq'], qual=r['qual']))
									
		
		
class FastqReader:
	def __init__(self, filename):
		self.filename = filename
		assert self.filename.endswith('.gz')
		self.fp = Popen(['zcat', filename], stdout=PIPE)
		self.f = io_method(self.fp.communicate()[0])

	def __iter__(self):
		return self

	def next(self):
		try:
			line = self.f.next()
		except StopIteration:
			raise StopIteration
		assert line.startswith('@')
		id = line.strip()[1:] # remove the '@'
		seq = self.f.next().strip() # sequence
		self.f.next()
		qual = self.f.next().strip() # qual
		return {'ID': id, 'seq': seq, 'qual': qual }

class FastqReaderPaired:
	def __init__(self, filename1, filename2):
		self.reader1 = FastqReader(filename1)
		self.reader2 = FastqReader(filename2)

	def __iter__(self):
		return self

	def next(self):
		r1 = self.reader1.next()
		r2 = self.reader2.next()
		return r1, r2


class BowTieReader:
	bowtie_struct = ['ID', 'strand', 'ref', 'offset', 'seq', 'qual', 'mismatch_num', 'mismatch_info']

	def __init__(self, filename, is_paired):
		"""
		Only works with .gz files for now!!
		"""
		self.filename = filename
		assert self.filename.endswith('.gz')
		self.fp = Popen(['zcat', filename], stdout=PIPE)
		self.f = io_method(self.fp.communicate()[0])
		self.is_paired = is_paired

	def get_base_frequency(self):
		"""
		reads through the whole file to get base freq, then rewind
		ignores ambiguous nucleotides
		"""
		self.fp = Popen(['zcat', self.filename], stdout=PIPE)
		self.f = io_method(self.fp.communicate()[0])
		base_freq = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
		while True:
			r = self.read()
			if r is None:
				break
			base_freq['A'] += r['seq'].count('A')
			base_freq['T'] += r['seq'].count('T') + r['seq'].count('U')
			base_freq['C'] += r['seq'].count('C')
			base_freq['G'] += r['seq'].count('G')
		self.fp = Popen(['zcat', self.filename], stdout=PIPE)
		self.f = io_method(self.fp.communicate()[0])
		_sum = sum(base_freq.values()) * 1.
		for k in base_freq:
			base_freq[k] /= _sum
		return base_freq

	def __iter__(self):
		return self

	def next(self):
		if self.is_paired: # paired end reads
			r1 = self.read()
			if r1 is None:
				raise StopIteration
			r2 = self.read()
			# sanity check that the IDs should have the same header!
			# ex: HWI-ST700693:106:D05FJACXX:1:1101:1273:1988 2:N:0:ATCACG/2 and
			# ex: HWI-ST700693:106:D05FJACXX:1:1101:1273:1988 1:N:0:ATCACG/1
			if len(r1['ID'].split(None)) == 2:
				assert r1['ID'].split()[0] == r2['ID'].split()[0]
			elif len(r1['ID'].split(None)) == 1:
			# the 2nd batch Illumina ID is different:
			# ex: HWI-ST700693:182:D0MGFACXX:7:1101:1472:1997/1  and
			# ex: HWI-ST700693:182:D0MGFACXX:7:1101:1472:1997/2
				#print r1['ID'][:-1]
				#print r2['ID'][:-1]
				assert r1['ID'][:r1['ID'].find('/')] == r2['ID'][:r2['ID'].find('/')]
			else:
				raise AssertionError, "ID format unknown! {0} and {1}".format(r1['ID'], r2['ID'])
			return r1, r2
				
		else: # single reads
			r = self.read()
			if r is None:
				raise StopIteration
			return r

	def read(self):
		"""
		Read a single line and return the results as dict
		A single bowtie output should be:
		0) ID
		1) strand
		2) ref
		3) (0-based) offset
		4) seq (as string)
		5) ASCII-33 qual (as string)
		"""
		try:	
			line = self.f.next()
		except StopIteration:
			return None
		line = line.strip().split('\t')
		return {'ID': line[0], \
				'strand': line[1], \
				'ref': line[2], \
				'offset': int(line[3]), \
				'seq': line[4], \
				'qual': line[5], \
				'mm': line[7] if len(line)>=8 else ''}

class BowTieWriter:
	def __init__(self, filename, mode='w'):
		assert mode!='w' or not os.path.exists(filename)
		self.mode = mode
		self.f = open(filename, mode)
	
	def close(self):
		self.f.close()

	def write_composite(self, r1, r2, seq, qual, overlap):
		"""
		r1, r2 is the paired reads 
		seq, qual is the composite read with <overlap>
		the new ID will be:
				r1['ID'] COMPOSED/overlap
		and the line follows bowtie output format:
		ID, strand, ref, offset, seq, qual 
		"""
		self.f.write("{0} COMPOSED/{1}\t+\t{2}\t{3}\t{4}\t{5}\n".format(r1['ID'], overlap,\
				r1['ref'], r1['offset'], seq, qual))

	def write(self, r):
		"""
		Implement LATER (TODO)
		"""
		self.f.write(r['ID'] + '\t' + (r['strand'] if 'strand' in r else '+') + '\t' + r['ref'] + '\t' + str(r['offset']) + '\t' + r['seq'] + '\t' + r['qual'] + '\n')
		

