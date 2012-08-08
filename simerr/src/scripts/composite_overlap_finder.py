import os,sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from PEassembly.miscBowTie import FastqReader, FastqWriter, BowTieWriter, FastqReaderPaired

MM_EPSILON = 2
MIN_OVERLAP = 7

def find_overlap(r1, r2, matchf, unf1, unf2):
	"""
	Find the longest overlap of 
	         s1 ............
			        s2 ........
	(p.s. remember to rev-comp s2)
	"""
	s1 = r1['seq']
	q1 = [10**-(ord(x)/10.) for x in r1['qual']]
	s2 = Seq(r2['seq']).reverse_complement().tostring()
	q2 = [10**-(ord(x)/10.) for x in r2['qual']]
	q2.reverse()
	n1 = len(s1)
	n2 = len(s2)
	for offset in xrange(n1):
		experr = 0
		mm = 0
		fail = False
		for i in xrange(offset, n1):
			if i-offset >= n2: break
			if s1[i] != s2[i-offset] and s1[i]!='N' and s2[i-offset]!='N':
				experr += q1[i] + q2[i-offset]
				mm += 1
				if mm > experr + MM_EPSILON:
					fail = True
					break
		if not fail and n1-offset >= MIN_OVERLAP:
			r1['strand'] = '+'
			r2['strand'] = '-'
			r1['ref'] = 'NA'
			r2['ref'] = 'NA'
			r1['offset'] = 0
			r2['offset'] = offset
			r2['seq'] = s2
			r2['qual'] = r2['qual'][::-1]
			matchf.write(r1)
			matchf.write(r2)
			#print >> sys.stderr, "is a match!"
			#print s1
			#print " "*offset + s2
			#print mm, experr
			#print r1
			#print r2
			return True
	#print >> sys.stderr, "no overlap!"
	#print s1
	#print s2
	unf1.write(r1)
	unf2.write(r2)
	return False

def main(fq1, fq2, output_prefix, abundance_filename):
	abundance = {}
	if abundance_filename is None:
		abundance = defaultdict(lambda: 1)
	else:	
		with open(abundance_filename) as f:
			for line in f:
				_id, _count = line.strip().split('\t')
				abundance[_id] = int(_count)
			
	matchf = BowTieWriter(output_prefix + '.overlap.aligned')
	unf1 = FastqWriter(output_prefix + '.overlap.1.unaligned')
	unf2 = FastqWriter(output_prefix + '.overlap.2.unaligned')

	total = 0
	total_expanded = 0
	aligned = 0
	aligned_expanded = 0
	for r1, r2 in FastqReaderPaired(fq1, fq2):
		realid = r1['ID'][:r1['ID'].find('/')]
		total += 1
		total_expanded += abundance[realid]
		if find_overlap(r1, r2, matchf, unf1, unf2): #overlap found
			aligned += 1
			aligned_expanded += abundance[realid]
	
	with open(output_prefix + '.overlap.log', 'w') as f:
		p = aligned*100./total
		f.write("# reads processed: {0}\n".format(total))
		f.write("# reads with at least one reported alignment: {0} ({1:.2f}%)\n".format(aligned,p))
		f.write("# reads that failed to align: {0} ({1:.2f}%)\n".format(total-aligned,100-p))
		f.write("Reported {0} paired-end alignments to 1 output stream(s)\n".format(aligned))

	with open(output_prefix + '.overlap.log_expanded', 'w') as f:
		p = aligned_expanded*100./total_expanded
		f.write("# reads processed: {0}\n".format(total_expanded))
		f.write("# reads with at least one reported alignment: {0} ({1:.2f}%)\n".format(aligned_expanded,p))
		f.write("# reads that failed to align: {0} ({1:.2f}%)\n".format(total_expanded-aligned_expanded,100-p))
		f.write("Reported {0} paired-end alignments to 1 output stream(s)\n".format(aligned_expanded))


	matchf.close()
	unf1.close()
	unf2.close()

if __name__ == "__main__":
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("-1", dest="fq1", required=True, help="Forward fastq filename (.gz)")
	parser.add_argument("-2", dest="fq2", required=True, help="Reverse fastq filename (.gz)")
	parser.add_argument("-o", dest="output", required=True, help="Output prefix")
	parser.add_argument("-a", dest="abundance", help="Abundance filename")
	
	args = parser.parse_args()
	
	main(args.fq1, args.fq2, args.output, args.abundance)
