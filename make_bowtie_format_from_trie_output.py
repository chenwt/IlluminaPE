import os,sys
import pdb
import numpy as np
from math import ceil
from miscBowTie import BowTieReader, BowTieWriter
"""
DS19190.Rprimer_good.gz.experror_good.gz
picked_otus/DS19190.Rprimer_good.gz.experror_good.gz_otus.txt
"""
def writeout(clusters, output_filename):
	w = BowTieWriter(output_filename)
	for cid, c in clusters.iteritems():
		qual = "".join(chr(int(round(x+33))) for x in np.array(c['qual']) / np.array(c['abun']))
		r = {'ID': cid, 'seq': c['seq'], 'qual': qual, 'strand': '+', \
				'ref': str(c['cycle']), 'offset': max(c['abun'])}
		w.write(r)
	w.close()

def main(input_bowtie, otu_filename):
	otu = {} # seqid --> cid
	clusters = {} # cid --> {'seq', 'qual', 'len'}
	for line in open(otu_filename):
		cid, rest = line.strip().split('\t', 1)
		for seqid in rest.split('\t'):
			otu[seqid] = cid

	for r in BowTieReader(input_bowtie, False):
		rid = r['ID'].split()[0]
		cid = otu[rid]
		if not cid in clusters:
			clusters[cid] = {'seq': r['seq'], 'qual': [ord(x)-33 for x in r['qual']], \
					'abun': [1]*len(r['seq']), 'len': len(r['seq']), 'cycle':int(r['offset'])}
		else:
			c = clusters[cid]
			_len = len(r['seq'])
			offset = r['seq'].find(c['seq'])
			if offset >= 0: # the cluster's seq succeeds the current seq
				tail = _len - c['len'] - offset
				if tail > 0:
					papa = r
					#pdb.set_trace()
				c['qual'] = [0]*offset + c['qual'] + [0]*tail
				c['abun'] = [0]*offset + c['abun'] + [0]*tail
				for i,phred in enumerate(r['qual']):
					c['qual'][i] += ord(phred) - 33
					c['abun'][i] += 1
				c['seq'] = r['seq'] 
				if tail < 0:
					c['seq'] += c['seq'][-tail:]
				c['len'] = len(c['seq'])
				c['cycle'] = int(r['offset'])
				assert c['len'] == len(c['qual']) == len(c['abun'])
			else: # the cluster's seq precedes the current seq
				offset = c['seq'].find(r['seq'])
				assert offset >= 0
				tail = _len - c['len'] + offset
				if tail > 0:
					papa = r
					#pdb.set_trace()
				c['qual'] = [0]*offset + c['qual'] + [0]*tail
				c['abun'] = [0]*offset + c['abun'] + [0]*tail
				for i,phred in enumerate(r['qual']):
					c['qual'][i+offset] += ord(phred) - 33
					c['abun'][i+offset] += 1
				if tail > 0:
					c['seq'] += r['seq'][-tail:]
				c['len'] = len(c['seq'])
				assert c['len'] == len(c['qual']) == len(c['abun'])

	return otu, clusters

if __name__ == "__main__":
	otu, clusters = main(*sys.argv[1:3])
	writeout(clusters, sys.argv[3])
	