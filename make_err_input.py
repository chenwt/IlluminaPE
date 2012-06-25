import os,sys
from Bio.Seq import MutableSeq
import pdb
from miscBowTie import BowTieReader, FastqWriter
from BaseErrPredict import BaseErrPredict

rint = lambda x: int(round(x))

class ClusterCorrect:
	def __init__(self, input_bowtie, otu_txt, primer, Rdata, max_allowed_remain_mismatch, mm_threshold):
		"""
		input_bowtie --- in BowTie format (has seq, qual, and primer offset), gzipped
		otu_txt --- otu file from Qiime (90%)
		"""
		self.predictor = BaseErrPredict(Rdata)
		self.input_bowtie = input_bowtie
		self.otu_txt = otu_txt
		self.primer = primer
		self.max_allowed_remain_mismatch = max_allowed_remain_mismatch
		self.mm_threshold = mm_threshold
		
		self.otu_info = {} # cid (after trie) --> 90% OTU id
		self.cluster_by_otu = {} # OTU id --> cid --> qual,len,seq... 
		self.create_clusters_from_bowtie()
	
	def create_clusters_from_bowtie(self):
		"""
		The 'offset' field is actually 'abundance'
		The 'ref' field is actually 'cycle' offset
		"""
		with open(self.otu_txt) as f:
			for line in f:
				otuid, rest = line.strip().split(None, 1)
				for x in rest.split():
					self.otu_info[x] = otuid
				self.cluster_by_otu[otuid] = {}

		for r in BowTieReader(self.input_bowtie, False):
			cid = r['ID']
			otuid = self.otu_info[r['ID']]
			self.cluster_by_otu[otuid][cid] = {'dirty':True, 'cids':[cid], 'len':len(r['seq']), 'seq': MutableSeq(r['seq']), 'size':int(r['offset']), \
					'qual': [ord(x)-33 for x in r['qual']], 'cycle': range(int(r['ref']), int(r['ref'])+len(r['seq']))}		
		
	def iter_merge_otu(self, otuid, logf):
		cluster = self.cluster_by_otu[otuid]
		ckeys = cluster.keys()
		ckeys.sort(key=lambda x: cluster[x]['size'], reverse=True)
		print >> sys.stderr, "cluster sizes", [cluster[x]['size'] for x in ckeys]
		i = 0
		while i < len(ckeys):
			r1 = cluster[ckeys[i]]
			j = i + 1
			changed = False
			while j < len(ckeys):
				r2 = cluster[ckeys[j]]
				# we must check if one or both of r1, r2 has been modified (dirty=True)
				# otherwise, we have checked this combo before and no need to do again
				if not r1['dirty'] and not r2['dirty']:
					j += 1
					continue
				justmm = r1['size'] >= self.mm_threshold and r1['size']*1./r2['size'] >= self.mm_threshold
				#justmm = r1['size'] >= self.mm_threshold and r2['size'] < self.mm_threshold # TODO: hard code 10 or decide to param!
				flag, diff_poses, correct = self.predictor.compare_seq(\
					r1['seq'], r2['seq'], r1['qual'], r2['qual'],\
					r1['cycle'], r2['cycle'], self.primer, self.max_allowed_remain_mismatch, justmm)
				if flag: # can correct!
					changed = True
					#print "Correcting!", justmm, diff_poses, correct
					#print r1
					#print r2
					#pdb.set_trace()
					logf.write("{cids}\t{pos}\n".format(cids=",".join(r2['cids']),pos=",".join(map(str,diff_poses))))
					for pos, base, q, c in correct:
						r1['seq'][pos] = base
						r1['qual'][pos] = q
						r1['cycle'][pos] = c
					r1['cids'] += r2['cids']
					r1['size'] += r2['size']
					# now we avg. the qual & cycles of non-corrected poses
					if r1['len'] >= r2['len']:
						for pos in xrange(r2['len']):
							if pos not in diff_poses:
								r1['qual'][pos] = rint((r1['qual'][pos]*r1['size']+r2['qual'][pos]*r2['size'])*1./(r1['size']+r2['size']))
								r1['cycle'][pos] = rint((r1['cycle'][pos]*r1['size']+r2['cycle'][pos]*r2['size'])*1./(r1['size']+r2['size']))
					else: # r2 is longer
						for pos in xrange(r1['len']):
							if pos not in diff_poses:
								r1['qual'][pos] = rint((r1['qual'][pos]*r1['size']+r2['qual'][pos]*r2['size'])*1./(r1['size']+r2['size']))
								r1['cycle'][pos] = rint((r1['cycle'][pos]*r1['size']+r2['cycle'][pos]*r2['size'])*1./(r1['size']+r2['size']))
						# append end of r2 to r1
						for pos in xrange(r1['len'],r2['len']):
							r1['seq'].append(r2['seq'][pos])
							r1['qual'].append(r2['qual'][pos])
							r1['cycle'].append(r2['cycle'][pos])
						r1['len'] = r2['len']
						assert r1['len']==len(r1['qual'])==len(r1['seq'])==len(r1['cycle'])
					#print("after merging")
					#print r1
					#print r2
					#pdb.set_trace()
					#raw_input('wait')
					
					# carefully remove j now that we merged it to i
					del cluster[ckeys[j]]
					ckeys.pop(j)
					#pdb.set_trace()
				else:
					#print "Cannot correct!", diff_poses, correct
					#print r1
					#print r2
					#raw_input('wait')
					j += 1
					#pdb.set_trace()
			r1['dirty'] = changed
			#pdb.set_trace()
			i += 1
			# re-sort in case sizes have changed the order
			ckeys = cluster.keys()
			ckeys.sort(key=lambda x: cluster[x]['size'], reverse=True)

	def iter_merge(self):
		numclust = sum(len(x) for x in self.cluster_by_otu.itervalues())
		print >> sys.stderr, "initial # of clusters: ", numclust
		while True:
			kkeys = self.cluster_by_otu.keys()
			kkeys.sort(key=lambda x: max(p['size'] for p in self.cluster_by_otu[x].itervalues()),reverse=True)
			for k in kkeys:#self.cluster_by_otu:
				self.iter_merge_otu(k, sys.stderr)
			tmp = sum(len(x) for x in self.cluster_by_otu.itervalues())
			print >> sys.stderr, "current # of clusters: ", tmp
			if tmp == numclust:
				break
			numclust = tmp

	def write(self, output_prefix):
		"""
		Output to clusters to a fasta file <output_prefix>.fasta
		>{cluster_index}
		{sequence here}
		
		And to a otu-style file <output_prefix>.otu.txt
		<cluster_index> \t <tab delimited seq IDs>
		
		"""
		w = FastqWriter(output_prefix + '.fq')
		f = open(output_prefix + '.fasta', 'w')
		h = open(output_prefix + '.otu.txt', 'w')
		a = open(output_prefix + '.abundance.txt', 'w')
		for cluster in self.cluster_by_otu.itervalues():
			for o in cluster.itervalues():
				# need to massage qual for fq writing
				o['qual'] = "".join(chr(o['qual'][i]+33) for i in xrange(o['len']))
				o['ID'] = o['cids'][0]
				w.write(o)
				f.write(">{0}\n{1}\n".format(o['ID'], o['seq']))
				h.write("{0}\t{1}\n".format(o['ID'], "\t".join(o['cids'])))
				a.write("{0}\t{1}\n".format(o['ID'], o['size']))
		w.close()
		f.close()
		h.close()
		a.close()
		os.system("gzip " + w.f.name)	

if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description="correct errors")
	parser.add_argument("-i", dest="input_bowtie", required=True, help="input bowtie")
	parser.add_argument("-o", dest="otu_txt", required=True, help="trie otu txt")
	parser.add_argument("-p", dest="prefix", required=True, help="output prefix")
	parser.add_argument("--mm", dest="mm", required=True, type=int, help="max allowed remain mismatch")
	parser.add_argument("--mm-threshold", dest="mm_threshold", required=True, type=int, help="Threshold above which we don't accept errcor")
	parser.add_argument("--primer", dest="primer", choices=['F', 'R'])
	parser.add_argument("--primerfile", required=True, help="primer file")
	parser.add_argument("--Rdata", required=True, help="model Rdata file")#default="/shared/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/EBB_Illumina/data.htseq.org/FHCRC_Hullar/FCD05FJ/IlluminaPE//simerr/MPDR.training_size2000.RData", help="model Rdata file")
	args = parser.parse_args()
	
	#Fprimer = "ACTCCTACGGGAGGCAGCAGT"  # 5' to 3'
	#Rprimer = "GTATTACCGCGGCTGCTGGCAC" # 5' to 3'
	#Rdata = "/shared/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/EBB_Illumina/data.htseq.org/FHCRC_Hullar/FCD05FJ/IlluminaPE//simerr/MPDR.training_size2000.RData"
	
	#primer = Fprimer if args.primer=='F' else Rprimer
	primerseq = None
	with open(args.primerfile) as handle:
		for line in handle:
			tag, seq = line.strip().split('\t')
			if tag == args.primer:
				primerseq = seq
				break
	assert primerseq is not None	
	
	o = ClusterCorrect(args.input_bowtie, args.otu_txt, primerseq, args.Rdata, args.mm, args.mm_threshold)
	o.iter_merge()
	o.write(args.prefix)

