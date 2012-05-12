import os,sys
from csv import DictReader
from Bio import SeqIO
from Bio.Seq import MutableSeq
import pdb
from collections import namedtuple, defaultdict
from miscBowTie import BowTieReader, FastqWriter
from BaseErrPredict import BaseErrPredict

AccCount = namedtuple('AccCount', ['acc', 'count'])
rint = lambda x: int(round(x))

class ClusterCorrect:
	def __init__(self, input_bowtie, otu_txt, primer, Rdata, max_allowed_remain_mismatch):
		"""
		input_bowtie --- in BowTie format (has seq, qual, and primer offset)
		otu_txt --- otu file from Qiime (should be done trie, 100% or suffix)
		"""
		self.predictor = BaseErrPredict(Rdata)
		self.input_bowtie = input_bowtie
		self.otu_txt = otu_txt
		self.primer = primer
		self.max_allowed_remain_mismatch = max_allowed_remain_mismatch
		
		self.otu_info = None
		self.cluster = self.create_clusters()
		# sort the cluster by size
		self.cluster_keys = self.cluster.keys()
		self.cluster_keys.sort(key=lambda k: self.cluster[k]['size'], reverse=True)
		
	def create_clusters(self):
		"""
		For each OTU, we avg. the phred score and cycles
		NOTE: make sure the cycles are GLOBAL not after primer removal
		"""
		self.otu_info = {} # cluster index --> list of seqs
		# the otu file is probably smaller than BowTie file so read it first
		cluster = {} # cluster id --> {size, avg. qual per pos, avg. cycle per pos, long seq}
		otu = {} # seqID --> cluster index
		with open(self.otu_txt) as f:
			for line in f:
				raw = line.strip().split('\t')
				cid = int(raw[0])
				for seqid in raw[1:]: otu[seqid] = cid
				self.otu_info[cid] = raw[1:]
		
		# now read the bowtie file
		for r in BowTieReader(self.input_bowtie, False):
			rid = r['ID'].split()[0] # sometimes seqids have extra stuff like COMPOSED/14
			cid = otu[rid]
			offset = int(r['offset'])
			if cid not in cluster:
				cluster[cid] = {'size':0, 'qual':defaultdict(lambda: {'acc':0,'count':0}), \
							'seq':None, 'cycle':defaultdict(lambda: {'acc':0,'count':0}),\
							'cid':cid, 'cids':[cid], 'len':None}
			cluster[cid]['size'] += 1
			# assign the longest sequence to the cluster as representative
			if cluster[cid]['seq'] is None or len(r['seq']) > cluster[cid]['len']:
				cluster[cid]['seq'] = r['seq']
				cluster[cid]['len'] = len(r['seq'])
			for pos in xrange(len(r['seq'])):
				cluster[cid]['qual'][pos]['acc'] += ord(r['qual'][pos]) - 33
				cluster[cid]['qual'][pos]['count'] += 1
				cluster[cid]['cycle'][pos]['acc'] += pos + offset
				cluster[cid]['cycle'][pos]['count'] += 1
		
		# avg. the qual and cycles
		for x in cluster.itervalues():
			x['seq'] = MutableSeq(x['seq'])
			for pos in x['qual']: x['qual'][pos] = rint(x['qual'][pos]['acc']*1./x['qual'][pos]['count'])
			for pos in x['cycle']: x['cycle'][pos] = rint(x['cycle'][pos]['acc']*1./x['cycle'][pos]['count'])
				
		return cluster				
		
	def iter_merge(self):
		i = 0
		while i < len(self.cluster_keys)-1:
			r1 = self.cluster[self.cluster_keys[i]]
			j = i + 1
			while j < len(self.cluster_keys):
				r2 = self.cluster[self.cluster_keys[j]]
				justmm = r1['size'] >= 10 and r2['size'] >= 10 # TODO: hard code 10 or decide to param!
				flag, diff_poses, correct = self.predictor.compare_seq(\
					r1['seq'], r2['seq'], r1['qual'], r2['qual'],\
					r1['cycle'], r2['cycle'], self.primer, self.max_allowed_remain_mismatch, True)
				if flag: # can correct!
					print "Correcting!", correct
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
							r1['qual'][pos] = r2['qual'][pos]
							r1['cycle'][pos] = r2['cycle'][pos]
						r1['len'] = r2['len']
					#print("after merging")
					#print r1
					#print r2
					#raw_input('wait')
					
					# carefully remove j now that we merged it to i
					del self.cluster[self.cluster_keys[j]]
					self.cluster_keys.pop(j)
					#pdb.set_trace()
				else:
					print "Cannot correct!", correct
					j += 1
					#pdb.set_trace()
			#pdb.set_trace()
			i += 1

	def write(self, output_prefix):
		"""
		Output to clusters to a fasta file <output_prefix>.fasta
		>{cluster_index}
		{sequence here}
		
		And to a otu-style file <output_prefix>.otu.txt
		<cluster_index> \t <tab delimited seq IDs>
		
		Output is in decreasing order of cluster size
		"""
		# sort the cluster by size
		self.cluster_keys = self.cluster.keys()
		self.cluster_keys.sort(key=lambda k: self.cluster[k]['size'], reverse=True)
		w = FastqWriter(output_prefix + '.fq')
		f = open(output_prefix + '.fasta', 'w')
		h = open(output_prefix + '.otu.txt', 'w')
		for cid in self.cluster_keys:
			o = self.cluster[cid]
			# need to massage qual for fq writing
			o['qual'] = "".join(chr(o['qual'][i]+33) for i in xrange(len(o['seq'])))
			o['ID'] = cid
			w.write(o)
			f.write(">{0}\n{1}\n".format(cid, o['seq']))
			h.write("{0}".format(cid))
			for member_cid in o['cids']:
				h.write("\t" + "\t".join(self.otu_info[member_cid]))
			h.write("\n")
		w.close()
		f.close()
		h.close()
		os.system("gzip " + w.f.name)	

if __name__ == "__main__":
	from argparse import ArgumentParser
	parser = ArgumentParser(description="correct errors")
	parser.add_argument("-i", dest="input_bowtie", required=True, help="input bowtie")
	parser.add_argument("-o", dest="otu_txt", required=True, help="trie otu txt")
	parser.add_argument("-p", dest="prefix", required=True, help="output prefix")
	parser.add_argument("--mm", dest="mm", required=True, type=int, help="max allowed remain mismatch")
	parser.add_argument("--primer", dest="primer", choices=['F', 'R'])
	args = parser.parse_args()
	
	Fprimer = "ACTCCTACGGGAGGCAGCAGT"  # 5' to 3'
	Rprimer = "GTATTACCGCGGCTGCTGGCAC" # 5' to 3'
	Rdata = "/shared/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/EBB_Illumina/data.htseq.org/FHCRC_Hullar/FCD05FJ/IlluminaPE//simerr/MPDR.training_size2000.RData"
	
	primer = Fprimer if args.primer=='F' else Rprimer
	o = ClusterCorrect(args.input_bowtie, args.otu_txt, primer, Rdata, args.mm)
	o.iter_merge()
	o.write(args.prefix)

