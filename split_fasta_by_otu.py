import os,sys
import resource
resource.setrlimit(resource.RLIMIT_NOFILE, (8000,8000))
from miscBowTie import BowTieReader, BowTieWriter
from Bio import SeqIO

#class ManyFileDict:
#	def __init__(self, max_no_files):
#		self.d = {}
#		self.count = 0
#
#	def __getitem__(self, k):
#		if self.d[k].closed:
#			self.d[k] = open(self.d[k].name, 'a')
#		return self.d[k]
#
#	def __setitem__(self, k, v):
#		self.d[k] = v
#
#	def has_key(self, k):
#		return k in self.d

def split_fasta_by_otu(fasta_filename, bowtie_gz_filename, otu_filename, output_dir):
	"""
	For each OTU, create a subdir <output_dir>/<cluster_index> and 
	put it in the OTU's fasta and bowtie (gzipped)
	"""
	otu = {}
	fa_d = {}
	bw_d = {}
	cids = set()
	with open(otu_filename) as f:
		for line in f:
			raw = line.strip().split()
			cid = raw[0]
			cids.add(cid)
			if not os.path.exists(os.path.join(output_dir, cid)): os.mkdir(os.path.join(output_dir, cid))
			#os.mkdir(os.path.join(output_dir, cid))
			for seqid in raw[1:]: 
				otu[seqid] = cid
	print >> sys.stderr, "finished reading", otu_filename

	for cid in cids:
		fa_d[cid] = open(os.path.join(output_dir, cid, cid+'.fasta'), 'w')
		bw_d[cid] = BowTieWriter(os.path.join(output_dir, cid, cid+'.bowtie'), mode='w')

	for r in SeqIO.parse(open(fasta_filename), 'fasta'):
		if r.id not in otu:
			continue
		cid = otu[r.id]
		fa_d[cid].write(">{0}\n{1}\n".format(r.id, r.seq))
	
	for r in BowTieReader(bowtie_gz_filename, False):
		try:
			cid = otu[r['ID'].split()[0]]
		except KeyError:
			continue
		bw_d[cid].write(r)

	for handle in fa_d.itervalues():
		handle.close()

	for handle in bw_d.itervalues():
		handle.close()
		os.system("gzip " + handle.f.name)


if __name__ == "__main__":
	split_fasta_by_otu(*sys.argv[1:])
