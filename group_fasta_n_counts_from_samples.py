import os,sys
from Bio import SeqIO

def group_samples(in_prefix, add_to_prefix, sample):
	"""
	<prefix>.count and <prefix>.fasta must exist
	ex: DS21061.aligned.composite.gz.primer_good.gz.phred10_passed.unique.count
	"""
	fasta = open(add_to_prefix+'.fasta', 'a+')
	count = open(add_to_prefix+'.count', 'a+')
	for r in SeqIO.parse(open(in_prefix+'.fasta'), 'fasta'):
		fasta.write(">{sample}_{id}\n{seq}\n".format(sample=sample,id=r.id,seq=r.seq))
	with open(in_prefix+'.count') as f:
		for line in f:
			count.write("{sample}_{rest}\n".format(sample=sample,rest=line.strip()))
	fasta.close()
	count.close()

if __name__ == "__main__":
	group_samples(sys.argv[1], sys.argv[2], sys.argv[3])
