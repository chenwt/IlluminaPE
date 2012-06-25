import os,sys
from miscBowTie import FastqReader, FastqWriter
from argparse import ArgumentParser

def consolidate_corrected_clusters(dir, output_prefix):
	"""
	in each cluster <dir>/<cluster_index>
	find the corrected <cluster_index>.errcor.{fq|otu.txt}
	and consolate them into one file
	rename the new seq ids to <cluster_index>_<otu_index>
	"""
	fqw = FastqWriter(output_prefix+'.fq')
	otuw = open(output_prefix+'.otu.txt', 'w')

	for cid in os.listdir(dir):
		d2 = os.path.join(dir, cid)
		with open(os.path.join(d2, cid+'.errcor.otu.txt')) as f:
			for line in f:
				otuw.write("{cid}_{rest}".format(cid=cid, rest=line))
		for r in FastqReader(os.path.join(d2, cid+'.errcor.fq.gz')):
			r['ID'] = cid + '_' + r['ID']
			fqw.write(r)
	
	otuw.close()
	fqw.close()

if __name__ == "__main__":
	parser = ArgumentParser(description="Consolidate clusters after correction")
	parser.add_argument("-d", dest="dir", required=True, help="directory containing clusters")
	parser.add_argument("-o", dest="output_prefix", required=True, help="output prefix")
	
	args = parser.parse_args()
	consolidate_corrected_clusters(args.dir, args.output_prefix)

