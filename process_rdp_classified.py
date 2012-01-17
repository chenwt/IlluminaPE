import os,sys,gzip
from collections import defaultdict

def process_rdp_classifier(sample):
	"""
	Given <sample> which should be ex: DS19175.aligned.composite.gz, there should be:
	 - <sample>.chained_otu.txt.gz
	 - <sample>.chained_otu.rdp_classified 

	IDs from .count file are: HWI-ST700693:106:D05FJACXX:2:1107:14610:198816 2:N:0:ATCACG/2 COMPOSED/12
	IDs from .rdp_classified only have the first part (HWI-ST700693:106:D05FJACXX:2:1107:14610:198816)

	Process the RDP output to generate a R-friendly table w/ headers:
	<ID> <overlap> <count> <phylum> <conf for phylum> <class> <conf> <family> <conf> <order> <conf> <genus> <conf>
	
	Output file is: <sample>.chained_otu.rdp_classified.table
	"""
	overlap = {} # shortened ID --> overlap size
	counts  = {} # shortened ID --> counts
	with gzip.open(sample + '.chained_otu.txt.gz') as f:
		for line in f:
			cluster_id, members = line.strip().split('\t', 1)
			overlap[cluster_id] = defaultdict(lambda: 0)
			for m in members.split('\t'):
				overlap[cluster_id][int(m.split(None)[2][len('COMPOSED/'):])] += 1
			counts[cluster_id] = len(members.split('\t'))
	
	levels = ["phylum", "class", "order", "family", "genus"]
	h = open(sample + '.chained_otu.rdp_classified.table', 'w')
	h.write("ID,overlap,count,phylum,phylumConf,class,classConf,family,familyConf,order,orderConf,genus,genusConf\n")
	with open(sample + '.chained_otu.rdp_classified') as f:
		for line in f:
			raw = line.strip().split('\t')
			id = raw[0]
			h.write(id + ',')
			p = overlap[id].items()
			p.sort(key=lambda x: x[0])
			for _overlap,_count in p: h.write("{0}:{1}|".format(_overlap,_count))
			h.seek(h.tell()-1)
			h.write(',')
			h.write(str(counts[id]) + ',')
			
			i = 9
			for lv in levels:
				if raw[i] == lv:
					h.write(raw[i-1] + ',' + raw[i+1] + ',')
				else:
					h.write('NA,NA,')
				i += 3
			h.seek(h.tell()-1)
			h.write('\n')
	h.close()

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Process <sample>.unique.fasta.rdp_classified output into R-friendly tables')
	parser.add_argument("-i", dest="sample", help="Sample name, ex: DS19198.aligned.composite.gz")

	args = parser.parse_args()
	process_rdp_classifier(args.sample)

