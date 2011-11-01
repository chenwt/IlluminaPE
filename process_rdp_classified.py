import os,sys,gzip

def process_rdp_classifier(sample):
	"""
	Given <sample>, there should be:
	 - <sample>.unique.count.gz
	 - <sample>.unique.fasta
	 - <sample>.unique.fasta.rdp_classified

	IDs from .count file are: HWI-ST700693:106:D05FJACXX:2:1107:14610:198816 1:N:0:ATCAGG/1 COMPOSED/25
	IDs from .rdp_classified only have the first part (HWI-ST700693:106:D05FJACXX:2:1107:14610:198816)

	Process the RDP output to generate a R-friendly table w/ headers:
	<ID> <overlap> <count> <phylum> <conf for phylum> <class> <conf> <family> <conf> <order> <conf> <genus> <conf>
	
	Output file is: <sample>.unique.fasta.rdp_classified.table
	"""
	overlap = {} # shortened ID --> overlap size
	counts  = {} # shortened ID --> counts
	with gzip.open(sample + '.unique.count.gz') as f:
		for line in f:
			long_id, count, ignore = line.split('\t')
			short_id, ignore, raw_overlap = long_id.split()
			overlap[short_id] = raw_overlap[len('COMPOSED/'):]
			counts[short_id] = count
	
	levels = ["phylum", "class", "order", "family", "genus"]
	h = open(sample + '.unique.fasta.rdp_classified.table', 'w')
	h.write("ID,overlap,count,phylum,phylumConf,class,classConf,family,familyConf,order,orderConf,genus,genusConf\n")
	with open(sample + '.unique.fasta.rdp_classified') as f:
		for line in f:
			raw = line.strip().split('\t')
			id = raw[0]
			h.write(id + ',' + overlap[id] + ',' + counts[id] + ',')
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
	os.system("gzip " + h.name)
	print >> sys.stderr, "Output written to", h.name, "then gzipped!"


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Process <sample>.unique.fasta.rdp_classified output into R-friendly tables')
	parser.add_argument("-i", dest="sample", help="Sample name, ex: DS19198.aligned.composite.gz")

	args = parser.parse_args()
	process_rdp_classifier(args.sample)

