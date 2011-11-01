import os,sys

def iter_qiime_otu(filename):
	"""
	Qiime's _otu.txt should have format:
	<cluster id> \t <tab-deilimited seq IDs>
	"""
	with open(filename) as f:
		for line in f:
			raw = line.strip().split('\t')
			yield raw[0], raw[1:]

def combine_chained_otus(from_filename, to_filename, out_filename):
	"""
	The from and to files should both be Qiime's _otu.txt output format
	which can be read by iter_qiime_otu()

	Outputs the combined chained otu output to <out_filename>
	"""
	d = dict((k,v) for (k,v) in iter_qiime_otu(from_filename))
	f = open(out_filename, 'w')
	for k,v in iter_qiime_otu(to_filename):
		members = []
		for c in v:
			members += d[c]
		f.write("{0}\t{1}\n".format(k, "\t".join(members)))
	f.close()

if __name__=="__main__":
	combine_chained_otus(sys.argv[1], sys.argv[2], sys.argv[3])


