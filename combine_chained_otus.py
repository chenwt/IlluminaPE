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

def combine_chained_otus(from_filename, to_filename, out_filename, joinby):
	"""
	The from and to files should both be Qiime's _otu.txt output format
	which can be read by iter_qiime_otu()

	Outputs the combined chained otu output to <out_filename>

	NOTE: ids are NOT allowed to have spaces! if there are, they will be removed!
	      only the first part of the ID will be retained as a result.
	"""
	d = dict((k,v) for (k,v) in iter_qiime_otu(from_filename))
	f = open(out_filename, 'w')
	for k,v in iter_qiime_otu(to_filename):
		members = []
		for c in v:
			members += d[c]
		for i in xrange(len(members)):
			members[i] = members[i].split(None,1)[0]
		f.write("{0}\t{1}\n".format(k, joinby.join(members)))
	f.close()

if __name__=="__main__":
	import argparse

	parser = argparse.ArgumentParser(description="Combined chained OTU lists")
	parser.add_argument("--from", dest="from_file", required=True, help="from filename")
	parser.add_argument("--to", dest="to_file", required=True, help="to filename")
	parser.add_argument("-o", dest="out_file", required=True, help="output filename")
	group = parser.add_mutually_exclusive_group()
	group.add_argument("--comma", action="store_true")
	group.add_argument("--tab", action="store_false")

	args = parser.parse_args()

	combine_chained_otus(args.from_file, args.to_file, args.out_file, "," if args.comma else "\t")


