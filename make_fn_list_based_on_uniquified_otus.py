import os,sys,gzip
from collections import defaultdict

def make_fn_list(otu_txt, out_filename, count_filename=None, cutoff='0.03'):
	"""
	Assume that otu_txt is the result from running Qiime's pick_otus.py using some cutoff
	Then we want to make a .fn.list in Mothur's style so that we can run single.rarefaction
	However remember that .unique_otus.txt is also uniquify-ed so must put the replicates back in!
	
	If otu_txt is DS19201.aligned.composite.gz.unique_otus.txt, then
	the unique.count file must be ../DS19201.aligned.composite.gz.unique.count.gz

	if count_filename is None, then that makes _otu.txt already takes it into account.

	Output is .list format that mothur uses for rarefaction n stuff, should be
	<cutoff> \t <# of clusters> \t <space-delimited clusters where IDs are comma-separated>
	ex: 0.03	3	a b,c,d e,f

	Since mothur's single.rarefaction doesn't actually pay attention to what the seq IDs are,
	we're just going to use 'a' to represent all seq IDs!
	"""
	#assert otu_txt.endswith('_otus.txt')
	d = {} # representative id --> count
	if count_filename is not None:
		for line in open(count_filename):
			rep_id, count, members = line.strip().split('\t')
			rep_id = rep_id.split(None)[0] # the ID had some clean up done here...
			d[rep_id] = int(count)
	else:
		d = defaultdict(lambda: 1)

	p = open(otu_txt).read().strip().split('\n')
	f = open(out_filename, 'w')
	f.write(cutoff + '\t')
	f.write(str(len(p)) + '\t')
	for x in p:
		# x is <cluster no.> \t <tab-delimited list of rep_ids in this cluster>
		# in this cluster, the represenative ids each represent 1 or more identical seqs
		_total = sum(d[rep_id]  for rep_id in x.split('\t')[1:])
		f.write(('a,'*_total)[:-1] + ' ')
		#raw_input("cluster {0} has {1} reps which sums to {2} total".format(x.split('\t')[0], len(x.split('\t'))-1, _total))
	f.close()

if __name__ == "__main__":
	make_fn_list(sys.argv[1], sys.argv[2])

