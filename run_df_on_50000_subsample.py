import os,sys
from DigitalFingerprint import Read, DF
refmap = Read.RefMap(os.environ['PCODE'] + '/src/Silva/SILVA100_justHumanCrap.1crap_masked.V3region.fna.gap_map.bz2', 520, os.environ['PCODE'] + '/src/Silva/SILVA100_justHumanCrap.1crap_masked.V3region_ungapped.fna')

from subprocess import Popen, PIPE
import cStringIO
import random
def g1(gz_filename, size, subsample_size):
	fp = Popen(["zcat", gz_filename], stdout=PIPE)
	f = cStringIO.StringIO(fp.communicate()[0])
	got = random.sample(xrange(size), subsample_size)
	got.sort()
	readdf = Read.ReadDF(gz_filename, refmap)
	for i in xrange(got[0]):
		f.next()
	line = f.next().strip().split('\t')
	r = Read.Read(line[0],line[-2],ref_seq_id=line[2],offset=int(line[3]))
	readdf.add_read_to_vec(r)
	for i in xrange(1, len(got)):
		for x in xrange(got[i]-got[i-1]): line = f.next()
		line = line.strip().split('\t')
		r = Read.Read(line[0],line[-2],ref_seq_id=line[2],offset=int(line[3]))
		readdf.add_read_to_vec(r)
	return readdf

if __name__ == "__main__":
	gz_filename = sys.argv[1]
	subsample_size = 5*10**5
	size = int(os.popen("zcat {0} | wc -l".format(gz_filename)).read())
	readdf = g1(gz_filename, size, subsample_size)
	with open(gz_filename + '.subsample' + str(subsample_size) + '.DF', 'w') as f:
		w = DF.DFWriter(f)
		w.write(readdf)

