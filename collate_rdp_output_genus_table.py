import os, sys, glob
import csv
from collections import defaultdict

d = {}
for r in csv.DictReader(open('sample.annotation.txt'), delimiter='\t'):
	d[r['Sample']] = r

files = glob.glob('../Sample_DS19*/filtered/*chained_otu.rdp_classified.table.order')
genera = set()

_sum_ = defaultdict(lambda: defaultdict(lambda: 0)) # sample -> genus --> count excluding Unclassified
# find the union genera
for x in files:
	f = open(x)
	f.readline()
	for line in f:
		name, count = line.strip().split('\t')
		name = name.replace('"', '') # remove the quotes
		if name!='Unclassified':
			_sum_[x][name] = int(count)
			genera.add(name)
genera = list(genera)
# remove all genera that results in less than 0.0001 fraction total
#genera = set()
#for v in _sum_.itervalues():
#	total = sum(v.itervalues())
#	for name,count in v.iteritems():
#		if count*1./total >= 0.0001:
#			genera.add(name)
#genera = list(genera) 
genera.sort(key=lambda x: _sum_[files[0]][x], reverse=True)

fout = sys.stdout

# ========= below is the R-table format output ============ #
fout.write("SAMPLE,WELL,FAT,GROUP")
for g in genera: fout.write(',' + g)
fout.write('\n')
for x in files:
	sample = os.path.basename(x)[2:7] # getting 19175
	fout.write(sample + ',')
	fout.write(d[sample]['Well'] + ',')
	fout.write(d[sample]['Fat'] + ',')
	fout.write(d[sample]['Name'])
	total = sum(_sum_[x].itervalues())
	for g in genera:
		fout.write("," + str(_sum_[x][g]*1./total))
	fout.write("\n")

# ========= below  is the PC-ord format output ============ #
#fout.write("{0} samples\n".format(len(files)))
#fout.write("{0} species\n".format(len(genera)))
#for i in xrange(len(genera)): fout.write("\tQ")
#fout.write("\n")
#for i in xrange(len(genera)): fout.write("\t" + genera[i])
#fout.write("\n")
#for x in files:
#	y = os.path.basename(x)
#	sample_name = d[y[:y.find('.')][2:]]['Notes']
#	fout.write(sample_name)
#	total = sum(_sum_[x].itervalues())
#	for g in genera:
#		fout.write("\t{0}".format(_sum_[x][g]*1./total))
#	fout.write("\n")
	
	

