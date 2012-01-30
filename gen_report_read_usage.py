import os,re,sys
import glob
from csv import DictReader

def countline_gz(file): return int(os.popen("zcat {0} | wc -l".format(file)).read())

def make_filtered_for_2nd_batch(dir):
	"""
	"""
	f = open(os.path.join(dir, 'filtered/filtered.txt'), 'w')
	afile = glob.glob(os.path.join(dir, '*L*.fastq.gz'))[0]
	original = countline_gz(afile) / 4
	afile = glob.glob(os.path.join(dir, 'filtered/*L*.fastq.gz'))[0]
	good = countline_gz(afile) / 4
	f.write("{0}:good {1} chastity {2} alignment 0\n".format(os.path.basename(afile),good,original-good))
	f.close()

rex_filtered = re.compile('(\S+):good (\d+) chastity (\d+) alignment (\d+)')
def gen_report_read_usage(dir):
	"""
	ex: <dir> is $EBB/Sample_DS19175
	"""
	report = {}
	# (1) get the number of reads that failed Illumina's chastity or aligned to contaminant
	# filtered.txt should look like: <file>:good 8352457 chastity 1409073 alignment 1109
	curdir = os.popen('pwd').read().strip()
	os.chdir(dir)
	os.chdir('filtered/')
	with open('filtered.txt') as f:
		m = rex_filtered.match(f.readline())
		file,good,chas,contam = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
		original = good+chas+contam
		report = {'original':original, 'chastity_failed':chas, 'contaminant':contam, 'good1':good}
	
	"""
	Read from bowtie:
	# reads processed: 7406026
	# reads with at least one reported alignment: 5419344 (73.17%)
	# reads that failed to align: 1986682 (26.83%)
	"""
	sample = os.path.basename(dir)[len('Sample_'):]
	report['sample'] = sample
	with open(sample + ".bowtie.log") as f:
		for line in f:
			if line.startswith('# reads processed: '):
				ori = int(line.strip().split(': ')[1])
				assert report['good1'] == ori
			if line.startswith('# reads with at least one reported alignment:'):
				good2 = int(line.strip().split(': ')[1].split()[0])
				report['good2'] = good2
				report['unaligned'] = report['good1'] - good2

	# (2) get number of bowtie-aligned composite reads
	# should be from file ex: DS19336.aligned.composite.gz
#	sample = os.path.basename(dir)[len('Sample_'):]
#	report['sample'] = sample
#	composed = countline_gz(sample + '.aligned.composite.gz')
#	report['unaligned'] = report['good1'] - composed
#	report['good2'] = composed

	# (3) get primer detection & removal
	"""
	# of original reads: 5419344
	# of reads removed: 358732 (0.07)
	# of reads remaining: 5060612 (0.93)
	"""
	with open(sample + '.aligned.composite.gz.primer.log') as f:
		for line in f:
			if line.startswith('# of reads removed: '):
				primer_removed = int(line.split(': ')[1].split()[0])
				report['good3'] = report['good2'] - primer_removed

	# (4) get number of low-qual composite reads
	# should be from file ex: DS19336.aligned.composite.gz.phred20_passed.unique.log
	with open(sample + '.aligned.composite.gz.primer_good.gz.phred10_passed.unique.log') as f:
		for line in f:
			if line.startswith('RemovedDueToLowQual: '):
				lowqual = int(line[len('RemovedDueToLowQual: '):])
			if line.startswith('RemainingTotal: '):
				remaintot = int(line[len('RemainingTotal: '):])
			if line.startswith('RemainingUnique: '):
				remainuni = int(line[len('RemainingUnique: '):])
	report['lowqual10'] = lowqual
	report['good4total'] = remaintot
	report['good4unique'] = remainuni

	# (4) get number of 97% OTUs and how much it takes to reach 80%,90%,100% OTU
	# read from: DS19336.aligned.composite.gz.chained_otu.rarefaction
	#
#	raw = open(sample + '.aligned.composite.gz.chained_otu.list').read().split('\t')
#	assert raw[0] == '0.03'
#	otu_num = int(raw[1])
#	get80, get90, get100 = None, None, None
#	for r in DictReader(open(sample + '.aligned.composite.gz.chained_otu.rarefaction'), delimiter='\t'):
#		if get80 is None and float(r['0.03']) >= .8*otu_num: get80 = r['numsampled']
#		if get90 is None and float(r['0.03']) >= .9*otu_num: get90 = r['numsampled']
#		if get100 is None and float(r['0.03']) >= otu_num: get100 = r['numsampled']
#	report['otu97'] = otu_num
#	report['otu97_get80'] = get80
#	report['otu97_get90'] = get90
#	report['otu97_get100'] = get100

	print >> sys.stderr, report
	os.chdir(curdir)
	return report


if __name__ == "__main__":
#	make_filtered_for_2nd_batch(sys.argv[1])
#	sys.exit(0)

	f = sys.stdout
	#f.write("SAMPLE,ORIGINAL,CHASTITY,CONTAM,PASS-FILTER,UNALIGNED,BOWTIE-ALIGNED,GOOD-PRIMER,LOWQUAL10,GOOD4TOTAL,GOOD4UNIQUE\n")
	if 1:
		d = sys.argv[1]
		#if not os.path.isdir(d): continue
		report = gen_report_read_usage(d)
		f.write(report['sample'] + ',')
		f.write(str(report['original']) + ',')
		f.write(str(report['chastity_failed']) + ',')
		f.write(str(report['contaminant']) + ',')
		f.write(str(report['good1']) + ',')
		f.write(str(report['unaligned']) + ',')
		f.write(str(report['good2']) + ',')
		f.write(str(report['good3']) + ',')
		f.write(str(report['lowqual10']) + ',')
		f.write(str(report['good4total']) + ',')
		f.write(str(report['good4unique']) + '\n')
#		f.write(str(report['otu97']) + ',')
#		f.write(str(report['otu97_get80']) + ',')
#		f.write(str(report['otu97_get90']) + ',')
#		f.write(str(report['otu97_get100']) + '\n')

