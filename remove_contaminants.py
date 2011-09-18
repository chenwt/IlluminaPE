import os,sys,gzip,glob

def remove_contaminants(fastq_gz_filename):
	f_d = os.path.join(os.path.dirname(fastq_gz_filename), 'filtered/')
	contaminant_filename = glob.glob(f_d+'*.contaminant')[0]
	print >> sys.stderr, "contaminant filename is", contaminant_filename
	out = gzip.open(os.path.join(f_d, os.path.basename(fastq_gz_filename)), 'w')
	cont = [line.split()[0] for line in open(contaminant_filename)]
	f = gzip.open(fastq_gz_filename)
	good, chastity_failed, aln_matched = 0, 0, 0
	while True:
		line = f.readline().strip()
		if len(line) == 0: 
			break
		if line.split(':',3)[3].split(None)[1].split(':')[1] == 'Y':
			# failed chastity filter
			f.readline()
			f.readline()
			f.readline()
			chastity_failed += 1
		elif line.split(':',3)[3].split(None)[0] in cont: # is contaminant!
			print >> sys.stderr, "removing....", line
			f.readline()
			f.readline()
			f.readline()
			aln_matched += 1
		else:
			good += 1
			out.write(line + '\n')
			out.write(f.readline())
			out.write(f.readline())
			out.write(f.readline())
	f.close()
	out.close()
	os.system("echo \"{0}:good {1} chastity {2} alignment {3}\" >> {4}/filtered.txt".format(fastq_gz_filename, good, chastity_failed, aln_matched, f_d))

if __name__ == "__main__":
	remove_contaminants(sys.argv[1])
