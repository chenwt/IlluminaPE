import os, sys, glob

refDB = 'Sample_DS19175/filtered/SILVA100_justHumanCrap.1crap_masked.V3region_ungapped_nonredundant'

for file in glob.iglob('Sample*/filtered/*fastq.gz'):
	if file.find('R1') > 0:
		file2 = file.replace('R1', 'R2')
		sample = file[:file.find('/')][len('Sample_'):]
		print("gunzip {0} {1}".format(file, file2))
		file = file[:-3]
		file2 = file2[:-3]
		print("bowtie --phred33-quals --un {5}/{3} {0} -1 {1} -2 {2} {5}/{4} 2> {5}/{6}".format(refDB, \
				file, file2, sample+'.unaligned', sample+'.aligned', os.path.dirname(file),\
				sample+'.bowtie.log'))
		print("gzip {0} {1} {2}/*aligned".format(file, file2, os.path.dirname(file)))
