import os, sys, glob

#refDB = '/shared/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/EBB_Illumina/data.htseq.org/FHCRC_Hullar/FCD05FJ/IlluminaPE/refDBs/SILVA100_justHumanCrap.1crap_masked.V3region_ungapped_nonredundant'
# below: complete & updated SILVA SSU db! but is HUGE, need lots of memory (--chunkmbs) 
refDB = '/shared/silo_researcher/Lampe_J/Gut_Bugs/SILVA/fastaExports/SSURef_108_tax_silva.DNA'

for file in glob.iglob('Sample*/filtered/*fastq.gz'): 
	if file.find('R1') > 0:
		file2 = file.replace('R1', 'R2')
		sample = file[:file.find('/')][len('Sample_'):]
		print("gunzip {0} {1}".format(file, file2))
		file = file[:-3]
		file2 = file2[:-3]
		msg = "bowtie --phred33-quals -p 12 --chunkmbs 2000 " # using 12 theads and 2GB first-chunk mem
		msg += "--un {5}/{3} {0} -1 {1} -2 {2} {5}/{4} 2> {5}/{6}".format(refDB, \
				file, file2, sample+'.unaligned', sample+'.aligned', os.path.dirname(file),\
				sample+'.bowtie.log')
		print(msg)
		print("gzip {0} {1} {2}/*aligned".format(file, file2, os.path.dirname(file)))
