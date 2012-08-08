import os, sys

if __name__ == "__main__":
	from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
	parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
	parser.add_argument("-1", dest="fq1", required=True, help="Forward fastq filename")
	parser.add_argument("-2", dest="fq2", required=True, help="Reverse fastq filename")
	parser.add_argument("-o", dest="output", required=True, help="Output prefix")
	parser.add_argument("-r", dest="refdb", default='/shared/silo_researcher/Lampe_J/Gut_Bugs/EBBforChangHao/SILVA/SSURef_108_tax_silva.DNA', help="BowTie reference database index")
	
	args = parser.parse_args()
	
	os.system("bowtie --chunkmbs 2000 -n 3 -y -l 10 -q --fr -1 {fq1} -2 {fq2} --un {output}.bowtie.unaligned {refdb} {output}.bowtie.aligned 2> {output}.bowtie.log".format(fq1=args.fq1,fq2=args.fq2, output=args.output, refdb=args.refdb))

