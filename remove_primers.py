import os,sys,time
import itertools
from cPickle import dump
from Bio.Seq import Seq
from miscBowTie import BowTieReader, BowTieWriter, FastqReader, FastqWriter
from c_PrimerMatch import PrimerMatch, match_primer_len

#def detect_primers_combined(input, f_primer, r_primer, max_mm_allowed=2, min_overlap=10):
#	"""
#	NOTE: this is for COMPOSITE sequences that have both F & R primers!!
#	"""
#	os.system("rm {0}.primer_*".format(input))
#	hgood = BowTieWriter(input + '.primer_good')
#	hbad  = BowTieWriter(input + '.primer_bad')
#	hlog = open(input + '.primer.log', 'w')
#	start_t = time.time()
#	good, bad = 0,0
#	for r in BowTieReader(input, False):
#		match_f_len = match_primer_len(r['seq'], f_primer, max_mm_allowed, min_overlap, False)
#		match_r_len = match_primer_len(r['seq'], r_primer, max_mm_allowed, min_overlap, True)
#		if match_f_len > 0 and match_r_len > 0:
#			r['seq'] = r['seq'][match_f_len:-match_r_len]
#			r['qual'] = r['qual'][match_f_len:-match_r_len]
#			r['offset'] = match_f_len
#			hgood.write(r)
#			good += 1
#		else:
#			hbad.write(r)
#			bad += 1
#	hlog.write("Primer detection and removal took {0} sec.\n".format(time.time()-start_t))
#	hlog.write("# of original reads: {0}\n".format(good+bad))
#	hlog.write("# of reads removed: {0} ({1:.2f})\n".format(bad,bad*1./(good+bad)))
#	hlog.write("# of reads remaining: {0} ({1:.2f})\n".format(good,good*1./(good+bad)))
#
#	hgood.close()
#	hbad.close()
#	hlog.close()
#	os.system("gzip " + hgood.f.name)
#	os.system("gzip " + hbad.f.name)

def detect_primers_PE(input1, input2, output_prefix, f_primer, r_primer, min_match_len, max_mm, max_de, max_in):
	"""
	NOTE: this is for paired end reads that comes in two separate files
	ex: DS19342_CTTGTA_L006_R1_001.fastq.gz and DS19342_CTTGTA_L006_R2_001.fastq.gz
	
	Given a pair of reads from input1, input2:
	1. Detect that F primer exists in one read and R primer in the other
	2. If both reads pass primer detection, output
	3. Otherwise, discard
	
	Output:  <output_prefix>.{F|R}primer_good
	         <output_prefix>.primer.bad
	         <output_prefix>.primer.log
	"""
	def process_primer(r, match_len, is_reverse):
		# get record into miscBowTie.BowTieReader format 
		# strip away primers from seq & qual, properly rev comp!
		r['offset'] = match_len
		r['seq'] = r['seq'][match_len:]
		r['qual'] = r['qual'][match_len:]
		r['ref'] = 'NA'
		if is_reverse:
			r['seq'] = Seq(r['seq']).reverse_complement().tostring()
			r['qual'] = r['qual'][::-1]
	
	os.system("rm {0}.*primer_*".format(output_prefix))
	Fgood = BowTieWriter(output_prefix + '.Fprimer_good')
	Rgood = BowTieWriter(output_prefix + '.Rprimer_good')
	hbad1 = FastqWriter(output_prefix + '.primer_bad.1')
	hbad2 = FastqWriter(output_prefix + '.primer_bad.2')
	hverbose = open(output_prefix + '.primer.verbose', 'w')
	hlog = open(output_prefix + '.primer.log', 'w')
	start_t = time.time()
	good, bad = 0,0
	
	pmF = PrimerMatch(f_primer)
	pmR = PrimerMatch(r_primer)

	for r1, r2 in itertools.izip(FastqReader(input1), FastqReader(input2)):
		# NOTE: in the case of PE reads
		#       regardless of whether we're matching for F or R primer
		#       they would all appear at the 5' end of the read
		#       which is why we call match_primer_len with is_reverse = False
		match_f_len1, mmf1 = match_primer_len(r1['seq'], f_primer, max_mm, min_match_len, False)
		match_r_len1, mmr1 = match_primer_len(r1['seq'], r_primer, max_mm, min_match_len, False)
		match_f_len2, mmf2 = match_primer_len(r2['seq'], f_primer, max_mm, min_match_len, False)
		match_r_len2, mmr2 = match_primer_len(r2['seq'], r_primer, max_mm, min_match_len, False)
		#match_f_len1 = match_f_len2 =match_r_len1=match_r_len2=0
		if match_f_len1 > 0 and match_r_len2 > 0:
			# case 1, read 1 is F, read 2 is R
			good += 1
			process_primer(r1, match_f_len1, False)
			Fgood.write(r1)
			process_primer(r2, match_r_len2, False)
			Rgood.write(r2)
		elif match_f_len2 > 0 and match_r_len1 > 0:
			# case 2, read 1 is R, case 2 is F
			good += 1
			process_primer(r2, match_f_len2, False)
			Fgood.write(r2)
			process_primer(r1, match_r_len1, False)
			Rgood.write(r1)
		else:
			pmF.make_suffix(r1['seq'])
			pmF.match(min_match_len, max_mm, max_in, max_de)
			if pmF.match_result is not None: 
				pmR.make_suffix(r2['seq'])
				pmR.match(min_match_len, max_mm, max_in, max_de)
				if pmR.match_result is not None:  # case 1, read 1 is F, read 2 is R
					good += 1
					process_primer(r1, pmF.match_result.match_len, False)
					Fgood.write(r1)
					hverbose.write("{0}\t{1}\t{2}\n".format(r1['ID'], pmF.match_result.match_len, pmF.match_result.miss))
					process_primer(r2, pmR.match_result.match_len, False)
					Rgood.write(r2)
					hverbose.write("{0}\t{1}\t{2}\n".format(r2['ID'], pmR.match_result.match_len, pmR.match_result.miss))
				else:
					hbad1.write(r1)
					hbad2.write(r2)
					bad += 1
			else:
				pmR.make_suffix(r1['seq'])
				pmR.match(min_match_len, max_mm, max_in, max_de)
				if pmR.match_result is not None:
					pmF.make_suffix(r2['seq'])
					pmF.match(min_match_len, max_mm, max_in, max_de)
					if pmF.match_result is not None:
						good += 1
						# case 2, read 1 is R, read 2 is F
						process_primer(r2, pmF.match_result.match_len, False)
						hverbose.write("{0}\t{1}\t{2}\n".format(r2['ID'], pmF.match_result.match_len, pmF.match_result.miss))
						Fgood.write(r2)
						process_primer(r1, pmR.match_result.match_len, False)
						Rgood.write(r1)
						hverbose.write("{0}\t{1}\t{2}\n".format(r1['ID'], pmR.match_result.match_len, pmR.match_result.miss))
					else:
						# case 3: unresolved, bad read pair
						hbad1.write(r1)
						hbad2.write(r2)
						bad += 1

	hlog.write("Input 1: {0}\nInput 2: {1}\n".format(input1, input2))
	hlog.write("F primer: {0}\nR primer: {1}\n".format(f_primer, r_primer))
	hlog.write("Min match len: {0}\n".format(min_match_len))
	hlog.write("Max mismatch: {0}\n".format(max_mm))
	hlog.write("Max deletion: {0}\n".format(max_de))
	hlog.write("Max insertion: {0}\n".format(max_in))
	hlog.write("Primer detection and removal took {0} sec.\n".format(time.time()-start_t))
	hlog.write("# of original reads: {0}\n".format(good+bad))
	hlog.write("# of reads removed: {0} ({1:.2f})\n".format(bad,bad*1./(good+bad)))
	hlog.write("# of reads remaining: {0} ({1:.2f})\n".format(good,good*1./(good+bad)))


	Fgood.close()
	Rgood.close()
	hbad1.close()
	hbad2.close()
	hlog.close()
	hverbose.close()
	os.system("gzip " + Fgood.f.name)
	os.system("gzip " + Rgood.f.name)
	os.system("gzip " + hbad1.f.name)
	os.system("gzip " + hbad2.f.name)
	os.system("gzip " + hverbose.name)
	
if __name__ == "__main__":
	import argparse
	
#	f_primer = 'ACTCCTACGGGAGGCAGCAGT'
#	r_primer = 'GTATTACCGCGGCTGCTGGCAC' # real reverse primer, use for paired ends
#	r_primer_revcomp = 'GTGCCAGCAGCCGCGGTAATAC' # reverse primer flipped to 5'-3', use for composite or aligned

	parser = argparse.ArgumentParser(description='Remove complete or partial forward & reverse primers',\
									formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("-i", dest="input", help="Composite/Aligned input file (BowTie format, gzipped)")
	parser.add_argument("-1", dest="input1", help="Paired end Fastq file 1")
	parser.add_argument("-2", dest="input2", help="Paired end Fastq file 2")
	parser.add_argument("-o", dest="output", help="Output prefix [used only for PE]")
	parser.add_argument("--mm", dest="mm", default=2, type=int, help="Max mismatch to primer allowed")
	parser.add_argument("--de", dest="de", default=1, type=int, help="Max deletion to primer allowed")
	parser.add_argument("--overlap", dest="overlap", default=10, type=int, help="Min primer overlap")
	parser.add_argument("--fprimer", type=str, required=True, help="forward primer")
	parser.add_argument("--rprimer", type=str, required=False, help="reverse primer")

	args = parser.parse_args()
	
	# sanity checking
	if args.input1 is not None or args.input2 is not None:
		assert args.input1 is not None and args.input2 is not None
		assert args.input is None
		assert args.output is not None
		detect_primers_PE(args.input1, args.input2, args.output, args.fprimer, args.rprimer, args.overlap, args.mm, args.de, 0)
	elif args.input is not None:
		assert args.input1 is None and args.input2 is None
		raise Exception, "Not implemented. Commented out@!"
		#detect_primers_combined(args.input, f_primer, r_primer_revcomp, args.mm, args.offset)
	else:
		print("NO INPUT GIVEN! ABORT!")

