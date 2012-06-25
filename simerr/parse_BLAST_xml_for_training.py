import os,sys
import pdb
from Bio.Blast import NCBIXML
from miscBowTie import BowTieReader

def parse_blast_xml_for_training(xml_filename, bowtie_filename, output_filename):
    """
    Parse the XML output, looking only at the 1st alignment for each query
    Write out in format:
        
    Phred   Cycle   B2  B1  B0  Class
    """
    fa_dict = dict((r['ID'], r) for r in BowTieReader(bowtie_filename))

    f = open(output_filename, 'w')
    f.write("Phred\tCycle\tB2\tB1\tB0\tClass\n")
    for blastout in NCBIXML.parse(open(xml_filename)):
        if len(blastout.alignments) == 0: # no match was found!
            continue
        hsp = blastout.alignments[0].hsps[0]
        record = fa_dict[blastout.query]
        primer_offset = int(record['offset'])
        for i in xrange(2, len(hsp.match)):# toDO: allow for i<2 and still get B2, B1
            # global position is i + (query_start-1) + primer_offset
            if hsp.match[i]==" " and hsp.query[i]!='-' and hsp.sbjct[i]!='-':
                pdb.set_trace()
                # is a mismatch!
                f.write(str(ord(record['qual'][i+hsp.query_start-1])-33) + '\t')
                f.write(str(i + hsp.query_start - 1 + primer_offset) + '\t')
                f.write(hsp.query[i-2] + '\t')
                f.write(hsp.query[i-1] + '\t')
                f.write(hsp.query[i] + '\t')
                f.write('-\n')

    f.close()
    
if __name__ == "__main__":
    parse_blast_xml_for_training(*sys.argv[1:])
