import os, sys
import random
import pdb
from Bio.Blast import NCBIXML
from miscBowTie import BowTieReader

def parse_blast_xml_for_training(xml_filename, bowtie_filename, output_filename):
    """
    Parse the XML output, looking only at the 1st alignment for each query
    Write out in format:
        
    Phred   Cycle   B2  B1  B0  Class
    """
    fa_dict = dict((r['ID'], r) for r in BowTieReader(bowtie_filename, False))

    f = open(output_filename, 'w')
    f.write("Phred\tCycle\tB2\tB1\tB0\tClass\n")
    for blastout in NCBIXML.parse(open(xml_filename)):
        if len(blastout.alignments) == 0: # no match was found!
            continue
        hsp = blastout.alignments[0].hsps[0]
        record = fa_dict[blastout.query]
        primer_offset = int(record['offset'])
        gap_offset = 0
        for i in xrange(2, len(hsp.match)):# toDO: allow for i<2 and still get B2, B1
            gap_offset += hsp.query[i]=='-'
            # global position is i + (query_start-1) + primer_offset - gap_offset
            doit = False
            if hsp.match[i]==" " and hsp.query[i]!='-' and hsp.sbjct[i]!='-':
                doit = True
                _class = '-'
            elif hsp.match[i]=='|' and hsp.query[i]==hsp.sbjct[i] and random.random() <= 1e-3:
                doit = True
                _class = '+'
            if doit:
                # is a mismatch!
                record_i = i + (hsp.query_start - 1) - gap_offset
                assert hsp.query[i] == record['seq'][record_i]
                f.write(str(ord(record['qual'][record_i])-33) + '\t')
                f.write(str(record_i + primer_offset) + '\t')
                b2, b1 = None, None
                j = i - 1
                while hsp.query[j]=='-': j -= 1
                b1 = hsp.query[j]
                j -= 1
                while hsp.query[j]=='-': j -= 1
                b2 = hsp.query[j]
                assert b1==record['seq'][record_i-1] and b2==record['seq'][record_i-2]
                f.write(b2 + '\t')
                f.write(b1 + '\t')
                f.write(hsp.query[i] + '\t')
                f.write(_class + '\n')

    f.close()
    
if __name__ == "__main__":
    parse_blast_xml_for_training(*sys.argv[1:])
