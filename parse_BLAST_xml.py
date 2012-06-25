import os,sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Blast import NCBIXML

def parse_blast_xml(xml_filename, query_filename, output_filename, abundance_filename=None):
    """
    Parse the XML output, looking only at the 1st alignment for each query
    Write out in format:
        
    ID \t COUNT \t LENGTH \t AMBIG \t QSTART \t QEND \t IDEN
    """
    if abundance_filename is None:
        abundance = defaultdict(lambda: 1)
    else:
        abundance = dict(line.strip().split('\t') for line in open(abundance_filename))
    handle = NCBIXML.parse(open(xml_filename))
    f = open(output_filename, 'w')
    f.write("ID\tCOUNT\tLENGTH\tAMBIG\tQSTART\tQEND\tIDEN\n")
    with open(query_filename) as h:
        for r in SeqIO.parse(h, 'fasta'):
            ambig = r.seq.count('N') + r.seq.count('?')
            blastout = handle.next()
            if len(blastout.alignments) == 0: # no match was found!
                f.write("{id}\t{count}\t{len}\t{ambig}\tNA\tNA\tNA\n".format(\
                id=r.id, count=abundance[r.id], len=len(r.seq), ambig=ambig))
            else:
                hsp = blastout.alignments[0].hsps[0]
                f.write("{id}\t{count}\t{len}\t{ambig}\t{qs}\t{qe}\t{iden}\n".format(\
                id=r.id, len=len(r.seq), qs=hsp.query_start, qe=hsp.query_end,\
                iden=hsp.identities, count=abundance[r.id], ambig=ambig))
    f.close()
    
if __name__ == "__main__":
    parse_blast_xml(*sys.argv[1:])
