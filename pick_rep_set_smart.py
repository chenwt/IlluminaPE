import os, sys
from Bio import SeqIO

def pick_rep_set(fasta_filename, otu_filename, abundance_filename, output_filename):
    abundance = {}
    with open(abundance_filename) as f:
        for line in f:
            _id, _count = line.strip().split('\t')
            abundance[_id] = int(_count)
            
    best_pick = {} # rep id --> cid
    with open(otu_filename) as f:
        for line in f:
            cid, rest = line.strip().split('\t', 1)
            cur_ab = 0
            cur_id = None
            for _id in rest.split('\t'):
                if abundance[_id] >= cur_ab:
                    cur_ab = abundance[_id]
                    cur_id = _id
            best_pick[cur_id] = cid
            
    f = open(output_filename, 'w')
    for r in SeqIO.parse(open(fasta_filename), 'fasta'):
        if r.id in best_pick:
            f.write(">{0} {1}\n{2}\n".format(best_pick[r.id], r.id, r.seq))
    f.close()
    
if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combined chained OTU lists")
    parser.add_argument("-i", dest="fasta_filename", required=True, help="fasta filename")
    parser.add_argument("-f", dest="otu_filename", required=True, help="otu filename")
    parser.add_argument("-a", dest="abundance_filename", required=True, help="abundance filename")
    parser.add_argument("-o", dest="output_filename", required=True, help="output filename")


    args = parser.parse_args()

    pick_rep_set(args.fasta_filename, args.otu_filename, args.abundance_filename, args.output_filename)