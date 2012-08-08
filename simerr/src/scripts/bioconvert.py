import os, sys
from PEassembly.miscBowTie import BowTieReader, BowTieWriter, FastqWriter, FastqReader
from Bio.Seq import Seq

def fastq2fasta(input):
    output = input + '.fasta'
    with open(output, 'w') as f:
        for r in FastqReader(input):
            f.write(">{id}\n{seq}\n".format(id=r['ID'], seq=r['seq']))

def bowtie2fasta(input, name, trimid):
    output = input + '.fasta'
    f = open(output, 'w')
    if name is not None:
        fmap = open(output + '.map', 'w')
        for r in BowTieReader(input, False):
            id = r['ID'].split()[0] if trimid else r['ID']
            f.write(">{name}_{id}\n{seq}\n".format(id=id, seq=r['seq'], name=name))
            fmap.write("{name}_{id}\t{name}\n".format(name=name, id=id))
    else:
        for r in BowTieReader(input, False):
            id = r['ID'].split()[0] if trimid else r['ID']
            f.write(">{id}\n{seq}\n".format(id=id, seq=r['seq']))
    f.close() 
    if name is not None:
        fmap.close()
    
def bowtie2fastq(input, revcomp=False):
    output = input + '.fq'
    f = FastqWriter(output)
    for r in BowTieReader(input, False):
        if revcomp:
            r['seq'] = Seq(r['seq']).reverse_complement().tostring()
            r['qual'] = r['qual'][::-1]
        f.write(r)
    f.close()       
    
if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description="""
                            converting various bio formats \
                            b2f: bowtie->fasta,
                            b2fq: bowtie->fastq""")
    parser.add_argument('type', choices=["b2f", "b2fq", "fq2f"])
    parser.add_argument('-i', dest='input', required=True, type=str, help="input filename")
    parser.add_argument("--trimid", action="store_true", default=False, help="trim ID anything after 1st space. For b2f only.")
    parser.add_argument("--name", required=False, type=str, help="extra <name> to put in front of IDs. For b2f only. Will also output a .fasta.map file with <name>_<oldid> \t <name>")
    parser.add_argument("--revcomp", default=False, action="store_true", help="rev comp. For b2fq only")
    args = parser.parse_args()
    
    if args.type == 'b2f':
        bowtie2fasta(args.input, args.name, args.trimid)
    elif args.type == 'b2fq':
        bowtie2fastq(args.input, args.revcomp)
    elif args.type == 'fq2f':
        fastq2fasta(args.input)
        
