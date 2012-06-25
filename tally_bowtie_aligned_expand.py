from miscBowTie import BowTieReader

def main(input_filename, abundance_filename, output):
    abundance = {}
    with open(abundance_filename) as f:
        for line in f:
            _id, _count = line.strip().split('\t')
            abundance[_id] = int(_count)
            
    total = sum(abundance.itervalues())
    aligned = 0
    for r1, r2 in BowTieReader(input_filename, True):
        realid = r1['ID'][:r1['ID'].find('/')]
        aligned += abundance[realid]

    with open(output, 'w') as f:
        p = aligned*100./total
        f.write("# reads processed: {0}\n".format(total))
        f.write("# reads with at least one reported alignment: {0} ({1:.2f}%)\n".format(aligned,p))
        f.write("# reads that failed to align: {0} ({1:.2f}%)\n".format(total-aligned,100-p))
        f.write("Reported {0} paired-end alignments to 1 output stream(s)\n".format(aligned))
        
        
if __name__ == "__main__":
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", required=True, help="Bowtie aligned output (.gz)")
    parser.add_argument("-a", "--abundance", required=True, help="Abundance filename")
    parser.add_argument("-o", "--output", required=True, help="Output filename")
    
    args = parser.parse_args()
    main(args.input, args.abundance, args.output)