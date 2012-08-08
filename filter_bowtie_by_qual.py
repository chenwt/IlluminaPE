import time
from miscBowTie import BowTieReader, BowTieWriter

def filter_low_qual_seqs(gz_filename, phred_offset, phred_cutoff):
    """
    Takes a BowTie-style gzipped file (ex: .aligned.composite.gz)
    and retain only seqs that have every base phred >= <cutoff>

    Outputs: .phred<cutoff>_passed for both files
    """
    assert phred_offset >= 0
    assert phred_cutoff >= 0
    bad = 0
    good = 0
    start_t = time.time()
    f = BowTieWriter(gz_filename + ".phred{0}_passed".format(phred_cutoff), 'w')
    for r in BowTieReader(gz_filename, False):
        if all(ord(x)-phred_offset >= phred_cutoff for x in r['qual']):
            good += 1
            f.write(r)
        else:
            bad += 1

    with open(gz_filename + ".phred{0}_passed.log".format(phred_cutoff), 'w') as f:
        f.write("Running filter_low_qual_seq took {0} secs\n".format(time.time()-start_t))
        f.write("Input: " + gz_filename + '\n')
        f.write("PhredCutoff: " + str(phred_cutoff) + '\n')
        f.write("RemovedDueToLowQual: " + str(bad) + '\n')
        f.write("RemainingTotal: " + str(good) + '\n')


if __name__ == "__main__":
#    filter_low_count_low_qual_seqs(sys.argv[1])
    import argparse

    parser = argparse.ArgumentParser(description='Discard reads with 1 or more base with phred score < T')
    parser.add_argument("-i", dest="input", required=True, help="Input file (BowTie format, gzipped)")
    parser.add_argument("-t", dest="threshold", required=True, type=int,\
            help="Phred score cutoff")
    parser.add_argument("-o", dest="offset", default=33, type=int, help="Phred score offset (default 33)")
    args = parser.parse_args()

    filter_low_qual_seqs(args.input, args.offset, args.threshold)
