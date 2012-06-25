# DSall.alien.summary.rdpgg_good.err_per_position.txt

# SAMPLE  FoR MATCH   POS REFPOS  REAL    OBS
from cPickle import dump
from csv import DictReader
from collections import defaultdict

VALIDNT = ('A', 'T', 'C', 'G')

def main(input_filename, output_filename):
    mutate_freq = defaultdict(lambda: defaultdict(lambda: 0))
    for r in DictReader(open(input_filename), delimiter='\t'):
        if r['OBS'] not in VALIDNT: continue
        mutate_freq[r['REAL']][r['OBS']] += 1
    
    for real in mutate_freq:
        _sum = sum(mutate_freq[real].itervalues()) * 1.
        for mut in mutate_freq[real]:
            mutate_freq[real][mut] /= _sum
        mutate_freq[real] = dict(mutate_freq[real]) 
    
    mutate_freq = dict(mutate_freq)
      
    with open(output_filename, 'w') as f:
        dump(mutate_freq, f)  
         
if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Calculate base mutation frequencies")
    parser.add_argument("-i", dest="input", required=True, help="Input summary text")
    parser.add_argument("-o", dest="output", required=True, help="Output pickle")

    args = parser.parse_args()

    main(args.input, args.output)
    
    
    
