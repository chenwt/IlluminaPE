
def iter_qiime_otu(filename):
    """
    Qiime's _otu.txt should have format:
    <cluster id> \t <tab-deilimited seq IDs>
    """
    with open(filename) as f:
        for line in f:
            raw = line.strip().split('\t')
            yield raw[0], raw[1:]
            
def main(from_filename, to_filename, output_filename, joinby):
    from collections import defaultdict
    if from_filename == "NONE":
        abundance = defaultdict(lambda:1)
    else:
        abundance = {}
        with open(from_filename) as f:
            for line in f:
                _id, _count = line.strip().split('\t')
                abundance[_id] = int(_count)
        
    otu_count = {}    
    for otu_id,v in iter_qiime_otu(to_filename):
        if otu_id not in otu_count:
            otu_count[otu_id] = 0
        for _id in v:
            otu_count[otu_id] += abundance[_id]
    
    with open(output_filename, 'w') as f:
        stuff = otu_count.items()
        stuff.sort(key=lambda x:x[1], reverse=True) # sort by abundance!
        for _id, _count in stuff:
            f.write("{0}\t{1}\n".format(_id, _count))

if __name__=="__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combined chained OTU lists")
    parser.add_argument("--from", dest="from_file", required=True, help="abundance filename")
    parser.add_argument("--to", dest="to_file", required=True, help="otu filename")
    parser.add_argument("-o", dest="out_file", required=True, help="output filename")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--comma", action="store_true")
    group.add_argument("--tab", action="store_false")

    args = parser.parse_args()

    main(args.from_file, args.to_file, args.out_file, "," if args.comma else "\t")
