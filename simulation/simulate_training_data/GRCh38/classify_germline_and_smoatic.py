import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="Demo of argparse")
    parser.add_argument("--summary", type=str, dest="summary", default="", help="summary file")
    parser.add_argument("--input", type=str, dest="input", default="", help="input BED file")
    parser.add_argument("--output", type=str, dest="output", default="", help="output BED file")
    parser.add_argument("--freq", type=float, dest="freq", default="", help="minimum frequency")
    return parser


def load_germ(summary_fn, min_freq):
    germs = {}
    for l in open(summary_fn, 'r'):
        if l:
            l = l.split()
            chrom = l[0]
            id = l[3]
            freq = float(l[6])
            if freq < min_freq:
                continue
            if chrom in germs:
                germs[chrom][id] = str(freq)
            else:
                germs[chrom] = {}
                germs[chrom][id] = str(freq)
    
    return germs


def filter_germ(input_fn, output_fn, germs, min_freq):
    out = open(output_fn, 'a')
    for l in open(input_fn, 'r'):
        if l:
            l = l.split()
            try:
                freq = germs[l[0]][l[3]]
                l[4] = freq
                out.write("\t".join(l) + '\n')
            except:
                l[4] = str(min_freq)
                out.write("\t".join(l) + '\n')
    
    out.close()


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    germs = load_germ(args.summary, args.freq)
    filter_germ(args.input, args.output, germs, args.freq)
