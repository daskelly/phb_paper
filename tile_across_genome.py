#!/usr/bin/env python

"""Write a BED file of windows across the genome"""
import sys, argparse, re, seqs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("genome", metavar="file.fa")
    parser.add_argument("--windowSize", default=1000, type=int, metavar="N", \
        help="window size for features (default %(default)d)")
    parser.add_argument("--addnames", action="store_true", \
        help="Add names for each BED feature (default chrom_start_end)")
    parser.add_argument("--nameprefix", default=None, metavar="prefix", \
        help="prefix for naming of BED features (default None)")
    parser.add_argument("--NignorePct", default=0.5, type=float, metavar="0.5",\
        help="Ignore the window/bin/tile if the percent N is greater than this")
    parser.add_argument("--center_bins", action="store_true", \
        help="miss same amount off of beginning & end of chrom")
    parser.add_argument("--chromname_sub_regex", nargs=2, metavar=('regexA', \
        'stringB'), help="substitute all matches to pattern A with string " + \
        "B in each chromosome name")
    args = parser.parse_args()

    for faObj in seqs.readFastaByChrom(args.genome):
        l = len(faObj.seq)
        n_full_bins = l/args.windowSize
        # center the full bins along the chrom
        # (i.e. leave the same number of bases off the beginning & end of chrom)
        if args.center_bins:
            s = (l - n_full_bins*args.windowSize)/2
        else:
            s = 0
        for i in range(s, l, args.windowSize):
            start, end = i, i + args.windowSize
            if end > l: continue
            binseq = faObj.seq[start:end]
            if binseq.count("N")/float(len(binseq)) > args.NignorePct:
                continue

            if args.chromname_sub_regex is not None:
                pattern, repl = args.chromname_sub_regex
                faObj.name = re.sub(pattern, repl, faObj.name)
            if args.addnames:
                name = ""
                if args.nameprefix is not None:
                    name += "{}_".format(args.nameprefix)
                name += "{}_{}_{}".format(faObj.name, start, end)
                print "{}\t{}\t{}\t{}".format(faObj.name, start, end, name)
            else:
                print "{}\t{}\t{}".format(faObj.name, start, end)


if __name__ == "__main__":
    sys.exit(main())
