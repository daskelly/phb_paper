#!/usr/bin/env python

"""Write a BED file of windows across the genome"""
import sys, argparse, re

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

    for faObj in readFastaByChrom(args.genome):
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

class Fasta:
    def __init__(self, name, seq=None):
        self.name = name
        self.seq = seq
    def __repr__(self):
        return ">%s\n%s".format(self.name, self.seq)

def readFastaByChrom(filename):
    """Read a file in fasta format. Return an iterator of Fasta objects where
    each object is a fasta sequence in the file. Uses a generator for
    memory-efficient reading in of the sequences"""
    f = open(filename)
    file_pos = f.tell()
    line = f.readline()
    while line.startswith('#'):
        file_pos = f.tell()
        line = f.readline()
    f.seek(file_pos)
    while True:
        yield parseFastaMatch(f)
    f.close()

def parseFastaMatch(filehandle):
    "Read in lines representing a single fasta sequence"
    readName = filehandle.readline().rstrip('\n')[1:]
    if readName == '':
        raise StopIteration
    fastaObj = Fasta(name=readName)
    filePos = filehandle.tell()
    line = filehandle.readline().rstrip('\n')
    seq = [line]
    while True:
        filePos = filehandle.tell()
        line = filehandle.readline().rstrip('\n')
        if line == '': break
        elif line.startswith('>'):
            filehandle.seek(filePos)
            break
        else:
            seq.append(line)
    fastaObj.seq = ''.join(seq)
    return fastaObj

if __name__ == "__main__":
    sys.exit(main())
