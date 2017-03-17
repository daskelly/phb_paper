#!/usr/bin/env python

"""Combine a set of individual fasta files into one file.
   The fasta file should be named strainName.fa.
"""
import sys, argparse, os, seqs

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafiles", nargs="+", metavar="file.fa", \
        help="fasta files to combine into a single large file")
    args = parser.parse_args()

    for filename in args.fastafiles:
        print >> sys.stderr, "\rworking on %s                           " \
            % os.path.basename(filename),
        strain = '.'.join(os.path.basename(filename).split('.')[:-1])
        for faObj in seqs.readFastaByChrom(filename):
            if faObj.name.startswith(strain):
                name = faObj.name
            else:
                name = "%s.%s" % (strain, faObj.name)
            print ">%s" % name
            seqs.pretty_print(faObj.seq)
    print >> sys.stderr, "\rDone!                                       "

if __name__ == "__main__":
    sys.exit(main())
