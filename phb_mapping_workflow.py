#!/usr/bin/env python

"""Starting with raw reads (fastq), map them and produce a BAM file
that can be input into the PHB HMM.
Outline:
1) Map reads
2) Sort BAM file
3) Mark duplicates in BAM file
4) Filter out reads that are duplicates, non-primary hits, unmapped, or fail QC
"""
import sys, argparse, os, re, subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq1")
    parser.add_argument("fastq2")
    parser.add_argument("ref", \
        help="ref genome in fasta format. Chromosomes be named strain.chrom")
    parser.add_argument("--skip_mark_dups", action="store_true")
    parser.add_argument("--max_ram", type=int, default=10, help="RAM in Gb")
    parser.add_argument("--min_mapq", type=int, default=10, metavar="N", \
        help="filter out reads that map with less than this mapping quality")
    parser.add_argument("--outdir", default=".")
    parser.add_argument("--outprefix", help="default: up to 1st [.-] in fastq")
    parser.add_argument("--save_intermediate", action="store_true")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    for filename in [args.fastq1, args.fastq2, args.ref]:
        if not os.path.exists(filename):
            print >> sys.stderr, "Could not find %s" % filename
            sys.exit(1)

    # Map the paired-end reads to the reference genome. Store in BAM
    home = os.environ["HOME"]
    if args.outprefix is not None:
        outprefix = "%s/%s" % (args.outdir, args.outprefix)
    else:
        outprefix = "%s/%s" % (args.outdir, re.split("[\.-]", os.path.basename(args.fastq1))[0])
    # Check that the outprefix will be valid:
    outprefixdir = dirname(os.path.expanduser(args.outdir))
    assert os.access(outprefixdir, os.W_OK), "could not write to {}".format(outprefixdir)
    if not os.path.exists("%s.bwt" % args.ref):
        print >> sys.stderr, "Please run 'bwa index' on your genome (%s) first" % args.ref
        sys.exit(1)
    cmd = "bwa mem -t %d %s %s %s | samtools view -bS - > %s.bam" % \
        (args.threads, args.ref, args.fastq1, args.fastq2, outprefix)
    system(cmd, message="bwa mapping", verbose=args.verbose)

    # Sort the BAM file. Delete the unsorted BAM
    cmd = "java -Xmx%dg -XX:ParallelGCThreads=%d" % (args.max_ram, args.threads) + \
        " -jar $PICARD_DIR/SortSam.jar INPUT=%s.bam" % outprefix + \
        " OUTPUT=%s.sorted.bam SORT_ORDER=coordinate" % outprefix + \
        " VALIDATION_STRINGENCY=LENIENT"
    system(cmd, message="sort bam", verbose=args.verbose)
    if not args.save_intermediate:
        os.remove("%s.bam" % outprefix)

    # Mark dups in the sorted BAM. Delete the unmarked BAM
    cmd = "java -Xmx%dg -XX:ParallelGCThreads=%d" % (args.max_ram, args.threads) + \
        " -jar $PICARD_DIR/MarkDuplicates.jar INPUT=%s.sorted.bam" % outprefix + \
        " OUTPUT=%s.sorted.markDups.bam METRICS_FILE=/dev/null" % outprefix + \
        " VALIDATION_STRINGENCY=LENIENT"
    if not args.skip_mark_dups:
        system(cmd, message="mark dups", verbose=args.verbose)
        if not args.save_intermediate:
            os.remove("%s.sorted.bam" % outprefix)
    else:
        os.rename("%s.sorted.bam" % outprefix, "%s.sorted.markDups.bam" % outprefix)

    # Filter out duplicates, non-primary hits, unmapped, or fail QC
    # (see https://broadinstitute.github.io/picard/explain-flags.html)
    cmd = "samtools view -hb -F 1796 -q %d %s.sorted.markDups.bam > %s_filtered.bam" \
        % (args.min_mapq, outprefix, outprefix)
    system(cmd, message="filtering", verbose=args.verbose)
    if not args.save_intermediate:
        os.remove("%s.sorted.markDups.bam" % outprefix)

def dirname(path):
    dir = os.path.dirname(path)
    if dir == '':
        dir = '.'
    return dir

def system(call, message, verbose=False):
    if verbose:
        print >> sys.stderr, message
    try:
        retcode = subprocess.call(call, shell=True)
        if retcode != 0:
            sys.exit("ERROR on %s. In %s. Call was %s" % (message, os.getcwd(), call))
    except OSError:
        print >> sys.stderr, "ERROR on %s" % message
        raise

if __name__ == "__main__":
    sys.exit(main())
