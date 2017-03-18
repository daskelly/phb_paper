#!/usr/bin/env python

"""Start with a file of ancestry probabilities for each spore.
This can be created by run_hmm_on_coverage.py.
Call breakpoints for where recombination events occurred.

The ancestry probabilities for each bin are based on tiling across
the genome and omitting any bins with >half N.
In this script I'll fill in those missing bins, and fill in the probabilities
with the mean of the probabilities before and after each missing bin.
"""
import sys, argparse, math, os
from itertools import groupby
import numpy as np
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("ancestry_probs", metavar="ancestry_probs.bed")
    parser.add_argument("--bp_filter", default=15001, type=int, metavar="N", \
        help="Filter out any switches in ancestry <= this many base pairs in size")
    parser.add_argument("--outdir", default='.', metavar='.')
    parser.add_argument("--print_all_breakpoints", action="store_true")
    parser.add_argument("--print_spore_genomic_segments", action="store_true", \
        help="print a BED file with BED intervals giving ancestry of each locus in each spore")
    parser.add_argument("--print_spore_breakpoints", action="store_true", \
        help="print a BED file giving the locations of spore breakpoints")
    parser.add_argument('--print_spore_breakpoints_with_ancestry', action="store_true", \
        help="print a BED file with locations of spore breakpoints and left/right ancestry")
    args = parser.parse_args()
    assert args.ancestry_probs.endswith('.ancestry_probs.bed')
    args.bp_filter = np.array([args.bp_filter])   # needed for boolean array in filter_small()

    dat = pd.read_table(args.ancestry_probs, comment='#', header=None,
        names=['chrom', 'start', 'end', 'spore1', 'spore2', 'spore3', 'spore4'])
    dat.set_index(['chrom', 'start', 'end'], inplace=True)

    spore_bkpts = dict((num, list()) for num in [1, 2, 3, 4])
    for chrom, df in dat.groupby(level='chrom'):
        if '2micron' in chrom: continue
        if 'chrM' in chrom: continue
        start = df.index.get_level_values(1).values
        end = df.index.get_level_values(2).values
        if np.all(start[1:] == end[:-1]):
            df2 = df
        else:
            df2 = fill_missing_bins(df)
            start = df2.index.get_level_values(1).values
            end = df2.index.get_level_values(2).values
            assert np.all(start[1:] == end[:-1])
        left, right, ancestry = probs2intervals(df2)
        for num in [1, 2, 3, 4]:
            # get BED intervals.
            j = num - 1
            leftf, rightf, ancestryf = filter_small_intervals(left[j], right[j], ancestry[j], args)
            assert len(leftf) == len(rightf) == len(ancestryf)
            assert np.all(ancestryf >= 0), "ancestryf was {}".format(ancestryf)
            n = len(leftf)
            spore_bkpts[num].append(pd.DataFrame({'chrom':np.repeat(chrom, n),
                'start':leftf, 'end':rightf, 'ancestry':ancestryf}))

    for sporenum in [1, 2, 3, 4]:
        df = pd.concat(spore_bkpts[sporenum])
        df = df[['chrom', 'start', 'end', 'ancestry']]
        df['call'] = 'ref'
        df.ix[df['ancestry'] == 1, 'call'] = 'nonref'
        df.drop('ancestry', axis=1, inplace=True)
        if args.print_spore_genomic_segments:
            outfile = "%s/%s.spore%d_segments.bed" % (args.outdir,
                os.path.basename(args.ancestry_probs)[:-19], sporenum)
            print >> sys.stderr, "writing to %s" % outfile
            df.to_csv(outfile, sep="\t", header=False, index=False)

    if args.print_spore_breakpoints:
        for sporenum in [1, 2, 3, 4]:
            outfile = "%s/%s.spore%d_breakpoints.bed" % (args.outdir,
                os.path.basename(args.ancestry_probs)[:-19], sporenum)
            f = open(outfile, 'w')
            df = pd.concat(spore_bkpts[sporenum])
            for chrom, df2 in df.groupby('chrom'):
                assert np.all(df2['start'].values[1:] == df2['end'].values[:-1])
                nums = df2['end'].values[:-1]
                newdf = pd.DataFrame({'chrom':np.repeat(chrom, len(nums)), 'start':nums,
                    'end':nums+1, 'y':np.repeat(1, len(nums))})
                newdf = newdf[['chrom', 'start', 'end', 'y']]
                newdf.to_csv(f, sep="\t", header=False, index=False, mode='a')
            f.close()

    # print a file giving each spore breakpoint but also including the
    # left/right ancestry on either side of the breakpoint.
    if args.print_spore_breakpoints_with_ancestry:
        for sporenum in [1, 2, 3, 4]:
            outfile = "%s/%s.spore%d_breakpoints_ancestry.bed" % (args.outdir,
                os.path.basename(args.ancestry_probs)[:-19], sporenum)
            f = open(outfile, 'w')
            f.write('#chrom\tstart\tend\tancestry_left\tancestry_right\n')
            df = pd.concat(spore_bkpts[sporenum])
            df['call'] = 'ref'
            df.ix[df['ancestry'] == 1, 'call'] = 'nonref'
            for chrom, df2 in df.groupby('chrom'):
                df3 = df2.sort_values(by=['start', 'end'])
                assert is_df_sorted(df3, 'start')  # sanity check
                assert np.all(df3['start'].values[1:] == df3['end'].values[:-1])
                if df3.shape[0] == 1:
                    # no breakpoints on this chrom
                    newdf = pd.DataFrame({'chrom':np.array([chrom]),
                        'start':np.array([-17]), 'end':np.array([-17]),
                        'ancestryLeft':df3['call'].values,
                        'ancestryRight':df3['call'].values})
                else:
                    # at least one breakpoint
                    nums = df3['end'].values[:-1]
                    newdf = pd.DataFrame({'chrom':np.repeat(chrom, len(nums)), 'start':nums,
                        'end':nums+1, 'ancestryLeft':df3['call'].values[:-1],
                        'ancestryRight':df3['call'].values[1:]})
                newdf = newdf[['chrom', 'start', 'end', 'ancestryLeft', 'ancestryRight']]
                newdf.to_csv(f, sep="\t", header=False, index=False, mode='a')
            f.close()

    # finally, print a file with ALL breakpoints:
    if args.print_all_breakpoints:
        alldat = pd.concat(spore_bkpts[1] + spore_bkpts[2] + spore_bkpts[3] + spore_bkpts[4])
        f = open("%s/%s.all_breakpoints.bed" % (args.outdir,
            os.path.basename(args.ancestry_probs)[:-19]), 'w')
        for chrom, df in alldat.groupby('chrom'):
            nums = np.unique(df['start'].values)
            nums.sort()
            assert 0 in nums
            if len(nums) == 1:
                continue
            assert len(nums) > 1
            nums = nums[1:]
            newdf = pd.DataFrame({'chrom':np.repeat(chrom, len(nums)), 'start':nums, 'end':nums+1,
                'y':np.repeat(1, len(nums))})
            newdf = newdf[['chrom', 'start', 'end', 'y']]
            newdf.to_csv(f, sep="\t", header=False, index=False, mode='a')
        f.close()

    # There is still one small problem with this script.
    # It filters by bin size going spore-by-spore, but this really should be done
    # considering all four spores simultaneously. For example,
    # an ancestry switch from spore 4 to spore 2 at, and a very nearby
    # ancestry switch from spore 2 to spore 3 at 380000, has a the really small
    # genomic segment in spore 2 that should be filtered out. However, we don't
    # see it if we look at spores individually.

def fill_missing_bins(dat):
    """Some bins are missing due to a high percentage of Ns.
    Fill them in by averaging the probs in the bins before and after."""
    probs = dat.reset_index()
    start = probs['start'].values
    end = probs['end'].values
    sizes = end - start
    windowsize = sizes[0]
    assert np.all(sizes == windowsize)
    gaps = start[1:] - start[:-1]
    wrong_gap = np.where(gaps != windowsize)[0]
    spores = ['spore%d' % i for i in range(1, 5)]
    for gap in wrong_gap:
        missing_start = probs.ix[gap, 'end']
        missing_end = probs.ix[gap+1, 'start']
        assert missing_end > missing_start
        new_starts = np.array(range(missing_start, missing_end, windowsize))
        new_ends = np.array([min(missing_end, i+windowsize) for i in new_starts])
        new_probs = probs.ix[gap:gap+1, spores].mean().values
        n = len(new_starts)
        missing_rows_dict = {'chrom':np.repeat(probs.ix[0, 'chrom'], n),
            'start':new_starts, 'end':new_ends, 'spore1':np.repeat(new_probs[0], n),
            'spore2':np.repeat(new_probs[1], n), 'spore3':np.repeat(new_probs[2], n),
            'spore4':np.repeat(new_probs[3], n)}
        missing_rows = pd.DataFrame(missing_rows_dict,
            columns=['chrom', 'start', 'end'] + spores)
        probs = probs.append(missing_rows, ignore_index=True)
    probs.sort_values(by=['start', 'end'], inplace=True)
    probs.set_index(['chrom', 'start', 'end'], inplace=True)
    return probs

def is_df_sorted(df, colname):
    return (np.diff(df[colname]) > 0).all()

def probs2intervals(probs):
    """Given ancestry probabilities, convert to segments of 0/1 ancestry.
    This is complicated because of places where the ancestry prob is close to 0.5,
    for stretches of the genome. Need to make sure the two complementary segregants
    match in terms of where I am calling the breakpoint.
    Returns (left coord of bin, right coord of bin, ancestry call)
    """
    probs = probs.reset_index()
    #
    ints = np.round(probs)
    sporecols = ['spore1', 'spore2', 'spore3', 'spore4']
    sums = ints[sporecols].sum(axis=1)
    problem_rows = sums != 2
    if problem_rows.sum() > 0:
        assert np.all((sums[problem_rows] == 3) | (sums[problem_rows] == 1))
        lens, problem = rle(problem_rows)
        if problem[0] == True:
            # start with a problematic rows. Just set the two max probs to 1, the other two to 0
            for i in range(lens[0]):
                current_probs = probs.ix[i, sporecols].values
                oo = np.argsort(current_probs)
                ints.ix[i, sporecols] = np.array([0, 0, 1, 1])[oo]
        if problem[-1] == True:
            # end with problematic rows. Again set the two max probs to 1, the other two to 0
            cslen = np.cumsum(lens)
            for i in range(cslen[-2], cslen[-1]):
                current_probs = probs.ix[i, sporecols].values
                oo = np.argsort(current_probs)
                ints.ix[i, sporecols] = np.array([0, 0, 1, 1])[oo]
        sums = ints[sporecols].sum(axis=1)
        problem_rows = sums != 2
    #
    if problem_rows.sum() > 0:
        lens, problem = rle(problem_rows)
        assert problem[0] == False and problem[-1] == False
        lens2, _ = rle(sums)
        ## Iterate through problem rows in reverse order
        group_idx_start = np.cumsum(lens) - lens
        group_idx_end = np.cumsum(lens)
        for i, j, prob in zip(group_idx_start, group_idx_end, problem)[::-1]:
            # iterate through in reverse because I'm expanding the data.frame as I go
            if prob == False: continue  # not a problem row
            ## retrieve the (rounded) ints before/after the problem row
            before = ints.iloc[i-1]
            this = ints.iloc[i:j]
            after = ints.iloc[j]
            start_coord, end_coord = before['end'], after['start']
            assert this['start'].values[0] == start_coord
            assert this['end'].values[-1] == end_coord
            ## the problem segregants are the ones where *before ancestry* != *after ancestry*
            #  should be only two of the four ...
            beforevals = before[sporecols].values
            aftervals = after[sporecols].values
            if np.all(beforevals == aftervals):
                # essentially a rounding error
                ints.ix[i:j, sporecols] = beforevals # (same as aftervals)
                continue
            prob_segregants = np.where(beforevals != aftervals)[0]
            ok_segregants = np.where(beforevals == aftervals)[0]
            if len(prob_segregants) == 4 and len(ok_segregants) == 0:
                # two very close-by breakpoints each in 2/4 spores
                # I'm going to cop out and not attempt a sophisticated
                # system for resolving this. I'll just take the max two
                # probabilities as 1 and min two probs as 0, for each bin.
                current_probs = probs.ix[i:j, sporecols].values
                for rnum in range(len(current_probs)):
                    oo = np.argsort(current_probs[rnum])
                    ints.ix[i+rnum, sporecols] = np.array([0, 0, 1, 1])[oo]
                continue
            assert len(prob_segregants) == 2 and len(ok_segregants) == 2
            ## if this.shape[0] % 2 == 0:
            ##     # even number of rows in the problematic region. Just split in half
            ## else:
            ##     # odd number of rows in this region. I'll need to split individual bins
            # Don't need to use the above strategy. Just split this bin down the middle.
            before_df = ints.iloc[0:i]
            after_df = ints.iloc[j:]
            this_df = pd.concat((this.iloc[0:1], this.iloc[0:1]),
                ignore_index=True).reset_index(drop=True)
            midpt = np.int(np.round((start_coord + end_coord)/2.0))
            this_df.ix[0, 'start'] = start_coord
            this_df.ix[0, 'end'] = midpt
            this_df.ix[1, 'start'] = midpt
            this_df.ix[1, 'end'] = end_coord
            this_df.ix[0, sporecols] = beforevals
            this_df.ix[1, sporecols] = aftervals
            assert np.all(this_df[sporecols].sum(axis=1) == 2)
            ints = pd.concat((before_df, this_df, after_df), ignore_index=True)
    #
    starts = ints['start'].values
    ends = ints['end'].values
    left, right, ancestry = [], [], []
    for i in [1, 2, 3, 4]:
        lens, anc = rle(ints['spore%d' % i])
        left.append(starts[np.cumsum(np.concatenate((np.array([0]), lens[:-1])))])
        right.append(ends[np.cumsum(lens)-1])
        ancestry.append(anc)
    return (left, right, ancestry)

def filter_small_intervals(left, right, ancestry, args):
    intervalsizes = right - left
    assert np.all(intervalsizes > 0)
    smallbins = intervalsizes < args.bp_filter
    if len(smallbins) == 1:
        assert smallbins[0] is not True
        return (left, right, ancestry)
    #
    # Filter out small intervals at end and start of chrom
    if smallbins[-1] == True:
        last_big_idx = np.max(np.where(~smallbins))
        left = left[:last_big_idx + 1]
        right = np.concatenate((right[:last_big_idx], np.array([right[-1]])))
        ancestry = ancestry[:last_big_idx + 1]
        intervalsizes = right - left
        smallbins = intervalsizes < args.bp_filter
    if smallbins[0] == True:
        first_big_idx = np.min(np.where(~smallbins))
        left = np.concatenate((left[:1], left[first_big_idx+1:]))
        right = right[first_big_idx:]
        ancestry = ancestry[first_big_idx:]
        intervalsizes = right - left
        smallbins = intervalsizes < args.bp_filter
    #
    # Filter out small intervals in middle of chrom:
    assert len(left) == len(right) == len(ancestry) == len(smallbins)
    assert np.all(left[1:] == right[:-1])
    if smallbins.sum() == 0:
        return (left, right, ancestry)
    # First, convert from RLE back to the sequence
    # RLE is (lens, vals)
    # (right - left, ancestry)
    # To reconstitute sequence use np.repeat(vals, lens)
    ancestry[smallbins] = -17
    seq = np.repeat(ancestry, right - left)
    # now back to RLE:
    lens, ancestry = rle(seq)
    positions = np.concatenate((np.array([0]), np.cumsum(lens)))
    left = positions[:-1]
    right = positions[1:]
    intervalsizes = right - left
    assert np.all(intervalsizes > 0)
    assert ancestry[0] >= 0 and ancestry[-1] >= 0   # i.e. not -17
    # small bins are no longer defined by interval size, because I could have combined
    # consecutive small bins that jointly exceed args.bp_filter in size.
    # Now we look for ancestry == -17 to indicate a "small" bin
    for i in range(len(ancestry)-2, 0, -1):     # can skip first and last bin
        # start from end so I can play with the lists
        if ancestry[i] != -17: continue
        before_ancestry = ancestry[i-1]
        after_ancestry = ancestry[i+1]
        if before_ancestry == after_ancestry:
            # combine bins i-1, i, and i+1 into one
            left = np.concatenate((left[:i], left[i+2:]))
            right = np.concatenate((right[:i-1], right[i+1:]))
            ancestry = np.concatenate((ancestry[:i-1], ancestry[i+1:]))
        else:
            binsize = right[i] - left[i]
            l1, l2 = splitbin(binsize)
            # change bins i-1, i, and i+1 from three to two bins
            # extend bin i-1 to the right by l1
            # extend bin i+1 to the left by l2
            # remove bin i
            left = np.concatenate((left[:i], np.array([left[i+1]-l2]), left[i+2:]))
            right = np.concatenate((right[:i-1], np.array([right[i-1]+l1]), right[i+1:]))
            ancestry = np.concatenate((ancestry[:i], ancestry[i+1:]))
    return (left, right, ancestry)

def splitbin(binsize):
    fb = float(binsize)
    return (int(math.floor(fb/2)), int(math.ceil(fb/2)))

def rle(sequence):
    """run-length encoding. To reconstitute sequence you can
       call np.repeat(vals, lens)"""
    lens, vals = list(), list()
    for key, val in groupby(sequence):
        lens.append(len(list(val)))
        vals.append(key)
    return (np.array(lens), np.array(vals))

if __name__ == "__main__":
    sys.exit(main())
