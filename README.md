## Private Haplotype Barcoding

This repository contains an example pipeline showing how to
implement private haplotype barcoding, as described in
Skelly et al. 2017 https://doi.org/10.1101/116582.
The data generated in the paper is available through NCBI's
Short Read Archive at
https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?study=SRP101897.

### Pipeline

The overall steps of the pipeline are:

* Create a *pool* genome consisting of all chromosomes of the
 sequenced individuals who have contributed ancestry to the pool.
 It is best to use *de novo*-assembled genomes for this step.
 However, an alternative approach if *de novo*-assembled genomes
 are not available would be to use a VCF of variation to create
 individualized genomes, and map to them. Examples of tools to do
 this include Seqnature http://dx.doi.org/10.1534/genetics.114.165886
 and g2gtools https://github.com/churchill-lab/g2gtools.
* Map short reads. You may use any method of mapping reads, but need to
 eventually end up with a BAM file of reads that map uniquely with high
 quality.
* Create a BED file of windows tiling across the genome of each
 strain contributing ancestry to your pool.
* For each strain, compute BAM read coverage across these windows.
* Run the HMM to obtain ancestry probabilities for each window.
* (If desired) call recombination breakpoints in the genomes of individual
 offspring.

### Requirements

* bedtools http://bedtools.readthedocs.io/en/latest/. We used version 2.25.0.
* bwa http://bio-bwa.sourceforge.net/. We used v0.7.12.
* Picard https://broadinstitute.github.io/picard/. We used version 1.101.
* samtools https://github.com/samtools/samtools. We used v0.1.19.
* python packages (We used python 2.7.13):
  * hmmlearn. We used version 0.2.0.
  * numpy. We used version 1.10.4.
  * pandas. We used version 0.19.2.
  * scipy. We used version 0.17.0.
* R packages (We used R 3.1.0):
  * data.table. We used version 1.10.0.
  * dplyr. We used version 0.4.1.
  * mmand. We used version 1.2.0.
  * assertthat. We used version 0.1.
* DWGSIM https://github.com/nh13/DWGSIM
* If you would like to convert strain-level breakpoints to reference
 coordinates, you will a tool like liftOver https://genome.ucsc.edu/util.html
 or CrossMap http://crossmap.sourceforge.net/.

### Notes

* The code here is one example of how this analysis could be carried
 out.
* If you are attempting to carry out this pipeline, set the environment
 variable $PICARD_DIR to the directory where Picard jarfiles are located.
* Chromosomes in the pooled genome should be named *strain*.*chrom*.
(This is done in the script `combineGenomes.fa` below).

### Example

An example of using this pipeline:

```bash
# Starting with individual genomes ref.fa ind1.fa ind2.fa ... indN.fa
python combineGenomes.fa ref.fa ind*.fa > pool_genome.fa
python phb_mapping_workflow.py poolA.R1.fastq poolA.R2.fastq pool_genome.fa
# --> maps paired-end reads, creating poolA.bam
# repeat for each pool.

python tile_across_genome.py --windowSize 7500 pool_genome.fa > windows.bed
for STRAIN in $(awk -F. '{ print $1; }' windows.bed | sort | uniq)
do
  grep $STRAIN windows.bed > ${STRAIN}_windows.bed
done
```

Using the single-sample HMM:
```bash
DIR=single_sample_HMM
for STRAIN in $(awk -F. '{ print $1; }' windows.bed | sort | uniq)
do
  bedtools coverage -a ${STRAIN}_windows.bed -b readsA.bam | cut -f-4 > ${STRAIN}_coverage.bed

  # smooth the data with a median filter
  Rscript $DIR/smooth_coverage.R ${STRAIN}_coverage.bed ${STRAIN}_coverage.filtered.bed

  # run the HMM:
  python $DIR/run_single_hmm_on_coverage.py ${STRAIN}_coverage.filtered.bed \
    ${STRAIN}_binAncestry_probs.bed ${STRAIN}_binNoAncestry_probs.bed \
    > ${STRAIN}_ancestry_probs.bed

  # filter out small bins
  Rscript $DIR/smooth_ancestry_probs.R ${STRAIN}_ancestry_probs.bed \
    ${STRAIN}_ancestry_probs.smoothed.bed

  # call breakpoints
  Rscript $DIR/call_breakpoints_from_smoothed.R ${STRAIN}_ancestry_probs.smoothed.bed \
    ${STRAIN}_breakpoints
done
```

Using the HMM that incorporates tetrad information:
```bash
# Let's say we have the four spores from a tetrad in each of pools A, B, C, and D
for STRAIN in $(awk -F. '{ print $1; }' windows.bed | sort | uniq)
do
  for POOL in A B C D
  do
    bedtools coverage -a ${STRAIN}_windows.bed -b reads${POOL}.bam | \
      cut -f-4 > ${STRAIN}_coverage_pool${POOL}.bed
  done

  # run tetrad HMM:
  python run_tetrad_hmm_on_coverage.py ${STRAIN}_coverage_pool[ABCD].bed \
    > ${STRAIN}_ancestry_probs.bed

  # filter out small bins and call breakpoints
  # run with -h to see all options
  python call_breakpoints.py ${STRAIN}_ancestry_probs.bed
done
```
