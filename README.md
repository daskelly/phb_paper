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
* hmmlearn python package. We used version 0.2.0.
* If you would like to convert strain-level breakpoints to reference
 coordinates, you will a tool like liftOver https://genome.ucsc.edu/util.html
 or CrossMap http://crossmap.sourceforge.net/.

### Notes

* The code here is one example of how this analysis could be carried
 out.
* If you are attempting to carry out this pipeline, set the environment
 variable $PICARD_DIR to the directory where Picard jarfiles are located.

### Example

An example of using this pipeline:

```bash
# Starting with individual genomes ref.fa ind1.fa ind2.fa ... indN.fa
cat *.fa > pool_genome.fa
python phb_mapping_workflow.py reads_1.fastq reads_2.fastq pool_genome.fa
# creates reads.bam
python tile_across_genome.py --windowSize 7500 pool_genome.fa > windows.bed
for STRAIN in $(awk -F. '{ print $1; }' windows.bed | sort | uniq)
do
  grep $STRAIN windows.bed > ${STRAIN}_windows.bed
  bedtools coverage -a ${STRAIN}_windows.bed -b reads.bam | cut -f-4 > ${STRAIN}_coverage.bed
  python run_hmm_on_coverage.py ${STRAIN}_coverage.bed > ${STRAIN}_ancestry_probs.bed
done
```
