# stacksPairwise

Calculate pairwise sequence divergence between samples from [Stacks](http://catchenlab.life.illinois.edu/stacks/) RAD genotyping output. Simply supply the path to a `samples.fa` file generated from a Stacks run and a text file with individuals to calculate divergence between. By default all pairwise comparisons are made.

Calculating pairwise sequence divergence across the genome from a vcf file output from Stacks (e.g., using vcftools --windowed-pi) will ignore invariant sites and lead to overestimation of the levels of segregating nucleotide diversity. Instead, it is more appropriate to calculate this on a "per-RADtag" basis. 

## Requirements

* Biopython

## Install

To install the program + biopython, run:
```
git clone https://github.com/gibsonMatt/stacks-pairwise.git
cd stacks-pairwise
python setup.py install
```

## Usage

```
usage: stacksPairwise [-h] [-v] [-o] individuals samples

Calculate pairwise sequence divergence (pairwise pi) from Stacks `samples.fa` output
file

positional arguments:
  individuals        Path to text file containing focal samples to compare
                     pairwise. One line per sample.
  samples            Path to `samples.fa` file (from Stacks output)

optional arguments:
  -h, --help         show this help message and exit
  -v, --verbose      Enable debugging messages to be displayed
  -o , --outputdir   Output directory/prefix
```


## Inputs

### `samples.fa` file
To generate a `samples.fa` file, you inclcude the `--fasta-samples` flag to the `populations` program in Stacks. For example: 
```
populations -P ./stacks_output/ --fasta-samples
```

> *Note that no internal filtering is applied, so any filtering should be done prior to calculating pairwise divergence.*


### Sample list
A text file containing the samples to be compared. One sample name per line.

## Output

The program will write a single tab delimited file (defaults to `stacksPairwise.out.tsv`) which contains loci as rows and sample comparisons as columns. See `stacksPairwise.out.tsv` for an example. By default each cell will  report sequence divergence as a fraction (# of diffs/# of sites where both samples have data). Alternatively, by including the `-s` flag, two files will be generated(`stacksPairwise.diffs.out.tsv` and `stacksPairwise.sites.out.tsv`).

|LocusID|Chr|StartPos|Sample1_Sample2|Sample2_Sample1|...|
|-|-|-|-|-|-|
|75|Chr1|303234|0/203|3/203|...|
|76|Chr1|305323|2/126|1/126|...|
|...|...|...|...|...|...|