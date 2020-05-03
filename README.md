# stacksPairwise

Calculate pairwise sequence divergence between samples from [Stacks](http://catchenlab.life.illinois.edu/stacks/) RAD genotyping output. Simply supply the path to a `samples.fa` file generated from a Stacks run and a text file with individuals to calculate divergence between.

Pairwise divergence is estimated seperately per-RAD tag. It is calculated as the average divergence between haplotypes from different samples. Within sample diversity (i.e., heterozygosity) is also calculated. 



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
usage: stacksPairwise [-h] [-names] [-o] samples

Calculate pairwise divergence (pairwise pi) from Stacks `samples.fa` output
fle

positional arguments:
  samples            Path to `samples.fa` file (from Stacks output)

optional arguments:
  -h, --help         show this help message and exit
  -names , --names   Names of samples to analyze. Either a text file or comma
                     seperated list.
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
A text file containing the samples to be compared. One name per line.

## Outputs

Four files will be written:

1. `.diffs.csv`, containing the average number of nucleotide differences for each RAD tag
2. `.sites.csv`, containing the total number of sites (per locus) for which each pairwise comparision has data
3. `.estimates.csv`, containing the estimated per-site divergence (per locus) for each pairwise comparison
4. `.summary.txt`, containing genome-wide summaries

