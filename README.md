# parallelLastz
## Lastz with multi-threads support.

[![Conda](https://anaconda.org/jnarayan81/parallellastz/badges/installer/conda.svg)](https://anaconda.org/jnarayan81/parallellastz)
[![Linux](https://anaconda.org/jnarayan81/parallellastz/badges/platforms.svg)](https://anaconda.org/jnarayan81/parallellastz)
<img src="https://img.shields.io/badge/Perl-Lang-informational?style=flat&logo=perl&logoColor=white&color=2bbc8a" />

Running Lastz (https://github.com/lastz/lastz) in parallel mode. This program is for single computer with multiple core processors.

When the query file format is fasta, you can specify many threads to process it. It can reduce run time linearly, and use almost equal memory as the original lastz program. This is useful when you lastz a big query file to a huge reference like human whole genome sequence.

The program is an extension on the original lastz program which was written by Bob Harris (the LASTZ guy).

parallelLastz can run on Linux and Mac OS.

It run lastz in parallel mode and generate <chr>.lz (tab file) file.

perl parallelLastz.pl -h for more help

```
--qfile|-q	query multifasta/fasta file
--tfile|-t	target genome file
--cfile|-c	config file
--speedup|-s	number of core to use
--length|-l	length below this is ignored
--help|-h	brief help message
```

## Conda
  
To install parallelLastz conda packages, in the terminal or an Anaconda Prompt, run:

```
conda install -c jnarayan81 parallellastz
```
The test data can be found at https://github.com/jnarayan81/parallelLastz/tree/master/testDATA, and the sample configuration file at https://github.com/jnarayan81/parallelLastz/blob/master/conf.

## Citation
Harris, R.S. (2007) Improved pairwise alignment of genomic DNA. Ph.D. Thesis, The Pennsylvania State University.

Please feel free to give this repository a few likes as encouragement. :+1: :pray: :clap: 

## Help
Contact me at jnarayan81@gmail.com or info@bioinformaticsonline.com
