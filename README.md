# parallelLastz v0.2
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
Usage: parallelLastz.pl --qfile <> --tfile <> --cfile <> --speedup <#>
Options:
   --qfile|-q      Query multifasta/fasta file
   --tfile|-t      Target genome file
   --cfile|-c      Config file
   --speedup|-s    Number of cores to use
   --length|-l     Minimum length of sequences to process
   --unmask|-u     Unmask lowercase in target and query files
   --wipe|-w       Wipe intermediate files
   --verbose|-v    Enable verbose logging
   --retry|-r      Number of retry attempts for failed jobs
   --output|-o     Output directory for saving results
   --help|-h       Show this help message

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
