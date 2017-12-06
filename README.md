# parallelLastz
## blastz with multi-threads support.

Running Lastz (https://github.com/lastz/lastz) in parallel mode. This program is for single computer with multiple core processors.

When the query file format is fasta, you can specify many threads to process it. It can reduce run time linearly, and use almost equal memory as the original blastz program. This is useful when you blastz a big query file to a huge reference like human whole genome sequence.

The program is an extension on the original blastz program which was written by Bob Harris.

parallelBlastz can run on Linux and Mac OS.

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

Contact me at jnarayan81@gmail.com
