# parallelLastz v0.3.3
## Lastz with multi-threads support.
[![Install with Conda](https://img.shields.io/badge/Install%20with-Conda-44A833?style=flat-square&logo=anaconda&logoColor=white)](https://anaconda.org/jnarayan81/parallellastz)
[![Conda Version](https://anaconda.org/jnarayan81/parallellastz/badges/version.svg)](https://anaconda.org/jnarayan81/parallellastz)
[![Conda Downloads](https://anaconda.org/jnarayan81/parallellastz/badges/downloads.svg)](https://anaconda.org/jnarayan81/parallellastz)
[![Platforms](https://anaconda.org/jnarayan81/parallellastz/badges/platforms.svg)](https://anaconda.org/jnarayan81/parallellastz)
[![Perl](https://img.shields.io/badge/Perl-5.x-2bbc8a?style=flat-square&logo=perl&logoColor=white)](https://www.perl.org/)
[![LASTZ](https://img.shields.io/badge/LASTZ-supported-blue?style=flat-square)](https://www.bx.psu.edu/~rsharris/lastz/)
[![License](https://anaconda.org/conda-forge/lcms2/badges/license.svg)](LICENSE)

Running Lastz (https://github.com/lastz/lastz) in parallel mode. This program is for single computer with multiple core processors.

# parallelLastz

**parallelLastz** is a Perl-based utility for running [LASTZ](https://www.bx.psu.edu/~rsharris/lastz/) sequence alignments in parallel. It is designed for large FASTA datasets where multiple target and query sequences need to be aligned efficiently while maintaining a resumable job manifest and detailed error reporting.

---

## Features

* 🚀 Run multiple LASTZ alignments in parallel
* 🧬 Supports multi-FASTA target and query files
* 📦 Automatically splits multi-FASTA targets into single-sequence temporary FASTA files
* ✂️ Supports query chunking based on a user-defined length
* ⚡ Control parallelism with `--jobs/-j`
* 🔄 Resume interrupted or incomplete runs
* 🔁 Retry failed LASTZ jobs
* 🧪 Dry-run mode for inspecting the planned jobs
* 📋 Maintains a job manifest for tracking execution
* 📝 Captures LASTZ stdout and stderr separately
* ❌ Reports the actual LASTZ error message instead of only an exit code
* 📊 Displays a compact progress indicator in non-verbose mode
* 🔍 Verbose mode for detailed execution information
* 🧹 Optional cleanup/wipe of previous results
* 🔤 Optional sequence unmasking
* 🛡️ Uses safe list-form command execution

---

## Requirements

### Software

The following software is required:

* Perl >= 5.20
* LASTZ
* Perl module `Parallel::ForkManager`
* Perl module `Bio::SeqIO` from BioPerl

### Install Perl dependencies

Using CPAN:

```bash
cpan Parallel::ForkManager Bio::SeqIO
```

Or using CPAN Minus:

```bash
cpanm Parallel::ForkManager Bio::SeqIO
```

Check LASTZ:

```bash
lastz --version
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/jnarayan81/parallelLastz.git
cd parallelLastz
```

Make the script executable:

```bash
chmod +x parallelLastz_v0.3.3.pl
```

Test the installation:

```bash
perl parallelLastz_v0.3.3.pl --help
```

---

# Basic Usage

```bash
perl parallelLastz_v0.3.3.pl \
    -q query.fasta \
    -t target.fasta \
    -c conf \
    -j 4 \
    -l 100
```

Where:

| Option            | Description                                   |
| ----------------- | --------------------------------------------- |
| `-q`, `--qfile`   | Query FASTA file                              |
| `-t`, `--tfile`   | Target FASTA file                             |
| `-c`, `--cfile`   | LASTZ configuration file                      |
| `-j`, `--jobs`    | Number of concurrent LASTZ jobs               |
| `-l`, `--length`  | Maximum query chunk length                    |
| `-w`, `--wipe`    | Remove previous output before starting        |
| `-r`, `--retry`   | Number of retries for failed jobs             |
| `-u`, `--unmask`  | Convert lowercase sequence bases to uppercase |
| `-v`, `--verbose` | Display detailed execution information        |
| `--resume`        | Resume an existing run                        |
| `--dry-run`       | Display planned jobs without executing LASTZ  |
| `-o`, `--output`  | Output directory                              |
| `-h`, `--help`    | Display help                                  |

---

# Parallelization

The `--jobs/-j` option controls the maximum number of LASTZ processes running simultaneously.

For example:

```bash
-j 4
```

allows up to four LASTZ jobs to run concurrently.

If there are eight alignment jobs:

```text
Job 1 ─┐
Job 2 ─┤
Job 3 ─┤  Running
Job 4 ─┘

Job 5
Job 6
Job 7
Job 8
       Waiting
```

As jobs finish, waiting jobs are started automatically.

The requested number of jobs cannot exceed the number of available CPUs.

---

# Multi-FASTA Handling

A major feature of `parallelLastz` is automatic handling of multi-FASTA target files.

LASTZ expects a single sequence in the target position in the standard invocation. Therefore, if the target contains multiple sequences, `parallelLastz` first creates temporary single-record FASTA files.

For example:

```text
target.fasta
├── chromosome_1
└── chromosome_2
```

is converted internally to:

```text
targets/
├── target_000001_chromosome_1.fa
└── target_000002_chromosome_2.fa
```

Each LASTZ invocation then receives exactly one target sequence.

---

# Query Chunking

Query sequences can be processed in chunks using:

```bash
-l 100000
```

For example:

```text
Query sequence
────────────────────────────────────────
       100 kb        100 kb        100 kb
          ↓             ↓             ↓
      chunk_1       chunk_2       chunk_3
```

This allows very large query sequences to be processed without loading the entire alignment workload into a single LASTZ process.

---

# Number of Alignment Jobs

For multiple target and query sequences, jobs are generated as:

```text
Number of jobs =
Number of target sequences × Number of query chunks
```

For example:

```text
Target:
  T1
  T2

Query:
  Q1
  Q2
  Q3
  Q4
```

produces:

```text
T1 × Q1
T1 × Q2
T1 × Q3
T1 × Q4

T2 × Q1
T2 × Q2
T2 × Q3
T2 × Q4
```

Total:

```text
2 × 4 = 8 alignment jobs
```

---

# LASTZ Configuration

The `-c/--cfile` option specifies a configuration file containing LASTZ arguments.

Example:

```text
--chain
--progress
--identity=90
--ambiguous=iupac
--format=general-
```

The configuration is parsed safely and passed to LASTZ as separate command-line arguments.

`parallelLastz` ensures that:

```text
--format=general-
```

is retained for the alignment output.

---

# Example Configuration

Create a file named:

```text
conf
```

with:

```text
--chain
--progress
--identity=90
--ambiguous=iupac
--format=general-
```

Then run:

```bash
perl parallelLastz_v0.3.3.pl \
    -q genomes/query.fna \
    -t genomes/target.fna \
    -c conf \
    -j 4 \
    -l 100000
```

---

# Progress Indicator

In normal mode, `parallelLastz` displays a compact progress indicator:

```text
LASTZ: [###############.........] 6/8 completed | 2 running | 0 failed
```

This avoids flooding the terminal with individual commands.

For detailed information, use:

```bash
-v
```

Example:

```bash
perl parallelLastz_v0.3.3.pl \
    -q query.fna \
    -t target.fna \
    -c conf \
    -j 4 \
    -v
```

Verbose mode reports information such as:

* FASTA records
* Job IDs
* LASTZ commands
* Retry attempts
* Output files
* Errors
* Job status

---

# Resume

A run can be resumed using:

```bash
--resume
```

Example:

```bash
perl parallelLastz_v0.3.3.pl \
    -q query.fna \
    -t target.fna \
    -c conf \
    -j 4 \
    --resume
```

The program reads:

```text
parallelLastz.manifest.tsv
```

and skips jobs that have already completed successfully.

This is particularly useful for large datasets where an alignment run may take several hours or days.

---

# Retry Failed Jobs

Failed LASTZ jobs can automatically be retried.

For example:

```bash
--retry 3
```

means that a failed job can be attempted again up to three times.

Example:

```bash
perl parallelLastz_v0.3.3.pl \
    -q query.fna \
    -t target.fna \
    -c conf \
    -j 8 \
    --retry 3
```

The manifest records the number of attempts for each job.

---

# Dry Run

Use:

```bash
--dry-run
```

to inspect the planned workload without running LASTZ.

Example:

```bash
perl parallelLastz_v0.3.3.pl \
    -q query.fna \
    -t target.fna \
    -c conf \
    -j 4 \
    --dry-run
```

This is useful for checking:

* Number of target sequences
* Number of query chunks
* Number of alignment jobs
* Output locations
* Parallelization settings

---

# Output Structure

A typical output directory contains:

```text
output/
├── chunks/
│   ├── targets/
│   │   ├── target_000001.fa
│   │   └── target_000002.fa
│   ├── job_000001.lz
│   ├── job_000002.lz
│   └── ...
│
├── logs/
│   ├── job_000001.stdout
│   ├── job_000001.stderr
│   ├── job_000002.stdout
│   └── job_000002.stderr
│
├── parallelLastz.manifest.tsv
└── finalAlign.tsv
```

### Alignment files

Individual LASTZ results are stored as:

```text
job_000001.lz
job_000002.lz
...
```

### Final alignment

Successful alignment outputs are combined into:

```text
finalAlign.tsv
```

### Manifest

The manifest tracks each job:

```text
job_id
target_file
target_id
query_chunk
query_id
query_start
query_end
output
status
attempts
exit_code
stdout_file
stderr_file
message
```

This provides a persistent record of the complete analysis.

---

# Error Handling

`parallelLastz` captures both:

```text
STDOUT
STDERR
```

from LASTZ.

For example, instead of reporting only:

```text
LASTZ failed with exit code 1
```

the program reports the actual diagnostic message from LASTZ.

This makes failures considerably easier to troubleshoot.

---

# Successful Empty Output

An important behavior is that an alignment producing **zero records** is still considered successful when LASTZ exits with code `0`.

For example:

```text
LASTZ exit code: 0
Alignment output: empty
Status: SUCCESS
```

This prevents biologically valid "no alignment found" results from being incorrectly classified as failed jobs.

---

# CPU Usage

By default, the program detects the number of available CPUs.

You can explicitly specify the number of parallel jobs:

```bash
-j 8
```

It is recommended to choose a value appropriate for both:

* CPU availability
* available RAM

More parallel processes do not necessarily produce better performance if the system becomes memory or I/O constrained.

---

# Example: Complete Run

```bash
perl parallelLastz_v0.3.3.pl \
    --qfile genomes/query.fna \
    --tfile genomes/target.fna \
    --cfile conf \
    --jobs 8 \
    --length 100000 \
    --retry 2 \
    --wipe
```

For detailed output:

```bash
perl parallelLastz_v0.3.3.pl \
    --qfile genomes/query.fna \
    --tfile genomes/target.fna \
    --cfile conf \
    --jobs 8 \
    --length 100000 \
    --retry 2 \
    --wipe \
    --verbose
```

---

# Resume an Interrupted Run

If the previous run was interrupted:

```bash
perl parallelLastz_v0.3.3.pl \
    --qfile genomes/query.fna \
    --tfile genomes/target.fna \
    --cfile conf \
    --jobs 8 \
    --length 100000 \
    --resume
```

Do **not** use `--wipe` when you want to resume an existing run.

---

# Recommended Workflow

For a large genome comparison:

```bash
# 1. Test the configuration
perl parallelLastz_v0.3.3.pl \
    -q query.fna \
    -t target.fna \
    -c conf \
    -j 8 \
    --dry-run

# 2. Start the alignment
perl parallelLastz_v0.3.3.pl \
    -q query.fna \
    -t target.fna \
    -c conf \
    -j 8 \
    -l 100000 \
    --retry 2

# 3. If interrupted, resume
perl parallelLastz_v0.3.3.pl \
    -q query.fna \
    -t target.fna \
    -c conf \
    -j 8 \
    -l 100000 \
    --resume
```

---

# Citation

If you use `parallelLastz` in your research, please cite the software repository and LASTZ.

### LASTZ

Harris RS.
**Improved pairwise alignment of genomic DNA.**
PhD thesis, Pennsylvania State University.

LASTZ:

https://www.bx.psu.edu/~rsharris/lastz/

---

# Author

**Jitendra Narayan** and chatGPT for this v0.3.3 version 

---

# License

Please add the license appropriate for your project.

For example:

```text
MIT License
```

See the `LICENSE` file for details.

---

# Version

**parallelLastz v0.3.3**

A parallel LASTZ workflow with multi-FASTA handling, resumable execution, retry support, job manifests, error capture, and terminal progress monitoring.


## Citation
Harris, R.S. (2007) Improved pairwise alignment of genomic DNA. Ph.D. Thesis, The Pennsylvania State University.

Please feel free to give this repository a few likes as encouragement. :+1: :pray: :clap: 

## Help
Contact me at jnarayan81@gmail.com or info@bioinformaticsonline.com
