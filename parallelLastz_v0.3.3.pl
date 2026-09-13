#!/usr/bin/env perl
#
# parallelLastz v0.3.3
#
# Parallel LASTZ runner with streaming FASTA chunking, resumable manifests,
# retry handling, safe list-form execution, and complete LASTZ diagnostics.
#
# v2-compatible core options:
#   --qfile/-q       query FASTA
#   --tfile/-t       target FASTA
#   --cfile/-c       LASTZ configuration file
#   --length/-l      chunk length in bp
#   --wipe/-w        remove/recreate working files
#   --unmask/-u      pass unmasked query/target sequence to LASTZ
#   --verbose/-v
#   --retry/-r
#   --output/-o
#
# v0.3 additions:
#   --jobs/-j        maximum concurrent LASTZ jobs
#   --resume         resume from parallelLastz.manifest.tsv
#   --dry-run        print planned jobs without running LASTZ
#
# LASTZ is always invoked with list-form exec(), never through a shell.
#

use strict;
use warnings;

use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Copy qw(copy);
use Cwd qw(abs_path getcwd);
use POSIX qw(strftime);
use Text::ParseWords qw(shellwords);
use Sys::Hostname qw(hostname);
use IO::Handle;

my $VERSION = '0.3.3';

my ($qfile, $tfile, $cfile, $length);
my $wipe = 0;
my $unmask = 0;
my $verbose = 0;
my $retry = 0;
my $output = 'parallelLastz_out';
my $jobs;
my $resume = 0;
my $dry_run = 0;
my $help = 0;
my $version = 0;

GetOptions(
    'qfile|q=s'   => \$qfile,
    'tfile|t=s'   => \$tfile,
    'cfile|c=s'   => \$cfile,
    'length|l=i'  => \$length,
    'wipe|w'      => \$wipe,
    'unmask|u'    => \$unmask,
    'verbose|v'   => \$verbose,
    'retry|r=i'   => \$retry,
    'output|o=s'  => \$output,
    'jobs|j=i'    => \$jobs,
    'resume'      => \$resume,
    'dry-run'     => \$dry_run,
    'help|h'      => \$help,
    'version'     => \$version,
) or die usage();

# -------------------------------------------------------------------------
# Welcome banner
# Print only when --help/-h or --verbose/-v is requested.
# -------------------------------------------------------------------------
print_welcome() if $help || $verbose;

if ($version) {
    print "parallelLastz v$VERSION\n";
    exit 0;
}

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined($qfile) && defined($tfile) &&
                    defined($cfile) && defined($length);

die "--length/-l must be > 0\n" unless $length > 0;
die "--retry/-r must be >= 0\n" unless $retry >= 0;
die "--jobs/-j must be > 0\n" if defined($jobs) && $jobs < 1;

$qfile = abs_path($qfile) // die "Cannot resolve query FASTA: $qfile\n";
$tfile = abs_path($tfile) // die "Cannot resolve target FASTA: $tfile\n";
$cfile = abs_path($cfile) // die "Cannot resolve config file: $cfile\n";

die "Query FASTA does not exist: $qfile\n" unless -f $qfile;
die "Target FASTA does not exist: $tfile\n" unless -f $tfile;
die "Config file does not exist: $cfile\n" unless -f $cfile;

# -------------------------------------------------------------------------
# CPU validation
# -------------------------------------------------------------------------
my $available_cpus = available_cpus();
$jobs = $available_cpus unless defined $jobs;

die "--jobs/-j ($jobs) exceeds available CPUs ($available_cpus). " .
    "Choose --jobs <= $available_cpus.\n"
    if $jobs > $available_cpus;

# -------------------------------------------------------------------------
# Output/work directories
# -------------------------------------------------------------------------
$output = File::Spec->rel2abs($output);

my $chunks_dir     = File::Spec->catdir($output, 'chunks');
my $logs_dir       = File::Spec->catdir($output, 'logs');
my $runs_dir       = File::Spec->catdir($output, 'runs');
my $target_tmp_dir = File::Spec->catdir($output, 'targets');

my $manifest   = File::Spec->catfile($output, 'parallelLastz.manifest.tsv');
my $final_align = File::Spec->catfile($output, 'finalAlign.tsv');

if ($wipe && -d $output) {
    print STDERR "Wiping output directory: $output\n" if $verbose;
    remove_tree($output) or die "Unable to wipe $output: $!\n";
}

make_path($output, $chunks_dir, $logs_dir, $runs_dir, $target_tmp_dir);

# A normal non-resume invocation starts a fresh manifest.
# This is deliberately an overwrite, not an append.
if (!$resume) {
    open(my $mf, '>', $manifest)
        or die "Cannot create $manifest: $!\n";

    print $mf join("\t",
        qw(job_id target_file target_id query_chunk query_id query_start query_end output status attempts
           exit_code stdout_file stderr_file message)
    ), "\n";

    close $mf;
}

# -------------------------------------------------------------------------
# LASTZ configuration
# -------------------------------------------------------------------------
my @config_args = read_config($cfile);

# v2 compatibility: preserve the historical configuration arguments while
# ensuring the required output format is explicitly present.
@config_args = grep {
    $_ ne '--format=general-' && $_ ne '--format'
} @config_args;

push @config_args, '--format=general-';

# --unmask/-u is retained as a compatibility switch.
# By default the FASTA records are written exactly as read.
# With --unmask, sequence masking is removed before writing chunks.

# -------------------------------------------------------------------------
# Streaming FASTA chunk generation
# -------------------------------------------------------------------------
my @jobs_manifest;

if ($resume && -s $manifest) {

    @jobs_manifest = read_manifest($manifest);

    die "Resume requested, but manifest contains no jobs: $manifest\n"
        unless @jobs_manifest;

    print STDERR "Resuming " . scalar(@jobs_manifest) .
                 " manifest jobs from $manifest\n"
        if $verbose;
}
else {

    # LASTZ requires a single target sequence. Split the target FASTA into
    # one temporary FASTA per record, then align every target record against
    # every query chunk. This avoids passing a multi-FASTA target to LASTZ.
    my @target_records = split_target_fasta(
        $tfile,
        $target_tmp_dir,
        $unmask,
        $verbose
    );

    my @query_jobs = create_chunks(
        $qfile,
        $chunks_dir,
        $length,
        $unmask,
        $verbose
    );

    @jobs_manifest = expand_target_query_jobs(
        \@target_records,
        \@query_jobs,
        $chunks_dir
    );

    write_manifest($manifest, \@jobs_manifest);
}

if ($dry_run) {
    print_plan(
        \@jobs_manifest,
        $tfile,
        \@config_args,
        $jobs
    );

    exit 0;
}

# -------------------------------------------------------------------------
# Execute jobs
# -------------------------------------------------------------------------
eval { require Parallel::ForkManager; 1 }
    or die "Parallel::ForkManager is required for parallel execution: $@\n";

my $pm = Parallel::ForkManager->new($jobs);

my $total_jobs     = scalar(@jobs_manifest);
my $completed_jobs = 0;
my $failed_jobs    = 0;
my $running_jobs   = 0;

if (!$verbose) {

    $pm->run_on_finish(sub {
        my (
            $pid,
            $exit_code,
            $ident,
            $exit_signal,
            $core_dump,
            $data
        ) = @_;

        ++$completed_jobs;

        if (defined($data) &&
            ref($data) eq 'HASH' &&
            !$data->{ok}) {

            ++$failed_jobs;
        }

        $running_jobs = $total_jobs - $completed_jobs;
        $running_jobs = 0 if $running_jobs < 0;

        print_progress(
            $completed_jobs,
            $total_jobs,
            $running_jobs,
            $failed_jobs
        );
    });

    print_progress(0, $total_jobs, 0, 0);
}

for my $job (@jobs_manifest) {

    # A completed successful job is skipped on resume.
    if ($resume &&
        $job->{status} eq 'success' &&
        -f $job->{output}) {

        print STDERR "RESUME: skipping completed job $job->{job_id}\n"
            if $verbose;

        next;
    }

    # Existing manifest jobs are retried from their current state.
    $job->{status}  = 'pending';
    $job->{message} = '';

    $pm->start and next;

    my $result = run_job(
        $job,
        \@config_args,
        $logs_dir,
        $runs_dir,
        $retry,
        $verbose
    );

    # Parent cannot receive lexical variables from a forked child.
    # Persist the complete result in a per-job status file and have the
    # parent read it after finish().
    #
    # This also fixes child attempt-count propagation.
    write_child_result(
        $runs_dir,
        $job->{job_id},
        $result
    );

    $pm->finish(
        $result->{ok} ? 0 : 1,
        $result
    );
}

$pm->wait_all_children;

print STDERR "\n" unless $verbose;

# -------------------------------------------------------------------------
# Collect child results and rewrite manifest atomically
# -------------------------------------------------------------------------
for my $job (@jobs_manifest) {

    my $result_file = File::Spec->catfile(
        $runs_dir,
        sprintf('%s.result.tsv', $job->{job_id})
    );

    if (-f $result_file) {

        my $r = read_child_result($result_file);

        for my $k (keys %$r) {
            $job->{$k} = $r->{$k};
        }
    }
    elsif ($job->{status} ne 'success') {

        $job->{status}    = 'failed';
        $job->{message}   = 'No child result file was produced';
        $job->{exit_code} = 255;
    }
}

write_manifest($manifest, \@jobs_manifest);

# -------------------------------------------------------------------------
# Build final alignment output from successful jobs
# -------------------------------------------------------------------------
my $successful = 0;
my $failed     = 0;

open(my $fa, '>', $final_align)
    or die "Cannot create $final_align: $!\n";

for my $job (@jobs_manifest) {

    if ($job->{status} eq 'success' &&
        -f $job->{output}) {

        open(my $in, '<', $job->{output})
            or die "Cannot read $job->{output}: $!\n";

        while (my $line = <$in>) {
            print $fa $line;
        }

        close $in;

        ++$successful;
    }
    else {
        ++$failed;
    }
}

close $fa;

print STDERR "\nparallelLastz v$VERSION completed\n";
print STDERR "  Alignment jobs: " . scalar(@jobs_manifest) . "\n";
print STDERR "  Successful   : $successful\n";
print STDERR "  Failed       : $failed\n";
print STDERR "  CPUs available: $available_cpus\n";
print STDERR "  Jobs         : $jobs\n";
print STDERR "  Manifest     : $manifest\n";
print STDERR "  Final output : $final_align\n";

exit($failed ? 1 : 0);

# =========================================================================
# SUBROUTINES
# =========================================================================

sub print_welcome {
    print <<'WELCOME';

       ╭─╮ ╭─╮ ╭─╮ ╭─╮
       │P│─│L│─│Z│─│▸│
       ╰─╯ ╰─╯ ╰─╯ ╰─╯
        parallelLastz
     Parallel LASTZ Runner
     ─────────────────────
       v0.3.3 | JitendraLab

WELCOME

    # Replace the placeholder with the actual version.
    # This keeps the banner definition readable while avoiding hard-coded
    # version strings.
}

sub usage {
    return <<"USAGE";
parallelLastz v$VERSION

Usage:
  parallelLastz.pl -q query.fa -t target.fa -c lastz.conf -l 100000 [options]

Required:
  -q, --qfile FILE       Query FASTA
  -t, --tfile FILE       Target FASTA
  -c, --cfile FILE       LASTZ configuration file
  -l, --length INT       Maximum query chunk length in bp

Compatibility:
  -u, --unmask            Remove FASTA lowercase masking before LASTZ
  -w, --wipe              Recreate output/work directory
  -v, --verbose           Verbose progress messages
  -r, --retry INT         Number of retries after an initial failure
  -o, --output DIR        Output directory (default: parallelLastz_out)

Parallel/resume:
  -j, --jobs INT          Concurrent LASTZ jobs; must not exceed CPUs
      --resume            Resume using parallelLastz.manifest.tsv
      --dry-run           Generate/inspect the plan but do not execute LASTZ
      --version            Print version
  -h, --help              Show this help

Outputs:
  parallelLastz.manifest.tsv
  finalAlign.tsv
  chunks/
  logs/
  runs/

LASTZ is invoked safely using list-form exec(), with:
  --format=general-

For diagnostics, both LASTZ stdout and stderr are captured per attempt.
A successful LASTZ invocation is considered successful even if stdout is empty.
USAGE
}

sub available_cpus {
    my $n;

    if ($^O =~ /MSWin32/i) {
        $n = $ENV{NUMBER_OF_PROCESSORS};
    }
    else {
        $n = `getconf _NPROCESSORS_ONLN 2>/dev/null`;
        chomp $n;
    }

    $n = int($n || 1);

    return $n > 0 ? $n : 1;
}

sub read_config {
    my ($file) = @_;

    open(my $fh, '<', $file)
        or die "Cannot read config $file: $!\n";

    my @args;

    while (my $line = <$fh>) {

        chomp $line;

        $line =~ s/^\s+//;
        $line =~ s/\s+$//;

        next if $line eq '';
        next if $line =~ /^#/;
        next if $line =~ /^;/;

        # Strip trailing comments only when the # begins a whitespace-
        # delimited comment, preserving # characters inside arguments.
        $line =~ s/\s+#.*$//;

        push @args, shellwords($line);
    }

    close $fh;

    return @args;
}

sub split_target_fasta {
    my ($file, $dir, $unmask, $verbose) = @_;

    eval { require Bio::SeqIO; 1 }
        or die "Bio::SeqIO is required: $@\n";

    my $in = Bio::SeqIO->new(
        -file   => $file,
        -format => 'fasta'
    ) or die "Cannot read target FASTA $file\n";

    my @records;
    my $n = 0;

    while (my $seq = $in->next_seq) {

        ++$n;

        my $id = $seq->display_id ||
                 $seq->id ||
                 "target_$n";

        my $safe = $id;

        $safe =~ s/[^A-Za-z0-9_.-]+/_/g;
        $safe = "target_$n" unless length $safe;

        my $path = File::Spec->catfile(
            $dir,
            sprintf(
                "target_%06d_%s.fa",
                $n,
                $safe
            )
        );

        my $s = $seq->seq;

        if ($unmask) {
            $s =~ s/[a-z]/uc($&)/ge;
        }

        open(my $fh, '>', $path)
            or die "Cannot create $path: $!\n";

        print $fh ">$id\n";

        for (
            my $i = 0;
            $i < length($s);
            $i += 80
        ) {
            print $fh substr($s, $i, 80), "\n";
        }

        close $fh
            or die "Cannot close $path: $!\n";

        push @records, {
            file   => $path,
            id     => $id,
            length => length($s)
        };

        print STDERR
            "Target FASTA: $id (" .
            length($s) .
            " bp) -> $path\n"
            if $verbose;
    }

    die "No FASTA records found in target file $file\n"
        unless @records;

    return @records;
}

sub expand_target_query_jobs {
    my ($targets, $queries, $chunks_dir) = @_;

    my @jobs;
    my $n = 0;

    for my $t (@$targets) {

        for my $q (@$queries) {

            ++$n;

            my $id = sprintf(
                'job_%06d',
                $n
            );

            push @jobs, {
                job_id      => $id,
                target_file => $t->{file},
                target_id   => $t->{id},
                query_chunk => $q->{query_chunk},
                query_id    => ($q->{query_id} // ''),
                query_start => $q->{query_start},
                query_end   => $q->{query_end},
                output      => File::Spec->catfile(
                    $chunks_dir,
                    "$id.lz"
                ),
                status      => 'pending',
                attempts    => 0,
                exit_code   => '',
                stdout_file => '',
                stderr_file => '',
                message     => ''
            };
        }
    }

    die "No target/query alignment jobs were created\n"
        unless @jobs;

    return @jobs;
}

sub create_chunks {
    my ($file, $dir, $max_bp, $unmask, $verbose) = @_;

    eval { require Bio::SeqIO; 1 }
        or die "Bio::SeqIO is required for streaming FASTA processing: $@\n";

    my $in = Bio::SeqIO->new(
        -file   => $file,
        -format => 'fasta'
    );

    my @jobs;

    my $job_no = 0;
    my $buffer = '';
    my $buffer_bp = 0;
    my $chunk_start = 1;

    while (my $seq = $in->next_seq) {

        my $id = $seq->display_id;
        my $s  = $seq->seq;

        if ($unmask) {
            $s =~ s/[a-z]/uc($&)/ge;
        }

        my $seq_len = length($s);

        # Preserve FASTA records. A single record longer than --length is
        # emitted as a single record rather than split in the middle.
        if ($buffer_bp > 0 &&
            $buffer_bp + $seq_len > $max_bp) {

            my $job_id = sprintf(
                'job_%06d',
                ++$job_no
            );

            my $path = File::Spec->catfile(
                $dir,
                "$job_id.fa"
            );

            write_fasta($path, $buffer);

            push @jobs, {
                job_id      => $job_id,
                query_chunk => $path,
                query_id    => "query_chunk_$job_no",
                query_start => $chunk_start,
                query_end   => $chunk_start + $buffer_bp - 1,
                output      => File::Spec->catfile(
                    $dir,
                    "$job_id.lz"
                ),
                status      => 'pending',
                attempts    => 0,
                exit_code   => '',
                stdout_file => '',
                stderr_file => '',
                message     => ''
            };

            $chunk_start += $buffer_bp;
            $buffer = '';
            $buffer_bp = 0;
        }

        $buffer .= ">$id\n$s\n";
        $buffer_bp += $seq_len;

        if ($verbose) {
            print STDERR
                "Streaming FASTA: $id ($seq_len bp)\n";
        }
    }

    if ($buffer_bp > 0) {

        my $job_id = sprintf(
            'job_%06d',
            ++$job_no
        );

        my $path = File::Spec->catfile(
            $dir,
            "$job_id.fa"
        );

        write_fasta($path, $buffer);

        push @jobs, {
            job_id      => $job_id,
            query_chunk => $path,
            query_start => $chunk_start,
            query_end   => $chunk_start + $buffer_bp - 1,
            output      => File::Spec->catfile(
                $dir,
                "$job_id.lz"
            ),
            status      => 'pending',
            attempts    => 0,
            exit_code   => '',
            stdout_file => '',
            stderr_file => '',
            message     => ''
        };
    }

    die "No FASTA records found in $file\n"
        unless @jobs;

    return @jobs;
}

sub write_fasta {
    my ($file, $text) = @_;

    open(my $fh, '>', $file)
        or die "Cannot create $file: $!\n";

    print $fh $text;

    close $fh
        or die "Cannot close $file: $!\n";
}

sub write_manifest {
    my ($file, $jobs) = @_;

    my $tmp = "$file.tmp.$$";

    open(my $fh, '>', $tmp)
        or die "Cannot write $tmp: $!\n";

    print $fh join("\t",
        qw(job_id target_file target_id query_chunk query_id query_start query_end output status attempts
           exit_code stdout_file stderr_file message)
    ), "\n";

    for my $j (@$jobs) {

        print $fh join("\t",
            map {
                tsv_escape($j->{$_} // '')
            }
            qw(job_id target_file target_id query_chunk query_id query_start query_end output status attempts
               exit_code stdout_file stderr_file message)
        ), "\n";
    }

    close $fh
        or die "Cannot close $tmp: $!\n";

    rename $tmp, $file
        or die "Cannot replace $file: $!\n";
}

sub read_manifest {
    my ($file) = @_;

    open(my $fh, '<', $file)
        or die "Cannot read $file: $!\n";

    my $header = <$fh>;

    die "Manifest is empty: $file\n"
        unless defined $header;

    chomp $header;

    my @head = split /\t/, $header, -1;

    my @jobs;

    while (my $line = <$fh>) {

        chomp $line;

        next if $line eq '';

        my @v = split /\t/, $line, -1;

        my %j;

        @j{@head} = map {
            tsv_unescape($_)
        } @v;

        # Resolve relative paths against the manifest directory.
        my ($vol, $dir, $name) =
            File::Spec->splitpath($file);

        my $base = File::Spec->catpath(
            $vol,
            $dir,
            ''
        );

        for my $key (
            qw(target_file query_chunk output stdout_file stderr_file)
        ) {

            next unless defined($j{$key}) &&
                        $j{$key} ne '';

            $j{$key} = File::Spec->rel2abs(
                $j{$key},
                $base
            ) unless File::Spec->file_name_is_absolute(
                $j{$key}
            );
        }

        push @jobs, \%j;
    }

    close $fh;

    return @jobs;
}

sub tsv_escape {
    my ($v) = @_;

    $v =~ s/\t/\\t/g;
    $v =~ s/\r/\\r/g;
    $v =~ s/\n/\\n/g;

    return $v;
}

sub tsv_unescape {
    my ($v) = @_;

    $v =~ s/\\t/\t/g;
    $v =~ s/\\r/\r/g;
    $v =~ s/\\n/\n/g;

    return $v;
}

sub print_progress {
    my ($done, $total, $running, $failed) = @_;

    return if $verbose;

    my $width = 24;

    my $filled = $total
        ? int(($done / $total) * $width)
        : 0;

    $filled = $width
        if $filled > $width;

    my $bar =
        ('#' x $filled) .
        ('.' x ($width - $filled));

    printf STDERR
        "\rLASTZ: [%s] %d/%d completed | %d running | %d failed",
        $bar,
        $done,
        $total,
        $running,
        $failed;

    STDERR->flush
        if STDERR->can('flush');
}

sub run_job {
    my (
        $job,
        $config,
        $logs_dir,
        $runs_dir,
        $max_retry,
        $verbose
    ) = @_;

    my $max_attempts = $max_retry + 1;
    my $attempt = 0;

    my $last_exit = 255;
    my $last_stdout = '';
    my $last_stderr = '';
    my $last_stdout_file = '';
    my $last_stderr_file = '';
    my $last_message = '';

    while ($attempt < $max_attempts) {

        ++$attempt;

        my $tag = sprintf(
            '%s.attempt_%02d',
            $job->{job_id},
            $attempt
        );

        my $stdout_file = File::Spec->catfile(
            $logs_dir,
            "$tag.stdout"
        );

        my $stderr_file = File::Spec->catfile(
            $logs_dir,
            "$tag.stderr"
        );

        my @cmd = (
            'lastz',
            $job->{target_file},
            $job->{query_chunk},
            @$config,
            "--output=$job->{output}",
        );

        # LASTZ output is requested in the command line; nevertheless stdout
        # is also redirected at the process level for diagnostics and for
        # LASTZ builds/configurations that write ordinary stdout.
        if ($verbose) {

            print STDERR
                "[$job->{job_id}] attempt " .
                "$attempt/$max_attempts\n";

            print STDERR
                "[$job->{job_id}] CMD: " .
                join(
                    ' ',
                    map {
                        quote_for_display($_)
                    } @cmd
                ) .
                "\n";
        }

        my $pid = fork();

        die "fork failed for $job->{job_id}: $!\n"
            unless defined $pid;

        if ($pid == 0) {

            open(STDOUT, '>', $stdout_file)
                or die "Cannot redirect stdout: $!\n";

            open(STDERR, '>', $stderr_file)
                or die "Cannot redirect stderr: $!\n";

            exec { $cmd[0] } @cmd
                or do {
                    print STDERR
                        "Cannot exec LASTZ: $!\n";

                    exit 127;
                };
        }

        waitpid($pid, 0);

        my $status = $?;

        my $exit_code;

        if ($status == -1) {
            $exit_code = 255;
        }
        elsif ($status & 127) {
            $exit_code = 128 + ($status & 127);
        }
        else {
            $exit_code = $status >> 8;
        }

        $last_exit = $exit_code;
        $last_stdout_file = $stdout_file;
        $last_stderr_file = $stderr_file;

        $last_stdout = slurp($stdout_file);
        $last_stderr = slurp($stderr_file);

        # LASTZ can legitimately produce an empty alignment.
        # Therefore output size is NOT used as a success criterion:
        # exit status is.
        if ($exit_code == 0) {

            return {
                status      => 'success',
                ok          => 1,
                attempts    => $attempt,
                exit_code   => 0,
                stdout_file => $stdout_file,
                stderr_file => $stderr_file,
                message     => '',
                output      => $job->{output},
            };
        }

        $last_message = diagnostic_message(
            $exit_code,
            $last_stdout,
            $last_stderr
        );

        print STDERR
            "[$job->{job_id}] LASTZ failed " .
            "(attempt $attempt): $last_message\n";

        last if $attempt >= $max_attempts;
    }

    return {
        status      => 'failed',
        ok          => 0,
        attempts    => $attempt,
        exit_code   => $last_exit,
        stdout_file => $last_stdout_file,
        stderr_file => $last_stderr_file,
        message     => $last_message,
        output      => $job->{output},
    };
}

sub diagnostic_message {
    my ($exit, $stdout, $stderr) = @_;

    my $detail = $stderr;

    $detail = $stdout
        if $detail eq '';

    $detail =~ s/\s+/ /g;
    $detail =~ s/^\s+|\s+$//g;

    if ($detail ne '') {
        return "LASTZ exit code $exit: $detail";
    }

    return
        "LASTZ exit code $exit " .
        "(no diagnostic output captured)";
}

sub slurp {
    my ($file) = @_;

    return ''
        unless -f $file;

    open(my $fh, '<', $file)
        or return '';

    local $/;

    my $x = <$fh>;

    close $fh;

    return defined($x) ? $x : '';
}

sub write_child_result {
    my ($dir, $job_id, $r) = @_;

    my $file = File::Spec->catfile(
        $dir,
        "$job_id.result.tsv"
    );

    my $tmp = "$file.tmp.$$";

    open(my $fh, '>', $tmp)
        or die "Cannot write child result $tmp: $!\n";

    for my $key (
        qw(status attempts exit_code stdout_file stderr_file message output)
    ) {

        print $fh
            tsv_escape($r->{$key} // ''),
            "\n";
    }

    close $fh
        or die "Cannot close child result $tmp: $!\n";

    rename $tmp, $file
        or die "Cannot install child result $file: $!\n";
}

sub read_child_result {
    my ($file) = @_;

    open(my $fh, '<', $file)
        or die "Cannot read $file: $!\n";

    my @keys = qw(
        status
        attempts
        exit_code
        stdout_file
        stderr_file
        message
        output
    );

    my %r;

    for my $key (@keys) {

        my $line = <$fh>;

        last unless defined $line;

        chomp $line;

        $r{$key} = tsv_unescape($line);
    }

    close $fh;

    $r{ok} =
        ($r{status} // '') eq 'success'
        ? 1
        : 0;

    return \%r;
}

sub print_plan {
    my (
        $jobs,
        $target,
        $config,
        $n_jobs
    ) = @_;

    print "parallelLastz v$VERSION dry-run\n";
    print "Target: $target\n";
    print "Jobs:   $n_jobs\n";
    print "LASTZ arguments:\n";

    for my $a (@$config) {
        print "  ",
            quote_for_display($a),
            "\n";
    }

    print "\nPlanned jobs:\n";

    for my $j (@$jobs) {

        print join(
            "\t",
            $j->{job_id},
            $j->{query_chunk},
            $j->{output},
            $j->{status}
        ), "\n";
    }
}

sub quote_for_display {
    my ($s) = @_;

    return "''"
        if $s eq '';

    return $s
        if $s =~ /^[A-Za-z0-9_\/.=:,+-]+$/;

    $s =~ s/'/'"'"'/g;

    return "'$s'";
}

