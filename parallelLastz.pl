#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use Parallel::ForkManager;
use Bio::SeqIO;
use File::Temp;
use File::Spec;
use Term::ProgressBar;
use Log::Dispatch;

# Author: Jitendra Narayan / jnarayan81@gmail.com
# Usage: perl parallelLastz.pl <multi_fasta_qfile> <tfile> <cfile> <thread> <length>
# perl parallelLastz.pl -q testDATA/qsample1.fa -t testDATA/tsample.fa -c conf -s 4 -w -l 10

print <<'WELCOME';
                       _ _      _   __           _       
 _ __   __ _ _ __ __ _| | | ___| | / /  __ _ ___| |_ ____
| '_ \ / _` | '__/ _` | | |/ _ \ |/ /  / _` / __| __|_  /
| |_) | (_| | | | (_| | | |  __/ / /__| (_| \__ \ |_ / / 
| .__/ \__,_|_|  \__,_|_|_|\___|_\____/\__,_|___/\__/___|v0.2
|_|                                                      
parallelLastz: Run lastz jobs in parallel
Contact: jnarayan81@gmail.com for support

WELCOME

# Variables declaration
my ($qfile, $tfile, $config, $thread, $length, $wipe, $help, $unmask, $verbose, $retry, $output_dir);
my $version = 0.2;

# Parse command line options
GetOptions(
    'qfile|q=s'        => \$qfile,
    'tfile|t=s'        => \$tfile,
    'cfile|c=s'        => \$config,
    'speedup|s=i'      => \$thread,
    'length|l=i'       => \$length,
    'wipe|w'           => \$wipe,
    'unmask|u'         => \$unmask,
    'verbose|v'        => \$verbose,
    'retry|r=i'        => \$retry,
    'output|o=s'       => \$output_dir,
    'help|h'           => \$help
) or die usage($version);
usage($version) if $help;

# Validate mandatory inputs
my @missing;
push @missing, 'query file (--qfile or -q)' unless $qfile;
push @missing, 'target file (--tfile or -t)' unless $tfile;
push @missing, 'config file (--cfile or -c)' unless $config;
push @missing, 'length (--length or -l)' unless $length;

if (@missing) {
    print "Error: Missing the following required option(s):\n";
    print " - $_\n" for @missing;
    exit;
}

$thread ||= `grep -c '^processor' /proc/cpuinfo` || 1;
$output_dir ||= '.';
mkdir $output_dir unless -d $output_dir;

$retry ||= 1;
$verbose ||= 0;

# Initialize logging
my $logger = Log::Dispatch->new(
    outputs => [
        ['Screen', min_level => 'info', newline => 1]
    ]
);

# Read configuration file
my $parameters = readConfig($config);
my $param = join(' ', @$parameters);

# Optionally convert query file to uppercase
if ($unmask) {
    my $qfile_corrected = File::Spec->catfile($output_dir, 'tmpUC');
    convertToUpperCase($qfile, $qfile_corrected);
    $qfile = $qfile_corrected;
}

# Locate lastz binary
my $lastZ_tool = locateLastz() or die "No lastz command found.\nTry installing it using conda: 'conda install -c bioconda lastz'\n";

# Load sequences from target file
my %sequences = loadSequences($tfile);

# Initialize parallel processing
my $pm = Parallel::ForkManager->new($thread);

# Initialize progress bar
my $progress = Term::ProgressBar->new({
    name   => 'Processing',
    count  => scalar(keys %sequences),
    fh     => \*STDOUT
});

$pm->run_on_finish(
    sub {
        my ($pid, $exit_code, $ident) = @_;
        $progress->update();  # Update progress bar
        if ($verbose) {
            $logger->info("Process $ident finished with exit code $exit_code.");
        }
    }
);

$pm->run_on_start(
    sub {
        my ($pid, $ident) = @_;
        $logger->info("Process $ident started with PID $pid.") if $verbose;
    }
);

# Process each sequence in parallel
foreach my $id (keys %sequences) {
    next if length($sequences{$id}) <= $length;
    
    my $pid = $pm->start($id) and next;
    my $attempts = 0;
    my $success = 0;
    
    while ($attempts < $retry && !$success) {
        eval {
            runLastz($id, $sequences{$id}, $qfile, $param, $unmask, $lastZ_tool, $output_dir);
            $success = 1;
        };
        if ($@) {
            $logger->error("Error running lastz for $id: $@");
            $attempts++;
            sleep 1;  # Wait before retrying
        }
    }
    
    $pm->finish if $success;
    if (!$success) {
        $logger->error("Failed to process $id after $retry attempts.");
    }
}

# Wait for all jobs to complete
$pm->wait_all_children;
$logger->info("All jobs completed.");

# Optional: Clean up intermediate files
if ($wipe) {
    cleanUpIntermediateFiles($output_dir);
    $logger->info("Alignment completed. Final results are in 'finalAlign.tsv'.");
} else {
    $logger->info("Alignment completed. Check the individual '.lz' files in $output_dir.");
}

# Subroutines

# Usage message
sub usage {
    my $ver = shift;
    print "\n parallelLastz v$ver\n";
    print "Usage: $0 --qfile <> --tfile <> --cfile <> --speedup <#>\n";
    print "Options:\n";
    print "   --qfile|-q      Query multifasta/fasta file\n";
    print "   --tfile|-t      Target genome file\n";
    print "   --cfile|-c      Config file\n";
    print "   --speedup|-s    Number of cores to use\n";
    print "   --length|-l     Minimum length of sequences to process\n";
    print "   --unmask|-u     Unmask lowercase in target and query files\n";
    print "   --wipe|-w       Wipe intermediate files\n";
    print "   --verbose|-v    Enable verbose logging\n";
    print "   --retry|-r      Number of retry attempts for failed jobs\n";
    print "   --output|-o     Output directory for saving results\n";
    print "   --help|-h       Show this help message\n";
    exit;
}

# Read sequences from a fasta file
sub loadSequences {
    my ($filename) = @_;
    my %sequences;
    my $seqio = Bio::SeqIO->new(-file => $filename, -format => "fasta");

    while (my $seqobj = $seqio->next_seq) {
        my $id  = $seqobj->display_id;
        my $seq = $seqobj->seq;
        
        # Replace spaces in the header with underscores
        $id =~ s/\s/_/g;
        
        $sequences{$id} = $seq;
    }

    return %sequences;
}


# Run lastz for a specific sequence
sub runLastz {
    my ($name, $seq, $qfile, $param, $unmask, $lastZ_tool, $output_dir) = @_;
    $logger->info("Processing $name...");

    # Create a temporary file for the sequence
    my $tmp_fh = File::Temp->new(UNLINK => 1);
    print $tmp_fh ">$name\n$seq\n";
    close $tmp_fh;

    # Build and execute lastz command
    my $output_file = File::Spec->catfile($output_dir, "seeALN_$name.lz");
    my $myLASTZ = "$lastZ_tool $tmp_fh->filename $qfile --output=$output_file $param";
    system($myLASTZ) == 0 or die "Error running lastz on $name: $!";
    $logger->info("lastz completed for $name, output saved to $output_file.");
}

# Locate lastz in the system path
sub locateLastz {
    my $lastZ_tool = "lastz";
    my $tool_path = '';
    my $found = 0;

    # Check if 'lastz' exists in the current directory (./)
    if (-f "./$lastZ_tool" && -x "./$lastZ_tool") {
        $logger->info("'$lastZ_tool' found in the current directory.");
        $tool_path = "./$lastZ_tool";
        $found = 1;
    } else {
        # If not found in the current directory, search in system's PATH
        for my $path (split /:/, $ENV{PATH}) {
            if (-f "$path/$lastZ_tool" && -x "$path/$lastZ_tool") {
                $logger->info("'$lastZ_tool' found in $path");
                $tool_path = "$path/$lastZ_tool";
                $found = 1;
                last;
            }
        }
    }

    return $tool_path if $found;
    die "No '$lastZ_tool' command found in the current directory or system's PATH.\nTry installing it using conda: 'conda install -c bioconda lastz'\n";
}

# Convert sequences to uppercase
sub convertToUpperCase {
    my ($infile, $outfile) = @_;
    open my $in_fh, '<', $infile or die "Cannot open input file $infile: $!";
    open my $out_fh, '>', $outfile or die "Cannot open output file $outfile: $!";

    while (<$in_fh>) {
        print $out_fh uc($_);
    }

    close $in_fh;
    close $out_fh;
}

# Read config file
sub readConfig {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot open config file $file: $!";
    my @lines;

    while (<$fh>) {
        chomp;
        next if /^#/ || /^\s*$/;
        push @lines, $_;
    }

    close $fh;
    return \@lines;
}

# Clean up intermediate files
sub cleanUpIntermediateFiles {
    my ($output_dir) = @_;
    my @files = glob(File::Spec->catfile($output_dir, '*.lz'));

    # Avoid using a large number of files with `cat` directly
    my $final_output = File::Spec->catfile($output_dir, 'finalAlign.tsv');

    open my $out_fh, '>', $final_output or die "Cannot open final output file $final_output: $!";

    for my $file (@files) {
        open my $in_fh, '<', $file or die "Cannot open file $file: $!";
        while (<$in_fh>) {
            print $out_fh $_;
        }
        close $in_fh;
        unlink $file; # Remove the intermediate file after processing
    }

    close $out_fh;
}

