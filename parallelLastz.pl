#!/usr/bin/perl

use strict;
use 5.010;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);
use Data::Dumper;
use Pod::Usage;
use Parallel::ForkManager;
use Bio::SeqIO;

#Author: Jitendra Narayan
#Usage: perl parallelLastz.pl <multi_fasta_qfile> <tfile> <cfile> <thread>
#Change the Lastz parameters at Line:81

my ($qfile, $tfile, $config, $thread, $debug, $help, $man);
my $version=0.1;
GetOptions(
    'qfile|q=s' => \$qfile,
    'tfile|t=s' => \$tfile,
    'cfile|c=s' => \$config,
    'speedup|s=i' => \$thread,
    'help|h' => \$help
) or die &help($version);
&help($version) if $help;
#pod2usage("$0: \nI am afraid, no files given.")  if ((@ARGV == 0) && (-t STDIN));

if (!$qfile or !$tfile or !$config) { help($version) }
if (!$thread) { $thread = `grep -c -P '^processor\\s+:' /proc/cpuinfo` }
my $parameters = readConfig ($config);
my $param = join (' ', @$parameters);

my %sequences;
my $seqio = Bio::SeqIO->new(-file => "$qfile", -format => "fasta");
while(my$seqobj = $seqio->next_seq) {
    my $id  = $seqobj->display_id;    # there's your key
    my $seq = $seqobj->seq;           # and there's your value
    $sequences{$id} = $seq;
}

  my $max_procs = $thread;
  my @names = keys %sequences;
  # hash to resolve PID's back to child specific information
  my $pm =  new Parallel::ForkManager($max_procs);
  # Setup a callback for when a child finishes up so we can
  # get it's exit code
  $pm->run_on_finish (
    sub { my ($pid, $exit_code, $ident) = @_;
      #print "** $ident just got out of the pool ". "with PID $pid and exit code: $exit_code\n";
    }
  );

  $pm->run_on_start(
    sub { my ($pid,$ident)=@_;
     #print "** $ident started, pid: $pid\n";
    }
  );

  $pm->run_on_wait(
    sub {
      #print "** Have to wait for one children ...\n"
    },
    0.5
  );

  NAMES:
  foreach my $child ( 0 .. $#names ) {
    my $pid = $pm->start($names[$child]) and next NAMES;
    runLastz($names[$child]);
    $pm->finish($child); # pass an exit code to finish
  }
  print "Waiting for all the jobs to complete...\n";
  $pm->wait_all_children;
  print "DONE ... Everybody is out of the computation pool!\n";

sub runLastz {
my $name=shift;
my $seq=$sequences{$name};
# remove the file when the reference goes away with the UNLINK option
	my $tmp_fh = new File::Temp( UNLINK => 1 );
	print $tmp_fh ">$name\n$seq\n";

        print "Working on $name sequence\n";
	my $myLASTZ="lastz $tmp_fh $tfile --output=seeALN_$name.lz $param";
        system ("$myLASTZ");

}

#Read config files
sub readConfig {
my ($file) = @_;
my $fh= read_fh($file);
my @lines;
while (<$fh>) {
    chomp;
    next if /^#/;
    next if /^$/;
    $_ =~ s/^\s+|\s+$//g;
    push @lines, $_;
}
close $fh or die "Cannot close $file: $!";
return \@lines;
}

#Open and Read a file
sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}

#Help section
sub help {
  my $ver = $_[0];
  print "\n parallelLastz $ver\n\n";

  print "Usage: $0 --qfile <> --tfile <> --core <#> \n\n";
  print	"Options:\n";
  print "	--qfile|-q	query multifasta/fasta file\n";
  print "	--tfile|-t	target genome file\n";
  print "	--cfile|-c	config file\n";
  print "	--speedup|-s	number of core to use\n";
  print "     	--help|-h	brief help message\n";

exit;
}
