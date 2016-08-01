#!/usr/bin/perl

my $version = "0.3";

## changes from v0.2
## changed group_name delimiters from | to # to work with nucmer_backbone_v0.3
## Will no longer by default align sequences to themselves. These are ignored by nucmer_backbone anyway and just take time and space. Will include a flag to include these if desired, but default will be to skip self vs self
## fixed a glitch in sequence grouping where sequences from different groups were being put in the same temporary sequence file.
## To save time and cycles, will also give the option to not perform back-alignments, i.e. instead of genA vs genB AND genB vs genA, just do genA_vs_genB. Will modify nucmer_backbone_v0.3 to accept this input and automatically include the reverse of all alignments. 
## output version

use strict;
use warnings;
use File::Which;
use File::Spec;

my $usage = "
nucmer_multi.pl - Performs all-vs-all nucmer alignments of multi-fasta files.

Required:
  -f    input fasta file
  
Optional:
  Load distribution options:
  -g	separate by groups. Fasta records with headers prefixed by \"#group_name#\"
        where \"group_name\" is the name of different strains, for example, will be
        separated. Records without headers in this format will be run individually.

    OR
	
  -d    maximum number of fasta records to run nucmer on at one time.
        WARNING: If this number is too small, may result in intense slowdown or
        crashing, too big and you may exceed nucmer's internal limits.
        (default: 1)
	
  Other options:
  -n    path to nucmer application (including nucmer). If not given, will try
        to find nucmer in your PATH
  -a    Additional arguments to pass to nucmer, surrounded by quotes,
        i.e. \"--maxmatch -b 200\".
  -t    Number of threads to use (Default: 15)
  -o    Output prefix (default: \"out\")
  -s    Include self-vs-self alignments (default: self-vs-self alignments are
        not done when -g is given)
  -b    Include back-alignments, i.e. will align genome A vs genome B and also
        align genome B vs genome A
        (default: each potential alignment will only be included once)
  -V    (uppercase) output version number and quit

";

use Getopt::Std;
our ($opt_f, $opt_d, $opt_n, $opt_a, $opt_t, $opt_o, $opt_g, $opt_s, $opt_b, $opt_V);
getopts('f:d:n:a:t:o:gsbV');
print "$version\n" and exit if $opt_V;
die $usage unless $opt_f;

my $fasin   = $opt_f if $opt_f;
my $dice    = $opt_d ? $opt_d : 1;
my $nucpath = $opt_n if $opt_n;
my $args    = $opt_a ? $opt_a : " ";
my $threads = $opt_t ? $opt_t : 15;
my $pref    = $opt_o ? $opt_o : "out";

if (!$nucpath){
    $nucpath = which('nucmer');
}

die "Nucmer could not be found in your PATH.\n" if (!$nucpath);

if (-d ($nucpath)){
    die "$nucpath is a directory. Please locate working nucmer executable.\n";
}
if (!-e ($nucpath)){
    die "$nucpath does not exist. Please locate working nucmer executable.\n";
}

while ($args =~ m/(-+p[^-]*)/){
    $args =~ s/$1//;  #removes any user-entered prefix commands to nucmer
}

my $abs_fa = File::Spec->rel2abs($fasin);

#read in fasta file and process into separate sequences
my ($id, $seq);
my @t_fa_files;
my $t_fa_count = 0;
my $d_count = 0;
my $out_string;
my $last_gen;
my @gens;
print STDERR "Preparing sequence files for alignment...\n";
open (my $infa, "<", $fasin) or die "Can't open $fasin: $!\n";
while (my $line = <$infa>){
    chomp $line;
    if ($line =~ m/^>/){
        if ($id){
            $out_string .= ">$id\n$seq\n";
            if ($opt_g){
                if (!$last_gen){
                    ($last_gen) = $id =~ m/^\#([^\#]+)\#/;
                }
                (my $gen) = $line =~ m/^>\#([^\#]+)\#/;
                if ($gen){
                    if (!$last_gen or ($gen ne $last_gen)){
                        $t_fa_count++;
                        open (my $out_fa, "> temp_$t_fa_count.fasta") or die "Can't open temporary output file: $!\n";
                        print $out_fa "$out_string";
                        close ($out_fa);
                        push @t_fa_files, "temp_$t_fa_count.fasta";
                        $out_string = "";
                        $d_count = 0;
                        push @gens, ([$last_gen, $t_fa_count]);
                    }
                    $last_gen = $gen;
                } else {
                    $t_fa_count++;
                    open (my $out_fa, "> temp_$t_fa_count.fasta") or die "Can't open temporary output file: $!\n";
                    print $out_fa "$out_string";
                    close ($out_fa);
                    push @t_fa_files, "temp_$t_fa_count.fasta";
                    $out_string = "";
                    $d_count = 0;
                    $last_gen = "";
                    push @gens, ([$id, $t_fa_count]);
                }
            } else {
                $d_count++;
                if ($d_count == $dice){
                    $t_fa_count++;
                    open (my $out_fa, "> temp_$t_fa_count.fasta") or die "Can't open temporary output file: $!\n";
                    print $out_fa "$out_string";
                    close ($out_fa);
                    push @t_fa_files, "temp_$t_fa_count.fasta";
                    $out_string = "";
                    $d_count = 0;
                }
            }
        }
        $id = substr($line, 1);
        $seq = "";
        next;
    }
    $line =~ s/\s.*$//;
    $seq .= $line;
}
if ($id){
    $out_string .= ">$id\n$seq\n";
    $t_fa_count++;
    open (my $out_fa, "> temp_$t_fa_count.fasta") or die "Can't open temporary output file: $!\n";
    print $out_fa "$out_string";
    close ($out_fa);
    push @t_fa_files, "temp_$t_fa_count.fasta";
    if ($opt_g){
        if ($last_gen){
            push @gens, ([$last_gen, $t_fa_count]);
        } else {
            push @gens, ([$id, $t_fa_count]);
        }
    }
    ($id, $seq) = ("") x 2;
} else {
    die "ERROR: $fasin does not appear to be a fasta file.  Please double-check.\n";
}

@gens = sort{$a->[0] cmp $b->[0]}@gens if @gens; #sort in the order that coords will eventually be output by show-coords

#make array of all possible combinations of alignments
my @combos;
if ($opt_b){ #include all back-alignments
    for my $i (0 .. $#t_fa_files){
        my $ref = $t_fa_files[$i];
        for my $j (0 .. $#t_fa_files){
            my $qry = $t_fa_files[$j];
            next if ($qry eq $ref and $opt_g and !$opt_s);
            push @combos, ([$ref, $qry]);
        }
    }
} else {
    if ($opt_g){
        for my $i (0 .. ($#gens - 1)){
            my $ref = $gens[$i][1];
            my $start = $i + 1;
            $start = $i if $opt_s;
            for my $j ($start .. $#gens){
                my $qry = $gens[$j][1];
                push @combos, (["temp_$ref.fasta", "temp_$qry.fasta"]);
            }
        }
    } else {
        for my $i (0 .. ($#t_fa_files - 1)){
            my $ref = $t_fa_files[$i];
            my $start = $i + 1;
            $start = $i if $opt_s;
            for my $j ($start .. $#t_fa_files){
                my $qry = $t_fa_files[$j];
                push @combos, ([$ref, $qry]);
            }
        }
    }
}
my $num_combos = scalar @combos;

#run alignments
my $num_threads_running = 0;
my $num_done = 0;
print STDERR "\rFinished $num_done of $num_combos alignments";
open (my $out, "> $pref.delta");
print $out "$abs_fa $abs_fa\nNUCMER\n";
while (@combos){
    if ($num_threads_running < $threads){
        my ($ref, $qry) = @{shift @combos};
        my $pid = fork;
        if (0 == $pid){
            my $tpref = "$$";
            `$nucpath -p $tpref $args $ref $qry > /dev/null 2>&1`;
            exit ($?);
        }
        $num_threads_running++;
    }
    if ($num_threads_running == $threads){
        my $pid = wait;
        my $status = $?;
        die "ERROR: Process $pid returned a status of $status\n" if ($status != 0);
        open (my $din, "< $pid.delta") or die "Can't open $pid.delta: $!\n";
        while (my $line = <$din>){
            chomp $line;
            next if $line =~ m/^\//;
            next if $line =~ m/^NUCMER/;
            next if $line =~ m/^\s*$/;
            print $out "$line\n";
        }
        unlink ("$pid.delta");
        $num_done++;
        print STDERR "\rFinished $num_done of $num_combos alignments";
        $num_threads_running--;
    }
}
while (my $pid = wait){
    last if $pid < 0;
    my $status = $?;
    die "ERROR: Process $pid returned a status of $status\n" if ($status != 0);
    open (my $din, "< $pid.delta") or die "Can't open $pid.delta: $!\n";
    while (my $line = <$din>){
        chomp $line;
        next if $line =~ m/^\//;
        next if $line =~ m/^NUCMER/;
        next if $line =~ m/^\s*$/;
        print $out "$line\n";
    }
    unlink ("$pid.delta");
    $num_done++;
    print STDERR "\rFinished $num_done of $num_combos alignments";
    $num_threads_running--;
}
print STDERR "\n";
close ($out);

#delete temporary files;
for my $i (0 .. $#t_fa_files){
    my $ref = $t_fa_files[$i];
    unlink ($ref);
}
