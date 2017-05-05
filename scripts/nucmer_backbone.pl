#!/usr/bin/perl

my $license = "
    nucmer_backbone.pl
    Copyright (C) 2016 Egon A. Ozer

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see [http://www.gnu.org/licenses/].
";

my $version = "0.3.2";

## Changes from v0.3.1
# Fixed bug where accessory regions composed entirely or mostly of ambiguous bases were causing an error when coordiates were being output

## Changes from v0.3
# Fixed bug where first CDS on each contig was not being output

## Changes from v0.2
# Output core and accessory loci from genomes with annotation data, regardless of whether they are multicontig or if non-annotated genomes are included
#  --> format of the output loci output files changed to include positions of loci on input as well as output core / accessory / backbone sequences
# Output both core and accessory element sequences from all genomes as well as the backbone sequence (i.e. species core from the reference(s)).
# remove genbank file processing. This will all be done by the Spine wrapper instead.
# Change the strain delimiter from "|" to "#" to avoid possible problems with NCBI sequence headers. Will also automatically restore the original contig names for the output
# Do not output all-N core, accessory, or backbone sequences. Trim leading and tailing ambiguous bases (N,X,-) from output sequences 
# Added annotation processing and post-processing to the forks to improve speed
# Several changes to remove the number of redundant loops to save time and removed reliance on loading large lists or per-coordinate hashes to greatly improve memory usage
# for annotation data, output locus IDs and product names of CDSs
# MAJOR CHANGE: output partial genes in the core and accessory lists, rather than assigning them either just to core or just to accesory based on a cutoff percentage. This will be to allow downstream programs (like ClustAGE) to identify variable portions of otherwise conserved core genes. 
# Option to output pangenome sequence as well

## Changes from v0.1
# Support for draft, multi-contig genomes
# Multithreading support
# Improved memory usage by not loading all CDS loci into a hash right from the beginning, but instead loading them genome by genome as needed.
# Fixed a glitch where a gene could potentially be output as both core and accessory if it was exactly 50% represented in either (or whatever was set by option -p)

use strict;
use warnings;

my $usage = "
nucmer_backbone.pl - Uses nucmer alignment to determine core and accessory
                     genome
                     
version = $version

Copyright (C) 2016  Egon A. Ozer
This program comes with ABSOLUTELY NO WARRANTY. This is free software, and you
are welcome to redistribute it under certain conditions;
see LICENSE.txt for details

Usage: perl nucmer_backbone.pl <options> [fasta_file_1] [fasta_file_2] ...

Required:
  -c		alignment [C]oordinates file. Produced by numcer alignment of
                all genomes vs. all genomes:
                    \$ nucmer -p out --maxmatch genomes.fasta genomes.fasta                
                followed by running show-coords and sorting by reference: 
                    \$ show-coords -rTH out.delta > out.coords.txt
                
  [fasta_file]	fasta files of sequences aligned by nucmer
			or
		one multi-fasta file of sequences aligned by nucmer
                
                Fasta records with the same header prefix will be grouped
                together, as in the case of multi-contig draft genomes. Make
                sure each header starts with \"#genome_name#\".
                i.e.
                If the genome \"strainA\"is contained in two contigs:
                >contig1
                AGCAAAG...
                >contig2
                AGAACCC...
                Rename the headers like so:
                >#strainA#contig1
                AGCAAAG...
                >#strainA#contig2
                AGAACCC...
  
Optional:  
  -a		number of genomes from which a section can be [A]bsent and still
		be included as core genome (default 0)
  -r		genome used as [R]eference for generating backbone fragment list
		with 1 being the first genome entered, 2 being the second, etc.
		(default is first genome entered, but you may try several
		runs altering this number as number and length of backbone
		segments output may be slightly different)
  -g		list of [G]enomes from which backbone should be created,
		separated by commas with no spaces between. Order of the list
		will determine priority in producing core genome. Genome given
		by -r will be given first priorty, with all other genomes moving
		down one step. Number of genomes entered must be at least (a+1),
                so if the value given by -a is 2, at least three genomes must be
                entered.
		Example: to create a backbone from the first, second, and fourth
		         genomes given and ignore the third, enter \"1,2,4\".
			 The first genome will be given highest priority,
			 second given next-highest, fourth given third-highest
		(default is to use genomes in the order their sequences were
                provided)
[deprecated]  -l		path to file listing paths to genbank files of genome
		annotations, given in the same order of the input fasta files.
		File locations should be separated into separate lines.
		    example:
		    /path/to/genbank_file_1.gbk
		    /path/to/genbank_file_2.gbk
		    etc.
		if this file is given, a list of locus IDs of genes present in
		the core genome (IDs from the reference genome) will be output,
		as will lists of accessory genome genes for each of the input
		sequences.
                *** This option currently only works for non-draft sequences
  -x		file listing all cds locus ids in the format:
		genome order number<tab>locus id<tab>contig_id<tab>start coord<tab>end coord<tab>strand<tab>product (optional)
		(one entry per line)
			example:
			1	gen1_00001	contig1 456	2176    +   exoU
			1	gen1_00002	contig2 3187	4599    +   spcU
			2	gen2_00001	contigA 679	3100    -   hypothetical protein
			etc.
		If an input is given for both -x and -l, only -x will be used as
		input. If the file given by -x does not exist, it will be
                created using the files from -l to create a locus id file.
[deprecated]  -p		minimum percent coverage of a locus for it to be output as core
		(default 50)	 
  -m		[M]aximum distance between backbone fragments.
		Distances less than this parameter between adjacent fragments
		will result in combination of fragments with interspersed N's
		(default 10)
  -h            minimum percent identity of a nucmer alignment to be considered
                homologous (default 85)
  -B		minimum [B]ackbone size to be output (default 10)
  -I		minimum [I]sland size to be output (default 10)
  -s		prefix of output files (default \"output\")
  -o		if given, will output coordinates of core and accessory genome
		segments for each genome. Takes longer (default: only output
                accessory statistics and coordinates for reference genome(s))
  -e            Output file of position values that can be used to calculate
                pangenome and core genome characteristics of the data set
                using core_and_pangenome.pl. If selected, will automatically
                trigger -o option.
  -n        Output a pangenome sequence from all input genomes. This sequence
            will be generated from sequences in the order given or from the
            order given to -g and -r above.
  -t            Number of threads (default 15)
  -v		verbose output
  -V (uppercase)    Output version number and quit.        

";

# command line processing
use Getopt::Std;
our ($opt_c, $opt_f, $opt_a, $opt_r, $opt_g, $opt_m, $opt_h, $opt_o, $opt_b, $opt_B, $opt_I, $opt_s, $opt_v, $opt_l, $opt_x, $opt_e, $opt_n, $opt_t, $opt_V, $opt_w);
getopts('c:f:a:r:g:s:m:h:b:B:I:l:x:t:oenvVw');
print "$version\n" and exit if $opt_V;
die $usage unless ($opt_c and @ARGV);

my ($coord, $fast, $abs, $ref, $minpct, $use, $stat, $max, $minhom, $outfile, $backfile, $backlen, $isllen, $gbks, $lid_file, $threads);

$coord	    = $opt_c if $opt_c;
$fast		= $opt_f if $opt_f;
$abs		= $opt_a ? $opt_a : 0;
$ref		= $opt_r ? $opt_r : 1;
$gbks		= $opt_l if $opt_l;
#$minpct	    = $opt_p ? $opt_p : 50;
$minpct     = 50;
$use		= $opt_g if $opt_g;
$stat		= $opt_s ? $opt_s : "output";
$max		= $opt_m ? $opt_m : 10;
$minhom     = $opt_h ? $opt_h : 85;
#$outfile	= $opt_o ? $opt_o : "core.fa";
$backfile	= $opt_b ? $opt_b : "backbone.fa";
$backlen	= $opt_B ? $opt_B : 10;
$isllen	    = $opt_I ? $opt_I : 10;
$lid_file	= $opt_x if $opt_x;
$threads    = $opt_t ? $opt_t : 15;

$opt_o = 1 if $opt_e;

my $make_lid_file;
if ($opt_x){
	$make_lid_file = 1 unless (-e $lid_file);
}

# Input sequence file(s) and store concatenated sequences in separate files
print STDERR "Reading in fasta sequences ... \n";
print STDERR "<br>\n" if $opt_w;
my $seq_count = 0;
my $g_count = 0;
my @order;
my (%seq_ids, %seq_counts, %g_size);
my (%contig_starts, %contig_ids, %contig_lengs);
for my $i (0 .. $#ARGV){
    my $infile;
    open ($infile, "<", $ARGV[$i]) or die "ERROR: Can't open input file $ARGV[$i]: $!\n$usage";
    my ($id, $seq, $c_id, $c_seq);
    while (my $line = <$infile>){
        chomp $line;
        #$line =~ s/\s.*$//; #removes everything after the first whitespace
        if ($line =~ /^>/){
            my $subid = substr($line, 1);
            if ($line =~ m/^>\#([^#]*)\#/){
                $subid = $1;
            }
            if ($id){
                unless ($subid eq $id){
                    open (my $out, ">tmp.seq.$seq_count.txt");
                    print $out "$seq\n";
                    close $out;
                    my $size = length($seq);
                    $g_size{$id} = $size;
                    $g_count++;
                    print STDERR "$id, code $seq_count, size $size\n" unless $opt_w;
                    push @order, $id;
                    
                }
                my $c_size = length $c_seq;
                push @{$contig_ids{$id}}, ([$c_id, $c_size]); #this is the only way to accurately capture genome/contig sizes. It means all of the sequence files have to be included.
                $contig_lengs{$c_id} = $c_size;
                $c_seq = "";
            }
            my $last_id = $id if $id;
            $id = $subid;
            $c_id = substr($line, 1);
            if (!$last_id or $id ne $last_id){
                $seq_count++;
                $seq_ids{$seq_count} = $id;
                $seq_counts{$id} = $seq_count;
                $seq = '';
                $contig_starts{$c_id} = 1;
            } else {
                $contig_starts{$c_id} = (length($seq) + 1);
            }
            next;
        }
        $line =~ s/\s//g;
        $seq .= $line;
        $c_seq .= $line;
    }
    open (my $out, ">tmp.seq.$seq_count.txt");
    print $out "$seq\n";
    close $out;
    my $size = length($seq);
    $g_size{$id} = $size;
    $g_count++;
    print STDERR "$id, code $seq_count, size $size\n" unless $opt_w;
    push @order, $id;
    my $c_size = length $c_seq;
    push @{$contig_ids{$id}}, ([$c_id, $c_size]); #this is the only way to accurately capture genome/contig sizes. It means all of the sequence files have to be included.
    $contig_lengs{$c_id} = $c_size;
    ($id, $seq, $c_seq) = ("") x 3;
}
die "ERROR: No sequences found in input fasta file(s)\n" if $seq_count == 0;

# Open Coordinates file and read into a hash
# I've been running into memory issues with loading a large number of coordinates into a single hash.
# Since nucmer_backbone should only be run through Spine now, I think I can safely presume that the coordinates file will come pre-sorted by genome, subsorted by contig.
# Using this assumption, I will instead load each genome's contigs, then output them to temporary files to be individually read in by the forks. Lots of I/O, but that's probably better than malloc errors.
# The only way this works is if nucmer_multi runs genomes in alphabetical order. Otherwise show-coords results won't be sorted correctly and this approach fails.

print STDERR "Reading alignments .";
print STDERR "..<br>\n" if $opt_w;
#my @order;
my (%hash);
my $d_start;
my $last_genome;
open (my $in, "<", $coord) or die "ERROR: Can't open $coord: $!\n";
while (my $line = <$in>){
    chomp $line;
    my ($r_start, $r_stop, $q_start, $q_stop, $x1, $x2, $pct_id, $r_id, $q_id) = split('\t', $line);
    next unless ($r_start and $r_start =~ m/\d+/); #skips unless the first value record is a number
    next if $pct_id < $minhom;
    my $subid = $r_id;
    if ($r_id =~ m/^\#([^#]*)\#/){
        $subid = $1;
    }
    my $q_subid = $q_id;
    if ($q_id =~ m/^\#([^#]*)\#/){
        $q_subid = $1;
    }
    if ($last_genome){
        if ($subid ne $last_genome){
            print STDERR "." unless $opt_w;
            if ($hash{$last_genome}){
                my $gcode = $seq_counts{$last_genome};
                open (my $out, ">tmp.crd.$gcode.txt") or die "ERROR: Can't open tmp.crd.$gcode.txt: $!\n";
                foreach my $ref (keys %{$hash{$last_genome}}){
                    my @array = @{$hash{$last_genome}{$ref}};
                    @array = sort{$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]}@array;
                    foreach my $slice (@array){
                        print $out "$ref\t", join("\t", @{$slice}), "\n";
                    }
                }
                close ($out);
                delete $hash{$last_genome};
            }
        }
    }
    
    if ($subid ne $q_subid){ #save some memory by not including self-hits (shouldn't be there anymore anyway with nucmer_multi_0.3)
        push @{$hash{$subid}{$r_id}}, ([$q_subid, $r_start, $r_stop]);
        ($q_start, $q_stop) = ($q_stop, $q_start) if ($q_start > $q_stop); #need to flip reversed query coordinates
        push @{$hash{$q_subid}{$q_id}}, ([$subid, $q_start, $q_stop]); #include the reverse alignment since this will no longer be performed by nucmer_multi_0.3
    }
    $last_genome = $subid;
}
close ($in);
if (%hash){
    die "Whoops, too many records left (".scalar(keys %hash).")\n" if (keys %hash) > 2; #for debugging. Should be able to remove.
    foreach my $gen (sort keys %hash) {
        print STDERR "." unless $opt_w;
        my $gcode = $seq_counts{$gen};
        open (my $out, ">tmp.crd.$gcode.txt") or die "ERROR: Can't open tmp.crd.$gcode.txt: $!\n";
        foreach my $ref (keys %{$hash{$gen}}){
            my @array = @{$hash{$gen}{$ref}};
            @array = sort{$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]}@array;
            foreach my $slice (@array){
                print $out "$ref\t", join("\t", @{$slice}), "\n";
            }
        }
        close ($out);
        delete $hash{$gen};
    }
}
print STDERR "\n";
print STDERR "<br>\n" if $opt_w;
%hash = ();

#checks whether -r input, if given, jibes with -c input characteristics
my $nog = scalar @order;
print STDERR "nog: $nog ref: $ref\n" unless $opt_w;
if (!$ref){
	$ref = $nog;
}
if ($ref < 1 or $ref > $nog or $ref=~/\D/){
	die "ERROR: Please enter a valid number for option -r (see below)***\n\n$usage\n";
}

#If locus id file was given, read locus IDs, coordinates, and product names to separate temporary files, one for each genome.
#  - Will do this instead of reading them to an array. With large numbers of genomes the array will take up too much memory when forking. 
my @loci_tmp_files;
if ($opt_x and !$make_lid_file){
    print STDERR "Reading locus IDs ... <br>\n" if $opt_w;
    open (my $fin, "<", $lid_file) or die "Can't open locus ID file $lid_file: $!\n";
    my $line_count;
    my $last_gbknum = -1;
    my $tmpcrdout;
    while (my $line = <$fin>){
        chomp ($line);
        next if $line =~ m/^\s*#/;  #skip any commented lines
        $line_count++;
        print STDERR "\rReading locus IDs from $lid_file: $line_count" unless $opt_w;
        my ($gbknum, $lid, $contig, $start, $stop, $dir, $prod) = split('\t', $line);
        if ($gbknum != $last_gbknum){
            close ($tmpcrdout) if $tmpcrdout;
            my $loci_tmp_file = "tmp_coordinates_$gbknum.txt";
            open ($tmpcrdout, ">$loci_tmp_file") or die "ERROR: Can't open $loci_tmp_file: $!\n";
            $loci_tmp_files[$gbknum] = "$loci_tmp_file";
        }
        print $tmpcrdout "$lid\t$contig\t$start\t$stop\t$dir";
        print $tmpcrdout "\t$prod" if $prod;
        print $tmpcrdout "\n";
        $last_gbknum = $gbknum;
    }
    close ($tmpcrdout) if $tmpcrdout;
    print STDERR " Done!\n" unless $opt_w;
    close ($fin);
}

$opt_l = 1 if ($opt_x); ### keep this for now. Will remove once I remove the -l option altogether.

## create array of genomes to use in the analysis (if given by -g option) 
my @to_use;
my $orig_ref = $ref;
if ($use){
    @to_use = split(/,/,$use);
    my ($bad_use, $ref_use);
    my %use_dup_buster;
    my @temp_use;
    my $new_ref;
    for my $i (0 .. $#to_use){
        my $use = $to_use[$i];
        if ($use =~ /\D/ or $use < 1 or $use > $nog){
            $bad_use = 1;
            last;
        }
        if ($use == $ref){
            $ref_use = 1;
            $use_dup_buster{$use}++;
            next;
        }
        next if ($use_dup_buster{$use}); #skips duplicate genome numbers
        push @temp_use, $use;
        $use_dup_buster{$use}++;
    }
    unshift @temp_use, $ref; #puts the ref given by -r first
    @to_use = @temp_use;
    die "***Please enter a valid list for option -g***\n\n" if ($bad_use);
    die "***Value entered for option -r ($ref) must be included in list given in option -g***\n\n" if (!$ref_use);
    die "***Not enough unique genomes given to option -g***\n\n" if (scalar @to_use < ($abs + 1));
    #add any non-included genomes to the end of the array
    for my $i (1 .. $nog){
        next if $use_dup_buster{$i}; #skips if the genome is already in the to_use array
        push @to_use, $i;
    }
} else {
    for my $i (1 .. $nog){
        push @to_use, $i unless $ref == $i;
    }
    unshift (@to_use, $ref);  #puts the ref given by -r first
}


if ($opt_v){
	print "-g list: ", join(",", @to_use), "\n";
}

# Initialize stats file
open (my $stats, ">$stat.statistics.txt");
select((select($stats), $|=1)[0]); #make the filehandle hot so it prints immediately (www.plover.com/FAQs/Buffering.html).
print $stats "Spine version: $version\n";
#print $stats "$0 version: $version\n";
print $stats "inputs: --pctcore $abs";
print $stats " --refs ", join(",", @to_use) unless $opt_w;
print $stats " --maxdist $max --pctid $minhom --minout $backlen";
print $stats " --pangenome" if $opt_n;

#print $stats "inputs: -a:$abs -r:$ref -g:";
#print $stats join(",", @to_use);
#print $stats " -m:$max -h:$minhom -B:$backlen -I:$isllen -c:$coord";
#print $stats " -l: $gbks" if ($gbks);
#print $stats " -x: $lid_file" if ($opt_x);
##print $stats " -p: $minpct" if ($opt_l);

print $stats "\n\n";
print $stats "gen_#\tgen_name\tgen_size\tsource\ttotal_bp\tgc_%\tnum_segs\tmin_seg\tmax_seg\tavg_leng\tmedian_leng";
print $stats "\tnum_cds" if ($opt_l);
print $stats "\n";

my $min_gen = $g_count - $abs;  # Minimum number of genomes in which a segement must be present to be "core"

print STDERR "<br>Processing...<br>\n" if $opt_w;

## Start the threads a-rollin'
my $finished = 0;
my ($part, $children);
my @child_order;
my $num_procs_running = 0;
my $eye = 0;
my @mt_to_use = @to_use;
while (@mt_to_use and $num_procs_running < $threads){
    last if $eye > $abs and !$opt_o; #skips if option to output core and accessory coords not selected
    my $gcode = shift @mt_to_use;
    start_next_process ($gcode, $eye);
    $eye++;
}

my $pos_count_out;
if ($opt_e){ 
    open ($pos_count_out, "> $stat.position_counts.txt");
    print $pos_count_out "ref_genome\tnum_genomes_sharing\tgenome_codes_not_sharing\tnum_shared_positions\n";
}

## Now round those threads on up
my @loci_order;
my %locuslengs;
my %loci_starts;
my %loci_stops;
my $done = 0;
my $total_children_wall_time;
my ($slowest, $fastest);
my @bb_gcodes;

do {
    $done = 1;
    my $pid = wait;
    my $status = $?;
    
    die "ERROR: wait returned pid: $pid\n" if (!$pid or $pid == -1);
    
    $children->{$pid}->{'done'} = 1;
    $children->{$pid}->{'exit'} = $status;
    $num_procs_running--;
    
    my $child_start_time = $children->{$pid}->{'start'};
    my $child_wall_time = time - $children->{$pid}->{'start'};
    $total_children_wall_time += $child_wall_time;
    if (!$slowest and !$fastest){
        $slowest = $child_wall_time;
        $fastest = $child_wall_time;
    } else {
        if ($child_wall_time < $fastest){
            $fastest = $child_wall_time;
        } elsif ($child_wall_time > $slowest){
            $slowest = $child_wall_time;
        }
    }
    $finished++;
    if ($opt_w){
        print STDERR "Finished processing $finished of $nog genomes<br>\n";
    } else {
        print STDERR "\nFinished $finished of $nog genomes.\n";
    }
    
    #set next process going
    if (@mt_to_use and $num_procs_running < $threads){
        unless ($eye > $abs and !$opt_o) { #skips if option to output core and accessory coords not selected
            my $gcode = shift @mt_to_use;
            start_next_process ($gcode, $eye);
            $eye++;
        }
    }
    
    while (@child_order > 0 and $children->{$child_order[0]}->{'done'}){
        $pid = shift @child_order;
        my $out = $children->{$pid}->{'out_st'};
        if ($children->{$pid}->{'exit'}){
            my $status = $children->{$pid}->{'exit'};
            die "\nChild $pid died with status $status\n" if $status;
        }
        if (not -e $out){
            die "\nERROR: no output file ($out) found for child: $pid\n";
        }
        #transfer the stats from the fork to the stats file
        open (my $in, "<", $out) or die "ERROR: Can't open $out: $!\n";
        while (my $line = <$in>){
            print $stats "$line";
        }
        close ($in);
        unlink ($out);
        
        #If the pangenome position file was requested, transfer those data to the universal position file
        if ($opt_e){
            if ($children->{$pid}->{'out_ph'}){
                my $phfile = $children->{$pid}->{'out_ph'};
                open (my $in, "<$phfile") or die "ERROR: Can't open $phfile: $!\n";
                while (my $line = <$in>){
                    print $pos_count_out "$line";
                }
                close ($in);
                unlink ($phfile);
            }
        }
        
        #we'll put all the gcodes of strains that belong to the backbone into an array for processing after all the threads are done. This way, we don't clog up the memory of subsequent forks as the bbone_seq string and bbone_lengs array are filled.
        my $gcode = $children->{$pid}->{'gcode'};
        my $proc_eye = $children->{$pid}->{'eye'};
        if ($proc_eye <= $abs){
            push @bb_gcodes, $gcode;
        }
    }
    
    my $dum;
    while ( ($pid, $dum) = each %$children){
        $done = 0 if ($children->{$pid}->{'done'} == 0);
    }
    $dum++;
    
} until ($done);

if ($opt_e){
    close ($pos_count_out);
}

#output backbone sequences and statistics
my @params = ("out");
push @params, @bb_gcodes;
process_final(\@params);

#output pangenome sequences and statistics (if requested)
if ($opt_n){
    @params = ("pan");
    push @params, @to_use;
    process_final(\@params);
}

close ($stats);


## output backbone multiple alignment (if requested);
my $ma; #dummy. This will be set by the user to request multiple alignment
if ($ma){
    #first, read the delta file to get insertion/deletion information for each of the alignments
    open (my $din, "<$stat.delta") or die "ERROR: Can't open $stat.delta: $!\n";
    <$din>; #skip the first line (paths to reference and query files)
    <$din>; #skip the second line ("NUCMER")
    my ($rid, $qid, $rstart, $rstop, $qstart, $qstop);
    #my @delta_array;
    my %delta_hash;
    my @deltas;
    while (my $line = <$din>){
        chomp $line;
        if ($line =~ m/^>/){
            if ($rid){
                (my $rsubid) = $rid =~ m/^#([^#]+)#/;
                (my $qsubid) = $qid =~ m/^#([^#]+)#/;
                push @{$delta_hash{$rid}{$qsubid}}, ([$qid, $rstart, $rstop, $qstart, $qstop, "N", [@deltas]]);
                #push @delta_array, ([$rid, $qid, $rstart, $rstop, $qstart, $qstop, "N", [@deltas]]);
                #change ref to query and query to ref to make sure we have all alignments
                my @inv_delta = map {$_ * -1} @deltas; #set all the distances to their inverse
                my ($n_rid, $n_qsubid, $n_qid, $n_rstart, $n_rstop, $n_qstart, $n_qstop) = ($qid, $rsubid, $rid, $qstart, $qstop, $rstart, $rstop);
                my $rev = "N"; #we'll keep track of whether the alignment was on the opposite strand so we can reverse the sequence later
                if ($qstart > $qstop){
                    ($n_rstart, $n_rstop, $rev) = ($qstop, $qstart, "Y");
                }
                push @{$delta_hash{$n_rid}{$n_qsubid}}, ([$n_qid, $n_rstart, $n_rstop, $n_qstart, $n_qstop, $rev, [@inv_delta]]);
                #push @delta_array, ([$n_rid, $n_qid, $n_rstart, $n_rstop, $n_qstart, $n_qstop, $rev, [@inv_delta]]);
            }
            @deltas = ();
            ($rid, $qid) = $line =~ m/^>(\S+)\s+(\S+)/;
            next;
        }
        if ($line =~ m/^(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/){
            ($rstart, $rstop, $qstart, $qstop) = ($1, $2, $3, $4);
        } else {
            next if $line =~ m/^\s*$/; #skip blank lines;
            push @deltas, $line;
        }
    }
    close ($din);
    if ($rid){
        (my $rsubid) = $rid =~ m/^#([^#]+)#/;
        (my $qsubid) = $qid =~ m/^#([^#]+)#/;
        push @{$delta_hash{$rid}{$qsubid}}, ([$qid, $rstart, $rstop, $qstart, $qstop, "N", [@deltas]]);
        #push @delta_array, ([$rid, $qid, $rstart, $rstop, $qstart, $qstop, "N", [@deltas]]);
        #change ref to query and query to ref to make sure we have all alignments
        my @inv_delta = map {$_ * -1} @deltas; #set all the distances to their inverse
        my ($n_rid, $n_qsubid, $n_qid, $n_rstart, $n_rstop, $n_qstart, $n_qstop) = ($qid, $rsubid, $rid, $qstart, $qstop, $rstart, $rstop);
        my $rev = "N";
        if ($qstart > $qstop){
            ($n_rstart, $n_rstop, $rev) = ($qstop, $qstart, "Y");
        }
        push @{$delta_hash{$n_rid}{$n_qsubid}}, ([$n_qid, $n_rstart, $n_rstop, $n_qstart, $n_qstop, $rev, [@inv_delta]]);
        #push @delta_array, ([$n_rid, $n_qid, $n_rstart, $n_rstop, $n_qstart, $n_qstop, $rev, [@inv_delta]]);
        undef @deltas;
    }
    
    #now, read the backbone coordinates back in
    my @bcoords;
    open (my $cin, "<$stat.backbone_coords.txt") or die "ERROR: Can't open $stat.backbone.coords.txt: $!\n";
    <$cin>; #skip the first line
    while (my $line = <$cin>){
        chomp $line;
        next if $line =~ m/^\s*$/; #skip blank lines
        my ($r_cont, $r_start, $r_stop, $r_gen) = split("\t", $line);
        my $r_full_id = "#".$r_gen."#".$r_cont;
        push @bcoords, ([$r_full_id, $r_gen, $r_start, $r_stop]);
    }
    close ($cin);
    
    #now let's go through one genome at a time, pull out the sequences, and add insertions and deletions
    my @outseqs;
    for my $i (0 .. $#order){
        my $gen = $order[$i];
        my $seq_count = $seq_counts{$gen};
        my $fullseq;
        open (my $sin, "<tmp.seq.$seq_count.txt") or die "ERROR: Can't open tmp.seq.$seq_count.txt: $!\n";
        chomp ($fullseq = <$sin>);
        close ($sin);
        my $last_gen;
        foreach my $slice (@bcoords){
            my ($rid, $rgen, $rstart, $rstop) = @{$slice};
            if ($gen eq $rgen){ #get the reference sequence
                my $offset = $rstart - $contig_starts{$rid};
                my $dist = ($rstop - $rstart) + 1;
                my $seq = substr($fullseq, $offset, $dist);
                push @{$outseqs[$i]}, ([$seq, $rstart, $rstop, "+"]);
            } else {
                if ($delta_hash{$rid}{$gen}){
                    my @tmp = @{$delta_hash{$rid}{$gen}};
                    foreach my $slice (@tmp){
                        my ($qid, $start, $stop, $qstart, $qstop, $rev) = @{$slice};
                        next unless ($start <= $rstart and $stop <= $rstop);
                        my @deltas = @{@{$slice}[6]};
                        #get the sequence
                        my $offset = $qstart - $contig_starts{$qid};
                        my $dist = $qstop - $qstart + 1;
                        my $seq = substr($fullseq, $offset, $dist);
                        my $dir = "+";
                        if (@deltas > 1){
                            my $newseq;
                            my $pos = 0;
                            foreach my $delta (@deltas){
                                if ($delta > 0){
                                    $newseq .= substr($seq, $pos, $delta - 1) unless $delta == 1;
                                    $newseq .= "-";
                                    $pos += ($delta - 1);
                                } elsif ($delta < 0){
                                    $newseq .= substr($seq, $pos, ($delta * -1) - 1 );
                                    $pos += ($delta * -1);
                                } else {
                                    my $fdist = (length($seq) - $pos);
                                    $newseq .= substr($seq, $pos, $fdist);
                                }
                            }
                            $seq = $newseq;
                        }
                        if ($rev eq "Y"){ 
                            $seq = reverse($seq);
                            $seq =~ tr/ACTGRYKMBVDHactgrykmbvdh/TGACYRMKVBHDtgacyrmkvbhd/;
                            $dir = "-";
                        }
                        push @{$outseqs[$i]}, ([$seq, $qstart, $qstop, $dir]);
                        last;
                    }
                }
            }
        }
    }
    

    
    
}

#clean up
for my $i (1 .. $nog){
    unlink ("tmp.seq.$i.txt") if -e "tmp.seq.$i.txt";
    unlink ("tmp.crd.$i.txt") if -e "tmp.crd.$i.txt";
    if ($loci_tmp_files[$i]){
        unlink("$loci_tmp_files[$i]") if -e "$loci_tmp_files[$i]";
    }
}
unlink ("$stat.delta") if (-e "$stat.delta" and $opt_w);

print "Done\n" unless $opt_w;

#--------------------------
sub start_next_process {
    my ($gcode, $proc_eye) = ($_[0], $_[1]);
    $part++;
    my $pid = fork;
    if (0 == $pid){
        my $this_part = $part;
        my $ref = $seq_ids{$gcode};
        my $gen_size = $g_size{$ref};
        my @all_gens;
        for my $i (1 .. $nog){
            push @all_gens, $i;
        }
        
        ## below is for progress indicator
        print STDERR sprintf("\rAnalysis # %3s is loading coords. ", $this_part) unless $opt_w;      
        
        #load the genome-specific alignment coordinates
        # will also group overlaps in a query (if they exist) at this time
        my %hash;
        open (my $in, "<tmp.crd.$gcode.txt") or exit(9);
        my ($lref, $lqry, $lstart, $lstop);
        while (my $line = <$in>){
            chomp $line;
            my ($ref, $qry, $start, $stop) = split("\t", $line);
            if (!$lref){
                ($lref, $lqry, $lstart, $lstop) = ($ref, $qry, $start, $stop);
                next;
            }
            #check for overlapping alignments in the same query
            if ($ref eq $lref){
                if ($qry eq $lqry){
                    if ($start <= $lstop){
                        if ($stop > $lstop){
                            $lstop = $stop;
                        }
                        next; #if the coordinate range is within the last coordinate range, will just go to the next
                    }
                }
            }
            push @{$hash{$lref}{$lqry}}, ([$lstart, $lstop]);
            ($lref, $lqry, $lstart, $lstop) = ($ref, $qry, $start, $stop); 
        }
        close ($in);
        push @{$hash{$lref}{$lqry}}, ([$lstart, $lstop]) if $lref;
        
        ## below is for progress indicator
        my ($last_pct, $pct_count) = (0) x 2;
        print STDERR sprintf("\rAnalysis # %3s is %3s pct complete.", $this_part, $last_pct) unless $opt_w; 
        
        my %proc_hash;
        my @core;
        my @out;
        my @acc;
        my @pan;
        my @contigs = @{$contig_ids{$ref}};
        my @ntbm = @to_use[0..($proc_eye - 1)]; #set the list of genome numbers that HAVE to be in %missh for this section to be included in the pangenome
        for my $l (0 .. $#contigs){
            my ($c_id, $c_size) = @{$contigs[$l]};
            my ($seg_start, $seg_stop);
            my ($core_start, $core_stop);
            my ($acc_start, $acc_stop);
            my ($pan_start, $pan_stop);
            
            my @ranges;
            ##load ranges array 
            for my $j (0 .. $#order){
                my @rng = (0,0);
                if (my $next = shift @{$hash{$c_id}{$order[$j]}}){
                    @rng = @{$next};
                }
                push @ranges, ([$rng[0], $rng[1]]);
            }
            my ($pos_counts, $pos_gens);
            
            for my $j (1 .. $c_size){
                #below is for progress indicator
                $pct_count++;
                my $pct = int(100* ($pct_count / $gen_size));
                if ($pct % 1 == 0 and $pct != $last_pct){
                    print STDERR sprintf("\rAnalysis # %3s is %3s pct complete.", $this_part, $pct) unless $opt_w;
                    $last_pct = $pct;
                }
                
                $pos_counts = 1;
                $pos_gens = ",$gcode,";
                my %missh;
                for my $k (0 .. $#ranges){
                    my $gen = $order[$k];
                    next if $gen eq $ref;
                    my $gencode = $seq_counts{$gen};
                    my ($start, $stop) = @{$ranges[$k]};
                    unless (!$hash{$c_id}{$gen} or @{$hash{$c_id}{$gen}} == 0){
                        while ($stop < $j and @{$hash{$c_id}{$gen}} > 0){
                            ($start, $stop) = @{shift @{$hash{$c_id}{$gen}}};
                            $ranges[$k] = [$start, $stop];
                        } 
                    }
                    if ($j >= $start and $j <= $stop){
                        $pos_counts++;
                        $pos_gens .= $gencode . ",";
                    } else {
                        $missh{$gencode}++;
                    }
                }
                if ($opt_e){
                    if ($pos_counts == $g_count){
                        $proc_hash{$gcode}{$pos_counts}{",0,"}++;
                    } else {
                        my @miss = sort{$a <=> $b} keys %missh;
                        my $missing = "," . join(",", @miss) . ",";
                        $proc_hash{$gcode}{$pos_counts}{$missing}++;
                    }
                }
                ## build pangenome regions (if requested)
                if ($opt_n){
                    if ($proc_eye > 0){ #will output entire first genome into pangenome sequence
                        my $miss_count = 0;
                        foreach my $test (@ntbm){
                            $miss_count++ if $missh{$test};
                        }
                        if ($miss_count == scalar @ntbm){ #all genomes after the first will contribute to pangenome only if there are in a subset of the genomes and are missing from all previously-queried genomes
                            if (!$pan_start){
                                ($pan_start, $pan_stop) = ($j, $j);
                            } elsif ($j > ($pan_stop + 1)) {
                                push @pan, ([$pan_start, $pan_stop, $c_id]);
                                ($pan_start, $pan_stop) = ($j, $j);
                            } else {
                                $pan_stop = $j;
                            }
                        }
                    }
                }
                ## build core and accessory regions
                if ($pos_counts < $min_gen){
                    if (!$acc_start){
                        ($acc_start, $acc_stop) = ($j, $j);
                    } elsif ($j > ($acc_stop + 1)) {
                        push @acc, ([$acc_start, $acc_stop, $c_id]);
                        ($acc_start, $acc_stop) = ($j, $j);
                    } else {
                        $acc_stop = $j
                    }
                    next;
                } else {
                    if (!$core_start){
                        ($core_start, $core_stop) = ($j, $j);
                    } elsif ($j > ($core_stop + 1)) {
                        push @core, ([$core_start, $core_stop, $c_id]);
                        ($core_start, $core_stop) = ($j, $j);
                    } else {
                        $core_stop = $j;
                    }
                }
                next if $pos_counts == $g_count;
                next if ($proc_eye == 0 or $proc_eye > $abs);  ## only run this analysis on genomes that aren't the first reference or those for which core genome sequence will be output
                my $miss_count = 0;
                for my $k (0 .. ($proc_eye - 1)){
                    my $missing = $to_use[$k];
                    $miss_count++ if ($pos_gens !~ m/,$missing,/);
                }
                if ($miss_count == $proc_eye){
                    if (!$seg_start){
                        ($seg_start, $seg_stop) = ($j, $j);
                        next;
                    }
                    if ($j > ($seg_stop + 1)){
                        push @out, ([$seg_start, $seg_stop, $c_id]);
                        ($seg_start, $seg_stop) = ($j, $j);
                        next;
                    }
                    $seg_stop = $j;
                } 
            }
            
            push @core, ([$core_start, $core_stop, $c_id]) if ($core_start);
            push @out, ([$seg_start, $seg_stop, $c_id]) if ($seg_start);
            push @acc, ([$acc_start, $acc_stop, $c_id]) if ($acc_start);
            if ($opt_n){
                if ($proc_eye == 0){ #if this is the reference genome, include all contig sequences in pangenome
                    ($pan_start, $pan_stop) = (1, $c_size);
                }
                push @pan, ([$pan_start, $pan_stop, $c_id]) if ($pan_start);
            }
        }
        if ($opt_e){
            my $outfile_ph = "$$.tmp_out.ph.txt";
            open (my $ph_out, ">$outfile_ph") or exit(1);
            while( my ($w, $x) = each %proc_hash){
                while (my ($y, $z) = each %$x){
                    while (my ($u, $v) = each %$z){
                        print $ph_out "$w\t$y\t$u\t$v\n";
                    }
                }
            }
            close ($ph_out); 
        }
        print STDERR sprintf("\rAnalysis # %3s is %3s pct complete.", $this_part, 100) unless $opt_w;
        
        #read in annotation information, create markers for start and stop coordinates
        @loci_order = ();
        %locuslengs = ();
        %loci_starts = ();
        %loci_stops = ();
        if (@loci_tmp_files){
            if ($loci_tmp_files[$gcode]){
                my $infile = $loci_tmp_files[$gcode];
                open (my $loc_in, "<$infile") or exit(2);
                while (my $line = <$loc_in>){
                    chomp $line;
                    my @tmparray = split("\t", $line);
                    my ($lid, $contig, $start, $stop) = @tmparray;
                    my $leng = ($stop - $start + 1);
                    $locuslengs{$lid}{$contig} = $leng;
                    push @loci_order, ([@tmparray]);
                }
                close ($loc_in);
                @loci_order = sort{$a->[1] cmp $b->[1] || $a->[2] <=> $b->[2] || $a->[3] <=> $b->[3]} @loci_order; #loci should be in order by start coordinates, but will sort just in case
                for my $i (0 .. $#loci_order){
                    my ($x1, $contig, $start, $stop) = @{$loci_order[$i]};
                    $loci_starts{$contig}{$start} = $i unless $loci_starts{$start};
                    $loci_stops{$contig}{$stop} = $i;
                }
            }
        }        
        
        #get characteristics of accessory sequences, output core genome sequence and genes and accessory coords and genes
        #will continue to output a "backbone" sequence, but will also output strain "core" and strain "accessory" sequences for each genome
        
        my $outfile_st = "$$.tmp_out.st.txt";
        open (my $out_st, ">$outfile_st") or exit(3);
        
        #first, process accessory sequences.
        my @params = ($gcode, "accessory", $proc_eye, $backlen, $this_part);
        push @params, @acc;
        my $stat_string = post_process(@params);
        exit (4) if $stat_string =~ m/^ERROR/;
        print $out_st "$stat_string";
        #next, process core sequences. If this is the reference genome, will simultaneously process backbone
        @params = ($gcode, "core", $proc_eye, $backlen, $this_part);
        push @params, @core;
        $stat_string = post_process(@params);
        exit (5) if $stat_string =~ m/^ERROR/;
        print $out_st "$stat_string";
        #next, process pangenome sequence (if requested)
        if ($opt_n){
            @params = ($gcode, "pan", $proc_eye, $backlen, $this_part);
            push @params, @pan;
            $stat_string = post_process(@params);
            exit (7) if $stat_string =~ m/^ERROR/;
        }
        #finally, if this is not a reference genome, but should be included in the backbone, output the backbone segments
        if ($proc_eye > 0 and $proc_eye <= $abs){
            @params = ($gcode, "out", $proc_eye, $backlen, $this_part);
            push @params, @out;
            $stat_string = post_process(@params);
            exit (6) if $stat_string =~ m/^ERROR/;
        }
        close ($out_st);
        exit(0);
    }
    my $outfile_ph = "$pid.tmp_out.ph.txt";
    my $outfile_st = "$pid.tmp_out.st.txt";
    $children->{$pid}->{'start'} = time;
    $children->{$pid}->{'out_ph'} = $outfile_ph;
    $children->{$pid}->{'out_st'} = $outfile_st;
    $children->{$pid}->{'num'}   = $part;
    $children->{$pid}->{'done'}  = 0;
    $children->{$pid}->{'gcode'}  = $gcode;
    $children->{$pid}->{'eye'}  = $proc_eye;
    push @child_order, $pid;
    $num_procs_running++;
    return ();
}

sub stats_old {
    my @vals = @{$_[0]};
    my $num = scalar @vals;
    return (0,0,0,0,0,0,0,0) if $num == 0;
    my ($tot, $num_gt100, $tot_gt100, $avg_gt100) = (0) x 4;
    @vals = sort {$a <=> $b} @vals;
    my ($low, $high) = ($vals[0], $vals[$#vals]);
    for my $i (0 .. $#vals){
        my $val = $vals[$i];
        $tot += $val;
        if ($val >= 100){
            $num_gt100 ++;
            $tot_gt100 += $val;
        }
    }
    my $avg = sprintf("%.1f", ($tot / $num));
    $avg_gt100 = sprintf("%.1f", ($tot_gt100 / $num_gt100)) if $num_gt100 > 0;
    return ($num, $tot, $avg, $num_gt100, $tot_gt100, $avg_gt100, $low, $high);
}
sub stats{
	my @lengths = @{$_[0]};
	my ($sum, $num, $min, $maxi, $mean, $median, $mode, $mode_freq);
        return (0,0,0,0,0,0,0,0) if (scalar @lengths == 0);
	my %seen;
	my @sorted_leng = sort {$a <=> $b} @lengths;
	for my $i (0 .. $#sorted_leng){
		$sum += $sorted_leng[$i];
		$seen{$sorted_leng[$i]}++;
	}
	$num = $#sorted_leng + 1;
	$min = $sorted_leng[0];
	$maxi = $sorted_leng[$#sorted_leng];
	$mean = $sum/$num;
	my $rounded_mean = sprintf("%.2f", $mean);
	my @modes;
	foreach my $leng (sort {$seen{$b} <=> $seen{$a}} keys %seen) {
		push @modes, ([$leng, $seen{$leng}]);
	}
	$mode = $modes[0][0];
	$mode_freq = $modes[0][1];
	my $mid = int @sorted_leng/2;
	if (@sorted_leng % 2){
		$median = $sorted_leng[$mid];
	} else {
		$median = ($sorted_leng[$mid-1] + $sorted_leng[$mid])/2;
	}
	return ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq);
}

sub gc_content {
    my $seq = $_[0];
    return ("0.0000") unless $seq;
    my $count = 0;
    while ($seq =~/G|C|g|c/g) {
        $count++;
    }
    my $n_count = 0;
    while ($seq =~/[^ACTGactg]/g) {
        $n_count++;
    }
    my $length = length($seq);
    my $leng_no_ns = $length - $n_count;
    $leng_no_ns = 1 if ($leng_no_ns == 0);
    my $gc_long = 100*($count/$leng_no_ns);
    my $gc = sprintf("%.4f", $gc_long);
    return ($gc);
}

#need to clear %locusids, %locuslengs, and @loci_order each time.
sub post_process {
    my $gcode = shift;
    my $type = shift;
    my $proc_eye = shift; 
    my $minlen = shift;
    my $this_part = shift; #for status update
    my @array = @_;
    my $ref = $seq_ids{$gcode};
    my $gen_size = $g_size{$ref};
    my $seq;
    open (my $in, "<tmp.seq.$gcode.txt") or return ("ERROR: Can't open tmp.seq.$gcode.txt: $!\n");
    while (my $line = <$in>){
        chomp $line;
        $seq .= $line;
    }
    close ($in);
    
    #status update
    my $stat_type = "Core";
    $stat_type = "Accessory" if $type eq "accessory";
    $stat_type = "Backbone" if $type eq "out";
    $stat_type = "Pangenome" if $type eq "pan";
    print STDERR sprintf("\r%9s # %3s is %3s pct complete.", $stat_type, $this_part, 0) unless $opt_w;
    
    my $fileid = "$stat.$ref.$type";
    my $seqid = "$ref\_$type\_";
    
    open (my $out_cor, ">$fileid\_coords.txt") or return "ERROR: Can't open $fileid\_coords.txt: $!\n";
    print $out_cor "contig_id\tcontig_length\tstart\tstop\tout_seq_id\n" unless ($type eq "out" or $type eq "pan");
    open (my $out_seq, ">$fileid.fasta") or return "ERROR: Can't open $fileid.fasta: $!\n";
    my ($bbone_cor, $bbone_seq); #to be used to avoid the time of double-processing core and backbone for the reference genome
    my $bbonefileid = "$stat.$ref.out";
    if ($type eq "core" and $proc_eye == 0){
        open ($bbone_cor, ">$bbonefileid\_coords.txt") or return "ERROR: Can't open $bbonefileid\_coords.txt: $!\n";
        open ($bbone_seq, ">$bbonefileid.fasta") or return "ERROR: Can't open $bbonefileid.fasta: $!\n";
    }
    my %loci;
    my @lengs;
    my $count = 0;
    my $tot_seq = "";
    my $last_pct_done = 0;
    for my $i (0 .. $#array){
        
        #status update
        my $pct_done = 0;
        $pct_done = int(100 * ($i / $#array)) if $#array > 0;
        if ($pct_done != $last_pct_done and $pct_done % 10 == 0){
            print STDERR sprintf("\r%9s # %3s is %3s pct complete.", $stat_type, $this_part, $pct_done) unless $opt_w;
        }
        $last_pct_done = $pct_done;
        
        my ($o_start, $o_stop, $c_id) = @{$array[$i]};
        my $c_size = $contig_lengs{$c_id};
        my $out_c_id = $c_id;
        $out_c_id =~ s/^\#[^#]*\#//;
        my $b_leng = ($o_stop - $o_start) + 1;
        #output sequences
        my $out_id = "-";
        my $bbone_id = "0,0";
        if ($b_leng >= $minlen){
            my $offset = $contig_starts{$c_id} + $o_start - 2;
            my $b_seq = substr($seq, ($offset), $b_leng);
            #trim ambiguous bases from front and back of sequence. If the size of the remaining sequence is smaller than the minimum length, drop it.
            my ($lead, $tail) = (0) x 2;
            (my $lead_seq) = $b_seq =~ m/^([NnXx-]*)/;
            if ($lead_seq){
                $lead = length($lead_seq);
            }
            if ($lead == $b_leng){ #skip if the entire sequence was N's
                #output coordinates
                if ($type eq "out" or $type eq "pan"){
                    print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\n"
                } else {
                    print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\n";
                    print $bbone_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\t$bbone_id\n" if ($type eq "core" and $proc_eye == 0);
                }
                next;
            }
            (my $tail_seq) = $b_seq =~ m/([NnXx-]*$)/;
            if ($tail_seq){
                $tail = length($tail_seq);
            }
            if ($lead > 0 or $tail > 0){
                my $dist = length($b_seq) - ($lead + $tail);
                if ($dist < $minlen){
                    #output coordinates
                    if ($type eq "out" or $type eq "pan"){
                        print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\n"
                    } else {
                        print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\n";
                        print $bbone_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\t$bbone_id\n" if ($type eq "core" and $proc_eye == 0);
                    }
                    next;
                }
                $b_seq = substr($b_seq, $lead, $dist);
                $o_start += $lead;
                $o_stop -= $tail;
            }
            push @lengs, length($b_seq);
            $count++;
            $out_id = $seqid.sprintf("%04d", $count)."_length\_$b_leng";
            $out_id = "$count,$b_leng" if ($type eq "out" or $type eq "pan");
            print $out_seq ">$out_id\n$b_seq\n";
            if ($type eq "core" and $proc_eye == 0){
                $bbone_id = "$count,$b_leng";
                print $bbone_seq ">$bbone_id\n$b_seq\n";
            }
            $tot_seq .= $b_seq;
            #output coordinates
            if ($type eq "out" or $type eq "pan"){
                print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\t$out_id\n"
            } else {
                print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$out_id\n";
                print $bbone_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\t$bbone_id\n" if ($type eq "core" and $proc_eye == 0);
            }
            #collect annotation information, if present
            my @locusids;
            if (@loci_order){ #first, will set a range of loci to screen so that we don't have to search the entire array each time
                my ($start_rec, $stop_rec);
                my $cleng = $contig_lengs{$c_id};
                for my $i (reverse 1 .. $o_start){
                    if (defined $loci_stops{$c_id}{$i}){
                        $start_rec = $loci_stops{$c_id}{$i};
                        last;
                    }
                }
                if (!defined $start_rec){
                    for my $i (reverse 1 .. $o_start){
                        if (defined $loci_starts{$c_id}{$i}){
                            $start_rec = $loci_starts{$c_id}{$i};
                            last;
                        }
                    }
                }
                if (!defined $start_rec){
                    for my $i ($o_start .. $o_stop){
                        if (defined $loci_starts{$c_id}{$i}){
                            $start_rec = $loci_starts{$c_id}{$i};
                            last;
                        }
                    }
                }
                if (defined $start_rec){ #no point in looking the other way if there's no start_rec
                    for my $i ($o_stop .. $cleng){
                        if (defined $loci_starts{$c_id}{$i}){
                            $stop_rec = $loci_starts{$c_id}{$i};
                            last;
                        }
                    }
                    if (!defined $stop_rec){
                        for my $i (reverse $o_start .. $o_stop){
                            if (defined $loci_starts{$c_id}{$i}){
                                $stop_rec = $loci_starts{$c_id}{$i};
                                last;
                            }
                        }
                    }
                    $stop_rec = $start_rec if !defined $stop_rec;
                    @locusids = @loci_order[$start_rec .. $stop_rec];
                }
            }
            if (@locusids){
                foreach my $slice (@locusids){
                    my ($lid, $c_id, $l_start, $l_stop) = @{$slice};
                    next if $l_stop < $o_start;
                    next if $l_start > $o_stop;
                    my ($over_front, $over_back) = (0) x 2; #keep track of how much a gene overhangs the end of the segment
                    if ($l_start < $o_start){
                        $over_front = $o_start - $l_start;
                        if ($l_stop <= $o_stop){
                            push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $o_start, $l_stop, $bbone_id, $over_front, $over_back]);
                        } else {
                            $over_back = $l_stop - $o_stop;
                            push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $o_start, $o_stop, $bbone_id, $over_front, $over_back]);
                        }
                        next;
                    }
                    if ($l_stop > $o_stop){
                        $over_back = $l_stop - $o_stop;
                        push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $l_start, $o_stop, $bbone_id, $over_front, $over_back]);
                        next;
                    }
                    push @{$loci{$lid}{$c_id}}, ([$out_id, $o_start, $l_start, $l_stop, $bbone_id, $over_front, $over_back]); #don't need an if statement. Everything left after the above ifs are loci with borders within o_star and o_stop, inclusive
                }
            }
        } else {
            #output coordinates
            if ($type eq "out" or $type eq "pan"){
                print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\n"
            } else {
                print $out_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\n";
                print $bbone_cor "$out_c_id\t$c_size\t$o_start\t$o_stop\t$ref\n" if ($type eq "core" and $proc_eye == 0);
            }
        }
    }
    print STDERR sprintf("\r%9s # %3s is %3s pct complete.", $stat_type, $this_part, 100) unless $opt_w; #status update
    close ($out_cor);
    close ($out_seq);
    if ($type eq "core" and $proc_eye == 0){
        close ($bbone_cor);
        close ($bbone_seq);
    }
    my $loci_count = 0;
    if (%loci){
        #sort and output the annotation information to file
        open (my $out_gen, ">$fileid\_loci.txt") or die "ERROR: Can't open $fileid\_loci.txt: $!\n";
        print $out_gen "locus_id\tgen_contig_id\tgen_contig_start\tgen_contig_stop\tstrand\tout_seq_id\tout_seq_start\tout_seq_stop\tpct_locus\toverhangs\tproduct\n" unless ($type eq "out" or $type eq "pan");
        my $bbone_gen;
        if ($type eq "core" and $proc_eye == 0){
            open ($bbone_gen, ">$bbonefileid\_loci.txt") or die "ERROR: Can't open $stat.backbone_loci.txt: $!\n";
        }
        foreach my $slice (@loci_order){
            my @slice_arr = @{$slice};
            #$lid, $contig, $start, $stop
            my ($locus, $contig) = @slice_arr;
            next unless $loci{$locus}{$contig};
            my $dir = $slice_arr[4];
            my $prod = "NoAnnotation";
            $prod = $slice_arr[5] if $slice_arr[5];
            my $out_contig = $contig;
            $out_contig =~ s/^\#[^#]*\#//;
            my $locusleng = $locuslengs{$locus}{$contig};
            my @hits = @{$loci{$locus}{$contig}};
            @hits = sort{$a->[0] cmp $b->[0]}@hits;
            #calculate total percentage of hit
            my $totleng = 0;
            foreach (@hits){
                $totleng += (${$_}[3] - ${$_}[2] + 1);
            }
            #if ($type eq "core" or $type eq "out"){
            #    next if ((100 * $totleng / $locusleng) < $minpct);
            #} else {
            #    next if ((100 * $totleng / $locusleng) <= $minpct);
            #}
            my $tot_pct = 0;
            foreach my $hit (@hits){
                my ($id, $hoffset, $hstart, $hstop, $bbid, $over_front, $over_back) = @{$hit};
                my $pct = sprintf("%.2f", 100 * (($hstop - $hstart + 1) / $locusleng));
                $tot_pct += $pct;
                my ($cstart, $cstop) = ($hstart - $hoffset + 1, $hstop - $hoffset + 1);
                print $out_gen "$locus\t$out_contig\t$hstart\t$hstop\t$dir\t$id\t$cstart\t$cstop\t$pct\t$over_front,$over_back";
                print $out_gen "\t$ref" if ($type eq "out" or $type eq "pan");
                print $out_gen "\t$prod" if $prod;
                print $out_gen "\n";
                if ($type eq "core" and $proc_eye == 0){
                    print $bbone_gen "$locus\t$out_contig\t$hstart\t$hstop\t$dir\t$bbid\t$cstart\t$cstop\t$pct\t$over_front,$over_back\t$ref";
                    print $bbone_gen "\t$prod" if $prod;
                    print $bbone_gen "\n";
                }
            }
            $loci_count++ if $tot_pct >= $minpct; #for the purposes of counting, will use a 50% cutoff (so that only if the gene is in the majority in either core of accessory it will be counted as belonging to that group)
        }
        close ($out_gen);
        close ($bbone_gen) if $bbone_gen;
    }
    #get stats (core or accessory. Backbone will be calculated later)
    my $stat_string = "x\n";
    if ($type eq "core" or $type eq "accessory"){
        my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq) = stats(\@lengs);
        my $gc = gc_content($tot_seq);
        $stat_string = "$gcode\t$ref\t$gen_size\t$type\t$sum\t$gc\t$num\t$min\t$maxi\t$rounded_mean\t$median";
        if ($opt_l){
            $stat_string .= "\t$loci_count";
        }
        $stat_string .= "\n";
    }
    return($stat_string);
}

sub process_final{
    my @gcodes = @{$_[0]};
    my $type = shift @gcodes;
    my $outtype = "backbone"; #if "out"
    $outtype = "pangenome" if $type eq "pan";
    #open files for writing
    open (my $out_cor, ">$stat.$outtype\_coords.txt") or die "ERROR: Can't open $stat.$outtype\_coords.txt: $!\n";
    print $out_cor "contig_id\tcontig_length\tstart\tstop\tsource_gen\tout_seq_id\n";
    open (my $out_seq, ">$stat.$outtype.fasta") or die "ERROR: Can't open $stat.$outtype.fasta: $!\n";
    my $out_gen;
    if (@loci_tmp_files){
        open ($out_gen, ">$stat.$outtype\_loci.txt") or die "ERROR: Can't open $stat.$outtype\_loci.txt: $!\n";
        print $out_gen "locus_id\tgen_contig_id\tgen_contig_start\tgen_contig_stop\tstrand\tout_seq_id\tout_seq_start\tout_seq_stop\tpct_locus\tsource_gen\toverhangs\tproduct\n" 
    }
    #now process backbone sequences and statistics
    my $out_count = 0;
    my $out_loci_count = 0;
    my $out_seqs;
    my @out_lengs;
    foreach my $gcode (@gcodes){
        my $ref = $seq_ids{$gcode};
        #first we'll process the backbone sequences
        my %out_ids;
        if (-e "$stat.$ref.$type.fasta"){
            open (my $in, "<$stat.$ref.$type.fasta") or die "ERROR: Can't open $stat.$ref.$type.fasta: $!\n";
            while (<$in>){
                #read the file two lines at a time
                chomp (my $id = $_);
                chomp (my $seq = <$in>);
                last unless $id =~ /^>/;
                my $subid = substr($id, 1);
                my ($cnt, $leng) = split(",", $subid);
                $out_count++;
                my $out_id = "$outtype\_".sprintf("%04d", $out_count)."_length\_$leng";
                print $out_seq ">$out_id\n$seq\n";
                $out_ids{$subid} = $out_id;
                $out_seqs .= $seq;
                push @out_lengs, $leng;
            }
            close ($in);
            unlink "$stat.$ref.$type.fasta";
        }
        #next we'll process the locus information, if it exists
        if (-e "$stat.$ref.$type\_loci.txt"){
            open (my $in, "<$stat.$ref.$type\_loci.txt") or die "ERROR: Can't open $stat.$ref.$type\_loci.txt: $!\n";
            my %unique_loci;
            while (my $line = <$in>){
                chomp $line;
                my ($locus,$out_contig,$hstart,$hstop,$dir,$id,$cstart,$cstop,$pct,$ref,$overhang,$prod) = split("\t", $line);
                if ($out_ids{$id}){
                    print $out_gen "$locus\t$out_contig\t$hstart\t$hstop\t$dir\t$out_ids{$id}\t$cstart\t$cstop\t$pct\t$overhang\t$ref";
                    print $out_gen "\t$prod" if $prod;
                    print $out_gen "\n";
                } else {
                    die "ERROR: The temporary backbone sequence ID $id was not seen in the sequences\n";
                }
                $unique_loci{$locus}+=$pct;
            }
            close ($in);
            foreach my $loc (keys %unique_loci){
                $out_loci_count++ if $unique_loci{$loc} >= $minpct;
            }
            #$bbone_loci_count += scalar(keys %unique_loci) if %unique_loci;
            unlink("$stat.$ref.$type\_loci.txt");
        }
        #last, we'll combine the coordinate information
        if (-e "$stat.$ref.$type\_coords.txt"){
            open (my $in, "<$stat.$ref.$type\_coords.txt") or die "ERROR: Can't open $stat.$ref.$type\_coords.txt: $!\n";
            while (my $line = <$in>){
                chomp $line;
                next if $line =~ m/^\s*$/; #skip blank lines
                my @outline = split("\t", $line);
                if ($outline[5]){
                    if (my $out_id = $out_ids{$outline[5]}){
                        $outline[5] = $out_id;
                    }
                }
                print $out_cor join("\t", @outline), "\n";
            }
            close ($in);
            unlink("$stat.$ref.$type\_coords.txt");
        }
    }
    close ($out_cor);
    close ($out_seq);
    close ($out_gen) if $out_gen;
    
    my ($sum, $num, $min, $maxi, $rounded_mean, $median, $mode, $mode_freq) = stats(\@out_lengs);
    my $gc = gc_content($out_seqs);
    if ($opt_w){
        print "<br><strong>Core Genome Statistics:</strong><br>\n" if $type eq "out";
        print "<strong>Pangenome Statistics:</strong><br>\n" if $type eq "pan";
    } else {
        print "\nCore Genome Statistics:\n" if $type eq "out";
        print "\nPangenome Statistics:\n" if $type eq "pan";
    }
    print "Number of segments >= $backlen bp: $num\n";
    print "<br>\n" if $opt_w;
    print "Total bp: $sum\n";
    print "<br>\n" if $opt_w;
    print "GC content: $gc%\n";
    print "<br>\n" if $opt_w;
    print "Shortest segment length: $min\n";
    print "<br>\n" if $opt_w;
    print "Longest segment length: $maxi\n";
    print "<br>\n" if $opt_w;
    print "Average segment length: $rounded_mean\n";
    print "<br>\n" if $opt_w;
    print "Median segment length: $median\n";
    print "<br>\n" if $opt_w;
    print "Mode segment length: $mode, Frequency: $mode_freq\n";
    print "<br>\n" if $opt_w;
    print $stats "-\t-\t-\t$outtype\t$sum\t$gc\t$num\t$min\t$maxi\t$rounded_mean\t$median";
    if ($opt_l){
            print "Number of CDS (if at least 50% by length): $out_loci_count\n";
            print "<br>\n" if $opt_w;
            print $stats "\t$out_loci_count";
    }
    print "\n";
    print "<br>\n" if $opt_w;
    print $stats "\n";
    return();
}
