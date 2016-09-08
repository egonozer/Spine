#!/usr/bin/perl

my $version = "0.2";

##Changes from v0.1.3 -> 0.2
# moved project to git
# Uses nucmer_backbone.pl version 0.3 (checks version before proceeding)
# Allows for multi-contig genbank files
# If mix of genbank and fasta files, will still output annotation information for genbank files
# Outputs sequence files for backbone (as before), but also sequences of accessory sequences and core sequences for each included genome
# Outputs annotation information for core, accessory, and backbone sequences (if annotations given)
# Keeps track of which strand genes are found on. Will help with figure drawing by ClustAGE
# change sequence name delimiters to "#" instead of "|". I think this will prevent potential problems with NCBI sequence and contig names. I'll also have nucmer_backbone.pl remove the delimiters automatically instead of outputting them.
# sequence headers with spaces in them will be truncated at the first space (since that's what nucmer is going to do anyway)
# in addition to locus IDs of coding sequences, will also output product names (if given)
# add in test for binary input files
# add in test for duplicate genome IDs
# Change to long options (Getopt::Long)
# Add Nucmer tuning options
# MAJOR CHANGE: output partial genes in the core and accessory lists, rather than assigning them either just to core or just to accesory based on a cutoff percentage. This will be to allow downstream programs (like ClustAGE) to identify variable portions of otherwise conserved core genes. 
# WISHLIST: Accept not just genbank files for annotation, but also other annotation types (ptt, gff, etc.) In this case both a sequence and annotation file would have to be given

##Changes from v0.1.2 -> 0.1.3
# Improve the gbk_convert subroutine (it was getting messed up on the Genome-Annotation-Data section)

##Changes from v0.1.1 -> 0.1.2
# Accept files with all types of line endings (Unix, Mac, PC)

##Changes from v0.1 -> v0.1.1
# Fixed bug entire reference genome list was being passed to nucmer_backbone.pl -r and -g instead of just the first genome
# Built in a little more error-checking for the -r option
# Fixed a little glitch in filetype guessing

my $license = "
    spine.pl
    Copyright (C) 2014 Egon A. Ozer

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

use strict;
use warnings;
use Cwd 'abs_path';
use File::Which;

$|++;

my $usage = "
spine.pl

This is a wrapper script to run the Spine algorithm

PREREQUISITES:
- Perl 5.10 or above
- MUMmer version 3.22 or above
- Mac OSX or Linux. No guarantees that this will work on Windows or other
  operating systems.

REQUIRED:
  -f or --file      file with list of input sequence files
                    This file should be formatted like so:
                    
                    path/to/file1<tab>unique_identifier<tab>fasta or gbk
                    path/to/file2<tab>unique_identifier<tab>fasta or gbk
                    
                    Example:
                    /home/seqs/PAO1.fasta   PAO1    fasta
                    /home/seqs/LESB58.gbk   LESB58  gbk
                    
                    The third column (fasta or gbk) is optional, but should
                    be given if your sequence files end with suffixes other
                    than \".fasta\" or \".gbk\".
                    
OPTIONS:
  -a or --pctcore   percentage of input genomes in which a region must be
                    found in order to be considered core
                    (default: 100)
  -g or --maxdist   maximum distance between core genome segments. Distances
                    less than this between adjacent segments will result in
                    combination of fragments with N's rather than separating
                    into two or more fragments.
                    (default: 10)
  -l  or --license  print license information and quit
  -m  or --nucpath  full path to folder containing MUMmer scripts and
                    executables, i.e. /home/applications/MUMmer/bin
                    (default: tries to find MUMmer in your PATH)
  -r or --refs      reference genome sequence(s) to use as primary output
                    source(s). This should be one or more integers corresponding
                    to the order of the genomes given in the file above, i.e.
                    \"1\" would use the first-listed sequence, \"3\" would use
                    the third-listed, etc. To prioritize multiple genome
                    sequences, separate the integers with commas, i.e. \"1,3\"
                    for giving sequence 1 the highest priority and sequence 3
                    the next-highest priority.
                    Reference sequences will serve as the source of backbone
                    sequences to be output, as well as the source of backbone
                    locus IDs, if applicable.
                    The number of reference genomes used will depend on the
                    definition of core genome given by option -a. For instance,
                    if core is determined from 10 input genomes and -a is given
                    as 100, then core sequence will only be taken from one
                    reference genome. If, for example, -a is given as 90 from
                    10 input genomes, then potentially two reference sequences
                    will be needed: The first for sequences present in all 10
                    genomes and for sequences present in 9 out of 10 genomes
                    including the first genome. The second reference sequence
                    would then be used as the source of all sequences present
                    in 9 out of 10 genomes, but not present in the first
                    reference genome.
                    (default: reference priority will be the same as the order
                    of genomes entered, with the first genome having the highest
                    priority and the last genome having the lowest priority)
  --pangenome       produce a pangenome sequence and characteristics from
                    sequences in the order given.
                    (default: no pangenome information will be output)
  -o or --prefix    Output prefix
                    (default: \"output\")        
  -p or --pctid     minimum percent identity for regions to be considered
                    homologous
                    (default: 85) 
  -s or --minout    minimum size of sequences to be output, in bases
                    (default: 10)
  -t or --threads   Number of parallel processes to run. Careful: This script
                    does not perform any verification of the number of
                    processers available. If you set this number higher than
                    the number of processors you have, performance is likely to
                    be significantly degraded.
                    (default: 4)
  -v or --version   print version information and quit

** Nucmer Options **
Advanced use only. Little reason to change defaults in most situations. See
MUMmer documentation for more information.
  --breaklen        Integer (default: 200)
  --mincluster      Integer (default: 65)
  --diagdiff        Integer (default: 5)
  --diagfactor      Float (default: 0.12)
  --minmatch        Integer (default: 20)
  --nosimplify      (default: simplify)

";

## command line processing.
use Getopt::Long;
Getopt::Long::Configure qw(gnu_getopt);

#set defaults
my $aval_pct        = 100;
my $nuc_pct         = 85;
my $min_out         = 10;
my $maxdist         = 10;
my $threads         = 4;
my $pref            = "output";
my $fof;
my $refs;
my $nucpath;
my $out_pan;
my $vers;
my $lic;
my $breaklen;
my $mincluster;
my $diagdiff;
my $diagfactor;
my $maxgap;
my $minmatch;
my $nosimplify;
my $web;

GetOptions(
    'file|f=s'      => \$fof,
    'pctcore|a=f'   => \$aval_pct,
    'refs|r=s'      => \$refs,
    'pctid|p=f'     => \$nuc_pct,
    'minout|s=i'    => \$min_out,
    'maxdist|g=i'   => \$maxdist,
    'threads|t=i'   => \$threads,
    'prefix|o=s'    => \$pref,
    'nucpath|m=s'   => \$nucpath,
    'pangenome'     => \$out_pan,
    'version|v'     => \$vers,
    'license|l'     => \$lic,
    'breaklen=i'    => \$breaklen,
    'mincluster=i'  => \$mincluster,
    'diagdiff=i'    => \$diagdiff,
    'diagfactor=f'  => \$diagfactor,
    'maxgap=i'      => \$maxgap,
    'minmatch=i'    => \$minmatch,
    'nosimplify'    => \$nosimplify,
    'web'           => \$web
) or die "$usage";
die "$usage" unless $fof;
die "version $version\n" if $vers;
die "$license\n" if $lic;

die ("ERROR: Input parameters must be non-negative numbers\n") if ($aval_pct =~ m/[^0123456789.]/);
die ("ERROR: Reference sequences must be entered as integers separated by commas\n") if ($refs and $refs =~m/[^0123456789,]/);
die ("ERROR: Input parameters must be non-negative numbers\n") if ($nuc_pct =~ m/[^0123456789.]/);
die ("ERROR: Input parameters must be non-negative numbers\n") if ($min_out =~ m/[^0123456789.]/);
die ("ERROR: Input parameters must be non-negative numbers\n") if ($maxdist =~ m/[^0123456789.]/);
die ("ERROR: Core genome percentage must be between 0 and 100\n") if ($aval_pct > 100);
die ("ERROR: Alignment percentage must be between 1 and 100\n") if ($nuc_pct > 100 or $nuc_pct < 1);
die ("ERROR: Minimum output size must be at least 1\n") if ($min_out < 1);
die ("ERROR: Maximum sequence distance must be at least 1\n") if ($maxdist < 1);

#if path to nucmer was given, check that the executables are present
my $nuc_loc = which("nucmer");
my $sc_loc = which("show-coords");
if ($nucpath){
    if (-e "$nucpath/nucmer"){
        $nuc_loc = "$nucpath/nucmer";
    } else {
        print STDERR "WARNING: Could not find nucmer at $nucpath. Searching PATH...\n";
        print STDERR "<br>\n" if $web;
    }
    if (-e "$nucpath/show-coords"){
        $sc_loc = "$nucpath/show-coords";
    }
    unless (-e "$nucpath/show-coords"){
        print STDERR "WARNING: Could not find show-coords at $nucpath. Searching PATH...\n";
        print STDERR "<br>\n" if $web;
    }
}
die "ERROR: Could not find nucmer in PATH. Make sure MUMmer is installed and executable.\n" if !$nuc_loc;
die "ERROR: Could not find show-coords in PATH. Make sure MUMmer is installed and executable.\n" if !$sc_loc;
print STDERR "nucmer found: $nuc_loc\n" unless $web;
print STDERR "show-coords found: $sc_loc\n" unless $web;

#check that nucmer-multi and nucmer_backbone are both present and accessible
unless ($web){ #skip most of this check if running the web version
    die "ERROR: Perl must be installed and in your PATH.\n" unless (which("perl"));
}
my $home_dir = abs_path($0); #get the absolute path to spine.pl
$home_dir =~ s/\/[^\/]*$//; #strip off "/spine.pl"
unless ($web){
    print STDERR "home_dir = $home_dir\n";
    die "ERROR: Can't find the required file \"nucmer_multi.pl\". Make sure it is in the \"scripts\" directory with spine.pl ($home_dir/scripts) and has not been renamed.\n" unless (-e "$home_dir/scripts/nucmer_multi.pl");
    die "ERROR: Can't find the required file \"nucmer_backbone.pl\". Make sure it is in the \"scripts\" directory with spine.pl ($home_dir/scripts) and has not been renamed.\n" unless (-e "$home_dir/scripts/nucmer_backbone.pl");
    my $nbb_vers = `perl $home_dir/scripts/nucmer_backbone.pl -V`;
    chomp $nbb_vers;
    my $min_nbb_vers = 0.3;
    my $wrong_vers = "unknown";
    if ($nbb_vers){
        $nbb_vers =~ m/^(\d+\.\d+)/;
        if ($1){
            if ($1 >= $min_nbb_vers){
                $wrong_vers = "";
            } else {
                $wrong_vers = "$1";
            }
        } else {
            $wrong_vers = $nbb_vers;
        }
    } 
    die "ERROR: Minimum version of scripts/nucmer_backbone.pl is $min_nbb_vers (detected version is $wrong_vers)\n" if $wrong_vers;
    print STDERR "nucmer_backbone.pl (version $nbb_vers) found: $home_dir/scripts/nucmer_backbone.pl\n";
    my $nm_vers = `perl $home_dir/scripts/nucmer_multi.pl -V`;
    chomp ($nm_vers);
    my $min_nm_vers = 0.3;
    $wrong_vers = "unknown";
    if ($nm_vers){
        $nm_vers =~ m/^(\d+\.\d+)/;
        if ($1){
            if ($1 >= $min_nm_vers){
                $wrong_vers = "";
            } else {
                $wrong_vers = "$1";
            }
        } else {
            $wrong_vers = $nm_vers;
        }
    } 
    die "ERROR: Minimum version of scripts/nucmer_multi.pl is $min_nm_vers (detected version is $wrong_vers)\n" if $wrong_vers;
    print STDERR "nucmer_multi.pl (version $nm_vers) found: $home_dir/scripts/nucmer_multi.pl\n";
}
    
#read in file of files
die "ERROR: Can't find file '$fof'. Please check path.\n" unless -e $fof;
die "ERROR: File '$fof' appears to be a binary file. Please check\n" if -B $fof;
open (my $fof_in, "<", $fof) or die "ERROR: Can't open file $fof: $!\n";
my @files;
my $nog = 0;
while (my $fline = <$fof_in>){
    $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
    my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
    while (@lines){
        my $line = shift @lines;
        unless ($line =~ m/^\s*$/){ #skip blank lines
            push @files, $line;
            $nog++;
        }
    }
}
close ($fof_in);
die "ERROR: At least 2 genome sequences must be given\n" if $nog < 2;

#rearrange file order such that the reference genome(s) are listed first
my $out_opt_r = 1;
my $out_opt_g = join(",", 1..$nog);
my @ref_order = 1 .. $nog;
if ($refs){
    #$out_opt_r = $refs;
    #if ($refs =~ m/^(\d+),/){
    #    $out_opt_r = $1;
    #}
    $refs =~ s/,+\s*$//; #remove any trailing commas, if they exist for some reason
    @ref_order = split(",", $refs);
    my %skip;
    foreach my $val (@ref_order){
        if ($val < 1 or $val > $nog or $val=~/\D/){
            die "ERROR: Please enter a valid number or numbers for option -r (only integers between 1 and $nog separated by commas).\n";
        }
        die "ERROR: Cannot use the same value more than once in option -r\n" if $skip{$val};
        $skip{$val}++;
    }
    for my $i (1 .. $nog){ #make sure the list of reference genomes includes all genomes if they weren't given
        next if $skip{$i};
        #$refs .= ",$i";
        push @ref_order, $i;
    }
    #$out_opt_g = $refs;
    
    ## sort the file list by the reference order
    my @idx1 = sort{$ref_order[$a] <=> $ref_order[$b]} 0 .. $#ref_order; 
    my @tmp = 1 .. $nog;
    @tmp = @tmp[@idx1]; #first need to sort a dummy list of input orders by the sorted indices of the reference order
    my @idx2 = sort{$tmp[$a] <=> $tmp[$b]} 0 .. $#tmp;
    @files = @files[@idx2]; #now sort the file list by the indices of the sorted dummy list
}

#my $no_genes_out = 1;
my $total_seqs = 0;
my %file_dup_check;
open (my $seqout, ">tmp_sequences.fasta") or die "ERROR: Can't open temporary file: $!\n";
open (my $crdout, ">tmp_coordinates.txt") or die "ERROR: Can't open temporary file: $!\n";
for my $i (0 .. $#files){
    my ($file, $fileid, $filetype) = split("\t", $files[$i]);
    die "ERROR: Can't find file '$file'. Please check path.\n" unless -e $file;
    die "ERROR: File '$file' appears to be a binary file. Please check.\n" if -B $file;
    #$nog++;
    $fileid =~ s/\s*$//; #removes any trailing space
    $fileid =~ s/^\s*//; #removes any leading space
    if (!$fileid){
        die "ERROR: All files must be given an ID\n";
    }
    $fileid =~ s/\s+/_/g; #changes any remaining spaces to underscores
    $fileid =~ s/[\/\\*\s]/_/g; #changes any unusual characters to underscores
    my $input_order = $ref_order[$i];
    print STDERR "ref order ". ($i+1) .", input order $input_order: $fileid ";
    die "ERROR: This genome ID is a duplicate of another ID in this same dataset (see above). Only unique genome IDs should be used.\n" if ($file_dup_check{$fileid});
    $file_dup_check{$fileid}++;
    if (!$filetype){
        $filetype = "fasta" if ($file =~ m/\.f[^.]*$/);
        $filetype = "gbk" if ($file =~ m/\.gb[^.]*$/);
        die "ERROR: Can't guess at file type for $fileid.\n" if !$filetype;
    }
    ##check whether the file matches the type given
    if ($filetype eq "fasta"){
        #$no_genes_out = 1;
        open (my $in, "<", $file) or die "ERROR: Can't open $file: $!\n";
        my $rec_count = 0;
        my @seqarray;
        while (my $fline = <$in>){
            $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
            my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
            while (@lines){
                my $line = shift @lines;
                $line =~ s/\s.*$//; #remove everything after the first space in the line (nucmer is going to do this to the headers anyway)
                push @seqarray, ("$line\n");
                if ($line =~ m/^>/){
                    $rec_count++;
                    next;
                }
                if ($line =~ m/[^ACTGNactgn]/){
                    #die "ERROR: File contains non-nucleotide letters. Please check file.\n";  #should I do this? Need to check how nucmer deals with ambiguous bases
                }
            }
        }
        close ($in);
        $total_seqs += $rec_count;
        if ($rec_count == 0){
            die "ERROR: File contains no records. Please check file.\n";
        } else {
            print STDERR "contains $rec_count record(s).\n";
            print STDERR "<br>\n" if $web;
            while (@seqarray){
                my $line = shift @seqarray;
                $line =~ s/>/>#$fileid#/;
                print $seqout "$line";
            }
        }
        #if ($rec_count == 1){
        #    print STDERR "contains 1 record.\n";
        #    while (@seqarray){
        #        my $line = shift @seqarray;
        #        $line =~ s/>.*/>$fileid/;
        #        print $seqout "$line";
        #    }
        #}
        #if ($rec_count > 1){
        #    print STDERR "contains $rec_count records.\n";
        #    while (@seqarray){
        #        my $line = shift @seqarray;
        #        $line =~ s/>/>|$fileid|/;
        #        print $seqout "$line";
        #    }
        #}
    } else {
        #$no_genes_out = "";
        my $g_status = gbk_convert($file, $fileid, $i+1);
        #nice_die ("File is not in Genbank format, CDS records do not have \"locus_tag\" tags, file does not contain DNA sequence, or sequence has non-nucleotide letters. Please check file.") if ($g_status != 0);
        die "ERROR: File does not contain DNA sequence. Please check file.\n" if ($g_status == 2);
        die "ERROR: CDS records missing \"locus_tag\" tags. Please check file and visit http://vfsmspineagent.fsm.northwestern.edu/gbk_reformat.cgi for conversion tool.\n" if ($g_status == 3);
        die "ERROR: DNA sequence has non-nucleotide letters. Please check file.\n" if ($g_status == 4);
        die "ERROR: File is not in Genbank format. Please check file.\n" if ($g_status != 0);
    }
}
close $seqout;
close $crdout;
%file_dup_check = ();

#print STDERR "\nFYI: Only fasta sequence files were given. No gene information will be output.\n" if ($no_genes_out);

#align sequences with nucmer
print STDERR "\nRunning nucmer with $threads processes...\n";
print STDERR "<br>\n" if $web;
my $a_params = "--maxmatch";
$a_params .= " -b $breaklen" if $breaklen;
$a_params .= " -c $mincluster" if $mincluster;
$a_params .= " -D $diagdiff" if $diagdiff;
$a_params .= " -d $diagfactor" if $diagfactor;
$a_params .= " -g $maxgap" if $maxgap;
$a_params .= " -l $minmatch" if $minmatch;
$a_params .= " --nosimplify" if $nosimplify;
print STDERR "\tnucmer options: $a_params\n";
print STDERR "<br>\n" if $web;
my $return;
{
    #local @ARGV = ("-ftmp_sequences.fasta", "-g", "-a\"--maxmatch\"", "-t$threads", "-o$pref", "-n$nuc_loc");
    our ($opt_f, $opt_g, $opt_a, $opt_t, $opt_o, $opt_n, $opt_w);
    local $opt_f = "tmp_sequences.fasta";
    local $opt_g = 1;
    local $opt_a = $a_params;
    local $opt_t = $threads;
    local $opt_o = $pref;
    local $opt_n = $nuc_loc;
    local $opt_w = 1 if $web;
    $return = do "$home_dir/scripts/nucmer_multi.pl";
}
unless ($return){
    die "ERROR: Couldn't run nucmer_multi.pl: $@\n" if $@;
    die "ERROR: Couldn't run nucmer_multi.pl: $!\n" unless defined $return;
}

#run show-coords on alignment delta file
print STDERR "\nRunning show-coords...\n";
print STDERR "<br>\n" if $web;
my @result = `$sc_loc -rTH $pref.delta > $pref.coords.txt 2>&1`;
my $error = $?;
die "ERROR: Show-coords failure (", join(",",@result), ")\n" if $error;

#run nucmer_backbone
my $rev_aval = roundup(($aval_pct / 100) * $nog);
my $aval = $nog - $rev_aval;
print STDERR "<br>\n" if $web;
print STDERR "\nGenerating core genome with core definition of at least $rev_aval of $nog genomes (>= $aval_pct%)\n";
print STDERR "<br>\n" if $web;
print STDERR "This can take a few minutes to a few hours depending on the sizes and number of the genomes.\n";
print STDERR "<br>\n" if $web;
$return = "";

print STDERR "out_pan = $out_pan\n" if $out_pan and !$web;

{
    local @ARGV = ("tmp_sequences.fasta");
    our ($opt_c, $opt_h, $opt_B, $opt_I, $opt_e, $opt_x, $opt_a, $opt_m, $opt_s, $opt_t, $opt_o, $opt_n, $opt_r, $opt_g, $opt_w);
    local $opt_c = "$pref.coords.txt";
    local $opt_a = $aval;
    local $opt_m = $maxdist;
    local $opt_h = $nuc_pct;
    local $opt_B = $min_out;
    local $opt_I = $min_out;
    local $opt_s = $pref;
    local $opt_t = $threads;
    local $opt_o = 1;
    local $opt_e = 1;
    local $opt_n = 1 if $out_pan;
    local $opt_x = "tmp_coordinates.txt";
    local $opt_r = $out_opt_r;
    local $opt_g = $out_opt_g;
    local $opt_w = 1 if $web;
    $return = do "$home_dir/scripts/nucmer_backbone.pl";
}
unless ($return){
    die "ERROR: Couldn't run nucmer_backbone.pl: $@\n" if $@;
    die "ERROR: Couldn't run nucmer_backbone.pl: $!\n" unless defined $return;
}
#unlink ("tmp_sequences.fasta");
#unlink ("tmp_coordinates.txt");
print STDERR "\nFinished!\n" unless $web;

#------------------------
sub gbk_convert{
    my $file = shift;
    my $filename = shift;
    my $filenum = shift;
    #$seqout, $crdout
    open (my $gbkin, "<", $file) or die "ERROR: Can't open $file: $!\n";
    my (@list);
    my @seqlist;
    my $loccount = 0;
    my $seqcount = 0;
    my ($c_id, $c_seq);
    #my ($start, $stop);
    my $is_origin;
    my $is_feature;
    my ($is_cds, $is_prod);
    my %tags;
    while (my $fline = <$gbkin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            if ($line =~ m/^LOCUS\s+([^\s]+)/){
                $is_origin = "";
                if ($c_id){
                    return (2) if !$c_seq;  #error if there is no sequence in the genbank file
                    push @seqlist, ([$c_id, $c_seq]);
                    $c_seq = "";
                    #($start, $stop) = ("") x 2;
                }
                $seqcount++;
                if ($line =~ m/^LOCUS\s+(\S+)\s+\d+ bp/){
                    $c_id = $1;
                } else {
                    $c_id = "rec$seqcount";
                }
                next;
            }
            if ($line =~ m/^FEATURES\s+Location\/Qualifiers/){
                $is_feature = 1;
                next;
            }
            next unless $is_feature;
            
            if ($line =~ m/^\s+(\S+)\s+(complement\()*<*(\d+)\.\.(\d+)>*\)*\s*$/){
                my ($type, $start, $stop) = ($1, $3, $4);
                my $dir = "+";
                $dir = "-" if $2;
                if (%tags){
                    return(3) unless $tags{'locus_tag'}; #no locus_tag was present on the last record
                    my ($o_id, $o_start, $o_stop, $o_dir) = ($tags{'locus_tag'}, $tags{'start'}, $tags{'stop'}, $tags{'dir'});
                    my $o_prod = $tags{'product'} if $tags{'product'};
                    print $crdout "$filenum\t$o_id\t#$filename#$c_id\t$o_start\t$o_stop\t$o_dir";
                    print $crdout "\t$o_prod" if $o_prod;
                    print $crdout "\n";
                    $loccount++;
                }
                undef %tags;
                $is_prod = "";
                $is_cds = "";
                if ($type eq "CDS"){
                    ($start, $stop) = ($stop, $start) if $start > $stop;
                    $tags{'start'} = $start;
                    $tags{'stop'} = $stop;
                    $tags{'dir'} = $dir;
                    $is_cds = 1;
                }
                next;
            }
            if ($is_cds and $line !~ m/^ORIGIN/){
                if ($line =~ m/^\s+\/(\S+)=\"*([^"]*)\"*/){
                    my ($key, $val) = ($1, $2);
                    if ($key eq "locus_tag" or $key eq "product"){
                        $val =~ s/\s*$//;
                        $tags{$key} = $val;
                        $is_prod = 1 if $key eq "product";
                        next;
                    } else {
                        $is_prod = "";
                    }
                } else {
                    if ($is_prod){
                        $line =~ s/^\s*//;
                        $line =~ s/\"//g;
                        $tags{'product'} .= " $line";
                    }
                }
                next;
            }
            if ($line =~ m/^ORIGIN/){
                if (%tags){
                    return(3) unless $tags{'locus_tag'}; #no locus_tag was present on the last record
                    my ($o_id, $o_start, $o_stop, $o_dir) = ($tags{'locus_tag'}, $tags{'start'}, $tags{'stop'}, $tags{'dir'});
                    my $o_prod = $tags{'product'} if $tags{'product'};
                    print $crdout "$filenum\t$o_id\t#$filename#$c_id\t$o_start\t$o_stop\t$o_dir";
                    print $crdout "\t$o_prod" if $o_prod;
                    print $crdout "\n";
                    $loccount++;
                }
                undef %tags;
                $is_prod = "";
                $is_cds = "";
                #if ($start){ #no locus_tag was present on the last record
                #    return (3);
                #}
                $is_origin = 1;
                next;
            }
            if ($is_origin and $line =~ m/^\s*\d+\s(.*)/){
                my $seqline = $1;
                $seqline =~ s/\s//g;
                if ($seqline =~ m/[^ACTGNactgn]/){
                    ###since sequences may contain some IUPAC ambiguity codes other than "N" I'll ignore this for now
                    #return (4);
                }
                $c_seq .= $seqline;
                next;
            }
        }
    }
    #$no_genes_out = 1 if $loccount == 0;
    if ($c_id){
        return (2) if !$c_seq;
        push @seqlist, ([$c_id, $c_seq]);
        $c_seq = "";
        #($start, $stop) = ("") x 2;
    }
    close ($gbkin);
    if ($seqcount == 0){
        return (2);
    } else {
        foreach my $slice (@seqlist){
            my ($a_id, $a_seq) = @{$slice};
            $a_id =~ s/\s.*$//; #remove everything after the first space in the sequence header (nucmer is going to do this to the headers anyway)
            print $seqout ">#$filename#$a_id\n$a_seq\n";
        }
        print STDERR "contains $seqcount record(s) with $loccount CDS.\n";
        print STDERR "<br>\n" if $web;
        return (0);
    }
    return (5);
}

sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}
