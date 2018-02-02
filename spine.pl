#!/usr/bin/perl

my $version = "0.3.1";
## Changes from v0.3 -> v0.3.1
# Can now accept fasta+gff3 annotations (i.e. Ensembl format)
# Improved gbk_convert subroutine to allow records where the locus id may be in the gene record, but not the CDS record

## Changes from v0.2.4 -> v0.3
# Can now accept multi-file genome sets (i.e. multiple chromosomes or plasmids)
# Requires version 0.4 of nucmer_backbone.pl
# Improves speed of --mini option by only outputting backone files and no core or accessory files

##Changes from v0.2.3 -> 0.2.4
# Removed File::Which dependency. Added subroutine to test for whether executable is in PATH that uses only core Perl modules
# Updated nucmer_backbone.pl to v0.3.3 to make it more memory-friendly to large data sets
# Added --mini option to allow spine to produce limited output, i.e. just backbone core genome sequence and coordinates.

##Changes from v0.2.2 -> 0.2.3
# Fixed bug in genbank file parsing where some genes that span the end of a contig might not appear in results

##Changes from v0.2.1 -> 0.2.2
# Fixed bug in nucmer_backbone.pl causing error messages when all-N regions were encountered

##Changes from v0.2 -> 0.2.1
# Fixed bug in nucmer_backbone.pl where first CDS on each contig was not being output

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
    Spine
    Copyright (C) 2016-2018 Egon A. Ozer

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
use File::Basename;
use File::Spec::Functions qw ( catfile path );

$|++;

my $usage = "
Spine (version $version);

This is a wrapper script to run the Spine algorithm

PREREQUISITES:
- Perl 5.10 or above
- MUMmer version 3.22 or above
- Mac OSX or Linux. No guarantees that this will work on Windows or other
  operating systems.

REQUIRED:
  -f or --file      file with list of input sequence files. Accepted file
                    formats include fasta sequence files, genbank sequence +
                    annotation files, or separate fasta sequence files with
                    corresponding gff3-formatted annotation files.
                    This file should be formatted like so:
                    
                    path/to/file1<tab>unique_identifier<tab>fasta or gbk or comb
                    path/to/file2<tab>unique_identifier<tab>fasta or gbk or comb
                    
                    Example:
                    /home/seqs/PAO1.fasta   PAO1    fasta
                    /home/seqs/LESB58.gbk   LESB58  gbk
                    /home/seqs/PA14.fasta,/home/seqs/PA14.gff3  PA14    comb
                    
                    The third column (fasta, gbk, or comb) is optional, but should
                    be given if your sequence files end with suffixes other
                    than \".fasta\" or \".gbk\" or if you are providing sequences
                    with gff3 annotation files, i.e. comb(ined).
                    
                    If you have genomes spread across multiple files (i.e.
                    chromosomes and/or plasmids), these can be combined by
                    either concatenating the files into one:
                    i.e. 'cat chrom_I.gbk chrom_II.gbk > combined.gbk'
                    or by including all the files in this input file,
                    separated by commmas:
                    Example:
                    /seqs/chrom_I.fasta,/seqs/chrom_II.fasta    mygenome    fasta
                    chrom_A.gbk,chrom_B.gbk,plasmid_X.gbk   myothergenome   gbk
                    seqA.fasta,seqB.fasta,seqA.gff3,seqB.gff3   genomeAB    comb
                    
                    IMPORTANT: When including multiple files for a strain or
                    joining multiple files within a strain, please ensure that
                    all chromosome/plasmid/contig IDs are unique across files
                    within a single genome. If sequence IDs are duplicated, the
                    results are likely to be wrong.
                    
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
  --mini            produce only limited output, i.e. just the backbone sequence
                    derived from the reference genome(s). This saves time on
                    large data sets, especially if you only need the backbone
                    sequence to get accessory sequences from AGEnt.
                    (default: core and accessory sequence sets will be output
                    for all included genomes)
  --pangenome       produce a pangenome sequence and characteristics from
                    sequences in the order given. This option will be ignored
                    if '--mini' option is given.
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
my $mini;
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
    'mini'          => \$mini,
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
die "version $version\n" if $vers;
die "$license\n" if $lic;
die "$usage" unless $fof;

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
my $nuc_loc = is_path("nucmer");
my $sc_loc = is_path("show-coords");
if ($nucpath){
    if (-x "$nucpath/nucmer"){
        $nuc_loc = "$nucpath/nucmer";
    } else {
        print STDERR "WARNING: Could not find nucmer at $nucpath. Searching PATH...\n";
        print STDERR "<br>\n" if $web;
    }
    if (-x "$nucpath/show-coords"){
        $sc_loc = "$nucpath/show-coords";
    } else {
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
    die "ERROR: Perl must be installed and in your PATH.\n" unless (is_path("perl"));
}
my $home_dir = abs_path($0); #get the absolute path to spine.pl
$home_dir =~ s/\/[^\/]*$//; #strip off "/spine.pl"
unless ($web){
    print STDERR "home_dir = $home_dir\n";
    die "ERROR: Can't find the required file \"nucmer_multi.pl\". Make sure it is in the \"scripts\" directory with spine.pl ($home_dir/scripts) and has not been renamed.\n" unless (-e "$home_dir/scripts/nucmer_multi.pl");
    die "ERROR: Can't find the required file \"nucmer_backbone.pl\". Make sure it is in the \"scripts\" directory with spine.pl ($home_dir/scripts) and has not been renamed.\n" unless (-e "$home_dir/scripts/nucmer_backbone.pl");
    my $nbb_vers = `perl $home_dir/scripts/nucmer_backbone.pl -V`;
    chomp $nbb_vers;
    my $min_nbb_vers = 0.4;
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
    my ($filestr, $fileid, $filetype) = split("\t", $files[$i]);
    $filestr =~ s/\s*$//;
    my @filelist = split(",", $filestr);
    #$nog++;
    $fileid =~ s/\s*$//; #removes any trailing space
    $fileid =~ s/^\s*//; #removes any leading space
    clean_exit ("ERROR: All files must be given an ID") if (!$fileid);
    $fileid =~ s/\s+/_/g; #changes any remaining spaces to underscores
    $fileid =~ s/[\/\\*\s]/_/g; #changes any unusual characters to underscores
    my $input_order = $ref_order[$i];
    print STDERR "ref order ". ($i+1) .", input order $input_order: $fileid";
    if (!$filetype){
        $filetype = "fasta" if ($filestr =~ m/\.f[^.]*$/i);
        $filetype = "gbk" if ($filestr =~ m/\.gb[^.]*$/i);
        $filetype = "comb" if ($filestr =~ m/\.f[^.]*(?:,|\Z)/i and $filestr =~ m/\.gf[^.]*(?:,|\Z)/i);
        clean_exit ("\nERROR: Can't guess at file type for $fileid. Please indicate 'fasta', 'gbk', or 'comb' on input form.") if !$filetype;
    }
    my $outtype = $filetype;
    $outtype = "fasta+gff3" if $filetype eq "comb";
    print STDERR " filetype: $outtype\n";
    print STDERR "<br>\n" if $web;
    clean_exit ("ERROR: This genome ID is a duplicate of another ID in this same dataset (see above). Only unique genome IDs should be used.") if ($file_dup_check{$fileid});
    $file_dup_check{$fileid}++;
    ##check whether the file matches the type given
    if ($filetype eq "comb"){
        my @ffiles;
        foreach my $file (@filelist){
            if ($file =~ m/\.f[^.]*\s*$/i){
                push @ffiles, $file;
            } elsif ($file =~ m/\.gf[^.]*\s*$/i){
                my $g_status = gff_convert($file, $fileid, $i+1);
                clean_exit ("ERROR: File '$file' cannot be opened. Please check file.") if ($g_status == 1);
                clean_exit ("ERROR: File '$file' appears to be binary. Please check file.") if ($g_status == 2);
                clean_exit ("ERROR: File '$file' is not in gff3 format. Please check file.") if ($g_status != 0);
            } else {
                #make an attempt to figure out if the file is a fasta file or a gff file
                clean_exit ("ERROR: Can't find file '$file'. Please check path.") unless -e $file;
                clean_exit ("ERROR: File '$file' appears to be a binary file. Please check.") if -B $file;
                open (my $test, "<$file") or die "ERROR: Can't open $file: $!\n";
                my $type;
                while (my $line = <$test>){
                    chomp $line;
                    next if $line =~ m/^\s*$/; #skip blank lines
                    if ($line =~ m/^##*gff-version/){
                        $type = "g";
                    } elsif ($line =~ m/^>/){
                        $type = "f";
                    }
                    last;
                }
                close ($test);
                unless ($type){
                    clean_exit ("ERROR: Can't determine whether file '$file' is fasta or gff3. Please verify the file format.");
                }
                if ($type eq "f"){
                    push @ffiles, $file;
                } else {
                    my $g_status = gff_convert($file, $fileid, $i+1);
                    clean_exit ("ERROR: File '$file' cannot be opened. Please check file.") if ($g_status == 1);
                    clean_exit ("ERROR: File '$file' appears to be binary. Please check file.") if ($g_status == 2);
                    clean_exit ("ERROR: File '$file' is not in gff3 format. Please check file.") if ($g_status != 0);
                }
            }
        }
        if (@ffiles){
            $filetype = "fasta";
            @filelist = @ffiles;
        } else {
            clean_exit ("ERROR: No fasta-formatted sequence files were identified for $fileid");
        }
    }
    if ($filetype eq "fasta"){
        #$no_genes_out = 1;
        foreach my $file (@filelist){
            clean_exit ("ERROR: Can't find file '$file'. Please check path.") unless -e $file;
            clean_exit ("ERROR: File '$file' appears to be a binary file. Please check.") if -B $file;
            open (my $in, "<", $file) or clean_exit ("ERROR: Can't open $file: $!");
            my $shortfile = basename($file);
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
                clean_exit ("ERROR: File $shortfile contains no records. Please check file.\n");
            } else {
                print STDERR "\t$shortfile contains $rec_count sequence record(s).\n";
                print STDERR "<br>\n" if $web;
                while (@seqarray){
                    my $line = shift @seqarray;
                    $line =~ s/>/>#$fileid#/;
                    print $seqout "$line";
                }
            }
        }
    } elsif ($filetype eq "gbk") {
        #$no_genes_out = "";
        foreach my $file (@filelist){
            my $g_status = gbk_convert($file, $fileid, $i+1);
            #nice_die ("File is not in Genbank format, CDS records do not have \"locus_tag\" tags, file does not contain DNA sequence, or sequence has non-nucleotide letters. Please check file.") if ($g_status != 0);
            clean_exit ("ERROR: File '$file' cannot be opened.") if ($g_status == 1);
            clean_exit ("ERROR: File '$file' does not contain DNA sequence. Please check file.") if ($g_status == 2);
            clean_exit ("ERROR: CDS records missing \"locus_tag\" tags. Please check file '$file' and visit http://vfsmspineagent.fsm.northwestern.edu/gbk_reformat.cgi for conversion tool.") if ($g_status == 3);
            clean_exit ("ERROR: DNA sequence has non-nucleotide letters. Please check file.") if ($g_status == 4);
            clean_exit ("ERROR: File '$file' may be binary instead of text. Please check file.") if ($g_status == 5);
            clean_exit ("ERROR: File '$file' is not in Genbank format. Please check file.") if ($g_status != 0);
        }
    }
}
close $seqout;
close $crdout;
%file_dup_check = ();

#print STDERR "\nFYI: Only fasta sequence files were given. No gene information will be output.\n" if ($no_genes_out);

#align sequences with nucmer
print STDERR "<h1>Running Nucmer</h1>\n" if $web;
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
print STDERR "<br>\n<h1>Generating Core Genome</h1><strong>" if $web;
print STDERR "\nGenerating core genome with core definition of at least $rev_aval of $nog genomes (>= $aval_pct%)\n";
print STDERR "</strong><br>\n" if $web;
print STDERR "This can take a few minutes to a few hours depending on the sizes and number of the genomes.\n";
print STDERR "<br>\n" if $web;
$return = "";

$out_pan = "" if $mini;
print STDERR "out_pan = $out_pan\n" if $out_pan and !$web;

{
    local @ARGV = ("tmp_sequences.fasta");
    our ($opt_c, $opt_h, $opt_B, $opt_I, $opt_e, $opt_x, $opt_a, $opt_m, $opt_s, $opt_t, $opt_o, $opt_n, $opt_r, $opt_g, $opt_w, $opt_z);
    local $opt_c = "$pref.coords.txt";
    local $opt_a = $aval;
    local $opt_m = $maxdist;
    local $opt_h = $nuc_pct;
    local $opt_B = $min_out;
    local $opt_I = $min_out;
    local $opt_s = $pref;
    local $opt_t = $threads;
    local $opt_o = 1 unless $mini;
    local $opt_e = 1 unless $mini;
    local $opt_n = 1 if $out_pan;
    local $opt_x = "tmp_coordinates.txt";
    local $opt_r = $out_opt_r;
    local $opt_g = $out_opt_g;
    local $opt_w = 1 if $web;
    local $opt_z = $version;
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
sub is_path {
    ## Subroutine based on StackOverflow post by Sinan Unur (https://stackoverflow.com/a/8243770)
    my $exe = shift;
    my @path = path;
    my @pathext = ( q{} );
    if ($^O eq 'MSWin32'){
        push @pathext, map { lc } split /;/, $ENV{PATHEXT};
    }
    for my $dir (@path){
        for my $ext (@pathext){
            my $f = catfile $dir, "$exe$ext";
            return ($f) if -x $f;
        }
    }
    return();
}

sub gbk_convert{
    my $file = shift;
    my $filename = shift;
    my $filenum = shift;
    my $return_status = 0;
    my $shortfile = basename($file);
    return(1) unless -e $file;
    return(5) if -B $file;
    open (my $gbkin, "<", $file) or return(1);
    my $loccount = 0;
    my $seqcount = 0;
    my ($c_id, $c_seq);
    my $is_prod;
    my @tags;
    my %crecs;
    my @ctg_order;
    my $reading = 1; # 1 = front material, 2 = annotations, 3 = sequence
    
    while (my $fline = <$gbkin>){
        $fline =~ s/\R/\012/g; #converts to UNIX-style line endings
        my @lines = split("\n", $fline); #need to split lines by line-ending character in the case of Mac-formatted files which only have CR line terminators, not both CR and LF like DOS
        while (@lines){
            my $line = shift @lines;
            next if $line =~ m/^\s*$/;
            if ($line =~ m/^LOCUS\s+\S*\s+\d+\sbp/){
                if ($reading == 2){ #no ORIGIN sequence record was found between LOCUS records
                    return (2);
                }
                if ($reading == 3){
                    if ($c_seq and $c_id){
                        print $seqout ">#$filename#$c_id\n$c_seq\n";
                        $c_seq = "";
                        $reading = 1;
                    } else {
                        return (2);
                    }
                }
            }
            if ($line =~ m/^\/\//){ #reached the end of the file (or record)
                if ($c_seq and $c_id){
                    print $seqout ">#$filename#$c_id\n$c_seq\n";
                    $c_seq = "";
                    $reading = 1;
                } else {
                    return (2);
                }
            }
            if ($reading == 1){
                if ($line =~ m/^LOCUS\s+([^\s]+)/){
                    $seqcount++;
                    if ($line =~ m/^LOCUS\s+(\S+)\s+\d+ bp/){
                        $c_id = $1;
                    } else {
                        $c_id = "rec$seqcount";
                    }
                    push @ctg_order, $c_id;
                    next;
                }
                if ($line =~ m/^FEATURES\s+Location\/Qualifiers/){
                    $reading = 2;
                    next;
                }
            } elsif ($reading == 2){
                if ($line =~ m/^\s+(\S+)\s+(complement\()*[<>]*(\d+)<*\.\.[<>]*(\d+)>*\)*\s*$/){
                    $is_prod = "";
                    my ($type, $start, $stop) = ($1, $3, $4);
                    my $dir = "+";
                    $dir = "-" if $2;
                    unless ($crecs{$c_id}{$start}{$stop}{$dir}){
                        @{$crecs{$c_id}{$start}{$stop}{$dir}} = (0);
                    }
                    if ($type eq "CDS"){
                        ${$crecs{$c_id}{$start}{$stop}{$dir}}[0] = 1;
                        $loccount++;
                    }
                    if (@tags){
                        my ($o_start, $o_stop, $o_dir) = @tags;
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[3] if $tags[3];
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[4] if $tags[4];
                        $loccount++;
                    }
                    @tags = ($start, $stop, $dir);
                    next;
                }
                if ($line =~ m/^ORIGIN\s*$/){
                    $is_prod = "";
                    if (@tags){
                        my ($o_start, $o_stop, $o_dir) = @tags;
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[1] = $tags[3] if $tags[3];
                        ${$crecs{$c_id}{$o_start}{$o_stop}{$o_dir}}[2] = $tags[4] if $tags[4];
                        $loccount++;
                    }
                    undef @tags;
                    $reading = 3;
                    next
                }
                if ($line =~ m/^\s+\/(\S+)=\"*([^"]*)\"*/){
                    $is_prod = "";
                    my ($key, $val) = ($1, $2);
                    if ($key eq "locus_tag"){
                        $tags[3] = $val;
                    }
                    if ($key eq "product"){
                        $tags[4] = $val;
                        $is_prod = 1;
                    }
                    next;
                }
                if ($is_prod){
                    $line =~ s/^\s*//;
                    $line =~ s/"*\s*$//;
                    $tags[4] .= " $line";
                }
            } elsif ($reading == 3){
                $line =~ s/\d//g;
                $line =~ s/\s//g;
                $c_seq .= $line;
                next;
            }
        }
    }
    if ($c_seq and $c_id){
        print $seqout ">#$filename#$c_id\n$c_seq\n";
        $c_seq = "";
        $reading = 1;
    }
    close ($gbkin);
    my $cds_count;
    foreach my $cid (@ctg_order){
        foreach my $start (sort{$a <=> $b} keys %{$crecs{$cid}}){
            foreach my $stop (sort{$a <=> $b} keys %{$crecs{$cid}{$start}}){
                foreach my $dir (sort keys %{$crecs{$cid}{$start}{$stop}}){
                    my ($is_cds, $lid, $prod) = @{$crecs{$cid}{$start}{$stop}{$dir}};
                    if ($is_cds){
                        $cds_count++;
                        unless ($lid){
                            print STDERR "ERROR: CDS at $start..$stop on contig $cid in file $file for strain $filename has no locus_id\n";
                            print STDERR "<br>\n" if $web;
                            return(3);
                        }
                        print $crdout "$filenum\t$lid\t#$filename#$cid\t$start\t$stop\t$dir\t";
                        print $crdout "$prod" if $prod;
                        print $crdout "\n";
                    }
                }
            }
        }
    }
    unless ($cds_count){
        print STDERR "\n<p>" if $web;
        print STDERR "FYI: No CDS annotations were found in genbank file $file. Only sequence information will be used from this file.\n";
        print STDERR "</p>\n" if $web;
    }
    print STDERR "\t$shortfile contains $seqcount sequence record(s) and $cds_count CDS annotations.\n";
    print STDERR "<br>\n" if $web;
    return (0);
}


sub old_gbk_convert{
    my $filestr = shift;
    my $filename = shift;
    my $filenum = shift;
    my @filelist = split(",", $filestr);
    my $return_status = 0;
    foreach my $file (@filelist){
        #$seqout, $crdout
        open (my $gbkin, "<", $file) or die "ERROR: Can't open $file: $!\n";
        my $shortfile = basename($file);
        my @seqlist;
        my $loccount = 0;
        my $seqcount = 0;
        my ($c_id, $c_seq);
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
                
                if ($line =~ m/^\s+(\S+)\s+(complement\()*<*(\d+)<*\.\.>*(\d+)>*\)*\s*$/){
                    my ($type, $start, $stop) = ($1, $3, $4);
                    my $dir = "+";
                    $dir = "-" if $2;
                    if (%tags){
                        return(3) unless $tags{'locus_tag'}; #no locus_tag was present on the last record
                        my ($o_id, $o_start, $o_stop, $o_dir) = ($tags{'locus_tag'}, $tags{'start'}, $tags{'stop'}, $tags{'dir'});
                        my $o_prod = $tags{'product'} if $tags{'product'};
                        print $crdout "$filenum\t$o_id\t#$filename#$c_id\t$o_start\t$o_stop\t$o_dir\t";
                        print $crdout "$o_prod" if $o_prod;
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
            print STDERR "\t$shortfile contains $seqcount record(s) with $loccount CDS.\n";
            print STDERR "<br>\n" if $web;
            #return (0);
        }
    }
    return(0);
}

sub gff_convert {
    my $file = shift;
    my $filename = shift;
    my $filenum = shift;
    my $shortfile = basename($file);
    my $return_status = 0;
    ## gff is tough because there seems to be very little standardization of tags for locus IDs and gene products
    ## Need to try to set some priorities.
    return(1) unless -e $file;
    return(2) if -B $file;
    open (my $in, "<$file") or return(1);
    my %ctg_order;
    my ($count, $ctg_num) = (0) x 2;
    my %crecs;
    while (my $line = <$in>){
        next if $line =~ m/^#/;
        next if $line =~ m/^\s*$/;
        my ($contig, $x1, $type, $start, $stop, $x2, $dir, $x3, $rest) = split("\t", $line);
        unless ($ctg_order{$contig}){
            $ctg_num++;
            $ctg_order{$contig} = $ctg_num;
        }
        unless ($crecs{$contig}{$start}{$stop}{$dir}){
            @{$crecs{$contig}{$start}{$stop}{$dir}} = (0); #initialize the record as non-cds
        }
        if ($type eq "gene"){ # in Ensembl gff3 files, records corresponding to locus_id and product in gbk files are found in the gene record, not the CDS record
            if ($rest =~ m/(?:\A|;)gene_id="*([^;"]+)/){
                ${$crecs{$contig}{$start}{$stop}{$dir}}[1] = $1;
            }
            if ($rest =~ m/(?:\A|;)description="*([^;"]+)/){
                ${$crecs{$contig}{$start}{$stop}{$dir}}[2] = $1;
            }
        }
        next unless $type eq "CDS";
        $count++;
        ${$crecs{$contig}{$start}{$stop}{$dir}}[0] = 1;
        #if the CDS record has "locus" or "name" records, will use these as the locus id and product names
        if ($rest =~ m/(?:\A|;)locus="*([^;"]+)/){
            ${$crecs{$contig}{$start}{$stop}{$dir}}[1] = $1;
        }
        if ($rest =~ m/(?:\A|;)name="*([^;"]+)/){
            ${$crecs{$contig}{$start}{$stop}{$dir}}[2] = $1;
        }
        #assign a locus ID if none was found
        unless (${$crecs{$contig}{$start}{$stop}{$dir}}[1]){
            my $a_count = sprintf("%05d", $count);
            $contig =~ s/^.*\|//;
            my $orfid = "$shortfile-$a_count";
            if ($rest =~ m/(?:\A|;)ID="*([^;"]+)/){
                $orfid = $1;
            }
            ${$crecs{$contig}{$start}{$stop}{$dir}}[1] = $orfid;
        }
    }
    close ($in);
    return(3) unless ($count > 0); # no CDS records identified
    foreach my $cid (sort{$ctg_order{$a} <=> $ctg_order{$b}} keys %ctg_order){
        foreach my $start (sort{$a <=> $b} keys %{$crecs{$cid}}){
            foreach my $stop (sort{$a <=> $b} keys %{$crecs{$cid}{$start}}){
                foreach my $dir (sort keys %{$crecs{$cid}{$start}{$stop}}){
                    my ($is_cds, $lid, $prod) = @{$crecs{$cid}{$start}{$stop}{$dir}};
                    if ($is_cds){
                        print $crdout "$filenum\t$lid\t#$filename#$cid\t$start\t$stop\t$dir\t";
                        print $crdout "$prod" if $prod;
                        print $crdout "\n";                        
                    }
                }
            }
        }
    }
    print STDERR "\t$shortfile contains $count CDS annotations.\n";
    print STDERR "<br>\n" if $web;
    return(0);
}

sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1));
}

sub clean_exit {
    my $message = shift;
    close ($seqout);
    unlink ("tmp_sequences.fasta");
    close ($crdout);
    unlink ("tmp_coordinates.txt");
    die "$message\n";
    return;
}
