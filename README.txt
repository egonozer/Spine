INTRODUCTION:
Nucmer_Backbone, along with alignment by the Nucmer function of
MUMmer is part of the Spine algorithm for identification of the
conserved core genome of bacteria and other small genome organisms. 
REQUIREMENTS:
Nucmer_Backbone was written and tested using Perl 5, version 12,
subversion 4. It has been tested on Mac OS X and Linux, but
theoretically should work on any system able to run Perl scripts.
Nucmer_Backbone requires output produced using the nucmer and
show-coords functions of the MUMmer software package
(http://mummer.sourceforge.net). Install MUMmer as directed by the
instructions included with the software.
INSTALLATION:
Simply copy the file nucmer_backbone.pl either to a scripts folder or
an application folder. Anywhere you can remember where it is.
USAGE:
Basic command: perl nucmer_backbone.pl -c alignment.coords.txt
genomes.fasta
For list of options, call the script without any inputs: perl
nucmer_backbone.pl
Required Inputs:

  -c	Alignment coordinates file. Produced by nucmer alignment of
all input genomes vs. all input genomes. To generate this file,
perform the following steps:**1. Concatenate all input sequence files, in fasta format, into a
single file. For example: "cat genome1.fasta genome2.fasta
genome3.fasta > all_genomes.fasta"**2. Produce alignment file using the MUMmer script nucmer with the
"--maxmatch" option to produce a .delta alignment file: "nucmer
--maxmatch -p out all_genomes.fasta all_genomes.fasta". If you are
using a large number of genomes (>10 or so) and you have several
processors to work with, consider downloading and using
nucmer_multi.pl from http://vfsmspineagent.fsm.northwestern.edu to
speed up the alignment step by using multiple processors. **3. Output coordinate file from the nucmer alignment using the MUMer
application show-coords: "show-coords -rTH out.delta >
out.coords.txt"
  fasta sequence file(s)	Can either be the path to the
multi-fasta file produced in step 1 above or paths to all of the
individual sequence files listed seqeuentially. Fasta records with
the same header prefix will be grouped together as one genome (i.e.
in the case of multi-contig draft genomes). In this case, make sure
each header starts with "|genome_name|", i.e ">contig1, >contig2,
>contig3, etc." should become ">|genome1|contig1, >|genome1|contig2,
>|genome1|contig3, etc." The script "rename_headers.pl" is included
with this script to automate this process for you.  Simply run "perl
rename_headers.pl input_file_name.fasta genome_name >
output_file_name.fasta"
Optional Inputs (need to be given before giving the path(s) to the
fasta sequence file(s)):
  -a	Number of genomes from which a section can be absent and
still be considered to be core genome (default 0)
  -r	Genome used as reference for generating backbone fragement
list with "1" being the first genome in the input fasta file(s), "2"
being the second genome in the input file(s), etc. (default 1)
  -g	List of genomes from which core sequence should be created,
separated by commas with no spaces between. Order of the list will
determine the priority in producing core genome. The genome given by
option -r will be given highest priority with all other genomes
moving down by one step. Number of genomes entered must be at least
(a+1) where "a" is the value given to option -a, os if the value
given to -a is 2, at least three genome numbers must be entered. As
with option -r above, "1" is the first genome in the input fasta
file(s), "2" is the second genome in the input file(s), etc.Example: To create a backbone sequence from the second, fourth, and
third genomes given and ignore the first, enter "-g 2,4,3". The
second genome will be given highest priority, the fourth genome
next-highest priority, and third genome next-highest priority.
(default is to prioritize genomes in the order their sequences were
provided, with the first genome having highest priority, second
genome having second-highest priority, etc).
  -l	Path to file listing paths to genbank files of genome
annotations, given in the same order as the input fasta files. File
locations should be one per line. If this file is given, a list of
locus IDs of genes present in the core genome (using IDs from the
reference genome(s)) will be output, as will lists of accessory
genome genes for each of the input sequences. ** WARNING: This option
currently only works if all of the input genomes are complete (i.e.
single scaffold or chromosome) and have a genbank file of
annotations.
  -x	File listing all CDS locus IDs for the input genomes in the
format "genome order number<tab>locus id<tab>start coord<tab>end
coord" (one entry per line). ** WARNING: This option currently only
works if all of the input genomes are complete (i.e. single scaffold
or chromosome) and have a genbank file of annotations. Example:
  1	gen1_0001	456	2176
  1	gen1_0002	3187	4599
  2	gen2_0001	679	3100
  etc.
If an input is given for both -x and -l, only the file given to -x
will be used as input. If an input is given for both -x and -l and
the file given to -x does not exist, it will be created using the
files from -l to create the locus ID file and saved under the name
given to -x.
  -p	Minimum percentage coverage of a locus for it to be output as
core (default 50)
  -m	Maximum distance between core fragments. Distances less than
this parameter will result in outupt of sequence of fragments
concatenated with N’s representing the non-backbone gaps. (default
10)
  -h	Minimum percent identity for nucmer alignments (default 85)

  -B	Minimum core genome fragment size to output, in bases
(default 10)
  -I	Minimum accessory genome fragment size to output, in bases
(default 10)
  -s	Prefix of output files (default: "output")

  -o	If given, wil output accessory genome coordinates for each
genome (default: only output accessory coordinates and statistics for
reference genome(s))
  -e	Output file of position values that can be used to calculate
pangenome and core genome characteristics of the entered genomes with
core_and_pangenome.pl. If selected, will force -o option. (default:
position values are not output)
  -t	Number of processes to run in parallel. Speeds up analysis on
multi-processor machines. This value should not exceed your number of
processors or you will likely experience massive slowdown. (default:
1)
  -v	Verbose output


OUTPUT FILES:

- statistics.txt
First line shows the current software version used.
Second line shows the input parameters given to the software.
Column descriptions:
* gen_#: Number assigned to the genome based on the order the
sequences were input to Spine.* gen_name: Name of the sequence, user-assigned
* gen_size: Total size of the sequence, in bases
* source: Shows whether the subsequent columns are refering to the
accessory genome or backbone (core) genome* total_bp: Size of the accessory or core genome component, in bases
* gc_%: Percent GC content of the accessory or core genome component
* num_segs: Number of separate sequence segements output
* min_seg: Smallest segment size, in bases
* max_seg: Largest segment size, in bases
* avg_leng: Average length of the output segments
* median_leng: Median length of the output segments
* num_cds (if annotation was provided): number of coding sequences
present found in the accessory or core genome component
- backbone_coords.txt
Core genome coordinates. Consists of three columns. 
Column 1 is the name of the source genome. If the source genome
consisted of multiple sequences, the sequence name will follow the
genome name, separated by "|". Example: if a portion of the core
genome was found in the sequence "contigA" in the genome "genomeA",
this will appear in Column 1 as "genomeA|contigA".Columns 2 and 3 are the start and stop coordinates of the core genome
sequence, respectively. Coordinates are 1-based.
- backbone.fasta
Nucleotide sequences of the core genome segments.

- backbone_loci.txt (if annotation was provided)
List of core genome coding sequences.
Column 1 is the locus ID of the coding sequence
Column 2 is the percent of the CDS length found within core genome
sequence
- accessory_coords.txt
Accessory genome coordinates for each input genome. Columns are as
described for backbone_coords.txt
- accessory_loci.txt (if annotation was provided)
List of accessory genome coding sequences. Columns are as described
for backbone_loci.txt
- position_counts.txt
This file should not be needed for routine use. Is meant to be used
as input for core_and_pangenome.pl to calculate core-, pan-, and new
genome sizes at permutations of the input genomic sequences.Column descriptions:
* ref_genome: ID number of the genome (based on the order of the
sequences input to Spine) used as a reference in alignment to the
other input genomes.* num_genomes_sharing: The number of genomes in which the sequences
were found, i.e "1" indicates that the bases indicated in column four
were only found in one genome (the "ref_genome"), "2" indicates that
the bases indicated in column four were found in the "ref_genome" and
one other genome, etc.* genome_codes_not_sharing: The ID numbers of those genomes in which
no alignment to the bases indicated in column four was found,
separated by commas. If there is a "0" in this column, the bases in
column four were found to have alignments in all of the input
genomes.* num_shared_positions: Total base count of all alignments with the
characteristics indicated in the preceding three columns.
LICENSE:
nucmer_backbone.pl is licensed under the GNU General Public License
version 3. See LICENSE.txt
CONTACT:
Contact Egon Ozer (e-ozer@northwestern.edu) with questions or
comments.

