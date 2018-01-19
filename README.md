# Spine

## INTRODUCTION:

Spine is a program for identification of the conserved core genome of bacteria and other small genome organisms. 

## REQUIREMENTS:

- Perl 5.10 or above
- [MUMmer](http://mummer.sourceforge.net) version 3.22 or above. Install MUMmer as directed by the instructions included with the software.
- Mac OSX or Linux. We provide no guarantees that this will work on Windows or other operating systems.

## INSTALLATION:

Simply move the Spine directory to the desired location. The "scripts" directory must remain in the same directory as spine.pl

## USAGE:

Basic command: `perl spine.pl -f genome_files.txt`

For list of options, call the script without any inputs: `perl spine.pl`

### Required Inputs:

`-f`	File with list of input sequence files. This file should beformatted like so:
```        
path/to/file1<tab>unique_identifier<tab>fasta or gbk
path/to/file2<tab>unique_identifier<tab>fasta or gbk
```        
Example:
```
/home/seqs/PAO1.fasta   PAO1    fasta
/home/seqs/LESB58.gbk   LESB58  gbk
```
        
The third column (fasta or gbk) is optional, but should be given if your sequence files end with suffixes other than ".fasta" or ".gbk".

If you have genomes spread across multiple files (i.e. chromosomes and/or plasmids), these can be combined by either concatenating the files into one:  
```
cat chrom_I.gbk chrom_II.gbk > combined.gbk
```
or by including all the files in this input file, separated by commmas:  
Example:
```
/seqs/chrom_I.fasta,/seqs/chrom_II.fasta    mygenome    fasta
chrom_A.gbk,chrom_B.gbk,plasmid_X.gbk   myothergenome   gbk
```

### Optional Inputs

`-a` or `--pctcore`  
Percentage of input genomes in which a region must be found in order to be considered core.   (default: 100)

  `-g` or `--maxdist`  
Maximum distance between core genome segments. Distances less than this between adjacent segments will result in combination of fragments with N's rather than separating into two or more fragments.  
(default: 10)

  `-l` or `--license`  
Print license information and quit  

  `-m` or `--nucpath`  
Full path to folder containing MUMmer scripts and executables, i.e. /home/applications/MUMmer/bin  
(default: tries to find MUMmer in your PATH)

  `-r` or `--refs`  
Reference genome sequence(s) to use as primary output source(s). This should be one or more integers corresponding to the order of the genomes given in the file above, i.e. "1" would use the first-listed sequence, "3" would use the third-listed, etc. To prioritize multiple genome sequences, separate the integers with commas, i.e. "1,3" for giving sequence 1 the highest priority and sequence 3 the next-highest priority. Reference sequences will serve as the source of backbone sequences to be output, as well as the source of backbone locus IDs, if applicable.

The number of reference genomes used will depend on the definition of core genome given by option -a. For instance, if core is determined from 10 input genomes and -a is given as 100, then core sequence will only be taken from one reference genome. If, for example, -a is given as 90 from 10 input genomes, then potentially two reference sequences will be needed: The first for sequences present in all 10 genomes and for sequences present in 9 out of 10 genomes including the first genome. The second reference sequence would then be used as the source of all sequences present in 9 out of 10 genomes, but not present in the first reference genome.

(default: reference priority will be the same as the order of genomes entered, with the first genome having the highest priority and the last genome having the lowest priority)

  `-o` or `--prefix`  
Output prefix.  
(default: "output")

  `-p` or `--pctid`  
Minimum percent identity for regions to be considered homologous.  
(default: 85)

  `-s` or `--minout`  
Minimum size of core region sequences to be output, in bases.  
(default: 10)

  `-t` or `--threads`  
Number of parallel processes to run.   
(default: 4)

_Careful:_ This script does not perform any verification of the number of processers available. If you set this number higher than the number of processors you have, performance is likely to be significantly degraded.
 
 `-v` or `--version`  
Print version information and quit.

** Nucmer Options **
Advanced use only. Little reason to change defaults in most situations. 
See MUMmer documentation for more information.  
  `--breaklen`        Integer (default: 200)  
  `--mincluster`      Integer (default: 65)  
  `--diagdiff`        Integer (default: 5)  
  `--diagfactor`      Float (default: 0.12)  
  `--minmatch`        Integer (default: 20)  
  `--nosimplify`      (default: simplify)  

## OUTPUT FILES:

__statistics.txt__  
First line shows the current software version used.  
Second line shows the input parameters given to the software.  
_Column headers and descriptions:_
* gen_#: Number assigned to the genome based on the order the sequences were input to Spine.
* gen_name: Name of the sequence, user-assigned
* gen_size: Total size of the sequence, in bases
* source: Indicates whether the row describes the strain's accessory or core genome. If this column is "backbone" it describes the characteristics of the sample core gneome. If this column is "pangenome" it describes the characteristics of the sample pangenome.
* total_bp: Sequence size, in bases
* gc_%: Percent GC content of the sequence
* num_segs: Number of separate sequence segements output
* min_seg: Smallest segment size, in bases
* max_seg: Largest segment size, in bases
* avg_leng: Average length of the output segments
* median_leng: Median length of the output segments
* num_cds (if annotation was provided): number of coding sequences present. A coding sequence is counted as present within either the core or the accessory genome if 50% or greater of the length of coding sequence is found in sequences within that genome fraction.

__coords.txt__  
Coordinates of genome sequences.   
"*.accessory_coords.txt": Accessory genome sequences for the indicated strain  
"*.core_coords.txt": Core genome sequences for the indicated strain  
"backbone_coords.txt": Core genome sequences for the group of strains  
"pangenome_coords.txt" (if requested): Pangeome sequences for the group of strains  
_Column headers and descriptions:_  
* contig_id: sequence ID of the source strain contig
* contig_length: length, in bases, of the source strain contig
* start: start coordinate of the genome segement on the source strain contig
* stop: stop coordinate of the genome segment on the source strain contig
* source_gen: (only for backbone or pangenome) genome name of the source strain
* out_seq_id: sequence ID of the segment as found in the corresponding sequence file output by Spine 

__*.fasta__  
Nucleotide sequences of the genome segments output by Spine. Original sources of the sequences can be determined by cross-referencing the sequence IDs with the cooresponding coords.txt file

__loci.txt__ (if annotated genbank file was provided for one or more genomes)  
List of coding sequences found in the core genome.  
"*.accessory_loci.txt": Accessory genome coding sequences for the indicated strain  
"*.core_loci.txt": Core genome coding sequences for the indicated strain  
"backbone_loci.txt": Core genome coding sequences for the group of strains  
"pangenome_loci.txt" (if requested): Pangeome coding sequences for the group of strains   
_Column headers and descriptions:_  
* locus_id: locusID of gene
* gen_contig_id: Source strain contig ID
* gen_contig_start: Gene start coordinate in source sequence (1-based)
* gen_contig_stop: Gene stop coordinate in source sequence (1-based)
* strand: Strand on which the gene is encoded
* out_seq_id: Output sequence ID (corresponds to sequence IDs in corresponding sequence fasta file above)
* out_seq_start: Gene start coordinate in output sequence (1-based)
* out_seq_stop: Gene stop coordinate in output sequence (1-based)
* pct_locus: Percentage of gene represented in the output sequence
* source_gen: (only for backbone or pangenome) Source genome name
* overhangs: Number of bases of the gene missing from the end(s) of the output segment. Values are separated by a comma. First value is the number of bases missing from the 5' end of the core segement, second value is the number of bases missing from the 3' end of the core segment.
* product: Gene product

__position_counts.txt__  
This file should not be needed for routine use. Is meant to be used as input for core_and_pangenome.pl to calculate core-, pan-, and new genome sizes at permutations of the input genomic sequences.  
_Column headers and descriptions:_  
* ref_genome: ID number of the genome (based on the order of the sequences input to Spine) used as a reference in alignment to the other input genomes.
* num_genomes_sharing: The number of genomes in which the sequences were found, i.e “1” indicates that the bases indicated in column four were only found in one genome (the “ref_genome”), “2” indicates that the bases indicated in column four were found in the “ref_genome” and one other genome, etc.
* genome_codes_not_sharing: The ID numbers of those genomes in which no alignment to the bases indicated in column four was found, separated by commas. If there is a “0” in this column, the bases in column four were found to have alignments in all of the input genomes.
* num_shared_positions: Total base count of all alignments with the characteristics indicated in the preceding three columns.

## LICENSE:

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
along with this program.  See LICENSE.txt

## CONTACT:

Contact [Egon Ozer](e-ozer@northwestern.edu) with questions or comments.


