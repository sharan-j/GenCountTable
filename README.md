# GenCountTable


This code generates a count table of single base substitutions from a .bam file.

Installation:
git clone https://github.com/sharan-j/GenCountTable

Pre-requisites:
  
  R version 3.5.2
  OS: Linux, Mac

Pre-processing:
The starting point of this code is a bam file. 
In order to generate a bam file, the sequenced fastq.gz file from a next-generation sequencing run is first trimmed with a suitable trimming tool, e.g. Cutadapt (version 2.2), upto the target site by trimming the adapters and scaffold sequence. The trimmed reads, e.g. 20 bases corresposnding to the length of the target site, are then mapped with an aligner, e.g. bowtie2 (version 2.3.5.1), to the reference. Bowtie2 is run with --no-unal or default parameters to generate a sam file. The sam file can be converted to a bam file using Samtools (version 1.9).

Usage:

Input:
1. Tab-seperated reference file containg the name of the guide, and the guide sequence. For example,
  
  ```
  CTRL_DOCK3NO2  TAAGACTGAACAAGAATGGT
  ```

2. bam file
  
  
Running the code:
On the terminal, run the following command

```
Rscript GenCountTableFromBAM_sub.R /path/to/reference.txt /path/to/bamfile /path/to/output.txt
```

Output:
A tab-seperated count table giving the total read count and percentage of substitutions for each position in the target sequence. For example, 
  
  ```
  ID	seq	allCounts	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20
  CTRL_DOCK3NO2	TAAGACTGAACAAGAATGGT	11722	0	0	0.001	0	0.523	0	0	0	0.001	0.001	0	0.001	0	0	0	0	0	0	0	0
  ```

