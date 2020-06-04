TEMP2
=====
   
## Overview   
TEMP2 is an algorithm for detecting transposon insertions using short-read whole-genome sequencing data. It can not only precisely detect germline transposon insertions, but also estimate the number of uninherited/somatic transposon insertions by removing artificial insertion introduced by chimeric reads.

If you use TEMP2 for transposon insertion detection, please cite:  
Yu T ... Weng Z. TEMP2: an algorithm for detecting germline and de novo transposon insertions using short-read whole-genome sequencing data. 2020.  
  
Current version

Author: Tianxiong Yu (yutianxiong@gmail.com) in Weng Lab   
If you have any questions or find any bugs please contact Tianxiong Yu through yutianxiong@gmail.com.
   
## Requirements   
* Linux x86_64 systems.
* Perl package [BioPerl](https://bioperl.org/) is needed for running **absence** module.
* Other third-party tools that have already pre-compiled and included with TEMP2. However, you may need later versions of [zlib](https://zlib.net/) and [glibc](https://www.gnu.org/software/libc/) installed, especially to get samtools and bedtools work. If TEMP2 output errors for any of the softwares below, you may need to install and have them in your PATH.
	* [bwa](http://bio-bwa.sourceforge.net/bwa.shtml)
	* [samtools](http://www.htslib.org/)
	* [bedtools](https://bedtools.readthedocs.io/en/latest/)
	* [bedops](https://bedops.readthedocs.io/en/latest/)
	* ParaFly from [trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)   
   
## Installation
With [git](https://git-scm.com/downloads) installed, simply type the following command to install **TEMP2**:
```
git clone https://github.com/weng-lab/TEMP2
ln -s $PWD/TEMP2/TEMP2 your_bin_path/
```

To avoid mixing the pre-compiled tools with your own versions, we do **not** recommend to add `/TEMP2/bin` to the `$PATH`.

Alternatively, you can also install TEMP2 after fetching [source code](http://publications.wenglab.org/TEMP2/Release):
```
tar -xzvf TEMP2.tar.gz
cd TEMP2
ln -s $PWD/TEMP2 your_bin_path/
```
  
## Testing
TEMP2 integrates a tested dataset in /TEMP2/test/  
To test if TEMP2 is successfully installed, you can run the command below:
```
cd TEMP2/test
TEMP2 insertion -l test.1.fastq.gz -r test.2.fastq.gz -I bwa_index/chr2L -g chr2L.fa -R transposon.fa -t rmsk.bed -o test_output -c 8
```

This command takes around 2 minutes with 8 avaiable CPU. A successul installation should outputs `test.insertion.bed` and `test.soma.summary.txt` in the folder `test_output`, which contains detected insertions and estimated *de novo* insertion number per genome respectively.
   
## Getting start
TEMP2 includes three modules:  
1. **insertion**. The most general module that detects both germline and *de novo* transposon insertions.
2. **insertion2**. A simplified module that only detects germline transposon insertions.
3. **absence**. The module annotates the absence of reference transposon copies. Same as `TEMP_absence.sh` in our previous version ––– **TEMP**.  
After installation, typing `TEMP2` in the terminal produces a list of usage options, which are:   

```
TEMP2 -- Version 0.1.1
 -- germline and de novo transposon insertion and deletion analysis

usage:
TEMP2 insertion -h	# for germline and somatic transposon insertion

TEMP2 insertion2 -h	# for only germline transposon insertion

TEMP2 absence -h	# for transposon absence or deletion
```

## Detect germline and *de novo* transposon insertions
### Arguments
To see the help information of **insertion** module, type `TEMP2 insertion` in the terminal.

```
Germline and somatic transposon insertion detection using short DNA-seq data.
Please send questions, suggestions and bug reports to yutianxiong@gmail.com.
Thank you for using it.
=============================================================================

usage:
TEMP2 insertion
	-l read1.fq	Read1 for paired-end DNA-seq, can be gzipped.
	-r read2.fq	Read2 for paired-end DNA-seq, can be gzipped.
	-i map.bam	Genome mapping file in sorted and indexed bam format. Use one of -i or -I
	-I bwa_index	You can also input bwa_index instead of input map.bam, TEMP2 can map DNA-seq to the genome.
	-g genome.fa	Genome sequences in fasta format.
	-R TE.fa	Transposon consensus sequences in fasta format.
	-t TE.bed6	Transposon copies annotated in genome in bed6 format.
			You can download it from UCSC or creat it by running RepeatMasker.
Options:
	-o out_path	Output directory for results. Default is assigned by the prefix of read1.fq.
	-p out_prefix	Prefix to output directory. Default is assigned by the prefix of read1.fq.
	-A 		Set this parameter to enable ALU mode, which allows the insertion between two concordantly mapped reads
	-M mismatch%	Percentage of mismatch allowed when anchor to genome. Default is 2.
	-m mismatch%	Percentage of mismatch allowed when mapping to TEs. Default is 5.
	-U ratio	The ratio between the second best alignment and the best alignment to judge if a read is uniquely mapped. Default is 0.8.
	-f frag_length	Fragment length of the library. Default is calculated based on the mapping result.
	-N reference_filter_window	window sizea (+-n) for filtering insertions overlapping reference insertions. Default is 300.
	-L		Set this parameter to use a looser criteria to filter reference annotated copy overlapped insertions; Default not allowed.
	-S		Set this parameter to skip insertion length checking; Default is to remove those insertions that are not full length of shorter than 500bp.
	-c cpu_number	Number of CPU used. Default is 1.
	-d		Set this parameter to delete tmp files. Default is moving them to folder tmpTEMP2.
	-h		Show this message.
	-v		Show version information.
```

The **insertion2** module accepts exactly the same arguments as **insertion**. However, you can also type `TEMP2 insertion2` to get the detailed help information.
### Outputs
**TEMP2 insertion** outputs several files. Assuming the output prefix is **test**, then the following files shall be present in the ouput directory:
1. **_test.insertion.bed_**: a bed format file with additional 8 columns (6+8).  

   **Column 1,2&3**: Reference genome position (chromosome, start, and end) of the transposon insertion.  
**Column 4**: Information of inserted transposon, including transposon name, start coordinate and end coordinate of the inserted transposon, and which strand it inserts. Separated by `:`  
**Column 5**: Frequency of the inserted transposon. It generally means what fraction of sequenced genome present this insertion.  
**Column 6**: Which strand this transposon inserts.  
**Column 7**: Type of the insertion. We typically separate insertions into three categories based on how many reads support them.  
*1p1*: Insertions supported by reads at both ends.  
*2p*: Insertions supported by multipe reads at only one end.  
*singleton*: Insertions supported by only one read.  
*Ps: We used to regard 1p1 insertions as confident germline insertions but are abondon this filtering creteria now. The reason is that our recent study found supporting reads showed bias in the two ends of some transposons, especially LINE elements in fruit fly.*  
**Column 8**: Number of reads supporting this insertion.  
*Ps: If you are annotate confident transposon insertions, make sure to filter insertions without enough supporting reads, as many of them are false positives.*  
**Column 9**: Number of reads that do not support this insertion, AKA reference reads.  
**Column 10**: Number of supporting reads at 5'end of the insertion.  
**Column 11**: Number of supporting reads at 3'end of the insertion.  
**Column 12**: Target site duplication (TSD) of the insertion. *unknown* is shown if not applicable.  
**Column 13**: Reliability of this insertion (0–100). 100 for *2p* and *1p1* insertions. For *singleton* insertions, TEMP2 already filtered out most of the false positives but not all of them. The reliability is a percentage stand for how many *singleton* insertions of a specific transposon is.  
**Column 14**: Number of supporting reads at 5'end of the insertion junction.  
**Column 15**: Number of supporting reads at 3'end of the insertion junction.  

2. **_test.soma.summary.txt_**: a tab delimited file includes 10 columns (not included when running **insertion2** module)

   **Column 1**: Name of transposon.  
   **Column 2**: Estimated number of *de novo* insertion of this transposon per genome.  
   **Column 3**: 95th percentile (lambda distribution) number of *de novo* insertion of this transposon per genome.   
   **Column 4**: Total estimated number of *de novo* insertion of this transposon.  
   **Column 5**: 95th percentile (lambda distribution) number of *de novo* insertion of this transposon.  
   **Column 6**: Number of singleton reads mapped to the end regions (positive regions) of this transposon.  
   **Column 7**: Number of singleton reads mapped to the center region (negative region) of this transposon.  
   **Column 8**: Number of reads mapped to the end regions (positive regions) of this transposon.  
   **Column 9**: Number of reads mapped to the center region (negative region) of this transposon.  
   **Column 10**: Status of the estimation.  

3. **_test.supportReadsUnfiltered.bb_**: bigBed6 file for all supporting reads.

4. **_test.supportingRead.dis.pdf_**: figures showing where the supporting reads mapped to each transposon.

5. **_tmpTEMP2_**: folder contains all the intermediate files. eg: you can find fragment length statistics in **_tmpTEMP2/test.fragL_**

## Detect transposon absense in the reference genome
This module renders the same code as **TEMP**. Typing `TEMP2 absence` shows the help information:
```

usage: ./TEMP2_absence.sh -i input_file.sorted.bam -s scripts_directory -o output_directory -r transposon_rpmk.bed -t reference.2bit -f fragment_size -c CPUs -h

TEMP is a software package for detecting transposable elements (TEs)
insertions and excisions from pooled high-throughput sequencing data.
Please send questions, suggestions and bug reports to:
jiali.zhuang@umassmed.edu

Options:
        -i     Input file in bam format with full path. Please sort and index the file before calling this program.
               Sorting and indexing can be done by 'samtools sort' and 'samtools index'
        -s     Directory where all the scripts are
        -o     Path to output directory. Default is current directory
        -r     Annotated transposon positions in the genome (e.g., repeakMask) in bed6 format with full path
        -t     2bit file for the reference genome (can be downloaded from UCSC Genome Browser)
        -f     An integer specifying the length of the fragments (inserts) of the library. Default is 500
        -x     The minimum score difference between the best hit and the second best hit for considering a read as uniquely mapped. For BWA MEM.
        -c     An integer specifying the number of CPUs used. Default is 4
        -h     Show help message
```

For transposon absence analysis, the summay output file remains exactly the same as TEMP.    
1. **_test.absence.refined.bp.summary_**: There are 9 columns in the summary file and their meanings are listed below:  
   **Column 1,2&3**: Reference genome position (chromosome, start, and end) of the transposon absence.  
**Column 4**: The TE family that the detected insertion belongs to.  
**Column 5**: Junctions at 5’ of the excised TE. The two numbers are the coordinates of the junctions on the two strands.  
**Column 6**: Junctions at 3’ of the excised TE. The two numbers are the coordinates of the junctions on the two strands.  
**Column 7**: The number of reads supporting the absence.  
**Column 8**: The number of reads supporting the reference (no absence).  
**Column 9**: Estimated population frequency of the detected absence event.  

