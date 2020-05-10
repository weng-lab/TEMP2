TEMP2
=====

## Overview
TEMP2 is an algorithm for detecting transposon insertions using short-read whole-genome sequencing data. It can not only precisely detect germline transposon insertions, but also estimate the number of uninherited/somatic transposon insertions by removing artificial insertion introduced by chimeric reads.

Current version v2.01

Author: Tianxiong Yu (yutianxiong@gmail.com) Weng Lab

## Requirements
TEMP runs on Linux x86_64 systems. 
Following softwares are required by TEMP and should be included in the path:
Samtools (http://samtools.sourceforge.net/),
bedtools (http://code.google.com/p/bedtools/),
bwa (http://sourceforge.net/projects/bio-bwa/),
bedops (https://bedops.readthedocs.io/)
Perl package BioPerl is also required for running TEMP2 absence (http://www.bioperl.org/wiki/Main_Page).

## Installation
To install TEMP2, please run the following commands:
`git clone https://github.com/weng-lab/TEMP2`
By link *TEMP2/TEMP2* to sytem path, you should be all set.
`ln -s $PWD/TEMP2 your_bin_path/TEMP2`

## Usage
After installation, simply run `TEMP2` for usage information.

To detect new germline transposon insertions, please run:
`TEMP2 insertion2`
If you also want to estimate the number inherited/somatic transposon insertions for each transposon, please run:
`TEMP2 insertion`
Please run the following command to retrieve the absence of reference annotated transposon copies:
`TEMP2 absence`
