#!/bin/bash
# TEMP2 (Transposable Element Movement present in a Population)
# 2019-04-04
# Tianxiong Yu(yutianxiong@gmail.com)
# Zhiping Weng Lab
# The School of Life Sciences and Technology, Tongji University
# Programs in Bioinformatics and Integrative Biology, University of Massachusetts Medical School

echo=echo0

#########
# usage #
#########
usage() {
	echo -e "usage"
}

###################
# read parameters #
###################
while getopts "hl:r:i:c:f:m:M:o:R:t:N:g:dI:vp:U:LSA" OPTION
do
        case $OPTION in
                h)	usage && exit 1;;
                i)	BAM=`readlink -f $OPTARG`;;
		I)	INDEX=`readlink -f $OPTARG`;;
                r)	RIGHT=`readlink -f $OPTARG`;;
                l)	LEFT=`readlink -f $OPTARG`;;
		p)      PREFIX=$OPTARG;;
	        f)	INSERT=$OPTARG;;
	        g)	GENOME=`readlink -f $OPTARG`;;
	        M)	MISMATCH=$OPTARG;;
		U)	UNIQ_RATIO=$OPTARG;;
	        m)	DIV=$OPTARG;;
                o)	OUTDIR=`readlink -f $OPTARG`;;
                c)	CPU=$OPTARG;;
	        R)	TESEQ=`readlink -f $OPTARG`;;
	        t)	RMSK=`readlink -f $OPTARG`;;
		N)	RMSK_WINDOW=$OPTARG;;
		d)	CLEAN=d;;
		L)	LOOSE_OVERLAP=1;;
		S)	SKIP_SHORT=1;;
		A)	ALU_MODE=1;;
		v)	echo -e "\nTEMP2:\t\tVersion ${TEMP2_VERSION}\nConstruct time:\tApril 29, 2019\n" && exit 1;;
                ?)	usage
        esac
done

####################
# check parameters #
####################
[ $# -lt 1 ] && usage
([ -z ${LEFT} ] || [ -z ${RIGHT} ]) && $echo 0 "Error: read1.fq and read2.fq not specified. Exiting..." && exit 1 
([ ! -f ${LEFT} ] || [ ! -f ${RIGHT} ]) && $echo 0 "Error: cannot access read1.fq (${LEFT}) or read2.fq (${RIGHT}). Exiting..." && exit 1 
[ -z ${BAM} ] && [ -z ${INDEX} ] && $echo 0 "Error: map.bam or bwa_index not specified. Exiting..." && exit 1 
[ ! -z ${BAM} ] && [ ! -f ${BAM} ] && $echo 0 "Error: cannot access map.bam (${BAM}). Exiting..." && exit 1 
[ ! -z ${BAM} ] && [ ! -f ${BAM}.bai ] && $echo 0 "Error: map.bam (${BAM}) is not indexed, please index it using samtools index. Exiting..." && exit 1 
[ ! -z ${INDEX} ] && [ ! -f ${INDEX}.amb ] && $echo 4 "Warning: cannot access bwa_index (${INDEX}). Will build it later..."
[ -z ${GENOME} ] && $echo 0 "Error: genome.fa not specified. Exiting..." && exit 1
[ ! -z ${GENOME} ] && [ ! -f ${GENOME} ] && $echo 0 "Error: cannot access genome.fa (${GENOME}). Exiting..." && exit 1 
[ -z ${TESEQ} ] && $echo 0 "Error: TE.fa not specified. Exiting..." && exit 1
[ ! -z ${TESEQ} ] && [ ! -f ${TESEQ} ] && $echo 0 "Error: cannot access TE.fa (${TESEQ}). Exiting..." && exit 1 
[ -z ${RMSK} ] && $echo 0 "Error: TE.bed6 not specified. Exiting..." && exit 1
[ ! -z ${RMSK} ] && [ ! -f ${RMSK} ] && $echo 0 "Error: cannot access TE.bed6 (${RMSK}). Exiting..." && exit 1 
[ -z ${PREFIX} ] && [ ! -z ${LEFT} ] && PREFIX=`basename ${LEFT}` && PREFIX=${PREFIX%[._]1.f*q*}
[ -z ${PREFIX} ] && [ ! -z ${BAM} ] && PREFIX=`basename ${BAM}` && PREFIX=${PREFIX%.sort*.bam}
[ -z ${OUTDIR} ] && OUTDIR=./${PREFIX}
[ ${CPU} -gt 0 ] 2>/dev/null || CPU=1
[ ${INSERT} -gt 0 ] 2>/dev/null || INSERT=Y
[ -z ${DIV%%*[!0-9.]*} ] && DIV=5
[ -z ${MISMATCH%%*[!0-9.]*} ] && MISMATCH=2
[ -z ${UNIQ_RATIO} ] && UNIQ_RATIO=0.8
[ -z ${RMSK_WINDOW%%*[!0-9]*} ] && RMSK_WINDOW=300
[ -z ${CLEAN} ] && CLEAN=p

######################
# check dependencies #
######################
function checkTools() {
        $echo 6 "${1}: \c" && which "$1" || { $echo 0 "Error: software ${1} not installed! Exiting..." && exit 1; }
}
$echo 6 "Testing required softwares:"
checkTools "bwa"
checkTools "samtools"
checkTools "bedtools"
checkTools "bedops"

##################
# Start pipeline #
##################
[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
cd ${OUTDIR} || { $echo 0 "cannot access out_path (${OUTDIR}). Exiting..." && exit 1; }
[ -z ${RMSK} ] && echo -e "chrNA\t1\t2\tNA\t0\t+" > ${PREFIX}.tmp.rmsk.bed && RMSK=${PREFIX}.tmp.rmsk.bed
$echo 4 "------ Start pipeline ------"
# Genome mapping if not done yet
if [ -z ${BAM} ];then
	$echo 2 "bam file not specified, map raw reads tp genome via bwa mem"
	bwa mem -t ${CPU} ${INDEX} ${LEFT} ${RIGHT} > ${PREFIX}.sam 2>${PREFIX}.bwamem.log || \
		{ $echo 0 "Error: bwa mem failed, please check ${OUTDIR}/${PREFIX}.bwamem.log. Exiting..." && exit 1; }
	$echo 2 "transform sam to sorted bam and index it"
	samtools view -bhS -@ ${CPU} ${PREFIX}.sam > ${PREFIX}.bam
	samtools sort -@ ${CPU} -o ${PREFIX}.sorted.bam ${PREFIX}.bam
	rm ${PREFIX}.sam ${PREFIX}.bam
	samtools index -@ ${CPU} ${PREFIX}.sorted.bam
	BAM=${PREFIX}.sorted.bam
else
	[ ! -f ${BAM}.bai ] && $echo 2 "index bam file" && samtools index -@ ${CPU} ${BAM}
fi

