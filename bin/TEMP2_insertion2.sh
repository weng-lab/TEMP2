#!/bin/bash
# TEMP2 (Transposable Element Movement present in a Population)
# 2019-04-04
# Tianxiong Yu(yutianxiong@gmail.com)
# Zhiping Weng Lab
# The School of Life Sciences and Technology, Tongji University
# Programs in Bioinformatics and Integrative Biology, University of Massachusetts Medical School

#########
# usage #
#########
usage() {
$echo 5 "Germline transposon insertion detection using short DNA-seq data."
$echo 5 "Please send questions, suggestions and bug reports to yutianxiong@gmail.com."
$echo 5 "Thank you for using it."
$echo 5 "=============================================================================\n"
$echo 4 "usage:"
$echo 4 "TEMP2 insertion"
$echo 4 "\t-l read1.fq\tRead1 for paired-end DNA-seq, can be gzipped."
$echo 4 "\t-r read2.fq\tRead2 for paired-end DNA-seq, can be gzipped."
$echo 4 "\t-i map.bam\tGenome mapping file in sorted and indexed bam format. Use one of -i or (-I + -l + -r)"
$echo 4 "\t-I bwa_index\tYou can also input bwa_index instead of input map.bam, TEMP2 can map DNA-seq to the genome."
$echo 4 "\t-g genome.fa\tGenome sequences in fasta format."
$echo 4 "\t-R TE.fa\tTransposon consensus sequences in fasta format."
$echo 4 "\t-t TE.bed6\tTransposon copies annotated in genome in bed6 format."
$echo 4 "\t\t\tYou can download it from UCSC or creat it by running RepeatMasker."
$echo 6 "Options:"
$echo 6 "\t-o out_path\tOutput directory for results. Default is assigned by the prefix of read1.fq."
$echo 6 "\t-p out_prefix\tPrefix to output directory. Default is assigned by the prefix of read1.fq."
$echo 6 "\t-A \t\tSet this parameter to enable ALU mode, which allows the insertion between two concordantly mapped reads"
$echo 6 "\t-M mismatch%\tPercentage of mismatch allowed when anchor to genome. Default is 2."
$echo 6 "\t-m mismatch%\tPercentage of mismatch allowed when mapping to TEs. Default is 5."
$echo 6 "\t-U ratio\tThe ratio between the second best alignment and the best alignment to judge if a read is uniquely mapped. Default is 0.8."
$echo 6 "\t-f frag_length\tFragment length of the library. Default is calculated based on the mapping result."
$echo 6 "\t-N reference_filter_window\twindow sizea (+-n) for filtering insertions overlapping reference insertions. Default is 300."
$echo 6 "\t-L\t\tSet this parameter to use a looser criteria to filter reference annotated copy overlapped insertions; Default not allowed."
$echo 6 "\t-S\t\tSet this parameter to skip insertion length checking; Default is to remove those insertions that are not full length of shorter than 500bp."
$echo 6 "\t-c cpu_number\tNumber of CPU used. Default is 1."
$echo 6 "\t-d\t\tSet this parameter to delete tmp files. Default is moving them to folder tmpTEMP2."
$echo 6 "\t-h\t\tShow this message."
$echo 6 "\t-v\t\tShow version information."
exit 1
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
		v)	$echo 5 "\nTEMP2:\t\tVersion ${TEMP2_VERSION}\nConstruct time:\tApril 29, 2019\n" && exit 1;;
                ?)	usage
        esac
done

####################
# check parameters #
####################
[ $# -lt 1 ] && usage
#([ -z ${LEFT} ] || [ -z ${RIGHT} ]) && $echo 0 "Error: read1.fq and read2.fq not specified. Exiting..." && exit 1 
#([ ! -f ${LEFT} ] || [ ! -f ${RIGHT} ]) && $echo 0 "Error: cannot access read1.fq (${LEFT}) or read2.fq (${RIGHT}). Exiting..." && exit 1 
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
[ -z ${RMSK_WINDOW%%*[!0-9]*} ] && RMSK_WINDOW=300
[ -z ${CLEAN} ] && CLEAN=p

######################
# check dependencies #
######################
export PATH=${PATH}:${BINDIR} # export TEMP2/bin as temporary path
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

# Get concordant-uniq-split mapped reads
samtools view -H ${BAM} > ${PREFIX}.tmp.header
if [ -f ${PREFIX}.pair.uniq.split.fastq ] && [ -f ${PREFIX}.pair.uniq.split.bed ];then
	$echo 2 "concordant-uniq-split reads exists, skip"
else
	$echo 2 "get concordant-uniq-edge-split reads"
	samtools view -@ ${CPU} -h -f 0X2 -F 0X800 ${BAM} | awk -v uniq_ratio=${UNIQ_RATIO} 'BEGIN{FS=OFS="\t"} {for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(xs/1<as*uniq_ratio){print $0}}' - | grep "SA:Z:" | awk '(and($2,16)==0 && $6~/^[0-9]+S/) || (and($2,16)==16 && $6~/S$/)' - | cat ${PREFIX}.tmp.header - | samtools view -@ ${CPU} -bhS - | samtools sort -@ ${CPU} -o ${PREFIX}.pair.uniq.split.bam -
	samtools fastq -@ ${CPU} ${PREFIX}.pair.uniq.split.bam | awk '{if(NR%4==1){print substr($1,1,length($1)-2)"_"substr($1,length($1))}else{print $0}}' - > ${PREFIX}.pair.uniq.split.fastq 
	bedtools bamtobed -i ${PREFIX}.pair.uniq.split.bam -tag NM | awk -v div=${MISMATCH} 'BEGIN{FS=OFS="\t"} {$4=substr($4,1,length($4)-2)"_"substr($4,length($4));$6=$6=="+"?"-":"+";if($5<=(int(($3-$2-1)*div/100)+1)){print $0}}' - > ${PREFIX}.pair.uniq.split.bed
	if [ ${ALU_MODE} ];then
		$echo 2 "ALU mode enabled, also get concordant-uniq-internal-split reads"
		samtools view -@ ${CPU} -h -f 0X2 -F 0X800 ${BAM} | awk -v uniq_ratio=${UNIQ_RATIO} 'BEGIN{FS=OFS="\t"} {for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(xs/1<as*uniq_ratio){print $0}}' - | grep "SA:Z:" | awk '(and($2,16)==16 && $6~/^[0-9]+S/) || (and($2,16)==0 && $6~/S$/)' - | cat ${PREFIX}.tmp.header - | samtools view -@ ${CPU} -bhS - | samtools sort -@ ${CPU} -o ${PREFIX}.pair.uniq.split.internal.bam -
		samtools fastq -@ ${CPU} ${PREFIX}.pair.uniq.split.internal.bam | awk 'BEGIN{FS=OFS="\t"} {if(NR%4==1){print substr($1,1,length($1)-2)"_"substr($1,length($1))}else{print $0}}' - > ${PREFIX}.tmp
		${BINDIR}/revertFq.sh ${PREFIX}.tmp >> ${PREFIX}.pair.uniq.split.fastq 
		bedtools bamtobed -i ${PREFIX}.pair.uniq.split.internal.bam -tag NM | awk -v div=${MISMATCH} 'BEGIN{FS=OFS="\t"} {$4=substr($4,1,length($4)-2)"_"substr($4,length($4));if($5<=(int(($3-$2-1)*div/100)+1)){print $0}}' - >> ${PREFIX}.pair.uniq.split.bed
	fi
fi

# check fragment length
$echo 2 "check fragment length"
samtools view -@ ${CPU} -F 0X4 -f 0X2 -F 0X800 ${BAM} | head -1000000 | awk '{if(and($2,128)==0 && $9!=0){print $9<0?-$9:$9}}' > ${PREFIX}.tmp
fragL=(`Rscript ${BINDIR}/checkFragL.R ${PREFIX}.tmp`)
echo -e "95_quantile\t${fragL[3]}\nstandard_deviation\t${fragL[1]}\naverage\t${fragL[5]}" > ${PREFIX}.fragL
[ "${INSERT}" == "Y" ] && INSERT=${fragL[3]} && $echo 1 "insert size set to 95 quantile: ${INSERT}"
[ ${fragL[1]} -gt 100 ] && $echo 4 "WARNING: standard deviation of insert size is higher than 100 (${fragL[1]})"

# Get the mate seq of the uniq-unpaired reads
if [ -f ${PREFIX}.unpair.uniq.1.fastq ] && [ -f ${PREFIX}.unpair.uniq.2.fastq ] && [ -f ${PREFIX}.unpair.uniq.bed ];then
	$echo 2 "mate seq of the uniq-unpaired reads exists, skip"
else
	$echo 2 "get mate seq of the uniq-unpaired"
	$BINDIR/pickUniqPairFastq.sh ${BAM} ${PREFIX} $CPU ${INSERT} ${UNIQ_RATIO}
fi

# Map to transposons
$echo 2 "map paired split uniqMappers and unpaired uniqMappers to transposons"
bwa index ${TESEQ} -p ${PREFIX}.tmp.te.index > /dev/null 2>&1
${BINDIR}/faToChromSize ${TESEQ} > ${PREFIX}.tmp.te.size
bwa mem -T 20 -Y -t ${CPU} ${PREFIX}.tmp.te.index ${PREFIX}.pair.uniq.split.fastq > ${PREFIX}.pair.uniq.split.transposon.sam 2>/dev/null # split reads
samtools view -@ ${CPU} -hSF 0X2 -F 0X800 ${PREFIX}.pair.uniq.split.transposon.sam > ${PREFIX}.t && mv ${PREFIX}.t ${PREFIX}.pair.uniq.split.transposon.sam 
bwa mem -T 20 -Y -t ${CPU} ${PREFIX}.tmp.te.index ${PREFIX}.unpair.uniq.1.fastq ${PREFIX}.unpair.uniq.2.fastq > ${PREFIX}.unpair.uniq.transposon.sam 2>/dev/null # unpaired reads
samtools view -@ ${CPU} -hSF 0x4 -F 0X800 ${PREFIX}.unpair.uniq.transposon.sam > ${PREFIX}.t && mv ${PREFIX}.t ${PREFIX}.unpair.uniq.transposon.sam
READ_LENGTH=`awk '{if($6~/M$/ && $6!~/S|I|D/){print substr($6,1,length($6)-1);exit}}' ${PREFIX}.unpair.uniq.transposon.sam`

# Merge fragments in genome and transposon; also remove PCR duplicates
$echo 2 "merge fragments in genome and transposon"
bash $BINDIR/mergeGenomeTESplit.sh ${PREFIX}.pair.uniq.split.transposon.sam ${PREFIX}.pair.uniq.split.bed ${DIV} ${PREFIX}.tmp.te.size ${READ_LENGTH} | awk 'BEGIN{FS=OFS="\t"} {if($5~/^;/){$5=substr($5,2)};print $0}' | awk '!a[$1"."$2"."$3"."$5]++' > ${PREFIX}.pair.uniq.split.transposon.bed
bash $BINDIR/mergeGenomeTEUnpair.sh ${PREFIX}.unpair.uniq.transposon.sam ${PREFIX}.unpair.uniq.bed ${DIV} ${PREFIX}.tmp.te.size ${READ_LENGTH} ${MISMATCH} ${INSERT} | awk 'BEGIN{FS=OFS="\t"} {if($5~/^;/){$5=substr($5,2)};print $0}' | awk '!a[$1"."$2"."$3"."$5]++' > ${PREFIX}.unpair.uniq.transposon.bed

# Fix multipleMappers due to LTRs in both sides of transposons
${BINDIR}/fixLTR.sh ${PREFIX}.pair.uniq.split.transposon.bed > ${PREFIX}.pair.uniq.split.transposon.fixLTR.bed
${BINDIR}/fixLTR.sh ${PREFIX}.unpair.uniq.transposon.bed > ${PREFIX}.unpair.uniq.transposon.fixLTR.bed
${BINDIR}/faToChromSize ${GENOME} > ${PREFIX}.tmp.chr.size
cat ${PREFIX}.unpair.uniq.transposon.fixLTR.bed ${PREFIX}.pair.uniq.split.transposon.fixLTR.bed | awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$5,0,$6}' | sort -k1,1 -k2,2n > ${PREFIX}.t && bedToBigBed ${PREFIX}.t ${PREFIX}.tmp.chr.size ${PREFIX}.supportReadsUnfiltered.bb && rm ${PREFIX}.t

# Merge supporting reads within insert size - read length and in the same direction
$echo 2 "merge support reads in the same direction within ${INSERT} - ${READ_LENGTH}"
mkdir ${PREFIX}.supportReads
awk -v prefix=${PREFIX} 'BEGIN{FS=OFS="\t"} {split($5,a,";");for(i in a){split(a[i],b,",");if(b[4]==$6){ts="-"}else{ts="+"};print $1,$2,$3,1/length(a),a[i],$6>>prefix".supportReads/"b[1]"."ts".bed"}}' ${PREFIX}.unpair.uniq.transposon.fixLTR.bed ${PREFIX}.pair.uniq.split.transposon.fixLTR.bed
MAXGAP=`expr ${INSERT} - ${READ_LENGTH}`
[ $MAXGAP -lt 0 ] && MAXGAP=0
for i in ${PREFIX}.supportReads/*.bed; do sort -k1,1 -k2,2n $i | awk 'BEGIN{FS=OFS="\t"} {$5=$5","$4;print $0}' | bedtools merge -i - -d 375 -s -c 4,5,6,2,3 -o collapse,collapse,first,collapse,collapse -delim "|" > ${i/.bed/.merged.bed}; done
for i in ${PREFIX}.supportReads/*merged.bed; do
	tel=`awk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){a[$1]=$2}else{if(FNR==1){split(FILENAME,fn,".merged.bed");split(fn[1],b,"/");print a[substr(b[2],1,length(b[2])-2)]}}}' ${PREFIX}.tmp.te.size $i`
	${BINDIR}/processMergedBed.py $i $tel $INSERT > ${i/merged.bed/processed.bed}
done

# Merge supporting reads within 2 * insert size - read length and in direction + then -
$echo 2 "merge support reads in different direction within 2 X ${INSERT} - ${READ_LENGTH}"
for i in ${PREFIX}.supportReads/*.processed.bed; do ${BINDIR}/mergeProcessedBed.sh $i $INSERT $READ_LENGTH ${i%.processed.bed}; done
cat ${PREFIX}.supportReads/*.final.bed > ${PREFIX}.final.bed

# filter false positive insertions which overlap with the same transposon insertion annotation or in high depth region
$echo 2 "filter candidate insertions which overlap with the same transposon insertion or in high depth region"
${BINDIR}/fixCoor.sh ${PREFIX}.final.bed | awk 'BEGIN{FS=OFS="\t"} {$1=$1"__"$4;print $0}' | sort -k1,1 -k2,2n | awk '$2>=0' > ${PREFIX}.tmp
RAWINS_UNFILTER=(`wc -l ${PREFIX}.tmp`)
awk 'BEGIN{FS=OFS="\t"} {$1=$1"__"$4;print $0}' ${RMSK} > ${PREFIX}.tmp.rmsk.bed && RMSK=${PREFIX}.tmp.rmsk.bed
if [ -z ${LOOSE_OVERLAP} ];then
	bedtools window -w ${RMSK_WINDOW} -a ${PREFIX}.tmp -b ${RMSK} -v | awk 'BEGIN{FS=OFS="\t"} {split($1,a,"__");$1=a[1];$4=a[2];print $0}' > ${PREFIX}.insertion.raw.bed
	bedtools window -w ${RMSK_WINDOW} -a ${PREFIX}.tmp -b ${RMSK} | awk 'BEGIN{FS=OFS="\t"} {split($1,a,"__");print a[1],$2,$3,0,0,$6}' > ${PREFIX}.removed.bed
	bedtools window -w ${RMSK_WINDOW} -a ${PREFIX}.tmp -b ${RMSK} | awk 'BEGIN{FS=OFS="\t"} {if($7=="1p1"){split($1,a,"__");print a[1],$2,$3,0,0,$6}}' > ${PREFIX}.removed.1p1.bed
else
	bedtools window -w ${RMSK_WINDOW} -a ${PREFIX}.tmp -b ${RMSK} -v -sm | awk 'BEGIN{FS=OFS="\t"} {split($1,a,"__");$1=a[1];$4=a[2];print $0}' > ${PREFIX}.insertion.raw.bed
	bedtools window -w ${RMSK_WINDOW} -a ${PREFIX}.tmp -b ${RMSK} -sm | awk 'BEGIN{FS=OFS="\t"} {split($1,a,"__");print a[1],$2,$3,0,0,$6}' > ${PREFIX}.removed.bed
	bedtools window -w ${RMSK_WINDOW} -a ${PREFIX}.tmp -b ${RMSK} -sm | awk 'BEGIN{FS=OFS="\t"} {if($7=="1p1"){split($1,a,"__");print a[1],$2,$3,0,0,$6}}' > ${PREFIX}.removed.1p1.bed
fi
awk '$7!="1p1"' ${PREFIX}.insertion.raw.bed | bedtools window -w 50 -a - -b ${PREFIX}.removed.bed -v > ${PREFIX}.t
awk '$7=="1p1"' ${PREFIX}.insertion.raw.bed | bedtools window -w 50 -a - -b ${PREFIX}.removed.1p1.bed -v >> ${PREFIX}.t
RAWINS_FILTERRMSK=(`wc -l ${PREFIX}.t`) && FILTERRMSK=`expr ${RAWINS_UNFILTER} - ${RAWINS_FILTERRMSK}`
# filter false positive 1p1 insertions which are not full length but short than 500bp 
if [ -z ${SKIP_SHORT} ];then
	awk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){tel[$1]=$2}else{if($7=="1p1"){if(tel[$4]<=500 || ($9-$8)>=500){print $0}}else{print $0}}}' ${PREFIX}.tmp.te.size ${PREFIX}.t | sort -k1,1 -k2,2n > ${PREFIX}.insertion.raw.bed
else
	cat ${PREFIX}.t | sort -k1,1 -k2,2n > ${PREFIX}.insertion.raw.bed
fi
RAWINS_FILTER1P1=(`wc -l ${PREFIX}.insertion.raw.bed`) && FILTER1P1=`expr ${RAWINS_FILTERRMSK} - ${RAWINS_FILTER1P1}`
# filter singleton insertions
awk '$7!="singleton"' ${PREFIX}.insertion.raw.bed > ${PREFIX}.t && mv ${PREFIX}.t ${PREFIX}.insertion.raw.bed 
# filter false positive insertions in high depth regions
$echo 2 "filter candidate insertions in high depth region"
bedtools random -g ${PREFIX}.tmp.chr.size -l 200 -n 1000 >${PREFIX}.tmp.random.bed
REGIONS=`awk '{r=r" "$1":"$2"-"$3} END{print r}' ${PREFIX}.tmp.random.bed`
AVE_DEPTH=`samtools view -bh -F 0X4 -@ ${CPU} ${BAM} ${REGIONS} | bedtools bamtobed -i - | intersectBed -a - -b ${PREFIX}.tmp.random.bed -wo | awk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){a[$10]++}else{s+=a[$4]}} END{print s/1000}' - ${PREFIX}.tmp.random.bed`
CTOF_DEPTH=`awk -v t=${AVE_DEPTH} 'BEGIN{print t*5}'` && $echo 3 "average read number for 200bp bins is ${AVE_DEPTH}, set read number cutoff to ${CTOF_DEPTH}"
awk 'BEGIN{FS=OFS="\t"} {s=int(($2+$3)/2)-100;e=int(($2+$3)/2)+100;if(s<0){s=0};print $1,s,e,$1","$2","$3","$4","$6,0,"+"}' ${PREFIX}.insertion.raw.bed > ${PREFIX}.tmp
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' ${PREFIX}.tmp | sort -k1,1 -k2,2n | bedtools merge -i - > ${PREFIX}.t
TN=(`wc -l ${PREFIX}.t`)
if [ $TN -lt 50000 ];then
	REGIONS=`awk 'BEGIN{FS=OFS="\t"} {printf $1":"$2"-"$3" "}' ${PREFIX}.t`
	samtools view -bh -F 0X4 -@ ${CPU} ${BAM} ${REGIONS} | bedtools bamtobed -i - | intersectBed -a - -b ${PREFIX}.tmp -wo | awk -v adp=${AVE_DEPTH} 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){a[$10]++}else{if(a[$1","$2","$3","$4","$6]/1<adp*5){print $0}}}' - ${PREFIX}.insertion.raw.bed > ${PREFIX}.t && mv ${PREFIX}.t ${PREFIX}.insertion.raw.bed 
else
	awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' ${PREFIX}.t > ${PREFIX}.tmp.region
	samtools view -bh -F 0X4 -@ ${CPU} -L ${PREFIX}.tmp.region ${BAM} | bedtools bamtobed -i - | intersectBed -a - -b ${PREFIX}.tmp -wo | awk -v adp=${AVE_DEPTH} 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){a[$10]++}else{if(a[$1","$2","$3","$4","$6]/1<adp*5){print $0}}}' - ${PREFIX}.insertion.raw.bed > ${PREFIX}.t && mv ${PREFIX}.t ${PREFIX}.insertion.raw.bed 
fi
RAWINS_FILTERDEPTH=(`wc -l ${PREFIX}.insertion.raw.bed`) && FILTERDEPTH=`expr ${RAWINS_FILTER1P1} - ${RAWINS_FILTERDEPTH}`
$echo 3 "Filtered insertion number: ${RAWINS_UNFILTER} - ${FILTERRMSK} (overlap rmsk) ${FILTER1P1} (short insertion) - ${FILTERDEPTH} (high depth) = ${RAWINS_FILTERDEPTH}"

# get TPregions from 1p1 insertions which have at least 3 support reads in both sides
awk -v ins=${INSERT} 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){tel[$1]=$2}else{if(int(substr($10,2))>=3 && int(substr($11,2))>=3){s1=$8;e1=$8+ins-1;if(e1>$9 && s1<$9-50){e1=$9}else if(e1>tel[$4]){e1=tel[$4]};e2=$9;s2=$9-ins+1;if(s2<$8 && e2>$8+50){s2=$8}else if(s2<1){s2=1};print $4,s1,e1,0,0,"-";print $4,s2,e2,0,0,"+"}}}' ${PREFIX}.tmp.te.size ${PREFIX}.insertion.raw.bed | sort -k1,1 -k2,2n - | bedtools merge -i - -c 4,5,6 -s -o count,first,first > ${PREFIX}.tmp
awk -v ins=${INSERT} 'BEGIN{FS=OFS="\t"} {s1=1;e1=ins;if(e1>$2){e1=$2};e2=$2;s2=$2-ins+1;if(s2<1){s2=1};print $1,s1,e1,0,0,"-";print $1,s2,e2,0,0,"+"}' ${PREFIX}.tmp.te.size | cat - ${PREFIX}.tmp | sort -k1,1 -k2,2n - | bedtools merge -i - -c 4,5,6 -s -o sum,first,first > ${PREFIX}.TPregion.bed

# process raw insertion to filtered insertion
$echo 2 "filter unreliable 2p insertions overlapped with similar reference transposon copies"
awk 'BEGIN{FS=OFS="\t"} {print $0,1}' ${PREFIX}.TPregion.bed > ${PREFIX}.tmp
if [ -z ${LOOSE_OVERLAP} ];then
	awk 'BEGIN{FS=OFS="\t"} {if($7=="2p"){split($13,a,",");print a[1],a[2],a[3],0,0,a[4],$0}}' ${PREFIX}.insertion.raw.bed | intersectBed -a - -b ${PREFIX}.tmp -s -f 1 -wo | cut -f 7-23,30 | awk 'BEGIN{FS=OFS="\t"} {if($7=="singleton"){$13=int($18*10000)/100"%"}else{$13="100%"};print $0}' | cut -f 1-17 > ${PREFIX}.tmp2
	awk 'BEGIN{FS=OFS="\t"} {if($7=="1p1"){$13="100%";print $0}}' ${PREFIX}.insertion.raw.bed | cat - ${PREFIX}.tmp2 > ${PREFIX}.insertion.filtered.bed
else
	awk 'BEGIN{FS=OFS="\t"} {if($7=="1p1" || $7=="2p"){$13="100%";print $0}}' ${PREFIX}.insertion.raw.bed > ${PREFIX}.insertion.filtered.bed
fi

# Calculate frequency of each transposon insertion
$echo 2 "Calculate frequency of each transposon insertion"
awk -v rl=60 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){cs[$1]=$2}else{n=$1","$2","$3","$4","$6;if(int(($14+$15)/2)-rl/2<0){s=0}else{s=int(($14+$15)/2)-rl/2};if(int(($14+$15)/2)+rl/2>cs[$1]){e=cs[$1]}else{e=int(($14+$15)/2)+rl/2};if(s>=e){s=e-1};print $1,s,e,n,0,"."}}' ${PREFIX}.tmp.chr.size ${PREFIX}.insertion.filtered.bed > ${PREFIX}.tmp
awk -v ins=${INSERT} 'BEGIN{FS=OFS="\t"} {s=$2-ins;e=$3+ins;if(s<0){s=0};print $1,s,e}' ${PREFIX}.tmp | sort -k1,1 -k2,2n | bedtools merge -i - > ${PREFIX}.t
TN=(`wc -l ${PREFIX}.t`)
if [ $TN -lt 50000 ];then
	REGIONS=`awk 'BEGIN{FS=OFS="\t"} {printf $1":"$2"-"$3" "}' ${PREFIX}.t`
	samtools view -@ ${CPU} -bh -F 0X800 -f 0X2 ${BAM} ${REGIONS} | samtools sort -@ ${CPU} -n - | bedtools bamtobed -bedpe -i - 2>/dev/null | awk '$1!="." && $2<$6 && $2>=0' | cut -f 1-2,6-9 > ${PREFIX}.tmp.bed
else
	awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' ${PREFIX}.t > ${PREFIX}.tmp.region
	samtools view -@ ${CPU} -bh -F 0X800 -f 0X2 -L ${PREFIX}.tmp.region ${BAM} | samtools sort -@ ${CPU} -n - | bedtools bamtobed -bedpe -i - 2>/dev/null | awk '$1!="." && $2<$6 && $2>=0' | cut -f 1-2,6-9 > ${PREFIX}.tmp.bed
fi
intersectBed -a ${PREFIX}.tmp -b ${PREFIX}.tmp.bed -f 1 -c | awk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){p[$4]=$7}else{n=$1","$2","$3","$4","$6;print $1,$2,$3,$4,int(1000*$5/(p[n]*2+$5))/1000,$6,$7,$8,$9,$5,p[n],substr($10,2),substr($11,2),$12,$13,$16,$17}}' - ${PREFIX}.insertion.filtered.bed > ${PREFIX}.insertion.bed

# Get TSD, remove redundant insertions and recalculate soma rate
$echo 2 "get TSD, remove redundant insertions and recalculate somatic insertion rate"
awk 'BEGIN{FS=OFS="\t"} {if($14!="unknown"){split($14,a,":");split(a[2],b,"-");print a[1],b[1],b[2],NR,0,"."}}' ${PREFIX}.insertion.bed | bedtools getfasta -fi ${GENOME} -bed - -fo ${PREFIX}.tmp -name
awk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){if(NR%2==1){i=substr($1,2)}else{seq[i]=$0}}else{if(seq[FNR]){$14=toupper(seq[FNR])};print $0}}' ${PREFIX}.tmp ${PREFIX}.insertion.bed  | awk 'BEGIN{FS=OFS="\t"} {if($2==$3){$3=$3+1};print $0}' | ${BINDIR}/removeRedundantIns.sh - | awk 'BEGIN{FS=OFS="\t";print "#Chr\tStart\tEnd\tTransposon:Start:End:Strand\tFrequency\tStrand\tType\tSupportReads\tUnspportReads\t5primeSupportReads\t3primeSupportReads\tTSD\tConfidenceForSomaticInsertion\t5splicSiteSupportReads\t3spiceSiteSupportReads"} {print $0}'> ${PREFIX}.t && mv ${PREFIX}.t ${PREFIX}.insertion.bed 
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4";"$5";"$7,0,$6}' ${PREFIX}.insertion.bed | sort -k1,1 -k2,2n > ${PREFIX}.t && bedToBigBed ${PREFIX}.t ${PREFIX}.tmp.chr.size ${PREFIX}.insertion.bb && rm ${PREFIX}.t

# Fix 2p insertions to full length
awk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){te[$1]=$2}else{if($7=="2p"){split($4,a,",");$4="";for(i=1;i<=length(a);i++){split(a[i],b,":");$4=$4","b[1]":1:"te[b[1]]":"b[4]};$4=substr($4,2)};print $0}}' ${PREFIX}.tmp.te.size ${PREFIX}.insertion.bed > ${PREFIX}.t && mv ${PREFIX}.t ${PREFIX}.insertion.bed

# Clean tmp files
$echo 2 "clean tmp files"
${BINDIR}/cleanTempFiles.sh ${PREFIX} ${CLEAN}
$echo 5 "Done, Congras!!!🍺🍺🍺"
