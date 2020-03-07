#!/bin/bash
# TEMP2 (Transposable Element Movement present in a Population)
# 2019-04-04
# Tianxiong Yu(yutianxiong@gmail.com)
# Zhiping Weng Lab
# The Schoold of Life Sciences and Technology, Tongji University
# Programs in Bioinformatics and Integrative Biology, University of Massachusetts Medical School

#usage function
usage() {
echo -en "\e[1;36m"
cat <<EOF

usage: $0 -i input_file.sorted.bam -o output_directory -r transposon_database.fa -t annotated_TEs.bed -m MISMATCH -f fragment_size -S transposon.size -c CPUs -h -F

TEMP2 is a software package for detecting transposable elements (TEs) 
insertions and excisions from pooled high-throughput sequencing data. 
Please send questions, suggestions and bug reports to:
yutianxiong@gmail.com

Options:
	-l     Read1.fastq
	-r     Read2.fastq
        -i     Input file in bam format with full path. Please sort and index the file before calling this program. 
               Sorting and indexing can be done by 'samtools sort' and 'samtools index'
        -o     Path to output directory. Default is current directory
        -R     Transposon sequence database in fasta format with full path
        -t     Annotated TEs in BED6 format with full path. Detected insertions that overlap with annoated TEs will be filtered. 
        -u     TE families annotations. If supplied detected insertions overlap with annotated TE of the same family will be filtered. Only use with -t.
	-m     Percentage of mismatch (AKA divergence; 0-40) allowed when mapping to TE concensus sequences. Default is 5
        -f     An integer specifying the length of the fragments (inserts) of the library. Default is 500
        -c     An integer specifying the number of CPUs used. Default is 8
        -h     Show help message

EOF
echo -en "\e[0m"
}

# taking options
while getopts "hl:r:i:c:f:m:o:R:s:t:u:F" OPTION
do
        case $OPTION in
                h)	usage && exit 1;;
                i)	BAM=`readlink -f $OPTARG`;;
                r)	LEFT=`readlink -f $OPTARG`;;
                l)	RIGHT=`readlink -f $OPTARG`;;
	        f)	INSERT=$OPTARG;;
	        m)	DIV=$OPTARG;;
                o)	OUTDIR=$OPTARG;;
                c)	CPU=$OPTARG;;
	        R)	TESEQ=$OPTARG;;
	        t)	ANNO=$OPTARG;;
                u)	FAMI=$OPTARG;;
                ?)	usage && exit 1;;
        esac
done

if [[ -z $BAM ]] || [[ -z $TESEQ ]]
then
        usage && exit 1
fi
BINDIR=`dirname $0`
[ ! -z "${CPU##*[!0-9]*}" ] || CPU=8
[ ! -z "${INSERT##*[!0-9]*}" ] || INSERT=500
[ ! -z "${DIV##*[!0-9]*}" ] || DIV=5
[ ! -z $OUTDIR ]  || OUTDIR=$PWD
[ ! -z $FILTER ]  || FILTER=0
mkdir -p "${OUTDIR}" || echo -e "\e[1;31mWarning: Cannot create directory ${OUTDIR}. Using the direcory of input fastq file\e[0m"
cd ${OUTDIR} || echo -e "\e[1;31mError: Cannot access directory ${OUTDIR}... Exiting...\e[0m" || exit 1
touch ${OUTDIR}/.writting_permission && rm -rf ${OUTDIR}/.writting_permission || echo -e "\e[1;31mError: Cannot write in directory ${OUTDIR}... Exiting...\e[0m" || exit 1

function checkExist {
        echo -ne "\e[1;32m\"${1}\" is using: \e[0m" && which "$1"
        [[ $? != 0 ]] && echo -e "\e[1;36mError: cannot find software/function ${1}! Please make sure that you have installed the pipeline correctly.\nExiting...\e[0m" && exit 1
}
echo -e "\e[1;35mTesting required softwares/scripts:\e[0m"
checkExist "echo"
checkExist "rm"
checkExist "mkdir"
checkExist "date"
checkExist "mv"
checkExist "sort"
checkExist "touch"
checkExist "awk"
checkExist "grep"
checkExist "bwa"
checkExist "samtools"
checkExist "bedtools"
echo -e "\e[1;35mDone with testing required softwares/scripts, starting pipeline...\e[0m"

cp $TESEQ ./
name=`basename $BAM`
te=`basename $TESEQ`
i=${name/.sorted.bam/}
echo $name
echo $i
if [[ ! -s $name ]]
then
    cp $BAM ./
fi
if [[ ! -s $name.bai ]]
then cp $BAM.bai ./
fi

# Get concordant-uniq-split mapped reads
samtools view -H $name > temp.header
if [ -f $i.pair.uniq.split.fastq ] && [ -f $i.pair.uniq.split.bed ];then
	$BINDIR/echo0 2 "concordant-uniq-split reads exists, skip"
else
	$BINDIR/echo0 2 "get concordant-uniq-split reads"
	samtools view -@ ${CPU} -h -f 0X2 -F 0X800 $name | awk 'BEGIN{FS=OFS="\t"} {for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(as-xs>as*0.8){print $0}}' - | grep "SA:Z:" | awk '(and($2,16)==0 && $6~/^[0-9]+S/) || (and($2,16)==16 && $6~/S$/)' - | cat temp.header - | samtools view -@ ${CPU} -bhS - | samtools sort -@ ${CPU} -o $i.pair.uniq.split.bam -
	samtools fastq -@ ${CPU} $i.pair.uniq.split.bam | awk '{if(NR%4==1){print substr($1,1,length($1)-2)"_"substr($1,length($1))}else{print $0}}' - > $i.pair.uniq.split.fastq 
	bedtools bamtobed -i $i.pair.uniq.split.bam | awk 'BEGIN{FS=OFS="\t"} {$4=substr($4,1,length($4)-2)"_"substr($4,length($4));print $0}' - > $i.pair.uniq.split.bed 
fi

# check fragment length
head -1000000 $i.pair.uniq.sam | awk '{if(and($2,128)==0){print $9<0?-$9:$9}}' > tmp
fragL=(`Rscript ${BINDIR}/checkFragL.R tmp`)
echo -e "95_quantile\t${fragL[1]}\nstandard_deviation\t${fragL[3]}\naverage\t${fragL[5]}" > $i.fragL
[ "${INSERT}" == "N" ] && INSERT=${fragL[1]} && ${BINDIR}/echo0 1 "insert size set to 95 quantile: ${INSERT}"
[ ${fragL[3]} -gt 100 ] && ${BINDIR}/echo0 4 "WARNING: standard deviation of insert size is higher than 100 (${fragL[3]})"

# Get the mate seq of the uniq-unpaired reads
if [ -f $i.unpair.uniq.1.fastq ] && [ -f $i.unpair.uniq.2.fastq ] && [ -f $i.unpair.uniq.bed ];then
	$BINDIR/echo0 2 "mate seq of the uniq-unpaired reads exists, skip"
else
	$BINDIR/echo0 2 "get mate seq of the uniq-unpaired"
	$BINDIR/pickUniqPairFastq.sh $name $i $CPU
fi

# Map to transposons
$BINDIR/echo0 2 "map unpaired uniqMappers to transposons"
bwa index -a is $te
${BINDIR}/faToChromSize $te > temp.te.size
bwa mem -T 0 -Y -t ${CPU} $te $i.pair.uniq.split.fastq > ${i}.pair.uniq.split.transposon.sam 2>/dev/null # split reads
samtools view -@ ${CPU} -hSF 0X2 -F 0X800 ${i}.pair.uniq.split.transposon.sam > t && mv t ${i}.pair.uniq.split.transposon.sam 
bwa mem -T 0 -Y -t ${CPU} $te $i.unpair.uniq.1.fastq $i.unpair.uniq.2.fastq > ${i}.unpair.uniq.transposon.sam 2>/dev/null # unpaired reads
samtools view -@ ${CPU} -hSF 0x2 -F 0X800 $i.unpair.uniq.transposon.sam > t && mv t $i.unpair.uniq.transposon.sam
READ_LENGTH=`awk '{if($6~/M$/ && $6!~/S|I|D/){print substr($6,1,length($6)-1);exit}}' $i.unpair.uniq.transposon.sam`

# Merge fragments in genome and transposon; also remove PCR duplicates
$BINDIR/echo0 2 "merge fragments in genome and transposon"
bash $BINDIR/mergeGenomeTESpit.sh $i.pair.uniq.split.transposon.sam $i.pair.uniq.split.bed ${DIV} temp.te.size | awk '!a[$1"."$2"."$3"."$5]++' > $i.pair.uniq.split.transposon.bed
bash $BINDIR/mergeGenomeTEUnpair.sh $i.unpair.uniq.transposon.sam $i.unpair.uniq.bed ${DIV} temp.te.size | awk '!a[$1"."$2"."$3"."$5]++' > $i.unpair.uniq.transposon.bed

# Fix multipleMappers due to LTRs in both sides of transposons
${BINDIR}/fixLTR.sh $i.pair.uniq.split.transposon.bed > $i.pair.uniq.split.transposon.fixLTR.bed
${BINDIR}/fixLTR.sh $i.unpair.uniq.transposon.bed > $i.unpair.uniq.transposon.fixLTR.bed

# Merge supporting reads within insert size - read length and in the same direction
mkdir temp
awk 'BEGIN{FS=OFS="\t"} {split($5,a,";");for(i in a){split(a[i],b,",");if(b[4]==$6){ts="-"}else{ts="+"};print $1,$2,$3,1/length(a),a[i],$6>>"temp/"b[1]"."ts".bed"}}' ${i}.unpair.uniq.transposon.fixLTR.bed ${i}.pair.uniq.split.transposon.fixLTR.bed
MAXGAP=`expr ${INSERT} - ${READ_LENGTH}`
[ $MAXGAP -lt 0 ] && MAXGAP=0
for i in temp/*.bed; do sort -k1,1 -k2,2n $i | bedtools merge -i - -d 375 -s -c 4,5,6,2,3 -o collapse,collapse,first,collapse,collapse -delim "|" > ${i/.bed/.merged.bed}; done
for i in temp/*merged.bed; do
	tel=`awk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){a[$1]=$2}else{if(FNR==1){split(FILENAME,fn,".merged.bed");split(fn[1],b,"/");print a[substr(b[2],1,length(b[2])-2)]}}}' temp.te.size $i`
	${BINDIR}/processMergedBed.py $i $tel $INSERT > ${i/merged.bed/processed.bed}
done

# Merge supporting reads within 2 * insert size - read length and in direction + then -
for i in temp/*.processed.bed; do ${BINDIR}/mergeProcessedBed.sh $i $INSERT $READ_LENGTH ${i%.processed.bed}; done
cat temp/*.TPregion.bed | sort -k1,1 -k2,2n | bedtools merge -i - -s -c 4,6 -o count,first > ${i}.TPregion.bed
cat temp/*.final.bed > ${i}.final.bed

# Estimate all transposon read distribution in each transposon
bwa mem -T 0 -Y -t ${CPU} $te ${LEFT} ${RIGHT} > ${i}.transposon.sam 2>${i}.transposon.bwamem.log
#samtools view -@ ${CPU} -hS -F 0X4 -F 0X800 ${i}.transposoddn.sam | awk -v div=$DIV 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){tel[$1]=$2}else{if($1~/^@/){print $0}else{for(i=12;i<=100;i++){if($i~/NM:i:/){mm=substr($i,6);break}};split($6,cigar1,"M|D|I|S");split($6,cigar2,"[0-9]+");aln=0;for(i=2;i<=length(cigar2);i++){if(cigar2[i]=="M" || cigar2[i]=="D"){aln+=cigar1[i-1]}};if(int(mm)<=int(aln*div/100) && ($4<50 || $4+aln+50>tel[$3] || aln>35)){print $0}}}}' temp.te.size - | samtools view -@ ${CPU} -bhS - > ${i}.transposon.bam
samtools view -@ ${CPU} -hS -F 0X4 -F 0X800 ${i}.transposon.sam > t && mv t ${i}.transposon.sam 
${BINDIR}/samToBdg.sh ${i}.transposon.sam ${DIV} temp.te.size > ${i}.transposon.bdg 

# 





