#!/bin/bash

if [ $# -lt 2 ];then
	echo -e "$0 in.bam|in.sam out.prefix cpu insert_size uniq_ratio"
	echo -e "extract uniquely mapped unpaired reads from bam or sam file"
	exit 1
fi

samtools view -@ $3 -h -F 0X2 $1 > $2.tmp
awk -v ins=$4 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){if(and($2,2048)==0){a=$9<0?-$9:$9;if(a>10 && a<ins+200){k[$1]=1}}}else{if($1~/^@/){print $0}else{if(!k[$1]){print $0}}}}' $2.tmp $2.tmp > $2.unpair.sam
awk -v uniq_ratio=$5 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(xs/1<as*uniq_ratio){a[$1]=1}}else{if(a[$1] || $1~/^@/){print $0}}}' $2.unpair.sam $2.unpair.sam | samtools view -@ $3 -bhS -  | samtools sort -@ $3 -n - -o $2.unpair.uniq.sortByName.bam
samtools fastq -@ $3 -1 $2.unpair.uniq.1.fastq -2 $2.unpair.uniq.2.fastq $2.unpair.uniq.sortByName.bam && rm $2.unpair.uniq.sortByName.bam 
awk -v uniq_ratio=$5 'BEGIN{FS=OFS="\t"} {if($1~/^@/){print $0}else{for(i=12;i<=100;i++){if($i~/AS:i:/){as=substr($i,6)};if($i~/XS:i:/){xs=substr($i,6);break}};if(xs/1<as*uniq_ratio){print $6"____"$0}}}' $2.unpair.sam | samtools view -@ $3 -bhS - | bedtools bamtobed -i - -tag NM | awk '
BEGIN{FS=OFS="\t"}
{
	split($4,a,"____");
	$4=a[2];
	cigar=a[1];
	split(cigar,ml,"M|I|D|S");
	gsub("[0-9]","",cigar);
	nm=$5;
	if(cigar~/S$/ && cigar!~/^S/ && $6=="+"){$5=1}
	else if(cigar~/^S/ && cigar!~/S$/ && $6=="-"){$5=1}
	else{$5=0};
	print $0,nm
}' > $2.unpair.uniq.bed
