#!/bin/bash

if [ $# -lt 4 ];then
	echo -e "usage:"
	echo -e "$0 in.1.fq in.2.fq chimeric_ratio out.prefix"
	exit 1
fi


awk -v cr=$3 -v p=$4 'BEGIN{FS=OFS="\t";srand()}
{
	if(ARGIND==1){
		i=int((FNR-1)/4+1);
		if(i%100000==0){printf i/1000000"M read1 processed\r"};
		if(FNR%4==1){split($0,a,"/");name1[i]=a[1]}
		else if(FNR%4==2){seq1[i]=$0}
		else if(FNR%4==0){qlt1[i]=$0}
	}
	else{
		i=int((FNR-1)/4+1);
		if(i%100000==0){if(FNR==400000-3){print ""};printf i/1000000"M read2 processed\r"};
		if(FNR%4==1){split($0,a,"/");name2[i]=a[1]}
		else if(FNR%4==2){seq2[i]=$0}
		else if(FNR%4==0){qlt2[i]=$0}
	}
}
END{
	for(i=1;i<=length(seq1);i++){
		if(i%100000==0){if(i==100000){print ""};printf i/1000000"M read pair written\r"};
		if(rand()<=cr){
			print name1[i]"\n"seq1[i]"\n+\n"qlt1[i] >> p".1.fastq"
			i++;
			print name2[i-1]"\n"seq2[i]"\n+\n"qlt2[i] >> p".2.fastq"
		}else{
			print name1[i]"\n"seq1[i]"\n+\n"qlt1[i] >> p".1.fastq"
			print name2[i]"\n"seq2[i]"\n+\n"qlt2[i] >> p".2.fastq"
		}	
	}
	print ""
}' $1 $2
