#!/bin/bash

if [ $# -lt 2 ];then
	echo -e "$0 in.sam out.prefix"
	exit 1
fi

[ -f $2.fastq ] && rm $2.fastq
[ -f $2.bed ] && rm $2.bed
gawk -v p=$2 'BEGIN{FS=OFS="\t";nul["A"]="T";nul["T"]="A";nul["C"]="G";nul["G"]="C";nul["N"]="N"}
{
	split($6,cigar1,"M|D|I|S");split($6,cigar2,"[0-9]+");aln=0;
	for(i=2;i<=length(cigar2);i++){
		if(cigar2[i]=="M" || cigar2[i]=="D"){aln+=cigar1[i-1]}
	};
	if(and($2,64)==64){rnum=1}else{rnum=2};
	if(and($2,16)==0){
		if(cigar2[2]=="S" && cigar1[1]>=20){
			seq=substr($10,1,cigar1[1]);qua=substr($11,1,cigar1[1]);
			print "@"$1"/"rnum"\n"seq"\n+"$1"/"rnum"\n"qua>>p".fastq";
			print $3,$4-1,$4-1+aln,$1"/"rnum,$5,"+">>p".bed"
		}
	}else{
		if(cigar2[length(cigar2)]=="S" && cigar1[length(cigar1)-1]>=20){
			t=substr($10,length($10)-cigar1[length(cigar1)-1]+1);seq="";
			for(i=length(t);i>=1;i--){seq=seq""nul[substr(t,i,1)]};
			qua=substr($11,length($11)-cigar1[length(cigar1)-1]+1);
			print "@"$1"/"rnum"\n"seq"\n+"$1"/"rnum"\n"qua>>p".fastq";
			print $3,$4-1,$4-1+aln,$1"/"rnum,$5,"-">>p".bed"
		}
	}
}
' $1
