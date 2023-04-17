#!/bin/bash

PATH_PRO=$(dirname `readlink -f $0`)

if [ $# -lt 2 ];then
	bash ${PATH_PRO}/echo0 5 "$0 in.sam in.bed divergence te.size read_length"
	exit 1
fi

#samtools view -S -F 0X4 $1 | gawk -v mm=$3 'BEGIN{FS=OFS="\t"}
gawk -v div=$3 -v rl=$5 'BEGIN{FS=OFS="\t"}
{if(ARGIND==1){
	tel[$1]=$2
}else if(ARGIND==2){
if($1~/^@/ || and($2,4)==4){next};
	for(i=12;i<=100;i++){
		if($i~/^NM:i:/){mm=int(substr($i,6))};
		XA="";if($i~/^XA:Z:/){XA=substr($i,6);break};
	};
	split($6,cigar1,"M|D|I|S");split($6,cigar2,"[0-9]+");aln=0;
	for(i=2;i<=length(cigar2);i++){
		if(cigar2[i]=="M" || cigar2[i]=="D"){aln+=cigar1[i-1]}
	};
	if(and($2,16)==16){strand="-"}else{strand="+"};
	if($4+aln>tel[$3]){aln=tel[$3]-$4+1};
	if(mm<=int(aln*div/100) && ($6!~/S[0-9a-zA-Z]+S/ || aln>=50 || aln>=0.75*rl) && aln>=25){
		record[$1]=$3","$4","$4+aln-1","strand",1";
	}
	split(XA,t,";");
	for(i=1;i<=length(t)-1;i++){
		txa=t[i];split(txa,trec,",");
		split(trec[3],cigar1,"M|D|I|S");split(trec[3],cigar2,"[0-9]+");aln=0;
		for(j=2;j<=length(cigar2);j++){
			if(cigar2[j]=="M" || cigar2[j]=="D"){aln+=cigar1[j-1]}
		};
		if(substr(trec[2],2)+aln>tel[trec[1]]){ttl=tel[trec[1]]-substr(trec[2],2)+1;if(ttl<25){continue}else{aln=ttl}};
		if(int(trec[4])<int(aln*div/100) && (trec[3]!~/S[0-9a-zA-Z]+S/ || aln>=50 || aln>=0.75*rl)){
			record[$1]=record[$1]";"trec[1]","substr(trec[2],2)","substr(trec[2],2)+aln-1","substr(trec[2],1,1)",1"
		}
	};
}else{
	bed[$4,1]=$1;bed[$4,2]=$2;bed[$4,3]=$3;bed[$4,4]=$4;bed[$4,5]=$6;
}}
END{
for(i in record){
	if(record[i] && bed[i,1]){
		print bed[i,1],bed[i,2],bed[i,3],bed[i,4],record[i],bed[i,5]
	}
}}' $4 $1 $2
