#!/bin/bash

BINDIR=$(dirname `readlink -f $0`)

if [ $# -lt 2 ];then
	bash ${BINDIR}/echo0 5 "$0 in.sam divergence te.size out.prefix read_length"
	bash ${BINDIR}/echo0 5 "filter sam record with more than (divergence*map_length) mismatch and potential tRNA bias; then make bedGraph"
	exit 1
fi

gawk -v div=$2 -v rl=$5 'BEGIN{FS=OFS="\t"}
{if(ARGIND==1){
	tel[$1]=$2
}else{
	if($1~/^@/){next};
	for(i=12;i<=100;i++){
		if($i~/^NM:i:/){mm=int(substr($i,6))};
		XA="";if($i~/^XA:Z:/){XA=substr($i,6);break};
	};
	split(XA,t,";");
	split($6,cigar1,"M|D|I|S");split($6,cigar2,"[0-9]+");aln=0;
	for(i=2;i<=length(cigar2);i++){
		if(cigar2[i]=="M" || cigar2[i]=="D"){aln+=cigar1[i-1]}
	};
	if(and($2,16)==0){strand="+"}else{strand="-"};
	if(mm<=int(aln*div/100) && ($6!~/S[0-9a-zA-Z]+S/ || aln>=50 || aln>=0.75*rl) && aln>=25){
		record=$3","$4","$4+aln-1","strand;
	}else{record=""};
	for(i=1;i<=length(t)-1;i++){
		txa=t[i];split(txa,trec,",");
		split(trec[3],cigar1,"M|D|I|S");split(trec[3],cigar2,"[0-9]+");aln=0;
		for(j=2;j<=length(cigar2);j++){
			if(cigar2[j]=="M" || cigar2[j]=="D"){aln+=cigar1[j-1]}
		};
		if(substr(trec[2],2)+aln>tel[trec[1]]){ttl=tel[trec[1]]-substr(trec[2],2)+1;if(ttl<25){continue}else{aln=ttl}};
		if(int(trec[4])<int(aln*div/100) && (trec[3]!~/S[0-9a-zA-Z]+S/ || aln>=50 || aln>=0.75*rl)){
			record=record";"trec[1]","substr(trec[2],2)","substr(trec[2],2)+aln-1","substr(trec[2],1,1)
		}
	};
	if(record){
		print 0,0,0,0,record,0
	}
}}' $3 $1 | gawk 'BEGIN{FS=OFS="\t"} {if($5~/^;/){$5=substr($5,2)};print $0}' - | ${BINDIR}/fixLTR.sh - > $1.tmp
gawk 'BEGIN{FS=OFS="\t"} {split($5,a,";");for(i in a){split(a[i],b,",");if(b[4]=="+"){print b[1],b[2],b[3],1/length(a)}}}' $1.tmp | sort -k1,1 -k2,2n > $1.tmp.sense && ${BINDIR}/mergeOverlappedBdg $1.tmp.sense > $4.sense.bdg
gawk 'BEGIN{FS=OFS="\t"} {split($5,a,";");for(i in a){split(a[i],b,",");if(b[4]=="-"){print b[1],b[2],b[3],1/length(a)}}}' $1.tmp | sort -k1,1 -k2,2n > $1.tmp.anti && ${BINDIR}/mergeOverlappedBdg $1.tmp.anti > $4.anti.bdg
gawk 'BEGIN{FS=OFS="\t"} {split($5,a,";");for(i in a){split(a[i],b,",");print b[1],b[2],b[3],1/length(a),0,b[4]}}' $1.tmp | sort -k1,1 -k2,2n > $4.bed
rm $1.tmp $1.tmp.sense $1.tmp.anti
