#!/bin/bash

if [ $# -lt 2 ];then
	echo -e "usage:"
	echo -e "$0 genome.fa in.insertion.bed insertion.tab out.prefix"
	echo -e "the script will randomly insert all the sequences in insertion.tab(1st column name, 2nd column insert times, 3rd column sequence) to genome.fa\n"
	exit 1
fi
#N=`awk '{s+=$2} END{print s}' $3`
#if [ $# -lt 4 ];then
#	bedtools random -l 1 -n $N -g $2 | awk 'BEGIN{FS=OFS="\t";tn=0} {if(ARGIND==1){a[FNR]=$0}else{for(i=1;i<=$2;i++){tn++;print a[tn],$1}}}' - $3 | sort -k1,1 -k2,2n | intersectBed -a - -b $5 -v > $4.bed
#else
#	bedtools random -l 1 -n $N -g $2 | awk 'BEGIN{FS=OFS="\t";tn=0} {if(ARGIND==1){a[FNR]=$0}else{for(i=1;i<=$2;i++){tn++;print a[tn],$1}}}' - $3 | sort -k1,1 -k2,2n > $4.bed
#fi
awk 'BEGIN{FS=OFS="\t";nc["A"]="T";nc["T"]="A";nc["C"]="G";nc["G"]="C"} 
{
	if(ARGIND==1){
		ts[$1]=$3
	}
	else if(ARGIND==2){
		chrN[$1]++
		if($6=="+"){seq=ts[$4]}
		else{seq="";for(i=length(ts[$4]);i>=1;i--){seq=seq""nc[substr(ts[$4],i,1)]}};
		insR[$1][chrN[$1]]=int($2/80)+1;
		insC[$1][chrN[$1]]=int($2%80);
		insSeq[$1][chrN[$1]]=seq;
	}
	else{
		if($1~/^>/){
			if(!chr){print $0}else{print "\n"$0};
			chr=substr($1,2);ti=1;rnum=0;
		}else{
			rnum++;
			if(rnum==insR[chr][ti]){
				printf substr($1,1,insC[chr][ti]);
				printf insSeq[chr][ti];
				printf substr($1,insC[chr][ti]+1);
				ti++;
			}else{
				printf $1
			}
		}

	}
}
' $3 $2 $1 > $4.fa
