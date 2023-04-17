#!/bin/bash

if [ $# -lt 1 ];then
	echo -e "$0 in.fq"
	echo -e "if the nucleotide is not A, T, C, G, will automatically transform to N"
	exit 1
fi

gawk '
BEGIN{
	FS=OFS="\t";
	a["A"]="T";a["T"]="A";a["C"]="G";a["G"]="C"
}
{
	if(NR%4==2){
		for(i=length($1);i>=1;i--){
			nucl=substr($1,i,1);
			if(a[nucl]){
				printf a[nucl]
			}else{
				printf "N"
			}
		}
		print ""
	}else if(NR%4==0){
		for(i=length($1);i>=1;i--){
			printf substr($1,i,1)
		}
		print ""	
	}else{
		print $0
	}
}
' $1
