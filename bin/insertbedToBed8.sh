#!/bin/bash

if [ $# -lt 1 ];then
	echo -e "$0 in.insertion.bed"
	echo -e "transform insertion.bed file to a standard bed format with penetrance in the 5th column, type in the 7th column and confidence in the 8th column"
	exit 1
fi

awk 'BEGIN{FS=OFS="\t"} {if($1!~/^#/){split($4,a,",");for(i in a){split(a[i],b,":");print $1,$2,$3,b[1],$5/length(a),b[4],$7,$13}}}' $1
