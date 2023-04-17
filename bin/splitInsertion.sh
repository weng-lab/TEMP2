#!/bin/bash

if [ $# -lt 1 ];then
	echo -e "$0 in.insertion.bed te.size"
	echo -e "This script is to separate transposons inserted into the same position, and make a easy read result."
	exit 1
fi

# NW_018343952.1	1864	1880	Ko.L1.3:12119:12174:+	0.0625	+	1p1	2	15	1	1	TATATATATATATATA	100	1	1
gawk 'BEGIN{FS=OFS="\t"} {if(ARGIND==1){te[$1]=$2}else{if(FNR>1){split($4,a,",");l=length(a);for(i in a){split(a[i],b,":");print $1,$2,$3,b[1],$5/l,b[4],$7,$8/l,$12,b[2],b[3],te[b[1]]}}}}' $2 $1
