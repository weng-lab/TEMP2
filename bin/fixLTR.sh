#!/bin/bash

if [ $# -lt 1 ];then
	echo -e "$0 in.fragment.bed"
	exit 1
fi

gawk 'BEGIN{FS=OFS="\t"}
{
	split($5,tes,";");delete o;
	if(length(tes)==1){print $0}
	else{
		for(i in tes){
			split(tes[i],t,",");
			if(!o[t[1]","t[4]]){o[t[1]","t[4]]=tes[i]}
			else{
				split(o[t[1]","t[4]],k,",");
				if(t[4]=="+"){
					if(t[2]>k[2]){o[t[1]","t[4]]=tes[i]}
				}else if(t[4]=="-"){
					if(t[2]<k[2]){o[t[1]","t[4]]=tes[i]}
				}
			}
		}
		oo="";for(i in o){oo=oo";"o[i]}
		print $1,$2,$3,$4,substr(oo,2),$6
	}
}
' $1
