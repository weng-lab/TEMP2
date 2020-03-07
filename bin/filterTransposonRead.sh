#!bin/bash

if [ $# -lt 3 ];then
	echo -e "$0 in.transposonRead.bed out.transposonRead.bed te.size"
	exit 1
fi

awk -v inst=${INSERT} 'BEGIN{FS=OFS="\t"} \
	{if(ARGIND==1){t[$1]=$2}\
	else{\
		rl=$3-$2;split($5,n,";");\
		on="";\
		for(i in n){\
			if(!n[i]){continue};
			split(n[i],m,",");trn=m[1];\
			strand=m[4];\
			pos=m[2];\
			if(strand=="+"){\
				if(int(t[trn]-pos)<int(inst-rl)){\
					on=on";"n[i]\
				}\
			}else{\
				if(int(pos)<int(inst-2*rl)){\
					on=on";"n[i]\
				}\
			}\
		}\
		if(on!=""){\
			$5=substr(on,2);print $0\
		}\
	}}' $3 $1 > $2
