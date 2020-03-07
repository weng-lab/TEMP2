#!/bin/bash

if [ $# -lt 3 ];then
	echo -e "$0 in.processed.bed insert_size read_length out.prefix"
	exit 1
fi

awk -v ins=$2 -v rl=$3 'BEGIN{FS=OFS="\t"}
{
	if($3-$2<30){next};
	if(ins*2>rl){ext=ins-rl}else{ext=0};
	if($6=="+"){
		if(p){
			split(p,a,"\t");
			print a[1],a[3],a[3]+ext,a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[4],0,a[5],0,a[11],a[13],0
		};
		p=$0
	}else{
		if(p){
			split(p,a,"\t");
			if(($2-a[3])<2*ext && $1==a[1]){
				if(a[8]=="-"){ts=a[9];te=$10}else{ts=$9;te=a[10]};
				print $1,a[3],$2,$4+a[4],$5+a[5],a[6]""$6,$7,a[8]""$8,ts,te,a[4],$4,a[5],$5,a[11]"|"$11,a[13],$12
			}else{
				print a[1],a[3],a[3]+ext,a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[4],0,a[5],0,a[11],a[13],0;
				print $1,$2-ext,$2,$4,$5,$6,$7,$8,$9,$10,0,$4,0,$5,$11,0,$12
			};
			p=0
		}else{
			print $1,$2-ext,$2,$4,$5,$6,$7,$8,$9,$10,0,$4,0,$5,$11,0,$12
		}
	}
}
END{
	if(p){
		split(p,a,"\t");
		print a[1],a[3],a[3]+ext,a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[4],0,a[5],0,a[11],a[13],0
	}
}
' $1 > $4".final.bed"
