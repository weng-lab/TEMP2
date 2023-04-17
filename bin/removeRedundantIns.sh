#!/bin/bash

if  [ $# -lt 1 ];then
	echo -e "usage:"
	echo -e "$0 in.ins.bed(17 columns bed file)"
	echo -e "this script will remove redundant insertions if several insertions are overlapped and the support read is two folder lesser than the max support reds in these insertions"
	exit 1
fi

sort -k1,1 -k2,2n $1 | bedtools merge -d 50 -i - -c 2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17 -delim "|" -o collapse | gawk 'BEGIN{
FS=OFS="\t";
"#Chr\tStart\tEnd\tTransposon:Start:End:Strand\tFrequency\tStrand\tType\tSupportReads\tUnspportReads\t5primeSupportReads\t3primeSupportReads\tTSD\tConfidenceForSomaticInsertion\t5splicSiteSupportReads\t3spiceSiteSupportReads"
} 
{
	if($4!~/|/){
		print $1,$4,$5,$6":"$10":"$11":"$8,$7,$8,$9,$12,$13,$14,$15,$16,$17,$18,$19
	}else{
		split($4,corS,"|");
		split($5,corE,"|");
		split($6,te,"|");
		split($7,freq,"|");
		split($8,strand,"|");
		split($9,type,"|");
		split($10,teS,"|");
		split($11,teE,"|");
		split($12,spR,"|");
		split($13,unspR,"|");
		split($14,lspR,"|");
		split($15,rspR,"|");
		split($16,TSD,"|");
		split($17,cfd,"|");
		split($18,lslR,"|");
		split($19,rslR,"|");
		maxspR=0;
		for(i in spR){maxspR=maxspR<spR[i]?spR[i]:maxspR};
		delete trueIdx;delete falseIdx;
		for(i in spR){
			if(spR[i]>maxspR/2){trueIdx[i]=1}else{falseIdx[i]=1};
			if(spR[i]==maxspR){mainIdx=i};
		};
		fspR=0;flspR=0;frspR=0;flslR=0;frslR=0;
		for(i in spR){
			fspR+=spR[i];
			flspR+=lspR[i];
			frspR+=rspR[i];
			flslR+=lslR[i];
			frslR+=rslR[i];
		};
		printf $1"\t"corS[mainIdx]"\t"corE[mainIdx]"\t";
		final_te="";final_strand=""
		for(i in trueIdx){
			final_te=final_te","te[i]":"teS[i]":"teE[i]":"strand[i];
			if(!final_strand){final_strand=strand[i]}
			else if(final_strand!=strand[i]){final_strand="."}
		};
		if(length(trueIdx)==2){
			split(final_te,ta,",");
			split(ta[2],tb1,":");
			split(ta[3],tb2,":");
			if(tb1[1]==tb2[1] && flspR>0 && frspR>0){
				type[mainIdx]="1p1"
			};
		};
		printf substr(final_te,2)"\t";
		printf fspR/(fspR+2*unspR[mainIdx])"\t";
		printf final_strand"\t"type[mainIdx]"\t"fspR"\t";
		printf unspR[mainIdx]"\t"flspR"\t"frspR"\t";
		printf TSD[mainIdx]"\t"cfd[mainIdx]"\t"flslR"\t"frslR"\n"
	}
}'
