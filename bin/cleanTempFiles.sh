#!/bin/bash

if [ $# -lt 1 ];then
	echo -e "$0 prefix d(delete)|p(preserve)"
	echo -e "clean all intermediate files produced by TEMP2"
	exit 1
fi 

if [ "$2" == "d" -o "$2" == "delete" ];then
	rm -rf $1.pair* $1.unpair* $1.transposon* $1.tmp* $1.transposonMapping $1.supportReads \
		$1.soma.rate.bed $1.parafile* $1.TPregion.bed \
		$1.final.bed $1.fragL $1.insertion.raw.bed $1.singleton* $1.1p1* $1.2p* \
		$1.insertion.filtered.bed $1.*.log $1.removed.*
else
	[ ! -d tmpTEMP2 ] && mkdir tmpTEMP2
	mv $1.pair* $1.unpair* $1.transposon* $1.tmp* $1.transposonMapping $1.supportReads \
		$1.soma.rate.bed $1.parafile* $1.TPregion.bed \
		$1.final.bed $1.fragL $1.insertion.raw.bed $1.singleton* $1.1p1* $1.2p* \
		$1.insertion.filtered.bed $1.*.log $1.removed.* tmpTEMP2 >/dev/null 2>&1
fi
