#!/usr/bin/env Rscript

### arguments ----
Args=commandArgs()
if(length(Args)<6){
	stop(paste(substr(Args[4],8,100),"in.fragL",sep=" "))
}

### main program ----
fl=read.table(Args[6],header=F,row.names=NULL)
print(as.integer(sd(fl[,1])))
print(as.integer(quantile(fl[,1],.95)))
print(as.integer(mean(fl[,1])))
