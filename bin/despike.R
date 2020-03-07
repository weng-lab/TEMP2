#!/usr/bin/env Rscript

### arguments ----
Args=commandArgs()
if(length(Args)<7){
	stop(paste(substr(Args[4],8,100),"in.cov out.spike.bed",sep=" "))
}

### main program ----
cov = read.table(Args[6], header=F, row.names=NULL)
spikes = c()
fun_spikeRange=function(tdf, center, cn){
  sdc = 0.5
  start = end1 = center # get spike start
  while(start>1){
    ws1 = mean(tdf[(max(1,start-10):(start-1)),cn])
    ws2 = mean(tdf[start:end1,cn])
    start = max(1, start-10)
    end1 = min(end1, start+10)
    if((ws2+sdc)/(ws1+sdc)>10){break}
  }
  start1 = end = center # get spike end
  while(end<nrow(tdf)){
    ws1 = mean(tdf[(end+1):min((end+10),nrow(tdf)),cn])
    ws2 = mean(tdf[start1:end,cn])
    start1 = max(start1, end-10)
    end = min(end+10, nrow(tdf))
    if((ws2+sdc)/(ws1+sdc)>10){break}
  }
  return(c(start, end))
}
for(te in levels(cov[,1])){
  tdf = cov[which(cov[,1]==te), 2:4]
  ms = mean(tdf[,2])
  ma = mean(tdf[,3])
  rn = 1:nrow(tdf)
  while(1){
    tdfC = tdf[rn,]
    max_sig = max(tdfC[,2])
    if(max_sig<=ms){break}
    center = tdfC[which(tdfC[,2]==max_sig)[1],1]
    spikeRange = fun_spikeRange(tdf, center, 2)
    if((spikeRange[2]-spikeRange[1])<=100){
      spikes = rbind(spikes, c(te, spikeRange[1], spikeRange[2], 0, 0, "+"))
    }
    rn = setdiff(rn, spikeRange[1]:spikeRange[2])
  }
  rn = 1:nrow(tdf)
  while(1){
    tdfC = tdf[rn,]
    max_sig = max(tdfC[,3])
    if(max_sig<=ms){break}
    center = tdfC[which(tdfC[,3]==max_sig)[1],1]
    spikeRange = fun_spikeRange(tdf, center, 3)
    if((spikeRange[2]-spikeRange[1])<=100){
      spikes = rbind(spikes, c(te, spikeRange[1], spikeRange[2], 0, 0, "-"))
    }
    rn = setdiff(rn, spikeRange[1]:spikeRange[2])
  }
}
if(length(spikes)>0){
	write.table(spikes, Args[7], row.names=F, col.names=F, sep="\t", quote=F)
}
