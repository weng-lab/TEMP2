#!/usr/bin/env Rscript

### arguments ----
Args=commandArgs()
if(length(Args)<6){
	stop(paste(substr(Args[4],8,100),"overall.sense.bdg overall.anti.bdg observed.sense.bdg observed.anti.bdg TPregion.bed transposon.size out.pdf",sep=" "))
}

### functions ----
fun_bdgS=function(x){
  s=0;for(i in 1:nrow(x)){s=s+(x[i,3]-x[i,2])*x[i,4]}
  return(s)
}

### main program ----

# read data frame
bdg_all_s=read.table(Args[6],header=F,row.names=NULL)
bdg_all_a=read.table(Args[7],header=F,row.names=NULL)
bdg_singleton_s=read.table(Args[8],header=F,row.names=NULL)
bdg_singleton_a=read.table(Args[9],header=F,row.names=NULL)
bdg_1p1_s=read.table(Args[10],header=F,row.names=NULL)
bdg_1p1_a=read.table(Args[11],header=F,row.names=NULL)
bdg_2p_s=read.table(Args[12],header=F,row.names=NULL)
bdg_2p_a=read.table(Args[13],header=F,row.names=NULL)
peak=read.table(Args[14],header=F,row.names=NULL)
te.size=read.table(Args[15],header=F,row.names=1)
csl=c("#E41A1C","#377EB8") # color

# plot distribution
pdf(Args[16], width=3.5,height=6,useDingbats=F)
par(mfrow=c(3,1),tcl=0.3,bty="n",cex=5/6)
for(tte in row.names(te.size)){ # for each transposon
  # get temp data frame for a specific transposon
  tdf_all_s=bdg_all_s[which(bdg_all_s[,1]==tte),]
  tdf_all_a=bdg_all_a[which(bdg_all_a[,1]==tte),]
  tdf_singleton_s=bdg_singleton_s[which(bdg_singleton_s[,1]==tte),]
  tdf_singleton_a=bdg_singleton_a[which(bdg_singleton_a[,1]==tte),]
  tdf_1p1_s=bdg_1p1_s[which(bdg_1p1_s[,1]==tte),]
  tdf_1p1_a=bdg_1p1_a[which(bdg_1p1_a[,1]==tte),]
  tdf_2p_s=bdg_2p_s[which(bdg_2p_s[,1]==tte),]
  tdf_2p_a=bdg_2p_a[which(bdg_2p_a[,1]==tte),]
  pk_s=peak[which(peak[,1]==tte & peak[,6]=="+"),]
  pk_a=peak[which(peak[,1]==tte & peak[,6]=="-"),]
  # plot 1p1
  par(mar=c(0,4,3,1))
  if(nrow(tdf_1p1_s)>0){ym1=max(tdf_1p1_s[,4])}else{ym1=0}
  if(nrow(tdf_1p1_a)>0){ym2=max(tdf_1p1_a[,4])}else{ym2=0}
  ym=max(ym1,ym2)
  plot(NA,xlim=c(1,te.size[tte,1]),ylim=c(-ym,ym),xlab="",ylab="1p1",
       xaxt="n",main=tte)
  if(nrow(pk_s)>0){for(i in 1:nrow(pk_s)){polygon(c(pk_s[i,2],pk_s[i,3],pk_s[i,3],pk_s[i,2]),c(0,0,ym,ym),col="#9ecae1",border="#9ecae1")}}
  if(nrow(pk_a)>0){for(i in 1:nrow(pk_a)){polygon(c(pk_a[i,2],pk_a[i,3],pk_a[i,3],pk_a[i,2]),c(0,0,-ym,-ym),col="#fdae6b",border="#fdae6b")}}
  xcor=c(1);ycor=c(0)
  if(nrow(tdf_1p1_s)>0){for(i in 1:nrow(tdf_1p1_s)){xcor=c(xcor,tdf_1p1_s[i,2],tdf_1p1_s[i,2],tdf_1p1_s[i,3],tdf_1p1_s[i,3]);ycor=c(ycor,0,tdf_1p1_s[i,4],tdf_1p1_s[i,4],0)}}
  polygon(xcor,ycor,col=csl[2],border=csl[2])
  xcor=c(1);ycor=c(0)
  if(nrow(tdf_1p1_a)>0){for(i in 1:nrow(tdf_1p1_a)){xcor=c(xcor,tdf_1p1_a[i,2],tdf_1p1_a[i,2],tdf_1p1_a[i,3],tdf_1p1_a[i,3]);ycor=c(ycor,0,-tdf_1p1_a[i,4],-tdf_1p1_a[i,4],0)}}
  xcor=c(xcor,xcor[length(xcor)]);ycor=c(ycor,0)
  polygon(xcor,ycor,col=csl[1],border=csl[1])
  # plot 2p
  par(mar=c(1.5,4,1.5,1))
  if(nrow(tdf_2p_s)>0){ym1=max(tdf_2p_s[,4])}else{ym1=0}
  if(nrow(tdf_2p_a)>0){ym2=max(tdf_2p_a[,4])}else{ym2=0}
  ym=max(ym1,ym2)
  plot(NA,xlim=c(1,te.size[tte,1]),ylim=c(-ym,ym),xlab="",ylab="2p",
       xaxt="n")
  if(nrow(pk_s)>0){for(i in 1:nrow(pk_s)){polygon(c(pk_s[i,2],pk_s[i,3],pk_s[i,3],pk_s[i,2]),c(0,0,ym,ym),col="#9ecae1",border="#9ecae1")}}
  if(nrow(pk_a)>0){for(i in 1:nrow(pk_a)){polygon(c(pk_a[i,2],pk_a[i,3],pk_a[i,3],pk_a[i,2]),c(0,0,-ym,-ym),col="#fdae6b",border="#fdae6b")}}
  xcor=c(1);ycor=c(0)
  if(nrow(tdf_2p_s)>0){for(i in 1:nrow(tdf_2p_s)){xcor=c(xcor,tdf_2p_s[i,2],tdf_2p_s[i,2],tdf_2p_s[i,3],tdf_2p_s[i,3]);ycor=c(ycor,0,tdf_2p_s[i,4],tdf_2p_s[i,4],0)}}
  polygon(xcor,ycor,col=csl[2],border=csl[2])
  xcor=c(1);ycor=c(0)
  if(nrow(tdf_2p_a)>0){for(i in 1:nrow(tdf_2p_a)){xcor=c(xcor,tdf_2p_a[i,2],tdf_2p_a[i,2],tdf_2p_a[i,3],tdf_2p_a[i,3]);ycor=c(ycor,0,-tdf_2p_a[i,4],-tdf_2p_a[i,4],0)}}
  xcor=c(xcor,xcor[length(xcor)]);ycor=c(ycor,0)
  polygon(xcor,ycor,col=csl[1],border=csl[1])
  # plot singleton
  par(mar=c(3,4,0,1))
  if(nrow(tdf_singleton_s)>0){ym1=max(tdf_singleton_s[,4])}else{ym1=0}
  if(nrow(tdf_singleton_a)>0){ym2=max(tdf_singleton_a[,4])}else{ym2=0}
  ym=max(ym1,ym2)
  plot(NA,xlim=c(1,te.size[tte,1]),ylim=c(-ym,ym),xlab="",ylab="singleton",
       xaxt="n")
  if(nrow(pk_s)>0){for(i in 1:nrow(pk_s)){polygon(c(pk_s[i,2],pk_s[i,3],pk_s[i,3],pk_s[i,2]),c(0,0,ym,ym),col="#9ecae1",border="#9ecae1")}}
  if(nrow(pk_a)>0){for(i in 1:nrow(pk_a)){polygon(c(pk_a[i,2],pk_a[i,3],pk_a[i,3],pk_a[i,2]),c(0,0,-ym,-ym),col="#fdae6b",border="#fdae6b")}}
  xcor=c(1);ycor=c(0)
  if(nrow(tdf_singleton_s)>0){for(i in 1:nrow(tdf_singleton_s)){xcor=c(xcor,tdf_singleton_s[i,2],tdf_singleton_s[i,2],tdf_singleton_s[i,3],tdf_singleton_s[i,3]);ycor=c(ycor,0,tdf_singleton_s[i,4],tdf_singleton_s[i,4],0)}}
  polygon(xcor,ycor,col=csl[2],border=csl[2])
  xcor=c(1);ycor=c(0)
  if(nrow(tdf_singleton_a)>0){for(i in 1:nrow(tdf_singleton_a)){xcor=c(xcor,tdf_singleton_a[i,2],tdf_singleton_a[i,2],tdf_singleton_a[i,3],tdf_singleton_a[i,3]);ycor=c(ycor,0,-tdf_singleton_a[i,4],-tdf_singleton_a[i,4],0)}}
  xcor=c(xcor,xcor[length(xcor)]);ycor=c(ycor,0)
  polygon(xcor,ycor,col=csl[1],border=csl[1])
  fa=fun_bdgS(tdf_all_s)/fun_bdgS(tdf_singleton_s)
  xcor=c();ycor=c()
  if(nrow(tdf_all_s)>0){for(i in 1:nrow(tdf_all_s)){xcor=c(xcor,tdf_all_s[i,3]);ycor=c(ycor,tdf_all_s[i,4])}}
  if(length(fa)>0){lines(xcor,ycor/fa,col="black")}
  fa=fun_bdgS(tdf_all_a)/fun_bdgS(tdf_singleton_a)
  xcor=c();ycor=c()
  if(nrow(tdf_all_a)>0){for(i in 1:nrow(tdf_all_a)){xcor=c(xcor,tdf_all_a[i,3]);ycor=c(ycor,-tdf_all_a[i,4])}}
  if(length(fa)>0){lines(xcor,ycor/fa,col="black")}
  axis(1,c(1,te.size[tte,1]),label=c(1,te.size[tte,1]),lwd=0)
}
dev.off()


