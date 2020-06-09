#!/usr/bin/env Rscript

### arguments ----
Args=commandArgs()
if(length(Args)<6){
	stop(paste(substr(Args[4],8,100),"in.soma.rate.bed",sep=" "))
}

### main program ----
soma = read.table(Args[6], header=F, row.names=NULL)
lamda_local = (soma[,3])*soma[,4]/(soma[,5])
soma_uniq = soma[!duplicated(soma[,1]),]
lamda_global = sum(soma_uniq[,3])*soma[,4]/sum(soma_uniq[,5])
lamda_local[which(is.na(lamda_local))]=0
lamda_global[which(is.na(lamda_global))]=0
lamda = lamda_global
for(i in 1:dim(soma)[1]){
	if((soma[i,2]+soma[i,3])>50 & soma[i,5]>50){
		lamda[i] = lamda_local[i]
	}
}
q95 = qpois(0.95, lamda)
fo = cbind(soma, round(soma[,2]-lamda,2), soma[,2]-q95, round(lamda,2), q95)
if(length(which(fo[,6]<0))>0){fo[which(fo[,6]<0),6] = 0}
if(length(which(fo[,7]<0))>0){fo[which(fo[,7]<0),7] = 0}
write.table(fo, Args[6], row.names=F, col.names=F, sep="\t", quote=F)
