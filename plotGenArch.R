library(RColorBrewer)
## SI surv, mweight, fweight, DT
gen1<-read.table("genarch_cowpea1.txt",header=F)
## CA surv, CA DT, SI surv., SI DT
gen2<-read.table("genarch_cowpea2.txt",header=F)
## total
gen3<-read.table("genarch_cowpea2_totals.txt",header=FALSE)
## SI early exits
gen4<-read.table("genarch_cowpea1_larv.txt",header=FALSE)
## CA and SI early exits
gen5<-read.table("genarch_cowpea2_larv.txt",header=FALSE)

## desired order
## experiment 1 surv, DT, mweight, fweight, exits
## experiment 2 SI surv, DT, eggs, adults, exits
## experiment 2 CA surv, DT, eggs, adults, exits

gen<-rbind(gen1[c(1,4,2,3),],gen4[2,],gen2[c(3,4),],gen3[c(4,3),],gen5[4,],
	gen2[c(1,2),],gen3[c(2,1),],gen5[2,])


cs<-brewer.pal(n=4,"Paired")

cl<-1.5
ca<-1.2
ph<-c("Surv","DT","MWg","FWg","EEx","Surv","DT","NEgg","NAdu","EEx","Surv","DT","Negg","NAdu","EEx")
pdf("fig_cowpea_pve.pdf",width=9,height=5)
par(mar=c(5,5,1,1))
x<-barplot(gen[,1]*100,ylim=c(0,100),col="white",names.arg=ph,xlab="Trait",ylab="Variance explained (%)",cex.lab=cl,cex.axis=ca,cex.names=ca,las=2)
x<-barplot(gen[,1]*gen[,4]*100,add=TRUE,col=cs[c(rep(3,10),rep(1,5))],axes=FALSE)
segments(x,gen[,2]*100,x,gen[,3]*100)
text(x,rep(95,15),gen[,7],cex.text=1.2)
box()
legend(.1,88,c("CA","SI"),fill=cs[c(1,3)],bty='n',cex=1.1)
dev.off()













############## sub plots ###################
ph<-c("Survival","DT","Survival","DT")
## CA survival, CA DT, SI survival, SI DT
gen<-read.table("genarch_cowpea2.txt",header=F)

cs<-brewer.pal(n=4,"Paired")

cl<-1.5
ca<-1.2

pdf("cowpea_pve_exp2.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
x<-barplot(gen[,1]*100,ylim=c(0,100),col="white",names.arg=ph,xlab="Trait",ylab="Variance explained (%)",cex.lab=cl,cex.axis=ca,cex.names=ca)
x<-barplot(gen[,1]*gen[,4]*100,add=TRUE,col=cs[c(1,1,3,3)],axes=FALSE)
segments(x,gen[,2]*100,x,gen[,3]*100)
box()
legend(.1,99,c("CA","SI"),fill=cs[c(1,3)],bty='n',cex=1.1)
dev.off()

ph<-c("Adults","Eggs","Adults","Eggs")
gen<-read.table("genarch_cowpea2_totals.txt",header=F)
pdf("cowpea_pve_exp2_totals.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
x<-barplot(gen[,1]*100,ylim=c(0,100),col="white",names.arg=ph,xlab="Trait",ylab="Variance explained (%)",cex.lab=cl,cex.axis=ca,cex.names=ca)
x<-barplot(gen[,1]*gen[,4]*100,add=TRUE,col=cs[c(1,1,3,3)],axes=FALSE)
segments(x,gen[,2]*100,x,gen[,3]*100)
box()
legend(.1,99,c("CA","SI"),fill=cs[c(1,3)],bty='n',cex=1.1)
dev.off()

ph<-c("Num","Prop","Num","Prop","Num","Prop")
gen<-read.table("genarch_cowpea2_larv.txt",header=F)
gen1<-read.table("genarch_cowpea1_larv.txt",header=FALSE)
gen<-rbind(gen1,gen)

pdf("cowpea_pve_larv.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
x<-barplot(gen[,1]*100,ylim=c(0,100),col="white",names.arg=ph,xlab="Trait",ylab="Variance explained (%)",cex.lab=cl,cex.axis=ca,cex.names=ca)
x<-barplot(gen[,1]*gen[,4]*100,add=TRUE,col=cs[c(3,3,1,1,3,3)],axes=FALSE)
segments(x,gen[,2]*100,x,gen[,3]*100)
box()
legend(.1,99,c("CA","SI"),fill=cs[c(1,3)],bty='n',cex=1.1)
dev.off()

ph<-c("Survival","MWght","FWght","DT")
## SI only
gen<-read.table("genarch_cowpea1.txt",header=F)

cs<-brewer.pal(n=4,"Paired")

cl<-1.5
ca<-1.2

pdf("cowpea_pve_exp1.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
x<-barplot(gen[,1]*100,ylim=c(0,100),col="white",names.arg=ph,xlab="Trait",ylab="Variance explained (%)",cex.lab=cl,cex.axis=ca,cex.names=ca)
x<-barplot(gen[,1]*gen[,4]*100,add=TRUE,col=cs[c(3,3,3,3)],axes=FALSE)
segments(x,gen[,2]*100,x,gen[,3]*100)
box()
dev.off()
