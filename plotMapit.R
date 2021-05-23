### Load in the R libraries ###
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(CompQuadForm)
library(scales)

### Load in functions to make QQ-plot plots ###
source("MAPIT/QQPlot.R")

#NOTE: This code assumes that the basic C++ functions are set up on the computer in use. If not, the MAPIT functions and Rcpp packages will not work properly. Mac users please refer to the homebrew applications and install the gcc commands listed in the README.md file before running the rest of the code [Warning: This step may take about an hour...].

### Load in the C++ MAPIT wrapper functions ###
source("MAPIT/Standard Version/MAPIT.R"); sourceCpp("MAPIT/Standard Version/MAPIT.cpp")

load("mapit.rdat")
#tits<-c("(A) Survival","(B) Male Wgt.","(C) Female Wgt.","(D) DT")
load("mapit2.rdat")
#tits<-c("(A) Survival, CA","(B) DT, CA","(C) Survival, SI","(D) DT, SI",
	#"(E) Adults, CA","(F) Eggs, CA","(G) Adults, SI","(H) Eggs, SI")

snps<-read.table("../gemma/output/snps_stan.txt",header=FALSE)
bnds<-which(snps[-1,1] != snps[-L,1])
mids<-tapply(X=1:L,INDEX=snps[,1],mean)    

cbg<-alpha("cadetblue",.3)

mapit_c<-list(mapit_ep1[[1]],mapit_ep1[[4]],mapit_ep2[[3]],mapit_ep2[[4]],mapit_ep2[[1]],mapit_ep2[[2]])

tits<-c("(A) SI survival","(B) SI development time",
	"(C) SI survival","(D) SI development time",
	"(E) CA survival","(F) CA development time")	

cm<-1.5
ca<-1.1
cl<-1.6
pdf("fig_cowpea_epi.pdf",width=8,height=9)
layout(matrix(c(1,1,3:4,2,2,5:8),nrow=5,ncol=2,byrow=TRUE),widths=c(4,4),heights=c(1,3.3,1,3.3,3.3))

par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
text(0.5,.5,"No competition",cex=2.8)
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
text(0.5,.5,"Competition",cex=2.8)


par(mar=c(4.5,5,3,1))
for(i in 1:6){
plot(-1 * log10(mapit_c[[i]]$pvalues),ylim=c(0,7),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,8,8),col=cbg,border=NA)
}   
points(-1 * log10(mapit_c[[i]]$pvalues),pch=20,col=alpha("gray10",.5))
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
abline(h=-1*log10(0.05/32130),lwd=2,lty=3,col="red")
}
dev.off()

######################################33


cm<-1.5
ca<-1.1
cl<-1.6
pdf("cowpea_epi2.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:8){
plot(-1 * log10(mapit_ep2[[i]]$pvalues),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,8,8),col=cbg,border=NA)
}   
points(-1 * log10(mapit_ep2[[i]]$pvalues),pch=20,col=alpha("black",.4))
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()

