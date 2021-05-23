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

## experiment 1
gdat<-read.table("../geno/exp1_geno",header=FALSE)
snps<-gdat[,1:3]
gdat<-as.matrix(gdat[,-c(1:3)]) ## get rid of snp dat
pdat<-read.table("../pheno/exp1_pheno",header=FALSE)

## standardize genotypes
X<-gdat
Xmn<-apply(X,1,mean,na.rm=TRUE)
Xsd<-apply(X,1,sd,na.rm=TRUE)
L<-dim(gdat)[1]
for(i in 1:L){
	X[i,]<-(X[i,]-Xmn[i])/Xsd[i]
}

## loop over traits, standardize and drop any with missing data
mapit_ep1<-vector("list",4)
for(j in 1:4){
	Y<-pdat[,j+1]
	mis<-which(is.na(Y)==TRUE)
	XX<-X
	if(length(mis) > 0){
		Y<-Y[-mis]
		XX<-X[,-mis]
	}
	Y<-(Y-mean(Y))/sd(Y)
	mapit_ep1[[j]]<-MAPIT(XX,as.matrix(Y))
}

save(list=ls(),file="mapit.rdat")

snps<-read.table("../gemma/output/snps.txt",header=FALSE)
bnds<-which(snps[-1,1] != snps[-L,1])
mids<-tapply(X=1:L,INDEX=snps[,1],mean)    

cbg<-alpha("cadetblue",.3)

tits<-c("(A) Survival","(B) Male Wgt.","(C) Female Wgt.","(D) DT")

cm<-1.5
ca<-1.1
cl<-1.6
pdf("cowpea_epi1.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
plot(-1 * log10(mapit_ep1[[i]]$pvalues),type='n',axes=FALSE,xlab="Chromosome",ylab="-log10(P)",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,6,6),col=cbg,border=NA)
}   
points(-1 * log10(mapit_ep1[[i]]$pvalues),pch=20,col=alpha("black",.4))
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()

pdat<-read.table("../pheno/exp1_larva_pheno",header=FALSE)
## loop over traits, standardize and drop any with missing data
mapit_ep1L<-vector("list",2)
for(j in 1:2){
	Y<-pdat[,j+1]
	mis<-which(is.na(Y)==TRUE)
	XX<-X
	if(length(mis) > 0){
		Y<-Y[-mis]
		XX<-X[,-mis]
	}
	Y<-(Y-mean(Y))/sd(Y)
	mapit_ep1L[[j]]<-MAPIT(XX,as.matrix(Y))
}

save(list=ls(),file="mapit1L.rdat")


## experiment 2
gdat<-read.table("../geno/exp2and3c_geno",header=FALSE)
snps<-gdat[,1:3]
gdat<-as.matrix(gdat[,-c(1:3)]) ## get rid of snp dat
pdat<-read.table("../pheno/exp2and3_pheno",header=FALSE)
pdatT<-read.table("../pheno/exp2and3_totals_pheno",header=FALSE)
pdat<-cbind(pdat,pdatT[,-1])

## standardize genotypes
X<-gdat
Xmn<-apply(X,1,mean,na.rm=TRUE)
Xsd<-apply(X,1,sd,na.rm=TRUE)
L<-dim(gdat)[1]
for(i in 1:L){
	X[i,]<-(X[i,]-Xmn[i])/Xsd[i]
}

## loop over traits, standardize and drop any with missing data
mapit_ep2<-vector("list",8)
for(j in 1:8){
	Y<-pdat[,j+1]
	mis<-which(is.na(Y)==TRUE)
	XX<-X
	if(length(mis) > 0){
		Y<-Y[-mis]
		XX<-X[,-mis]
	}
	Y<-(Y-mean(Y))/sd(Y)
	mapit_ep2[[j]]<-MAPIT(XX,as.matrix(Y))
}

save(list=ls(),file="mapit2.rdat")

tits<-c("(A) Survival, CA","(B) DT, CA","(C) Survival, SI","(D) DT, SI",
	"(E) Adults, CA","(F) Eggs, CA","(G) Adults, SI","(H) Eggs, SI")

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

pdat<-read.table("../pheno/exp2and3_larv_pheno",header=FALSE)
## loop over traits, standardize and drop any with missing data
K<-as.matrix(read.table("../gemma/output_stan/o_cowpea_fit_ph1_ch0.cXX.txt",header=FALSE))
mapit_ep2L<-vector("list",4)
for(j in 1:4){
        Y<-pdat[,j+1]
        mis<-which(is.na(Y)==TRUE)
        XX<-X
        if(length(mis) > 0){
                Y<-Y[-mis]
                XX<-X[,-mis]
        }
        Y<-(Y-mean(Y))/sd(Y)
        mapit_ep2L[[j]]<-MAPIT(XX,as.matrix(Y),C=K)
}

save(list=ls(),file="mapit2L.rdat")

## NOTE larval characters have p-values totally inflated relative to null, and adding K doesn't help at all.
