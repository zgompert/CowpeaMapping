library(RColorBrewer)
library(scales)
cs<-alpha("gray30",.7)
pf<-list.files(pattern="RIL_exp2and3_")


## desired order
## experiment 1 surv, DT, mweight, fweight, exits
## experiment 2 SI surv, DT, eggs, adults, exits
## experiment 2 CA surv, DT, eggs, adults, exits


######### exp. 1 ############

cl<-1.7
ca<-1.2
cm<-1.6

pdf("fig_cowpea_traits_exp1.pdf",width=7,height=10.5)
par(mfrow=c(3,2))
par(mar=c(4.5,5,3,2))
ph1<-read.table("RIL_exp1_clean.csv",header=TRUE,sep=",")
ph2<-read.table("RIL_exp1_larv_clean.csv",header=TRUE,sep=",")
magic<-grep("MAGIC",ph1$Line)
N<-dim(ph1)[1]
css<-rep(1,N)
css[magic]<-2
barplot(sort(ph1$Survival),col=c("red",cs)[css][order(ph1$Survival)],border=NA,ylab="Survival proportion",cex.lab=cl,xlab="Line")
title(main="(A) Survival",cex.main=cm)
barplot(sort(ph1$DevtTime),col=c("red",cs)[css][order(ph1$DevtTime)],border=NA,ylab="Development time (days)",cex.lab=cl,xlab="Line")
title(main="(B) Development time",cex.main=cm)
barplot(sort(ph1$MeanMaleWt),col=c("red",cs)[css][order(ph1$MeanMaleWt)],border=NA,ylab="Weight (mg)",cex.lab=cl,xlab="Line")
title(main="(C) Male weight",cex.main=cm)
barplot(sort(ph1$MeanFemaleWt),col=c("red",cs)[css][order(ph1$MeanFemaleWt)],border=NA,ylab="Weight (mg)",cex.lab=cl,xlab="Line")
title(main="(D) Female weight",cex.main=cm)
magic<-grep("MAGIC",ph2$Line)
N<-dim(ph2)[1]
css<-rep(1,N)
css[magic]<-2
y<-ph2$No.larvae/apply(ph2[,-1],1,sum)
y[y==0]<-0.001
barplot(sort(y),col=c("red",cs)[css][order(y)],border=NA,ylab="Proportion early exit larva",cex.lab=cl,xlab="Line")
title(main="(E) Early exiting larva",cex.main=cm)
plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
legend(0.1,0.9,c("MAGIC RILs","Parents"),fill=c(cs,"red"),cex=1.5,bty='n')

dev.off()





########## exp. 2 and 3 #################

cl<-1.7
ca<-1.2
cm<-1.6

pdf("fig_cowpea_traits_exp2.pdf",width=7,height=10.5)
par(mfrow=c(3,2))
par(mar=c(4.5,5,3,2))

ph<-read.table(pf[1],sep=",",header=TRUE)
magic<-grep("MAGIC",ph$MAGIC)
plot(ph$CA_SURVIVAL,ph$SI_SURVIVAL,pch=20,col=cs,xlab="Survival prop. CA",ylab="Survival prop. SI",cex.lab=cl,cex.axis=ca)
points(ph$CA_SURVIVAL[-magic],ph$SI_SURVIVAL[-magic],pch=19,col="red")
r<-cor(ph$CA_SURVIVAL,ph$SI_SURVIVAL)
mtext(round(r,2),3,adj=0.05,line=-2,cex=ca)
title(main="(A) Survival",cex.main=cm)

plot(ph$CA_DEV,ph$SI_DEV,pch=20,col=cs,xlab="Dev. time (days) CA",ylab="Dev. time (days) SI",cex.lab=cl,cex.axis=ca)
points(ph$CA_DEV[-magic],ph$SI_DEV[-magic],pch=19,col="red")
r<-cor.test(ph$CA_DEV,ph$SI_DEV,na.rm=TRUE)
mtext(round(r$estimate,2),3,adj=0.05,line=-2,cex=ca)
title(main="(B) Development time",cex.main=cm)

ph<-read.table(pf[3],sep=",",header=TRUE)
magic<-grep("MAGIC",ph$MAGIC)
plot(ph$CA_EGGS,ph$SI_EGGS,pch=20,col=cs,xlab="No. eggs CA",ylab="No. eggs SI",cex.lab=cl,cex.axis=ca)
points(ph$CA_EGGS[-magic],ph$SI_EGGS[-magic],pch=19,col="red")
r<-cor(ph$CA_EGGS,ph$SI_EGGS)
mtext(round(r,2),3,adj=0.05,line=-2,cex=ca)
title(main="(C) Number eggs",cex.main=cm)

plot(ph$CA_ADULTS,ph$SI_ADULTS,pch=20,col=cs,xlab="No. adults CA",ylab="No. adults SI",cex.lab=cl,cex.axis=ca)
points(ph$CA_ADULTS[-magic],ph$SI_ADULTS[-magic],pch=19,col="red")
r<-cor(ph$CA_ADULTS,ph$SI_ADULTS)
mtext(round(r,2),3,adj=0.05,line=-2,cex=ca)
title(main="(D) Number adults",cex.main=cm)

ph<-read.table(pf[2],sep=",",header=TRUE)
magic<-grep("MAGIC",ph$MAGIC)
plot(ph$CA.PROP.LARVAE,ph$SI.PROP.LARVAE,pch=20,col=cs,xlab="Early exit CA",ylab="Early exit SI",cex.lab=cl,cex.axis=ca)
points(ph$CA.PROP.LARVAE[-magic],ph$SI.PROP.LARVAE[-magic],pch=19,col="red")
r<-cor(ph$CA.PROP.LARVAE,ph$SI.PROP.LARVAE)
mtext(round(r,2),3,adj=0.05,line=-2,cex=ca)
title(main="(E) Prop. early exits",cex.main=cm)

plot(c(0,1),c(0,1),type='n',axes=FALSE,xlab="",ylab="")
legend(0.1,0.9,c("MAGIC RILs","Parents"),pch=c(20,19),col=c(cs,"red"),cex=1.5,bty='n')

dev.off()

########

#ph1<-read.table("RIL_exp1_clean.csv",header=TRUE,sep=",")
#ph2<-read.table("RIL_exp1_larv_clean.csv",header=TRUE,sep=",")
#magic1<-grep("MAGIC",ph1$Line)
#magic2<-grep("MAGIC",ph2$Line)
   

