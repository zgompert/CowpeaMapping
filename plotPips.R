library(scales)
library(RColorBrewer)


snps<-read.table("snps.txt",header=FALSE)

## pip plots
bnds<-which(snps[-1,1] != snps[-L,1])
mids<-tapply(X=1:L,INDEX=snps[,1],mean)

cbg<-alpha("cadetblue",.3)

## exp 1
pf<-list.files(pattern="pip_o_cowpea1_fit_ph")
L<-32130
pips_exp1<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips_exp1[,j]<-scan(pf[j])
    }
#tits<-c("(A) Survival","(B) Male Wgt.","(C) Female Wgt.","(D) DT")

## all larva early exits 
pf<-c(list.files(pattern="pip_o_cowpea1_larv_fit_ph"),list.files(pattern="pip_o_cowpea_larv_fit_ph"))
L<-32130
pips_l<-matrix(NA,nrow=L,ncol=6)
for(j in 1:6){
    pips_l[,j]<-scan(pf[j])
    }
#tits<-c("(A) No. early, SI1","(B) Prop. early, SI1","(C) No. early, CA","(D) Prop. early, CA","(E) No. early, SI","(F) Prop. early, SI")


## read exp 2 survival and DT
pf<-list.files(pattern="pip_o_cowpea_fit_ph")
L<-32130
pips_exp2<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips_exp2[,j]<-scan(pf[j])
    }

#### make first set, surv, dt, early exit
pips1<-cbind(pips_exp1[,c(1,4)],pips_l[,2],pips_exp2[,c(3,4)],pips_l[,6],pips_exp2[,c(1:2)],pips_l[,4])

tits<-c("(A) SI survival","(B) SI development time","(C) SI early-exiting larvae",
	"(D) SI survival","(E) SI development time","(F) SI early-exiting larvae",
	"(G) CA survival","(H) CA development time","(I) CA early-exiting larvae")

cm<-1.5
ca<-1.1
cl<-1.6
pdf("fig_cowpea_pips_main.pdf",width=12,height=11)
par(mfrow=c(3,3))
layout(matrix(c(1,1,1,3:5,2,2,2,6:11),nrow=5,ncol=3,byrow=TRUE),widths=c(4,4,4),heights=c(1,3,1,3,3))

par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
text(0.5,.5,"No competition",cex=3.4)
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
text(0.5,.5,"Competition",cex=3.4)

par(mar=c(4.5,5,3,1))
for(i in 1:9){
xx<-1:L
plot(xx,pips1[,i],type='n',axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,1,1),col=cbg,border=NA)
}
points(xx,pips1[,i],pch=20,col=alpha("gray10",.5))
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()


###############################################
## read exp 2 totals
pf<-list.files(pattern="pip_o_cowpea_totals_fit_ph")
L<-32130
pips_t2<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips_t2[,j]<-scan(pf[j])
    }
#tits<-c("(A) No. adults, CA","(B) No. eggs, CA","(C) No. adults, SI","(D) No. eggs, SI")

pips2<-cbind(pips_exp1[,c(2,3)],pips_t2[,c(4,3,2,1)])


tits<-c("(A) SI male weight","(B) SI female weight",
	"(C) SI no. eggs","(D) SI no. adults",
	"(E) CA no. eggs","(F) CA no. adults")

cm<-1.5
ca<-1.07
cl<-1.6
pdf("fig_cowpea_pips_supp.pdf",width=8,height=8)
layout(matrix(c(1,1,3:4,2,2,5:8),nrow=5,ncol=2,byrow=TRUE),widths=c(4,4),heights=c(1,3.3,1,3.3,3.3))

par(mar=c(0,0,0,0))
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
text(0.5,.5,"No competition",cex=2.8)
plot(c(0,1),c(0,1),type='n',xlab="",ylab="",axes=FALSE)
text(0.5,.5,"Competition",cex=2.8)

par(mar=c(4.5,5,3,1))
for(i in 1:6){
xx<-1:L
plot(xx,pips2[,i],type='n',axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,1,1),col=cbg,border=NA)
}
points(xx,pips2[,i],pch=20,col=alpha("gray10",.5))
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()

###### quantitative summary
pipsc<-cbind(pips1,pips2)

snps[apply(pips1,1,max) > 0.1,]
#      V1      V2
#6956   3 36.5060
#6960   3 36.5060
#7121   3 36.8662
#13453  5 13.7582
#13469  5 13.7582
#13473  5 13.7582
#13545  5 15.1234
#13548  5 15.1234
#13549  5 15.1234
#13553  5 15.1234
#13554  5 15.1234
#13559  5 15.1234
#13564  5 15.1234
#13568  5 15.1234
#13579  5 15.1234
#13586  5 15.1234
#13593  5 16.0963
#13596  5 16.0963
#13655  5 18.5364
#13657  5 18.5364
#13659  5 18.5364
#13661  5 18.5364
#16241  6  2.3486
#16280  6  2.7266
#19612  7 43.5156
#21223  7 96.6338
#23596  8 57.8903
#23602  8 57.8903
#23604  8 57.8903
#23664  8 62.4202
#23689  8 62.4202
#23692  8 62.4202
#24016  8 72.4500
#24105  8 75.8805

lg5reg<-which(snps[,1] == 5 & (snps[,2] >=13.7582 & snps[,2] <= 18.5364))

(1-apply(1-pipsc[lg5reg,],2,prod))[c(2,4:5,7:8,11,13)]
#[1] 0.8903607 0.6785720 0.9061177 0.6864878 0.9914207 0.8475235 0.5695464
 mean((1-apply(1-pipsc[lg5reg,],2,prod))[c(2,4:5,7:8,11,13)])
#[1] 0.7957184



ano<-read.table("/uufs/chpc.utah.edu/common/home/gompert-group3/projects/cowpea_mapping/annotation/snpChromPosGene.txt",header=TRUE)

unique(as.character(ano$Gene[lg5reg[(apply(pipsc[lg5reg,c(2,4:5,7:8,11,13)],1,max) > 0.2)]]) 
#[1] "Vigun05g046000" "Vigun05g047400" NA               "Vigun05g052600"
# [5] "Vigun05g052800" "Vigun05g052900" "Vigun05g053000" "Vigun05g053200"
# [9] "Vigun05g054900" "Vigun05g060200" "Vigun05g060500"

## get direction of effect
## exp 1
pf<-list.files(pattern="mav_o_cowpea1_fit_ph")
L<-32130
mav_exp1<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    mav_exp1[,j]<-scan(pf[j])
    }
#tits<-c("(A) Survival","(B) Male Wgt.","(C) Female Wgt.","(D) DT")

## all larva early exits 
pf<-c(list.files(pattern="mav_o_cowpea1_larv_fit_ph"),list.files(pattern="mav_o_cowpea_larv_fit_ph"))
L<-32130
mav_l<-matrix(NA,nrow=L,ncol=6)
for(j in 1:6){
    mav_l[,j]<-scan(pf[j])
    }
#tits<-c("(A) No. early, SI1","(B) Prop. early, SI1","(C) No. early, CA","(D) Prop. early, CA","(E) No. early, SI","(F) Prop. early, SI")


## read exp 2 survival and DT
pf<-list.files(pattern="mav_o_cowpea_fit_ph")
L<-32130
mav_exp2<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    mav_exp2[,j]<-scan(pf[j])
    }

## read exp 2 totals
pf<-list.files(pattern="mav_o_cowpea_totals_fit_ph")
L<-32130
mav_t2<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    mav_t2[,j]<-scan(pf[j])
    }
#tits<-c("(A) No. adults, CA","(B) No. eggs, CA","(C) No. adults, SI","(D) No. eggs, SI")


mav1<-cbind(mav_exp1[,c(1,4)],mav_l[,2],mav_exp2[,c(3,4)],mav_l[,6],mav_exp2[,c(1:2)],mav_l[,4])
mav2<-cbind(mav_exp1[,c(2,3)],mav_t2[,c(4,3,2,1)])
mavc<-cbind(mav1,mav2)

tits<-c("(A) SI survival","(B) SI development time","(C) SI early-exiting larvae",
        "(D) SI survival","(E) SI development time","(F) SI early-exiting larvae",
        "(G) CA survival","(H) CA development time","(I) CA early-exiting larvae",
        "(A) SI male weight","(B) SI female weight",
        "(C) SI no. eggs","(D) SI no. adults",
        "(E) CA no. eggs","(F) CA no. adults")

ano[lg5reg[(apply(pipsc[lg5reg,c(2,4:5,7:8,11,13)],1,max) > 0.2)],]

#          SNP OldChr   OldCm NewChr  NewPos           Gene
#13453 2_01390      5 13.7582      5 3829851 Vigun05g046000
#13469 2_07583      5 13.7582      5 3966407 Vigun05g047400
#13473 2_06553      5 13.7582      5 3994844           <NA>
#13549 2_00868      5 15.1234      5 4497448 Vigun05g052600
#13553 2_08113      5 15.1234      5 4517159 Vigun05g052800
#13554 2_01320      5 15.1234      5 4522988 Vigun05g052900
#13559 2_51427      5 15.1234      5 4551502 Vigun05g053000
#13564 2_13010      5 15.1234      5 4558350 Vigun05g053200
#13593 2_41936      5 16.0963      5 4686167           <NA>
#13596 2_22350      5 16.0963      5 4700276 Vigun05g054900
#13655 2_34617      5 18.5364      5 5156508 Vigun05g060200
#13657 2_30291      5 18.5364      5 5165555           <NA>
#13661 2_24125      5 18.5364      5 5172115 Vigun05g060500
wm<-apply(abs(mavc[lg5reg[(apply(pipsc[lg5reg,c(2,4:5,7:8,11,13)],1,max) > 0.2)],c(2,4:5,7:8> 1,13)]),1,which.max)
meff<-rep(NA,length(wm))
for(i in 1:length(wm)){
	meff[i]<-mavc[lg5reg[(apply(pipsc[lg5reg,c(2,4:5,7:8,11,13)],1,max) > 0.2)][i],c(2,4:5,7:8,11,13)[wm[i]]]
}

cbind(ano[lg5reg[(apply(pipsc[lg5reg,c(2,4:5,7:8,11,13)],1,max) > 0.2)],],meff)
#          SNP OldChr   OldCm NewChr  NewPos           Gene         meff
#13453 2_01390      5 13.7582      5 3829851 Vigun05g046000 -0.005026806
#13469 2_07583      5 13.7582      5 3966407 Vigun05g047400 -0.014884395
#13473 2_06553      5 13.7582      5 3994844           <NA> -0.012133986
#13549 2_00868      5 15.1234      5 4497448 Vigun05g052600 -0.170354123
#13553 2_08113      5 15.1234      5 4517159 Vigun05g052800  0.220501423
#13554 2_01320      5 15.1234      5 4522988 Vigun05g052900  0.169864156
#13559 2_51427      5 15.1234      5 4551502 Vigun05g053000  0.211375487
#13564 2_13010      5 15.1234      5 4558350 Vigun05g053200 -0.435266453
#13593 2_41936      5 16.0963      5 4686167           <NA>  0.185914087
#13596 2_22350      5 16.0963      5 4700276 Vigun05g054900  0.189039996
#13655 2_34617      5 18.5364      5 5156508 Vigun05g060200  0.076913912
#13657 2_30291      5 18.5364      5 5165555           <NA>  0.078598979
#13661 2_24125      5 18.5364      5 5172115 Vigun05g060500  0.080620373
cbind(ano[lg5reg[(apply(pipsc[lg5reg,c(2,4:5,7:8,11,13)],1,max) > 0.2)],],meff,tits[c(2,4:5,7:8,11,13)[wm]])
#1  -0.005026806             (B) SI female weight
#2  -0.014884395             (B) SI female weight
#3  -0.012133986             (B) SI female weight
#4  -0.170354123          (B) SI development time
#5   0.220501423          (E) SI development time
#6   0.169864156          (B) SI development time
#7   0.211375487          (E) SI development time
#8  -0.435266453          (E) SI development time
#9   0.185914087          (B) SI development time
#10  0.189039996          (B) SI development time
#11  0.076913912                  (D) SI survival
#12  0.078598979                  (D) SI survival
#13  0.080620373                  (D) SI survival



lg8reg<-which(snps[,1] == 8 & (snps[,2] >= 57.8903 & snps[,2] <= 62.4202))

(1-apply(1-pipsc[lg8reg,],2,prod))[c(3,6,9)]


unique(as.character(ano$Gene[lg8reg[(apply(pipsc[lg8reg,c(3,6,9)],1,max) > 0.2)]]))
#[1] "Vigun08g171700" "Vigun08g179000" "Vigun08g180500"

wm<-apply(abs(mavc[lg8reg[(apply(pipsc[lg8reg,c(3,6,9)],1,max) > 0.2)],c(3,6,9)]),1,which.max)
meff<-rep(NA,length(wm))
for(i in 1:length(wm)){
	meff[i]<-mavc[lg8reg[(apply(pipsc[lg8reg,c(3,6,9)],1,max) > 0.2)][i],c(3,6,9)[wm[i]]]
}

cbind(ano[lg8reg[(apply(pipsc[lg8reg,c(3,6,9)],1,max) > 0.2)],],meff)
#          SNP OldChr   OldCm NewChr   NewPos           Gene          meff
#23596 2_01925      8 57.8903      8 34260474 Vigun08g171700 -0.3211541521
#23664 2_28941      8 62.4202      8 34887678 Vigun08g179000  0.1338002696
#23689  1_0226      8 62.4202      8 35020285 Vigun08g180500  0.0005673578

cbind(ano[lg8reg[(apply(pipsc[lg8reg,c(3,6,9)],1,max) > 0.2)],],meff,tits[c(3,6,9)[wm]])
#23596 (C) SI early-exiting larvae
#23664 (C) SI early-exiting larvae
#23689 (I) CA early-exiting larvae
 









##### old/individual figures
pf<-list.files(pattern="pip_o_cowpea_fit_ph")
L<-32130
pips<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips[,j]<-scan(pf[j])
    }
    
snps<-read.table("snps.txt",header=FALSE)

## pip plots
bnds<-which(snps[-1,1] != snps[-L,1])
mids<-tapply(X=1:L,INDEX=snps[,1],mean)    

cbg<-alpha("cadetblue",.3)

tits<-c("(A) Survival, CA","(B) DT, CA","(C) Survival, SI","(D) DT, SI")


cm<-1.5
ca<-1.1
cl<-1.6
pdf("cowpea_pips2.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
plot(pips[,i],type='n',axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,1,1),col=cbg,border=NA)
}   
points(pips[,i],pch=20)
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()

cs<-brewer.pal(n=4,"Paired")[c(1,1,3,3)]


pdf("cowpea_qtl2.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
barplot(tapply(X=pips[,i],INDEX=snps[,1],sum),xlab="Chromosome",ylab="Number QTL",cex.lab=cl,
    cex.names=ca,cex.axis=ca,col=cs[i])
title(main=tits[i],cex.main=cm)
}
dev.off()

############# experiment 2 totals
pf<-list.files(pattern="pip_o_cowpea_totals_fit_ph")
L<-32130
pips<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips[,j]<-scan(pf[j])
    }
tits<-c("(A) No. adults, CA","(B) No. eggs, CA","(C) No. adults, SI","(D) No. eggs, SI")

cm<-1.5
ca<-1.1
cl<-1.6
pdf("cowpea_pips2_totals.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
plot(pips[,i],type='n',axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,1,1),col=cbg,border=NA)
}   
points(pips[,i],pch=20)
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()

cs<-brewer.pal(n=4,"Paired")[c(1,1,3,3)]


pdf("cowpea_qtl2_totals.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
barplot(tapply(X=pips[,i],INDEX=snps[,1],sum),xlab="Chromosome",ylab="Number QTL",cex.lab=cl,
    cex.names=ca,cex.axis=ca,col=cs[i])
title(main=tits[i],cex.main=cm)
}
dev.off()

############# larv early exit
pf<-c(list.files(pattern="pip_o_cowpea1_larv_fit_ph"),list.files(pattern="pip_o_cowpea_larv_fit_ph"))
L<-32130
pips<-matrix(NA,nrow=L,ncol=6)
for(j in 1:6){
    pips[,j]<-scan(pf[j])
    }
tits<-c("(A) No. early, SI1","(B) Prop. early, SI1","(C) No. early, CA","(D) Prop. early, CA","(E) No. early, SI","(F) Prop. early, SI")

cm<-1.5
ca<-1.1
cl<-1.6
pdf("cowpea_pips_larv.pdf",width=11,height=13.5)
par(mfrow=c(3,2))
par(mar=c(4.5,5,3,1))
for(i in 1:6){
plot(pips[,i],type='n',axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,1,1),col=cbg,border=NA)
}   
points(pips[,i],pch=20)
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()

cs<-brewer.pal(n=4,"Paired")[c(1,1,3,3)]


pdf("cowpea_qtl2_larv.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
barplot(tapply(X=pips[,i],INDEX=snps[,1],sum),xlab="Chromosome",ylab="Number QTL",cex.lab=cl,
    cex.names=ca,cex.axis=ca,col=cs[i])
title(main=tits[i],cex.main=cm)
}
dev.off()


############## experiment 1
library(scales)
library(RColorBrewer)

pf<-list.files(pattern="pip_o_cowpea1_fit_ph")
L<-32130
pips<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips[,j]<-scan(pf[j])
    }
    
snps<-read.table("snps.txt",header=FALSE)

## pip plots
bnds<-which(snps[-1,1] != snps[-L,1])
mids<-tapply(X=1:L,INDEX=snps[,1],mean)    

cbg<-alpha("cadetblue",.3)

tits<-c("(A) Survival","(B) Male Wgt.","(C) Female Wgt.","(D) DT")


cm<-1.5
ca<-1.1
cl<-1.6
pdf("cowpea_pips1.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
plot(pips[,i],type='n',axes=FALSE,xlab="Chromosome",ylab="Inclusion probability",cex.lab=cl)
for(j in seq(1,9,2)){
    polygon(c(bnds[c(j,j+1,j+1,j)]),c(-.1,-.1,1,1),col=cbg,border=NA)
}   
points(pips[,i],pch=20)
axis(2,cex.axis=ca)
axis(1,at=mids,c(1:11),cex.axis=ca)
box()
title(main=tits[i],cex.main=cm)
}
dev.off()

cs<-brewer.pal(n=4,"Paired")[c(3,3,3,3)]


pdf("cowpea_qtl1.pdf",width=11,height=9)
par(mfrow=c(2,2))
par(mar=c(4.5,5,3,1))
for(i in 1:4){
barplot(tapply(X=pips[,i],INDEX=snps[,1],sum),xlab="Chromosome",ylab="Number QTL",cex.lab=cl,
    cex.names=ca,cex.axis=ca,col=cs[i])
title(main=tits[i],cex.main=cm)
}
dev.off()
