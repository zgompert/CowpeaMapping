## calculate GEBVs for RG and GB with and w/o mapit epistasis SNPs 

G1<-as.matrix(read.table("../../geno/exp1_geno",header=FALSE)[,-c(1:3)])
for(k in 1:dim(G1)[1]){
	G1[k,]<-G1[k,]-mean(G1[k,])
}


## get list of effect files from GEMMA
f_eff_ph1<-list.files(pattern=glob2rx("o_cowpea1_fit_ph1_ch*param.txt"))
f_eff_ph2<-list.files(pattern=glob2rx("o_cowpea1_fit_ph2_ch*param.txt"))
f_eff_ph3<-list.files(pattern=glob2rx("o_cowpea1_fit_ph3_ch*param.txt"))
f_eff_ph4<-list.files(pattern=glob2rx("o_cowpea1_fit_ph4_ch*param.txt"))
f_eff_2ph1<-list.files(pattern=glob2rx("o_cowpea_fit_ph1_ch*param.txt"))
f_eff_2ph2<-list.files(pattern=glob2rx("o_cowpea_fit_ph2_ch*param.txt"))
f_eff_2ph3<-list.files(pattern=glob2rx("o_cowpea_fit_ph3_ch*param.txt"))
f_eff_2ph4<-list.files(pattern=glob2rx("o_cowpea_fit_ph4_ch*param.txt"))
f_eff_2tph1<-list.files(pattern=glob2rx("o_cowpea_totals_fit_ph1_ch*param.txt"))
f_eff_2tph2<-list.files(pattern=glob2rx("o_cowpea_totals_fit_ph2_ch*param.txt"))
f_eff_2tph3<-list.files(pattern=glob2rx("o_cowpea_totals_fit_ph3_ch*param.txt"))
f_eff_2tph4<-list.files(pattern=glob2rx("o_cowpea_totals_fit_ph4_ch*param.txt"))
f_eff_lph1<-list.files(pattern=glob2rx("o_cowpea1_larv_fit_ph1_ch*param.txt"))
f_eff_lph2<-list.files(pattern=glob2rx("o_cowpea1_larv_fit_ph2_ch*param.txt"))
f_eff_2lph1<-list.files(pattern=glob2rx("o_cowpea_larv_fit_ph1_ch*param.txt"))
f_eff_2lph2<-list.files(pattern=glob2rx("o_cowpea_larv_fit_ph2_ch*param.txt"))
f_eff_2lph3<-list.files(pattern=glob2rx("o_cowpea_larv_fit_ph3_ch*param.txt"))
f_eff_2lph4<-list.files(pattern=glob2rx("o_cowpea_larv_fit_ph4_ch*param.txt"))

## function to calcualte GEBVs
calcGEBV<-function(ff=NA,G=NA){
	N<-length(ff)
	gebv<-vector("list",N)
	for(j in 1:N){
		eff<-read.table(ff[j],header=TRUE)
		L<-dim(eff)[1]
		I<-dim(G)[2]
		gebv[[j]]<-rep(NA,I)
		mav<-eff$alpha + eff$beta * eff$gamma
		for(i in 1:I){
			gebv[[j]][i]<-sum(mav*G[,i])
		}
	}
	gebv<-matrix(unlist(gebv),nrow=I,ncol=N)
	return(gebv)
}

gebv_p1<-apply(calcGEBV(ff=f_eff_ph1,G=G1),1,mean)
gebv_p2<-apply(calcGEBV(ff=f_eff_ph2,G=G1),1,mean)
gebv_p3<-apply(calcGEBV(ff=f_eff_ph3,G=G1),1,mean)
gebv_p4<-apply(calcGEBV(ff=f_eff_ph4,G=G1),1,mean)

gebv_2p1<-apply(calcGEBV(ff=f_eff_2ph1,G=G1),1,mean)
gebv_2p2<-apply(calcGEBV(ff=f_eff_2ph2,G=G1),1,mean)
gebv_2p3<-apply(calcGEBV(ff=f_eff_2ph3,G=G1),1,mean)
gebv_2p4<-apply(calcGEBV(ff=f_eff_2ph4,G=G1),1,mean)

gebv_2tp1<-apply(calcGEBV(ff=f_eff_2tph1,G=G1),1,mean)
gebv_2tp2<-apply(calcGEBV(ff=f_eff_2tph2,G=G1),1,mean)
gebv_2tp3<-apply(calcGEBV(ff=f_eff_2tph3,G=G1),1,mean)
gebv_2tp4<-apply(calcGEBV(ff=f_eff_2tph4,G=G1),1,mean)

gebv_2lp1<-apply(calcGEBV(ff=f_eff_2lph1,G=G1),1,mean)
gebv_2lp2<-apply(calcGEBV(ff=f_eff_2lph2,G=G1),1,mean)
gebv_2lp3<-apply(calcGEBV(ff=f_eff_2lph3,G=G1),1,mean)
gebv_2lp4<-apply(calcGEBV(ff=f_eff_2lph4,G=G1),1,mean)
gebv_lp1<-apply(calcGEBV(ff=f_eff_lph1,G=G1),1,mean)
gebv_lp2<-apply(calcGEBV(ff=f_eff_lph2,G=G1),1,mean)


gebvM<-cbind(gebv_p1,gebv_p4,gebv_p2,gebv_p3,gebv_lp2,gebv_2p3,gebv_2p4,gebv_2tp4,gebv_2tp3,gebv_2lp4,gebv_2p1,gebv_2p2,gebv_2tp2,gebv_2tp1,gebv_2lp2)
nms<-c("Surv","D","MWg","FWg","EEx","Surv","DT","NEgg","Nadu","EEx","Surv","DT","Negg","Nadu","EEx")

colnames(gebvM)<-nms
library(RColorBrewer)
library(fields)
cs<-brewer.pal(11,"BrBG")
pdf("cowpea_gmatrix.pdf",width=6,height=6)
par(mar=c(6,6,6,6))
brks<-c(rev(seq(.1,1.1,0.2))*-1,seq(.1,1.1,.2))
image(1:15,1:15,cor(gebvM),col=cs,breaks=brks,axes=FALSE,xlab="",ylab="")
image.plot(cor(gebvM),col=cs,breaks=brks,legend.only=TRUE,horizontal=FALSE)
axis(1,1:15,nms,las=2,cex.axis=.9)
axis(2,1:15,nms,las=2,cex.axis=.9)
box()
dev.off()

############ OLD #################
gebvM<-cbind(gebv_p1,gebv_p2,gebv_p3,gebv_p4,gebv_2p1,gebv_2p2,gebv_2p3,gebv_2p4,gebv_2tp1,gebv_2tp2,gebv_2tp3,gebv_2tp4,gebv_lp1,gebv_lp2,gebv_2lp1,gebv_2lp2,gebv_2lp3,gebv_2lp4)
nms<-c("e1_SI_Su","e1_SI_MW","e1_SI_FW","e1_SI_DT","e2_CA_Su","e2_CA_DT","e2_SI_Su","e2_SI_DT","e2_CA_Ad","e2_CA_Eg","e2_SI_Ad","e2_SI_Eg","e1_SI_LN","e1_SI_LP","e2_CA_LN","e2_CA_LP","e2_SI_LN","e2_SI_LP")
colnames(gebvM)<-nms
library(RColorBrewer)
library(fields)
cs<-brewer.pal(11,"BrBG")
pdf("cowpea_gmatrix.pdf",width=6,height=6)
par(mar=c(6,6,6,6))
brks<-c(rev(seq(.1,1.1,0.2))*-1,seq(.1,1.1,.2))
image(1:18,1:18,cor(gebvM),col=cs,breaks=brks,axes=FALSE,xlab="",ylab="")
image.plot(cor(gebvM),col=cs,breaks=brks,legend.only=TRUE,horizontal=FALSE)
axis(1,1:18,nms,las=2,cex.axis=.9)
axis(2,1:18,nms,las=2,cex.axis=.9)
box()
dev.off()

