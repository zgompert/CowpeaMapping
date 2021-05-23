## wasn't able to install qvalue on the cluster, ran on my computer
library(qvalue)
load("../Scratch/mapit.rdat")
load("../Scratch/mapit2.rdat")
L<-length(mapit_ep2[[1]]$pvalues)
qq<-matrix(NA,nrow=L,ncol=12)

for(i in 1:4){
    o<-qvalue(mapit_ep1[[i]]$pvalues)
    qq[,i]<-o$qvalues
    }
for(i in 1:8){    
    o<-qvalue(mapit_ep2[[i]]$pvalues)
    qq[,i+4]<-o$qvalues
    }

table(apply(qq < 0.01,1,sum))

#    0     1     2     3 
#31636   273    82   139 
## same answer for just DT
table(apply(qq[,c(4,6,8)] < 0.01,1,sum))


## keep those with 3 = 139 = (139*138)/2 = 9591
kepi<-which(apply(qq < 0.01,1,sum)==3)

## add epi interactions, do for each data set
load("../Scratch/mapit.rdat")
ge<-gdat[kepi,]
for(j in 1:139){
    ge[j,]<-(ge[j,]-mean(ge[j,]))/sd(ge[j,])
    }

Le<-9591
gee<-matrix(NA,nrow=Le,ncol=dim(ge)[2])
x<-1
for(i in 1:138){for(j in (i+1):139){
    gee[x,]<-ge[i,]*ge[j,]
    x<-x+1
    }}
    
#hdr<-data.frame(paste("epi",1:Le,sep=""),rep("A",Le),rep("T",Le))   
#o<-data.frame(hdr,gee)
o<-gee
colnames(o)<-1:293
colnames(gdat)<-1:293
o<-rbind(gdat,o)
write.table(file="exp1_epi_geno",o,row.names=FALSE,col.names=FALSE,quote=FALSE)

## add epi interactions, do for each data set
load("../Scratch/mapit2.rdat")
ge<-gdat[kepi,]
for(j in 1:139){
    ge[j,]<-(ge[j,]-mean(ge[j,]))/sd(ge[j,])
    }

Le<-9591
gee<-matrix(NA,nrow=Le,ncol=dim(ge)[2])
x<-1
for(i in 1:138){for(j in (i+1):139){
    gee[x,]<-ge[i,]*ge[j,]
    x<-x+1
    }}
    
#hdr<-data.frame(paste("epi",1:Le,sep=""),rep("A",Le),rep("T",Le))   
#o<-data.frame(hdr,gee)
o<-gee
colnames(o)<-1:276
colnames(gdat)<-1:276
o<-rbind(gdat,o)
write.table(file="exp2_epi_geno",o,row.names=FALSE,col.names=FALSE,quote=FALSE)

write.table(kepi,file="epiIndexes.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)





