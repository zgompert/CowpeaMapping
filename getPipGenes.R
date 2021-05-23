## get gene annotations for top PIP genes
#### load annotation workspace
load("../../annotation/annotation.rdat")
## out is the annotation, includes columns for gene id and protein

## cowpea survival and DT from experiment 2 
pf<-list.files(pattern="pip_o_cowpea_fit_ph")
L<-32130
pips<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips[,j]<-scan(pf[j])
    }
    
## 14 with pips > 0.1
out[which(apply(pips,1,max) > 0.1),]
#          SNP OldChr   OldCm NewChr  NewPos           Gene
#13453 2_01390      5 13.7582      5 3829851 Vigun05g046000
#13469 2_07583      5 13.7582      5 3966407 Vigun05g047400
#13473 2_06553      5 13.7582      5 3994844           <NA>
#13545 2_00867      5 15.1234      5 4488638 Vigun05g052500
#13548 2_24797      5 15.1234      5 4496586 Vigun05g052600
#13549 2_00868      5 15.1234      5 4497448 Vigun05g052600
#13553 2_08113      5 15.1234      5 4517159 Vigun05g052800
#13559 2_51427      5 15.1234      5 4551502 Vigun05g053000
#13564 2_13010      5 15.1234      5 4558350 Vigun05g053200
#13579 2_08171      5 15.1234      5 4630314 Vigun05g054000
#13655 2_34617      5 18.5364      5 5156508 Vigun05g060200
#13657 2_30291      5 18.5364      5 5165555           <NA>
#13659 2_20595      5 18.5364      5 5167944 Vigun05g060400
#13661 2_24125      5 18.5364      5 5172115 Vigun05g060500

#                                                                        prot
#13453        Leucine-rich repeat receptor-like protein kinase family protein
#13469                          ethylene-responsive element binding factor 13
#13473                                                                   <NA>
#13545                                                     purine permease 11
#13548 signal recognition particle 19 kDa protein, putative / SRP19, putative
#13549 signal recognition particle 19 kDa protein, putative / SRP19, putative
#13553                               RNA-binding KH domain-containing protein
#13559                                  sister chromatid cohesion 1 protein 4
#13564                                       2-phosphoglycolate phosphatase 1
#13579                                                         calreticulin 3
#13655                                    Late embryogenesis abundant protein
#13657                                                                   <NA>
#13659                                   Diacylglycerol kinase family protein
#13661                              alpha/beta-Hydrolases superfamily protein

## 2 with pips > 0.5
#          SNP OldChr   OldCm NewChr  NewPos           Gene
#13559 2_51427      5 15.1234      5 4551502 Vigun05g053000
#13564 2_13010      5 15.1234      5 4558350 Vigun05g053200
                                       prot
#13559 sister chromatid cohesion 1 protein 4
#13564      2-phosphoglycolate phosphatase 1


############# experiment 2 totals
pf<-list.files(pattern="pip_o_cowpea_totals_fit_ph")
L<-32130
pips<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips[,j]<-scan(pf[j])
    }
tits<-c("(A) No. adults, CA","(B) No. eggs, CA","(C) No. adults, SI","(D) No. eggs, SI")
## 7 with pip > 0.1, none > 0.5
out[which(apply(pips,1,max) > 0.1),]
#          SNP OldChr   OldCm NewChr   NewPos           Gene
#13411 2_53900      5 11.4385      5  3089275 Vigun05g038400
#13453 2_01390      5 13.7582      5  3829851 Vigun05g046000
#13469 2_07583      5 13.7582      5  3966407 Vigun05g047400
#13473 2_06553      5 13.7582      5  3994844           <NA>
#23633 2_14488      8 59.5472      8 34544843 Vigun08g175400
#23636 2_12166      8 59.5472      8 34550091 Vigun08g175500
#23637 2_12162      8 59.5472      8 34553766           <NA>
#                                                                 prot
#13411                                          SNF1 kinase homolog 10
#13453 Leucine-rich repeat receptor-like protein kinase family protein
#13469                   ethylene-responsive element binding factor 13
#13473                                                            <NA>
#23633                                                            <NA>
#23636                                                            <NA>
#23637 


############# larv early exit
pf<-c(list.files(pattern="pip_o_cowpea1_larv_fit_ph"),list.files(pattern="pip_o_cowpea_larv_fit_ph"))
L<-32130
pips<-matrix(NA,nrow=L,ncol=6)
for(j in 1:6){
    pips[,j]<-scan(pf[j])
    }
tits<-c("(A) No. early, SI1","(B) Prop. early, SI1","(C) No. early, CA","(D) Prop. early, CA","(E) No. early, SI","(F) Prop. early, SI")

## only prop. earl exit
out[which(apply(pips[,c(2,4,6)],1,max) > 0.1),]

## 14 pip > 0.1
#          SNP OldChr   OldCm NewChr   NewPos           Gene
#6956  2_29933      3 36.5060      3 12436989           <NA>
#6960  2_42431      3 36.5060      3 12452088 Vigun03g128800
#7121  2_27843      3 36.8662      3 13405896 Vigun03g136900
#16241 2_02155      6  2.3486      6 12737996 Vigun06g028900
#16280 2_40955      6  2.7266      6 13518774           <NA>
#19612 2_02375      7 43.5156      7 22014760 Vigun07g118600
#21223 2_02849      7 96.6338      7 39270470 Vigun07g277700
#23596 2_01925      8 57.8903      8 34260474 Vigun08g171700
#23602 2_07173      8 57.8903      8 34298139 Vigun08g172000
#23604 2_55032      8 57.8903      8 34304994 Vigun08g172100
#23664 2_28941      8 62.4202      8 34887678 Vigun08g179000
#23689  1_0226      8 62.4202      8 35020285 Vigun08g180500
#23692 2_33467      8 62.4202      8 35038005 Vigun08g180800
#24016 2_09995      8 72.4500      8 37207239 Vigun08g209500
#                                                 prot
#6956                                             <NA>
#6960                                             <NA>
#7121                                             <NA>
#16241                                            <NA>
#16280                                            <NA>
#19612 SPX (SYG1/Pho81/XPR1) domain-containing protein
#21223                       hydroxypyruvate reductase
#23596                                            <NA>
#23602                                            <NA>
#23604                                            <NA>
#23664                                            <NA>
#23689                                            <NA>
#23692                                            <NA>
#24016                                            <NA>

## one > 0.5
out[which(apply(pips[,c(2,4,6)],1,max) > 0.5),]
#          SNP OldChr   OldCm NewChr   NewPos           Gene prot
#23596 2_01925      8 57.8903      8 34260474 Vigun08g171700 <NA>


############## experiment 1
pf<-list.files(pattern="pip_o_cowpea1_fit_ph")
L<-32130
pips<-matrix(NA,nrow=L,ncol=4)
for(j in 1:4){
    pips[,j]<-scan(pf[j])
    }
    

tits<-c("(A) Survival","(B) Male Wgt.","(C) Female Wgt.","(D) DT")
## 14 > 0.1, none > 0.5
out[which(apply(pips,1,max) > 0.1),]
#          SNP OldChr   OldCm NewChr   NewPos           Gene
#13545 2_00867      5 15.1234      5  4488638 Vigun05g052500
#13548 2_24797      5 15.1234      5  4496586 Vigun05g052600
#13549 2_00868      5 15.1234      5  4497448 Vigun05g052600
#13554 2_01320      5 15.1234      5  4522988 Vigun05g052900
#13568 2_19807      5 15.1234      5  4592401 Vigun05g053500
#13579 2_08171      5 15.1234      5  4630314 Vigun05g054000
#13586 2_16833      5 15.1234      5  4661053 Vigun05g054500
#13593 2_41936      5 16.0963      5  4686167           <NA>
#13596 2_22350      5 16.0963      5  4700276 Vigun05g054900
#13617 2_03774      5 18.2410      5  4861405 Vigun05g056800
#13624 2_17212      5 18.2410      5  4900723           <NA>
#13627 2_41956      5 18.2410      5  4911889 Vigun05g057300
#13652 2_07432      5 18.5364      5  5116382 Vigun05g059700
#24105 2_54350      8 75.8805      8 37701772 Vigun08g216600
#                                                                        prot
#13545                                                     purine permease 11
#13548 signal recognition particle 19 kDa protein, putative / SRP19, putative
#13549 signal recognition particle 19 kDa protein, putative / SRP19, putative
#13554                                                     LETM1-like protein
#13568                         Integrase-type DNA-binding superfamily protein
#13579                                                         calreticulin 3
#13586                           chloroplastic acetylcoenzyme A carboxylase 1
#13593                                                                   <NA>
#13596                            Eukaryotic aspartyl protease family protein
#13617                                   Protein of unknown function (DUF707)
#13624                                                                   <NA>
#13627                   Protein kinase family protein with ARM repeat domain
#13652                               jasmonic acid carboxyl methyltransferase
#24105                                                HEAT SHOCK PROTEIN 89.1

