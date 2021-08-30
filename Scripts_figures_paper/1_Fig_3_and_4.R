
###############################################################################################################
### File: 1_Fig_3_and_4.R
###
### Created: 25/03/2018
### Author: Tanja Miethe
###
### R 3.5.1
################################################################################################################

rm(list = ls())

## set working directory here if needed

## setwd("")



####################################################################################################################

source("./0_Functions_Lmax5refpointLBSPR.R")

library(R.utils)
library(LBSPR)

p1<-read.csv("./parameters_manuscript.csv",sep=",")                            # read life history parameters

outDir <- "./Figures/"                 # Directory for plots.

dir.create(outDir)


species<-rbind(c("LBE","HER"), c("RJN","ISO"))                                 # species code in table
name2<-rbind(c("H. gammarus","C. harengus"), c("L. naevus","I. oxyrinchus"))   # species name as title
sexes<-rbind(c("male","male"),c("female","female"))                            # sexes

colours<-c("black","grey20","grey60")   # colours for plots



for(i in 1:2){
  ii<-i+2
  tiff(paste(outDir,"Figure_",ii,".tiff",sep=""), bg="white", pointsize=11, res=300, width = 170, height = 140, units="mm")
  layout(matrix(c(1,3,2,4),2,2,byrow=T))

  for(j in 1:2){

    name<-name2[i,j]   # choose species
    sex<-sexes[i,j]

    p<-p1[p1$species==species[i,j],]
    Lmat<-p[p$parameter=="Lmat",sex]        # maturation size in mm
    Linf<-p[p$parameter=="Linf",sex]        # asymptotic length
    k<-p[p$parameter=="K",sex]              # growth parameter
    b<-p[p$parameter=="b",sex]              # length weight relationship parameter

    start<-1   # minimum length
    if(i==1)ls<-seq(start,(Linf*0.99),2)         # length range of reference point calculations,end point just below Linf
    if(i==2)ls<-seq(start,(Linf*0.99),10)
    if(i==2 & j==2)ls<-seq(start,(Linf*0.99),30)

    M<-4.118*k^0.73*(Linf/10)^(-0.33)       # estimate of natural mortality,  Linf in mm, Then at al (2015)

    lengthfinal0<-NULL
    lengthfinal01<-NULL
    lengthfinalSPR<-NULL
    lengthfinalSPR1<-NULL

    FSPR<-NULL
    FSPR1<-NULL


        for(l in ls){

          ## LB-SPR implementation

          LFSPR<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.0001, BinWidth=1,
                                    Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=l, Ls95=l*1.00001, SPRref = 0.4)
          LFSPR1<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.1, BinWidth=1,
                                     Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=l, Ls95=l*1.00001, SPRref = 0.4)

          LFzero<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.0001, BinWidth=1,
                                     Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=l, Ls95=l*1.00001, SPRref = NULL, FMref=0)
          LFzero1<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.1, BinWidth=1,
                                      Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=l, Ls95=l*1.00001, SPRref = NULL, FMref=0)

          lengthfinal0<-c(lengthfinal0,LFzero[[1]][1])                  # F=0
          lengthfinal01<-c(lengthfinal01,LFzero1[[1]][1])               # F=0


          lengthfinalSPR<-c(lengthfinalSPR,LFSPR[[1]][1])        # FSPR=40%
          lengthfinalSPR1<-c(lengthfinalSPR1,LFSPR1[[1]][1])      # FSPR=40%
          #F
          FSPR<-c(FSPR,LFSPR[[1]][4])
          FSPR1<-c(FSPR1,LFSPR1[[1]][4])

        }


        # Reference points at different values Lc
        lll0<-ceiling(lengthfinal0[which(round(ls)==round(Lmat))])      # F=0
        lllSPR<-ceiling(lengthfinalSPR[which(round(ls)==round(Lmat))])  # FSPR=0.4

        ls1<-seq(1,Linf*1.2,by=10)

        par(mar=c(4,5,2,1),cex=0.9, family="serif", lwd=1.1)

      plot(ls,lengthfinal0, ylab=expression(paste(L["max5%"], "  (mm)")),
           xlab=expression(paste(L["c"], "  (mm)")), type="l",bty="l", lwd=0.8, col="transparent",
           main=substitute(italic(xx)~(ss),list(xx=name,ss=sex)), xlim=c(start,Linf*1.2),
           yaxt="n", xaxt="n",ylim=c(0,Linf*1.2),lty="solid",cex.axis=1.2, cex.lab=1.3,cex.main=1.2)


          if(i==1 ){axis(side=1, at=seq(0,max(ls1),50), labels=seq(0,max(ls1),50),cex.axis=1.2)
            axis(side=2, at=seq(0,max(ls1),50), labels=seq(0,max(ls1),50),cex.axis=1.2)}
          if(i==2 & j==1){axis(side=1, at=seq(0,max(ls1),200), labels=seq(0,max(ls1),200),cex.axis=1.2)
            axis(side=2, at=seq(0,max(ls1),200), labels=seq(0,max(ls1),200),cex.axis=1.2)}
          if(i==2 & j==2){axis(side=1, at=seq(0,max(ls1),1000), labels=seq(0,max(ls1),1000),cex.axis=1.2)
            axis(side=2, at=seq(0,max(ls1),1000), labels=seq(0,max(ls1),1000),cex.axis=1.2)}

        lines(ls1,rep(Linf, length(ls1)),  type="l", lwd=0.6, col=colours[2],lty="dotted")
        lines(rep(Lmat, length(ls1)),ls1,  type="l", lwd=0.6, col=colours[2],lty="dashed")
        lines(ls,lengthfinal0,  type="l", lwd=0.8, col=colours[1])
        lines(ls,lengthfinal01,  type="l", lwd=0.8, col=colours[3])
        lines(ls[which(FSPR<=2)],lengthfinalSPR[which(FSPR<=2)],  type="l", lwd=2, col=colours[1], lty="solid")
        lines(ls[which(FSPR1<=2)],lengthfinalSPR1[which(FSPR1<=2)],  type="l", lwd=2, col=colours[3], lty="solid")
        lines(ls[which(FSPR>=2)],lengthfinalSPR[which(FSPR>=2)],  type="l", lwd=2, col=colours[1], lty="dashed")
        lines(ls[which(FSPR1>=2)],lengthfinalSPR1[which(FSPR1>=2)],  type="l", lwd=2, col=colours[3], lty="dashed")
        mtext(expression(L["mat"]),side=1,0.4,cex=1,at=Lmat, col=colours[1])
        if(j==1){mtext("(a)", side=3,0.1,at=Linf*1.2,cex=1.2, col=colours[1])}else{
          mtext("(c)", side=3,0.1,at=Linf*1.2,cex=1.2, col=colours[1])}

        if(j==2){legend("bottomleft", bty="n", legend=c(expression(L[infinity]),"F=0, CV=0","F=0, CV=0.1", expression(paste(F["40%SPR"],", CV=0")), expression(paste(F["40%SPR"],", CV=0.1"))) ,
               col=c(colours[1],colours[1],colours[3],colours[1],colours[3]),lty=c("dotted","solid", "solid","solid","solid"), lwd=c(0.6,0.8, 0.8, 2,2),ncol=1, y.intersp=1.3, cex=1.0)
          }

        LFSPR<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.0001, BinWidth=1, Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=Lmat, Ls95=Lmat*1.00001, SPRref = 0.4)
        LFSPR1<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.1, BinWidth=1, Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=Lmat, Ls95=Lmat*1.00001, SPRref = 0.4)

        LFzero<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.0001, BinWidth=1, Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=Lmat, Ls95=Lmat*1.00001, SPRref = NULL, FMref=0)
        LFzero1<-LmaxPrefpointLBSPR( M, k, Linf=Linf, Wbeta=b, FecB=b, CVLinf=0.1, BinWidth=1, Lm50=Lmat, Lm95=Lmat*1.00001, Ls50=Lmat, Ls95=Lmat*1.00001, SPRref = NULL, FMref=0)


       if(i==2 & j==2){ plot(LFSPR$simSizeDist@LMids,LFSPR$simSizeDist@pLCatch, ylab="Frequency", xlab="Length (mm)", type="l",bty="l", lwd=2, col="transparent",yaxt="n", xaxt="n",
             xlim=c(start,Linf*1.2),  ylim=c(0,1.3),lty="solid",cex.axis=1.1,cex.lab=1.3)}else{plot(LFSPR$simSizeDist@LMids,LFSPR$simSizeDist@pLCatch, ylab="Frequency", xlab="Length (mm)", type="l",bty="l", lwd=2, col="transparent",yaxt="n", xaxt="n",
                                                                                        xlim=c(start,Linf*1.2),  ylim=c(0,1.05),lty="solid", cex.axis=1.1, cex.lab=1.3)}
        if(i==1 & j==1){axis(side=1, at=seq(0,max(ls1),50), labels=seq(0,max(ls1),50),cex.axis=1.2)
          axis(side=2, at=seq(0,1,0.2), labels=seq(0,1,0.2),cex.axis=1.2)}
        if(i==1 & j==2){axis(side=1, at=seq(0,max(ls1),50), labels=seq(0,max(ls1),50),cex.axis=1.2)
          axis(side=2, at=seq(0,1,0.2), labels=seq(0,1,0.2),cex.axis=1.2)}
        if(i==2 & j==1){axis(side=1, at=seq(0,max(ls1),200), labels=seq(0,max(ls1),200),cex.axis=1.2)
          axis(side=2, at=seq(0,1.3,0.2), labels=seq(0,1.3,0.2),cex.axis=1.2)}
        if(i==2 & j==2){axis(side=1, at=seq(0,max(ls1),1000), labels=seq(0,max(ls1),1000),cex.axis=1.2)
          axis(side=2, at=seq(0,1.3,0.2), labels=seq(0,1.3,0.2),cex.axis=1.2)}

        #lines(rep(Linf, 11) ,seq(0,1,0.1), type="l", lwd=0.4, col=colours[2],lty="dotted")
        #lines(rep(Lmat, 11),seq(0,1,0.1),  type="l", lwd=0.4, col=colours[2],lty="dashed")
        lines(LFzero1$simSizeDist@LMids,LFzero1$simSizeDist@pLCatch/LFzero1$simSizeDist@pLCatch[LFzero1$simSizeDist@LMids==Lmat+0.5],  type="l", lwd=0.8, col=colours[3])
        lines(LFzero$simSizeDist@LMids,LFzero$simSizeDist@pLCatch/LFzero1$simSizeDist@pLCatch[LFzero1$simSizeDist@LMids==Lmat+0.5],  type="l", lwd=0.8, col=colours[1])

       # lines(LFSPR1$simSizeDist@LMids,LFSPR1$simSizeDist@pLCatch/max(LFSPR1$simSizeDist@pLCatch),  type="l", lwd=1, col=colours[3], lty="solid")
       # lines(LFSPR$simSizeDist@LMids,LFSPR$simSizeDist@pLCatch/max(LFSPR$simSizeDist@pLCatch),  type="l", lwd=1, col=colours[1], lty="solid")

        mtext(expression(L["mat"]),side=1,0.4,cex=1,at=Lmat, col=colours[1])
        if(j==1){mtext("(b)", side=3,0.1,at=Linf*1.2,cex=1.2, col=colours[1])}else{
          mtext("(d)", side=3,0.1,at=Linf*1.2,cex=1.2, col=colours[1])
        }

       if(j==2){ legend("topleft", bty="n", legend=c("F=0, CV=0", "F=0, CV=0.1") ,
               col=c(colours[1],colours[3]),lty=c("solid", "solid"), lwd=c(0.8, 0.8),ncol=1, y.intersp=1.3, cex=1)}



  }

  dev.off()

}

# example refererence points for cukoo ray (WKSHARK4)  for M/k=1.5 and Lmat=590, Linf 731, beta=3.1396 as used in WKSHARK4
mean_L<- mean(c(0.93,0.96,0.96,0.96,0.92,0.96,0.95)) # mean 2010-2016
mean_Lc_Lmat<-mean(c(0.39,0.88,0.34,0.41,0.43,0.56,0.69))

LFSPR_Cuckoray<-LmaxPrefpointLBSPR(M=0.29, k=0.29/1.5, Linf=731, Wbeta=3.1396, FecB=3.1396, CVLinf=0.1, BinWidth=1, Lm50=590, Lm95=590*1.00001, Ls50=590*mean_Lc_Lmat, Ls95=590*mean_Lc_Lmat*1.00001, SPRref = 0.4)

LF<-LFSPR_Cuckoray[[1]][1]/731



