###############################################################################################################
### File: 4_Fig1and2.R
###
### Created: 01/11/2018
### Author: Tanja Miethe
###
### R 3.5.1
################################################################################################################

rm(list = ls())

## set working directory here if needed
## setwd("")


###############################################################################################################
## function to determine Lmax5 at FSPR, Fmax and at F=0 analytically
###############################################################################################################

Lmax5FSPR <- function(Linf, M, k, b, Lmat, Lc, SPRref)  #Reference point at 0.4SPR
{
  lm <- Lmat/Linf          # stand. length at maturity
  lmls <- Lc/Linf          # stand. length at first capture

  tm <- -(log(1-lm))/k           # age at maturity
  tmls <- -(log(1-lmls))/k       # age at first capture
  tmax <- Inf                    # assuming infinite age for integral

  ## SPR=0.4 find F
  fa <- function(t)
  { ((1-exp(-k*t))^b)*(exp(-M*t))
  }  # ssb0
  ia <- integrate(fa,tm,tmax)
  lim1 <- ia[[1]] * SPRref

  if(lm>=lmls)
  { fr <- function(Fi)
    {fc <- function(t)
      {((1-exp(-k*t))^b)*(exp(Fi*tmls))*(exp(-(Fi+M)*t))
      }   # ssb fished

      i <-  integrate(fc,tm,tmax)[[1]]
      abs(i-lim1)
    }
  }

  if(lm<lmls)
  { fr <- function(Fi)
    { fb <- function(t)
      { ((1-exp(-k*t))^b)*(exp(-M*t))
      }
      fc <- function(t)
      { ((1-exp(-k*t))^b)*(exp(Fi*tmls))*(exp(-(Fi+M)*t))
      }   #ssb fished
      i <-  (integrate(fb,tm,tmls)[[1]]  + integrate(fc,tmls,tmax)[[1]] )
      return(abs(i - lim1))
    }
  }

  opt <- optimize(fr,lower=0,upper=7, tol=0.00001, maximum=F)   #find F (with upper limit F=2)
  Fi1 <- opt$minimum  # Fspr=40%

  ll1 <- 1-(0.05^(k/(M+Fi1)))*(1-lmls)   # lstar critical length with F for 5%
  tl1 <- -log(1-ll1)/k                   # tstar age
  lf1 <- (1-(exp(-k*tl1))*(Fi1+M)/(Fi1+M+k))*Linf  # reference point for mean length of top 5%

  ff<-function(t){Fi1*(exp(Fi1*tmls))*exp(-(Fi1+M)*t)*(1-exp(-(k*t)))^b } #yield
  if1<-integrate(ff,tmls,tmax)[[1]]

  return(c(Lmax5 = lf1, F = Fi1, yield=if1))
}


Lmax5Fmax <- function(Linf, M, k, b, Lmat, Lc){

  lm<-Lmat/Linf  # standardized length at maturity
  lmls<-l/Linf   # standardized at first capture

  tm<--(log(1-lm))/k           # age at maturity
  tmls<--(log(1-lmls))/k       # age at first capture
  tmax<-Inf                    # assuming infinite age for integral

  ## find Fmax
  fmax<-function(Fi){
    fd<-function(t){Fi*(exp(Fi*tmls))*exp(-(Fi+M)*t)*(1-exp(-(k*t)))^b }
    id<-integrate(fd,tmls,tmax)[[1]]
    abs(id)
  }
  opt<- optimize(fmax,lower=0,upper=7, tol=0.00001, maximum=T)  # find max
  Fimax<-opt$maximum

  ll2<-1-(0.05^(k/(M+Fimax)))*(1-lmls)   # lstar critical length with F
  tl2<--log(1-ll2)/k                     # tstar age
  lf2<-(1-(exp(-k*tl2))*(Fimax+M)/(Fimax+M+k))*Linf  # mean length of to 5%

  ff<-function(t){Fimax*(exp(Fimax*tmls))*exp(-(Fimax+M)*t)*(1-exp(-(k*t)))^b }
  if2<-integrate(ff,tmls,tmax)[[1]]

return(c(Lmax = lf2, F = Fimax, yield=if2))
  }


Lmax5Fzero <- function(Linf, M, k, b, Lmat, Lc){

  lmls<-Lc/Linf   # standardized at first capture
  ll0<-1-(0.05^(k/(M)))*(1-lmls)    # lstar critical length F=0
  tl0<--log(1-ll0)/k
  lf0<-(1-(exp(-k*tl0))*(M)/(M+k))*Linf

  return(c(Lmax5 = lf0, F = 0, yield=0))
  }


#############################################################################################################
## plots Figure 1 and 2
#############################################################################################################


p1<-read.csv("./parameters_manuscript.csv",sep=",")                            # read life history parameters

outDir <- "./Figures/"                 # Directory for plots.

dir.create(outDir)


species<-rbind(c("LBE","HER"), c("RJN","ISO"))                                 # species code in table
name2<-rbind(c("H. gammarus","C. harengus"), c("L. naevus","I. oxyrinchus"))   # species name as title
sexes<-rbind(c("male","male"),c("female","female"))                            # sexes

colours<-c("black","grey20","grey60")   # colours for plots



for(i in 1:2){

 tiff(paste(outDir,"Figure_",i,".tiff",sep=""), bg="white", pointsize=11, res=300, width = 170, height = 190, units="mm")
 layout(matrix(c(1,4,2,5,3,6),3,2,byrow=T))

for(j in 1:2){
 name<-name2[i,j]   # choose species
 sex<-sexes[i,j]


p<-p1[p1$species==species[i,j],]
Lmat<-p[p$parameter=="Lmat",sex]        # maturation size in mm
Linf<-p[p$parameter=="Linf",sex]        # aymptotic length
k<-p[p$parameter=="K",sex]              # growth parameter
b<-p[p$parameter=="b",sex]              # length weight relationship parameter

start<-1   # minimum length
ls<-c(start:(Linf-Linf*0.0001))         # length range of reference point calculations,end point just below Linf

M<-4.118*k^0.73*(Linf/10)^(-0.33)       # estimate of Natural mortality,  Linf in mm, Then at al (2015)

lengthfinal0<-NULL
lengthfinalSPR<-NULL
lengthfinalMAX<-NULL

yieldSPR<-NULL
yieldMAX<-NULL

FSPR<-NULL
FMAX<-NULL

for(l in ls){

  LFSPR<-Lmax5FSPR(Linf, M, k, b, Lmat, Lc = l, SPRref = 0.4)
  LFMAX<-Lmax5Fmax(Linf, M, k, b, Lmat, Lc = l)
  LFzero<-Lmax5Fzero(Linf, M, k, b, Lmat, Lc = l)

  lengthfinal0<-c(lengthfinal0,LFzero[[1]])      # F=0
  lengthfinalSPR<-c(lengthfinalSPR,LFSPR[[1]])   # FSPR=40%
  lengthfinalMAX<-c(lengthfinalMAX,LFMAX[[1]])   # Fmax

  ## F
  FSPR<-c(FSPR,LFSPR[[2]])
  FMAX<-c(FMAX,LFMAX[[2]])

  ## yield
  yieldSPR<-c(yieldSPR,LFSPR[[3]])
  yieldMAX<-c(yieldMAX,LFMAX[[3]])
}


## Reference points at different values Lc
lll0<-ceiling(lengthfinal0[which(round(ls)==round(Lmat))])      # F=0
lllSPR<-ceiling(lengthfinalSPR[which(round(ls)==round(Lmat))])  # FSPR=0.4
lllMAX<-ceiling(lengthfinalMAX[which(round(ls)==round(Lmat))])  # Fmax


par(mar=c(4,5,2,1),cex=0.9, family="serif", lwd=1.1)

 plot(ls,lengthfinal0, ylab=expression(paste(L["max5%"], "  (mm)")), xlab="", type="l",
      bty="l", lwd=0.8, col="transparent",
      main=substitute(italic(xx)~(ss),list(xx=name,ss=sex)), xlim=c(start,max(ls)),  yaxt="n",
      xaxt="n", ylim=c(0,Linf),lty="solid",cex.axis=1.2, cex.lab=1.3,cex.main=1.3)

if(i==1){axis(side=1, at=seq(0,max(ls),50), labels=seq(0,max(ls),50),cex.axis=1.2)
  axis(side=2, at=seq(0,max(ls),50), labels=seq(0,max(ls),50), cex.axis=1.2)}
if(i==2 & j==1){axis(side=1, at=seq(0,max(ls),200), labels=seq(0,max(ls),200),cex.axis=1.2)
  axis(side=2, at=seq(0,max(ls),200), labels=seq(0,max(ls),200), cex.axis=1.2)}
if(i==2 & j==2){axis(side=1, at=seq(0,max(ls),1000), labels=seq(0,max(ls),1000),cex.axis=1.2)
  axis(side=2, at=seq(0,max(ls),1000), labels=seq(0,max(ls),1000),cex.axis=1.2)}

lines(ls,rep(Linf, length(ls)),  type="l", lwd=0.6, col=colours[2],lty="dotted")
lines(rep(Lmat, length(ls)),ls,  type="l", lwd=0.6, col=colours[2],lty="dashed")
lines(ls,lengthfinal0,  type="l", lwd=0.8, col=colours[1])
lines(ls[which(FMAX<6)],lengthfinalMAX[which(FMAX<6)],  type="l", lwd=2, col=colours[3], lty="solid")
lines(ls[which(FMAX>=6)],lengthfinalMAX[which(FMAX>=6)],  type="l", lwd=2, col=colours[3], lty="dashed")
lines(ls[which(FSPR<6)],lengthfinalSPR[which(FSPR<6)],  type="l", lwd=2, col=colours[2], lty="solid")

legend("bottomleft", bty="n", legend=c(expression(L[infinity]),"F=0", expression(F["40%SPR"]),expression(F["MAX"])) ,
       col=c(colours[1],colours[1],colours[2],colours[3]),lty=c("dotted","solid", "solid", "solid"),
       lwd=c(0.6,0.8, 2, 2),ncol=1, y.intersp=1.3, cex=1.1)

mtext(expression(L["mat"]),side=1,0.4,cex=1,at=Lmat, col=colours[1])

v<-bquote(L["mat"]/L[infinity]==.(round(Lmat/Linf,2)))
text(Linf*0.98,((Linf)*0.33),v, pos=2, cex=1.2)
v1<-bquote(M/k==.(round(M/k,2)))
text(Linf*0.98,((Linf)*0.18),v1, pos=2, cex=1.2, lwd=0.3)
if(j==1){mtext("(a)", side=3,0.1,at=Linf,cex=1.2, col=colours[1])}else{
mtext("(d)", side=3,0.1,at=Linf,cex=1.3, col=colours[1])}

ls1<-seq(0,8,by=0.001)

 plot(ls,FSPR, ylab="F", xlab="", type="l",bty="l", lwd=2, col="transparent",
      xaxt="n", xlim=c(start,max(ls)),ylim=c(0,5), lty="solid",cex.axis=1.2, cex.lab=1.3)
if(i==1){axis(side=1, at=seq(0,max(ls),50), labels=seq(0,max(ls),50),cex.axis=1.2)}
if(i==2 & j==1){axis(side=1, at=seq(0,max(ls),200), labels=seq(0,max(ls),200),cex.axis=1.2)}
if(i==2 & j==2){axis(side=1, at=seq(0,max(ls),1000), labels=seq(0,max(ls),1000),cex.axis=1.2)}

lines(ls,FMAX,  type="l", lwd=2, col=colours[3],lty="solid")
lines(ls,FSPR,  type="l", lwd=2, col=colours[2],lty="solid")
lines(rep(Lmat, length(ls1)),ls1,  type="l", lwd=0.6, col=colours[2],lty="dashed")
if(j==1){mtext("(b)", side=3,0.1,at=Linf,cex=1.3, col=colours[1])}else{
  mtext("(e)", side=3,0.1,at=Linf,cex=1.3, col=colours[1])
}
legend("topleft", bty="n", legend=c(expression(F["40%SPR"]),expression(F["MAX"])) ,
                     col=c(colours[2],colours[3]),lty=c("solid", "solid"), lwd=c(2,2),ncol=1, cex=1.1, y.intersp=1.3)
mtext(expression(L["mat"]), side=1,0.4, cex=1, at=Lmat, col=colours[1])

 plot(ls,yieldMAX, ylab="Yield per recruit",xlab=expression(paste(L["c"], " (mm)")) ,
      yaxt="n",xaxt="n",  type="l",bty="l", lwd=2, col="transparent",  xlim=c(start,max(ls)),
      lty="solid",cex.axis=1.1, cex.lab=1.3)   # "survival", "weight","mat", "sel"
if(i==1 & j==1){axis(side=1, at=seq(0,max(ls),50), labels=seq(0,max(ls),50),cex.axis=1.2)
  axis(side=2, at=seq(0,max(yieldMAX*1.01),0.004), labels=seq(0,max(yieldMAX*1.01),0.004), cex.axis=1.2)}
if(i==1 & j==2){axis(side=1, at=seq(0,max(ls),50), labels=seq(0,max(ls),50),cex.axis=1.2)
  axis(side=2, at=seq(0,max(yieldMAX*1.01),0.01), labels=seq(0,max(yieldMAX*1.01),0.01), cex.axis=1.2)}
if(i==2 & j==1){axis(side=1, at=seq(0,max(ls),200), labels=seq(0,max(ls),200),cex.axis=1.2)
  axis(side=2, at=seq(0,max(yieldMAX*1.01),0.01), labels=seq(0,max(yieldMAX*1.01),0.01), cex.axis=1.2)}
if(i==2 & j==2){axis(side=1, at=seq(0,max(ls),1000), labels=seq(0,max(ls),1000),cex.axis=1.2)
  axis(side=2, at=seq(0,max(yieldMAX*1.01),0.02), labels=seq(0,max(yieldMAX*1.01),0.02),cex.axis=1.2)}

lines(ls,yieldMAX, type="l", lwd=2, col=colours[3],lty="solid")
lines(ls,yieldSPR, type="l", lwd=2, col=colours[2],lty="solid")
lines(rep(Lmat, length(ls1)),ls1,  type="l", lwd=0.6, col=colours[2],lty="dashed")
mtext(expression(L["mat"]),side=1,0.4,cex=1,at=Lmat, col=colours[1])
if(j==1){mtext("(c)", side=3,0.1,at=Linf,cex=1.3, col=colours[1])}else{
  mtext("(f)", side=3,0.1, at=Linf, cex=1.3, col=colours[1])
}

}
dev.off()

}
