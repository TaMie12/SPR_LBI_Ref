#########################################################################################################################
### File: SPE_LBI_Ref_fun1.R
### Created: 16/08/2021
### Author: Tanja Miethe
###
#########################################################################################################################
### Description:
###   Functions to calculate SPR-based reference points for LBIs: Lmax5% and Lmean based on equations from Miethe et
###   al. (2019) using CV=0, knife-edge selectivity and maturity ogive
###
### Miethe, T., Reecht, Y., and Dobby, H. 2019. Reference points for length-based indicator Lmax5% to support
### assessment of data-limited stocks and fisheries.  ICES Journal of Marine Science, 76:7, 2125-2139.
#########################################################################################################################


##' @title Lmax5_Lmean_LBSPR1: Lmax<5>% and Lmean reference point estimation using LBSPR.
##' @param k A length-one numeric vector for VBGF k parameter.
##' @param M A length-one numeric vector for natural mortality M.
##' @param Linf A length-one numeric vector for the assymptotic length.
##' @param Lm50 A length-one numeric vector for Length 50% maturity.
##' @param Ls50 A length-one numeric vector for Selectivity length 50%.
##' @param SPRref A length-one numeric vector for MSY proxy for SPR (default 0.4).
##' @param Wbeta A length-one numeric vector for Weight-length relationship beta param (W = a * L^b).
##' @return Lmax<5>% and Lmean value estimated using LBSPR, unfished values (UF), F/M, F at 40%SPR, M/k, SPR.
##' @author Tanja Miethe



Lmax5_Lmean_LBSPR1 <- function(Linf, M, k, Wbeta, Lm50, Ls50, SPRref=0.4)  #Reference point at 0.4SPR
{
  lm <- Lm50/Linf            # stand. length at maturity
  lmls <- Ls50/Linf          # stand. length at first capture

  tm <- -(log(1-lm))/k           # age at maturity
  tmls <- -(log(1-lmls))/k       # age at first capture
  tmax <- Inf                    # assuming max infinite age for integral

  ## SPR=0.4 find F, fished
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

  opt <- optimize(fr,lower=0, upper=7, tol=0.00001, maximum=F)   #find F (with upper limit F=7)
  Fi1 <- opt$minimum  # F(spr=40%)

  ll1 <- 1-(0.05^(k/(M+Fi1)))*(1-lmls)   # lstar critical length with F for 5% largest ind
  tl1 <- -log(1-ll1)/k                   # tstar age
  lf1 <- (1-(exp(-k*tl1))*(Fi1+M)/(Fi1+M+k))*Linf    # reference point for mean length of top 5%
  lf2 <- (1-(exp(-k*tmls))*(Fi1+M)/(Fi1+M+k))*Linf   # reference point for mean length above Lc

  ff<-function(t){Fi1*(exp(Fi1*tmls))*exp(-(Fi1+M)*t)*(1-exp(-(k*t)))^b } #yield
  if1<-integrate(ff,tmls,tmax)[[1]]

  ## unfished
  lmls<-Ls50/Linf   # standardized at first capture
  tmls <- -(log(1-lmls))/k
  ll0<-1-(0.05^(k/(M)))*(1-lmls)    # lstar critical length F=0
  tl0<--log(1-ll0)/k
  lf0<-(1-(exp(-k*tl0))*(M)/(M+k))*Linf
  lf3<-(1-(exp(-k*tmls))*(M)/(M+k))*Linf # reference point for mean length above Ls50

  return(round(c(Lmax5RP = lf1, Lmax5UF = lf0 , LmeanRP=lf2, LmeanUF=lf3, F_M=Fi1/M, F = Fi1, M_k=M/k, SPR=SPRref ),3))
}

### End:
