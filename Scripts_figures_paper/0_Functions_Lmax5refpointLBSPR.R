#-*- coding: latin-1 -*-

### File: Lmax5refpointLBSPR.R
### Time-stamp: <2019-01-17 08:32:07 yreecht>
###
### Created: 04/10/2017	13:40:28
### Author: Yves Reecht, Tanja Miethe
###
####################################################################################################
### Description:
###
###
####################################################################################################

##'
##' @title Lmax<P>% estimation from abundances/proportions at length.
##' @param Lmids A numeric vector of mid-values of length-bins.
##' @param NbyBin A numeric vector of abundance/proportion by length-bin.
##' @return Lmax<P>% value estimated from abundances/proportions at length (length-one numeric vector).
##' @author Yves Reecht
LmaxPcalc1 <- function(Lmids, NbyBin, P = 0.05)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 24 Oct 2017, 09:47

    if (P[1] > 1 & P[1] <= 100)
    {
        P[1] <- P[1] / 100
    }

    stopifnot(length(NbyBin) == length(Lmids),
              is.numeric(Lmids),
              is.numeric(NbyBin),
              P[1] >= 0,
              P[1] <= 100)

    ## Making sure it is ordered by size:
    NbyBin <- NbyBin[order(Lmids)]
    Lmids <- sort(Lmids)

    ## Cumulative proportion of individuals represented by bin (including this one):
    cumProp <- cumsum(NbyBin) / sum(NbyBin)
    idx <- which(cumProp >= (1 - P[1]))[1]

    ## If the prop is over 100 - P%, some of the individuals in the first bin must be kept
    ##   (0 otherwise):
    if (cumProp[idx] > (1 - P[1]))
    {
        ## Number of "individuals" from the first bin to keep:
        NbyBin[idx] <- (cumProp[idx] - (1 - P[1])) * sum(NbyBin)
    }else{
        ## Previous block could be also used, but adds some rounding error when exact.
        idx <- idx + 1
    }

    ## Weighted mean on the bin mid point values:
    LmaxP <- weighted.mean(x = Lmids[idx:length(NbyBin)], w = NbyBin[idx:length(NbyBin)])
    return(LmaxP)
}

##'
##' @title Lmax<P>% estimation from abundances/proportions at length.
##' @param Lmids A numeric vector of mid-values of length-bins.
##' @param NbyBin A numeric vector of abundance/proportion by length-bin.
##' @return Lmax<P>% value estimated from abundances/proportions at length (length-one numeric vector).
##' @author Yves Reecht
LmaxPcalc <- function(Lmids, NbyBin, P = 0.05)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 24 Oct 2017, 10:00

    sapply(P,
           function(x, Lmids, NbyBin)
    {
        LmaxPcalc1(P = x, Lmids = Lmids, NbyBin = NbyBin)
    }, Lmids = Lmids, NbyBin = NbyBin)
}


##'
##' @title Lmax5% estimation from abundances/proportions at length.
##' @param Lmids A numeric vector of mid-values of length-bins.
##' @param NbyBin A numeric vector of abundance/proportion by length-bin.
##' @return Lmax5% value estimated from abundances/proportions at length (length-one numeric vector).
##' @author Yves Reecht
Lmax5calc <- function(Lmids, NbyBin)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  5 Oct 2017, 12:00

    return(LmaxPcalc(Lmids = Lmids, NbyBin = NbyBin, P = 0.05))
}

##' @title Lmax5refpointLBSPR: Lmax<P>% reference point from Length-Based Spawning Potential Ratio.
##' @param M A length-one numeric vector for Natural mortality.
##' @param k A length-one numeric vector for VBGF k parameter.
##' @param MK A length-one numeric vector for the M/k ratio.
##' one of MK or M and k must be specified.
##' @param Linf A length-one numeric vector for the assymptotic length.
##' @param Lm50 A length-one numeric vector for Length 50% maturity.
##' @param Lm95 A length-one numeric vector for Length 95% maturity.
##' @param Ls5 A length-one numeric vector for Selectivity length 5%.
##' One of Ls5 and Ls50 must be provided.
##' @param Ls50 A length-one numeric vector for Selectivity length 50%.
##' @param Ls95 A length-one numeric vector for Selectivity length 95%.
##' @param P A length-one numeric vector for the proportion of largest individuals used.
##' Default 0.05 (Lmax5%).
##' @param SPRref A length-one numeric vector for MSY proxy for SPR (default 0.4).
##' @param FMref A length-one numeric vector for MSY proxy for F/M (override SPRref).
##' @param Mpow A length-one numeric vector for Size varying mortality parameter (default = 0, constant mortality).
##' @param Wbeta A length-one numeric vector for Weight-length relationship beta param (W = a * L^b).
##' @param Walpha A length-one numeric vector for Weight-length relationship alpha param.
##' @param FecB A length-one numeric vector for beta parameter of length-fecundity relationship
##' @param R0
##' @param CVLinf A length-one numeric vector for Coefficient of variation on Linf.
##' @param Steepness A length-one numeric vector for steepness of SPR.
##' @param MLL A length-one numeric vector for minimum landing size.
##' @param sdLegal A length-one numeric vector for the scale parameter of on-board selectivity
##' (Sigmoid centered on MLL. !!! not a sd: scale = sd * sqrt(3) / pi).
##' @param fDisc A length-one numeric vector for mortality of discards.
##' @param BinWidth A length-one numeric vector for width of the size bins.
##' @param Control Controls for the LBSPRsim function (see LBSPR package help).
##' @return List with:
##'    $res: names vector of main outputs: LmaxX% refpoint, LmaxX% with no exploitation,
##'                                        F/M, F (if M available), SPR.
##'    $simSizeDist: LBSPR simulation object (e.g. for plotting purpose).
##' @author Yves Reecht
LmaxPrefpointLBSPR <- function(M, k, MK = M/k,
                               Linf, Lm50, Lm95,
                               Ls5, Ls50, Ls95,
                               P = 0.05,
                               SPRref = 0.4,
                               FMref = NULL,
                               Mpow = 0, Wbeta = 3, Walpha = 1e-04,
                               FecB = 3, R0 = 10000,
                               CVLinf = 0.1, Steepness = 0.7,
                               MLL = numeric(0), sdLegal = numeric(0),
                               fDisc = numeric(0),
                               BinWidth = numeric(0),
                               Control = list(modtype = "GTG", ngtg = 101))
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 24 Oct 2017, 10:27

    require(LBSPR)

    ## Parameters initialisation:
    suppressMessages(pars <- new("LB_pars"))
    if (missing(MK) || ! is.numeric(MK))
    {
        if (missing(M) || (! is.numeric(M)) ||
            missing(k) || (! is.numeric(k)))
        {
            stop("One of MK or (M and k) must be given")
        }else{
            pars@MK <- M/k
        }
    }else{
        pars@MK <- MK
    }

    if (! missing(M))
    {
        pars@M <- M
    }
    pars@Linf <- Linf
    if (missing(Ls50))
    {
        if (missing(Ls5)) stop("One of Ls5 and Ls50 must be given.")
        pars@SL50 <- (Ls95 + Ls5) / 2
    }else{
        pars@SL50 <- Ls50
    }
    pars@SL95 <- Ls95
    pars@L50 <- Lm50
    pars@L95 <- Lm95
    pars@Mpow <- Mpow
    pars@Wbeta <- Wbeta
    pars@Walpha <- Walpha
    pars@FecB <- FecB
    pars@CVLinf <- CVLinf
    pars@Steepness <- Steepness
    pars@R0 <- R0
    pars@MLL <- MLL
    pars@sdLegal <- sdLegal
    pars@fDisc <- fDisc

    pars@BinMin <- 0
    pars@BinMax <- ifelse(pars@MK > 0.5, 1.3, 1.8) * pars@Linf

    ## Default bin width:
    if (missing(BinWidth) || length(BinWidth) == 0)
    {
        pars@BinWidth <- (pars@Linf - pars@BinMin) / 100
    }else{
        pars@BinWidth <- BinWidth
    }

    ## Force FM if provided, and ignore SPRref (only if FM refpoint calculated externaly):
    if (missing(FMref) || is.null(FMref) || ! is.numeric(FMref))
    {
        pars@SPR <- SPRref
    }else{
        pars@FM <- FMref
    }

    suppressMessages(simSizeDist <- LBSPRsim(pars, Control = Control, verbose = FALSE))

    ## Unfished LmaxP%:
    LmaxPUF <- LmaxPcalc(Lmids = simSizeDist@pLPop[ , "LMids"],
                         NbyBin=simSizeDist@pLPop[ , "VulnUF"],
                         P = P)

    ## LmaxP% at SPR = <SPRref>:
    LmaxPF <- LmaxPcalc(Lmids = simSizeDist@pLPop[ , "LMids"],
                        NbyBin=simSizeDist@pLPop[ , "VulnF"],
                        P = P)


    ## In case it is too selective, bring the SPR down to the reference point:
    if (simSizeDist@SPR > SPRref &&
        length(pars@SPR) > 0)            # Not relevant if FM is provided.
    {

        warning("SPR ref cannot be reached, too high length at first capture:",
                "\n  Lmax", round(P * ifelse(P > 1, 1, 100)),
                "% will be overestimated and the associated F underestimated")
    }

    if (LmaxPF/LmaxPUF > 0.9)
        warning("Lmax",
                round(P * ifelse(P > 1, 1, 100)),
                "% ref point (",
                round(LmaxPF, 1),") close to Lmax",
                round(P * ifelse(P > 1, 1, 100)),
                "% un-fished (",
                round(LmaxPUF, 1),"):",
                "\n  maybe too low a contrast for an HCR to perform decently.")

    return(list(res = c(LmaxPRP = LmaxPF, LmaxPUF = LmaxPUF,
                        FM = as.vector(simSizeDist@FM),
                        F = as.vector(simSizeDist@FM * simSizeDist@M),
                        SPR = simSizeDist@SPR),
                simSizeDist = simSizeDist))
}


##' @title Lmax5refpointLBSPR: Lmax5% reference point from Length-Based Spawning Potential Ratio.
##' @param M A length-one numeric vector for Natural mortality.
##' @param k A length-one numeric vector for VBGF k parameter.
##' @param MK
##' @param Linf A length-one numeric vector for the assymptotic length.
##' @param Lm50 A length-one numeric vector for Length 50% maturity.
##' @param Lm95 A length-one numeric vector for Length 95% maturity.
##' @param Ls5 A length-one numeric vector for Selectivity length 5%.
##' One of Ls5 and Ls50 must be provided.
##' @param Ls50 A length-one numeric vector for Selectivity length 5%.
##' @param Ls95 A length-one numeric vector for Selectivity length 95%.
##' @param SPRref A length-one numeric vector for MSY proxy for SPR (default 0.4).
##' @param FMref A length-one numeric vector for MSY proxy for F/M (override SPRref).
##' @param Mpow A length-one numeric vector for Size varying mortality parameter (default = 0, constant mortality).
##' @param Wbeta A length-one numeric vector for Weight-length relationship beta param (W = a * L^b).
##' @param Walpha A length-one numeric vector for Weight-length relationship alpha param.
##' @param FecB A length-one numeric vector for beta parameter of length-fecundity relationship
##' @param R0
##' @param CVLinf A length-one numeric vector for Coefficient of variation on Linf.
##' @param Steepness A length-one numeric vector for steepness of SRR
##' @param MLL
##' @param sdLegal
##' @param fDisc
##' @param BinWidth
##' @param Control Controls for the LBSPRsim function (see LBSPR package help).
##' @return List with:
##'    $res: names vector of main outputs: Lmax5% refpoint, Lmax5% no exploitation,
##'                                        F/M, F (if M available), SPR.
##'    $simSizeDist: LBSPR simulation object (e.g. for plotting purpose).
##' @author Yves Reecht
Lmax5refpointLBSPR <- function(M, k, MK = M/k,
                               Linf, Lm50, Lm95,
                               Ls5, Ls50, Ls95,
                               SPRref = 0.4,
                               FMref = NULL,
                               Mpow = 0, Wbeta = 3, Walpha = 1e-04,
                               FecB = 3, R0 = 10000,
                               CVLinf = 0.1, Steepness = 0.7,
                               MLL = numeric(0), sdLegal = numeric(0),
                               fDisc = numeric(0),
                               BinWidth = numeric(0),
                               Control = list(modtype = "GTG", ngtg = 101))
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date:  4 Oct 2017, 12:37

    require(LBSPR)

    ## Parameters initialisation:
    suppressMessages(pars <- new("LB_pars"))
    if (missing(MK) || ! is.numeric(MK))
    {
        if (missing(M) || (! is.numeric(M)) ||
            missing(k) || (! is.numeric(k)))
        {
            stop("One of MK or (M and k) must be given")
        }else{
            pars@MK <- M/k
        }
    }else{
        pars@MK <- MK
    }

    if (! missing(M))
    {
        pars@M <- M
    }
    pars@Linf <- Linf
    if (missing(Ls50))
    {
        if (missing(Ls5)) stop("One of Ls5 and Ls50 must be given.")
        pars@SL50 <- (Ls95 + Ls5) / 2
    }else{
        pars@SL50 <- Ls50
    }
    pars@SL95 <- Ls95
    pars@L50 <- Lm50
    pars@L95 <- Lm95
    pars@Mpow <- Mpow
    pars@Wbeta <- Wbeta
    pars@Walpha <- Walpha
    pars@FecB <- FecB
    pars@CVLinf <- CVLinf
    pars@Steepness <- Steepness
    pars@R0 <- R0
    pars@MLL <- MLL
    pars@sdLegal <- sdLegal
    pars@fDisc <- fDisc

    pars@BinMin <- 0
    pars@BinMax <- ifelse(pars@MK > 0.5, 1.3, 1.8) * pars@Linf

    ## Default bin width:
    if (missing(BinWidth) || length(BinWidth) == 0)
    {
        pars@BinWidth <- (pars@Linf - pars@BinMin) / 100
    }else{
        pars@BinWidth <- BinWidth
    }

    ## Force FM if provided, and ignore SPRref (only if FM refpoint calculated externaly):
    if (missing(FMref) || is.null(FMref) || ! is.numeric(FMref))
    {
        pars@SPR <- SPRref
    }else{
        pars@FM <- FMref
    }

    suppressMessages(simSizeDist <- LBSPRsim(pars, Control = Control, verbose = FALSE))

    ## Unfished Lmax5%:
    Lmax5UF <- Lmax5calc(Lmids = simSizeDist@pLPop[ , "LMids"],
                         NbyBin=simSizeDist@pLPop[ , "VulnUF"])

    ## Lmax5% at SPR = 0.4:
    Lmax5F <- Lmax5calc(Lmids = simSizeDist@pLPop[ , "LMids"],
                        NbyBin=simSizeDist@pLPop[ , "VulnF"])


    ## In case it is too selective to bring the SPR down to the reference point:
    if (simSizeDist@SPR > SPRref &&
        length(pars@SPR) > 0)            # Not relevant if FM is provided.
    {
        ## warning("SPR ref cannot be reached, too high length at first capture:",
        ##         "\n  setting selectivity ogive = maturity ogive",
        ##         "\n  (-> conservative reference point).")

        ## pars@SL50 <- pars@L50
        ## pars@SL95 <- pars@L95

        ## simSizeDist <- LBSPRsim(pars, Control = Control, verbose = FALSE)

        ## warning("SPR ref cannot be reached, too high length at first capture:",
        ##         "\n  Lmax5% ref point set to SL50.")

        ## Lmax5F <- pars@SL50

        warning("SPR ref cannot be reached, too high length at first capture:",
                "\n  Lmax5% will be overestimated and the associated F underestimated")
    }

    if (Lmax5F/Lmax5UF > 0.9)
        warning("Lmax5% ref point (",
                round(Lmax5F, 1),") close to Lmax5% un-fished (",
                round(Lmax5UF, 1),"):",
                "\n  maybe too low a contrast for an HCR to perform decently.")

    return(list(res = c(Lmax5RP = Lmax5F, Lmax5UF = Lmax5UF,
                        FM = as.vector(simSizeDist@FM),
                        F = as.vector(simSizeDist@FM * simSizeDist@M),
                        SPR = simSizeDist@SPR),
                simSizeDist = simSizeDist))
}


## suppressWarnings(Lmax5refpointLBSPR(M=0.2, k=0.8, Linf=100, Lm50=45, Lm95=50, Ls5=50, Ls95=55,
##                                     SPRref=0.4, Mpow=0.3))$res

## ####################################################################################################
## Tanja's implementation:
##
## Cutting edge selectivity/maturity, no CVLinf.
## SL = Lmat by default:


Lmax5refpointMiethe <- function(Linf, M, k, b, Lmat, Lc = Lmat, SPRref = 0.4)
{
    ## Purpose:
    ## ----------------------------------------------------------------------
    ## Arguments:
    ## ----------------------------------------------------------------------
    ## Author: Tanja Miethe, adapted by Yves Reecht, Date:  4 Oct 2017, 15:47

    lm <- Lmat/Linf          # stand. length at maturity
    lmls <- Lc/Linf          # stand. length at first capture

    tm <- -(log(1-lm))/k           # age at maturity
    tmls <- -(log(1-lmls))/k       # age at first capture
    tmax <- Inf                    # assuming infinite age for integral

    ## SPR=0.4 find F
    fa <- function(t)
    {
        ((1-exp(-k*t))^b)*(exp(-M*t))
    }  # ssb0
    ia <- integrate(fa,tm,tmax)

    lim1 <- ia[[1]] * SPRref

    ##
    if(lm>=lmls)
    {
        fr <- function(Fi)
        {
            fc <- function(t)
            {
                ((1-exp(-k*t))^b)*(exp(Fi*tmls))*(exp(-(Fi+M)*t))
            }   # ssb fished

            i <-  integrate(fc,tm,tmax)[[1]]
            abs(i-lim1)
        }
    }

    ##
    if(lm<lmls)
    {
        fr <- function(Fi)
        {
            fb <- function(t)
            {
                ((1-exp(-k*t))^b)*(exp(-M*t))
            }
            fc <- function(t)
            {
                ((1-exp(-k*t))^b)*(exp(Fi*tmls))*(exp(-(Fi+M)*t))
            }   #ssb fished
            i <-  (integrate(fb,tm,tmls)[[1]]  + integrate(fc,tmls,tmax)[[1]] )
            return(abs(i - lim1))
        }
    }

    opt <- optimize(fr,lower=0,upper=2, tol=0.00001, maximum=F)   #find F (with upper limit F=2)
    Fi1 <- opt$minimum  # Fspr=40%


    ll1 <- 1-(0.05^(k/(M+Fi1)))*(1-lmls)   # lstar critical length with F for 5%

    tl1 <- -log(1-ll1)/k   #tstar age

    SPR.ref.point <- (1-(exp(-k*tl1))*(Fi1+M)/(Fi1+M+k))*Linf  # reference point for mean length of top 5%

    return(c(Lmax5RP = SPR.ref.point, FM = Fi1 / M, F = Fi1))
}

## ####################################################################################################
## Tests / usage:

## Linf <- 839   # VBGF
## k <- 0.197     # VBGF
## b <- 3.147      # length-weight relationship exponent
## Lmat <- 562   # maturity size


## M <- 4.118*k^0.73*Linf^(-0.33)   # natural mortality Then et al 2015
## M/k

## Lmax5refpointMiethe(Linf = Linf, M = M, k = k, b = b, Lmat = Lmat,
##                     Lc = Lmat, SPRref = 0.4)

## ## Similar settings with the full simulation model:
## res <- Lmax5refpointLBSPR(MK = M/k, M = M, Linf = Linf, Lm50 = Lmat, Lm95 = Lmat * 1.0001,
##                           Ls5 = Lmat * 0.9999, Ls95 = Lmat * 1.0001, CVLinf = 0.0001,
##                           FecB = b, Wbeta = b, ##Walpha = 1,
##                           ## Control = list(modtype = "absel", Nage = 1000, P = 0.0001),
##                           BinWidth = 1,
##                           # FM = Lmax5refpointMiethe(Linf = Linf, M = M, k = k, b = b, Lmat = Lmat,
##                          #                          Lc = Lmat, SPRref = 0.4)["FM"],
##                           SPRref = 0.4)
## ##
# res$res

## plotSim(res$simSizeDist)


## library(cowplot)
## library(ggplot2)
## simDist <- plotSize(res$simSizeDist) # +

##     geom_vline(xintercept = c(res$res[1:2]),
##                colour="red", linetype = c("solid", "longdash"),
##                size = 0.8)

## X11()
## plot(res$simSizeDist@pLPop[ , "LMids"], res$simSizeDist@pLPop[ , "VulnF"], type = "l")
## lines(res$simSizeDist@pLPop[ , "LMids"], res$simSizeDist@pLPop[ , "PopF"], col = "red", lty = 2)


## Lmids <- res$simSizeDist@pLPop[ , "LMids"]
## NbyBin <- res$simSizeDist@pLPop[ , "VulnF"]

## Lmax5calc(Lmids=Lmids, NbyBin=NbyBin)


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
  
  ff<-function(t){Fi1*(exp(Fi1*tmls))*exp(-(Fi1+M)*t)*(1-exp(-(k*t)))^b } #ield 
  if1<-integrate(ff,tmls,tmax)[[1]]  
  
  return(c(Lmax5 = lf1, F = Fi1, yield=if1))
}


Lmax5Fmax <- function(Linf, M, k, b, Lmat, Lc){
  
  lm<-Lmat/Linf  # standardized length at maturity
  lmls<-l/Linf   # standardized at first capture
  
  tm<--(log(1-lm))/k           # age at maturity
  tmls<--(log(1-lmls))/k       # age at first capture
  tmax<-Inf                    # assuming infinite age for integral          
  
  # find Fmax
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
  
  lmls<-l/Linf   # standardized at first capture
  ll0<-1-(0.05^(k/(M)))*(1-lmls)    # lstar critical length F=0
  tl0<--log(1-ll0)/k
  lf0<-(1-(exp(-k*tl0))*(M)/(M+k))*Linf 
  
  return(c(Lmax5 = lf0, F = 0, yield=0))
}













### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
