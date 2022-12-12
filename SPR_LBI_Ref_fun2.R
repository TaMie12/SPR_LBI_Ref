
#########################################################################################################################
### File: SPR_LBI_Ref_fun2.R
### Created: 16/08/2021
### Author: Yves Reecht, Tanja Miethe
###
#########################################################################################################################
### Description:
###   Functions to calculate SPR-based reference points for LBIs: Lmax5% and Lmean using the simulation tool in
###   R-library LBSPR with setting CVLinf and define selectivity and maturity ogives
###
### Miethe, T., Reecht, Y., and Dobby, H. 2019. Reference points for length-based indicator Lmax5% to support
### assessment of data-limited stocks and fisheries.  ICES Journal of Marine Science, 76:7, 2125-2139.
#########################################################################################################################


##' @title Lmax<P>% Lmean estimation from abundances/proportions at length.
##' @param Lmids A numeric vector of mid-values of length-bins.
##' @param NbyBin A numeric vector of abundance/proportion by length-bin.
##' @return Lmax<P>% and Lmean value estimated from abundances/proportions at length (length-one numeric vector).
##' @author Yves Reecht
##'
LmaxPcalc1 <- function(Lmids, NbyBin, P = 0.05)
{

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



Lmeancalc1 <- function(Lmids, NbyBin)
{


    ## Weighted mean on the bin mid point values:
    Lmean <- weighted.mean(x = Lmids, w = NbyBin)
    return(Lmean)
}


LmaxPcalc <- function(Lmids, NbyBin, P = 0.05)
{


    sapply(P,
           function(x, Lmids, NbyBin)
    {
        LmaxPcalc1(P = x, Lmids = Lmids, NbyBin = NbyBin)
    }, Lmids = Lmids, NbyBin = NbyBin)
}




Lmeancalc <- function(Lmids, NbyBin)
{

  return(Lmeancalc1(Lmids = Lmids, NbyBin = NbyBin))
}


##' @title LmaxP_mean_ref_LBSPR2: Lmax<P>% and Lmean reference point from Length-Based Spawning Potential Ratio.
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
##' @author Yves Reecht, Tanja Miethe


LmaxP_Lmean_ref_LBSPR2 <- function(M, k, MK = M/k,
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
                         NbyBin = simSizeDist@pLPop[ , "VulnUF"],
                         P = P)

    ## LmaxP% at SPR = <SPRref>:
    LmaxPF <- LmaxPcalc(Lmids = simSizeDist@pLPop[ , "LMids"],
                        NbyBin = simSizeDist@pLPop[ , "VulnF"],
                        P = P)


    ## Unfished Lmean:
    LmeanUF <- Lmeancalc(Lmids = simSizeDist@pLPop[ , "LMids"],
                         NbyBin = simSizeDist@pLPop[ , "VulnUF"])

    ## Lmean at SPR = <SPRref>:
    LmeanF <- Lmeancalc(Lmids = simSizeDist@pLPop[ , "LMids"],
                        NbyBin = simSizeDist@pLPop[ , "VulnF"])



    ## In case it is too selective, bring the SPR down to the reference point:
    if (simSizeDist@SPR > SPRref &&
        length(pars@SPR) > 0)            # Not relevant if FM is provided.
    {

        warning("SPR ref cannot be reached, too high length at first capture:",
                "\n  Lmax", round(P * ifelse(P > 1, 1, 100)),
                "% will be overestimated and the associated F underestimated")
    }

    if (LmaxPF / LmaxPUF > 0.9)
        warning("Lmax",
                round(P * ifelse(P > 1, 1, 100)),
                "% ref point (",
                round(LmaxPF, 1),") close to Lmax",
                round(P * ifelse(P > 1, 1, 100)),
                "% un-fished (",
                round(LmaxPUF, 1),"):",
                "\n  maybe too low a contrast for an HCR to perform decently.")

    return(list(res = round(c(LmaxPRP = LmaxPF, LmaxPUF = LmaxPUF,
                              LmeanRP=LmeanF, LmeanUF = LmeanUF,
                              F_M = as.vector(simSizeDist@FM),
                              M_k = MK,
                              SPR = simSizeDist@SPR),3),
                simSizeDist = simSizeDist))
}




### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
