#-*- coding: latin-1 -*-

### File: 2_discussion_examples.R
### Time-stamp: <2021-08-30 11:16:23 a23579>
###
### Created: 30/04/2018	13:59:59
### Author: Yves Reecht
###
####################################################################################################
### Description:
###
###
####################################################################################################

## ####################################################################################################
##

rm(list = ls())

## setwd("")

source("./0_Functions_Lmax5refpointLBSPR.R")

library(R.utils)
library(LBSPR)

## Space to explore:
MK <- c(0.5, 1.5)
CVLinf <- c(0.0001, 0.1)
Linf <- 100
Mpow <- 0
Lm50 <- 50
Lm95 <- 51
Ls50 <- Lm50
Ls95 <- Lm95
BinWidth <- 1


## possible combinations of factors:
dataGridComp <- expand.grid(MK = MK, Linf = Linf, Mpow = Mpow,
                            Lm50 = Lm50, Lm95 = Lm95, Ls50 = Ls50, Ls95 = Ls95,
                            BinWidth = BinWidth, CVLinf = CVLinf, P = 0.05)

refpointsComp <- by(data = dataGridComp,
                    INDICES = as.list(dataGridComp),
                    FUN = function(x)
                    {
                        suppressWarnings(res <- do.call(what = LmaxPrefpointLBSPR, args = as.list(x)))
                        ## Appending results as supplementary columns:
                        cbind(x, as.matrix(t(res$res)))
                    }, simplify = FALSE)

refpointsComp <- do.call(rbind, refpointsComp)

## Metrics with Lmax5%
refpointsComp$Lmax5RP.pcUF <- 100 * refpointsComp$LmaxPRP / refpointsComp$LmaxPUF
refpointsComp$Lmax5RP.pcLinf <- 100 * refpointsComp$LmaxPRP / refpointsComp$Linf
refpointsComp$Lmax5RP.pcDelta <- 100 * (refpointsComp$LmaxPRP - refpointsComp$LmaxPUF) / refpointsComp$LmaxPUF
refpointsComp

## Increase in the RP value CVLinf 0 -> 0.1 at M/k = 0.5 and 1.5:
print(diff(refpointsComp[refpointsComp$MK == 0.5, "Lmax5RP.pcLinf"]))
print(diff(refpointsComp[refpointsComp$MK == 1.5, "Lmax5RP.pcLinf"]))


## Example misspecification of CVLinf (0, instead of 0.1):
objFun <- function(par, dat, refVal)
{
    ## Purpose: Objective function to retrieve SPR corresponding to a given
    ##          Lmax5%.
    ## ----------------------------------------------------------------------
    ## Arguments: par: SPR value (optimized parameter).
    ##            dat: one-row data.frame with named arguments for
    ##                 LmaxPrefpointLBSPR(...).
    ##            refVal: the Lmax5% ref point to reach.
    ## ----------------------------------------------------------------------
    ## Author: Yves Reecht, Date: 30 Apr 2018, 14:17

    set.seed(234)
    suppressWarnings(res <- do.call(LmaxPrefpointLBSPR, c(as.list(dat), SPRref = par[1])))
    SqE <- (res$res["LmaxPRP"] - refVal) ^ 2
    return(SqE)
}

SPRcorr0.5 <- optim(par = 0.3, fn = objFun, lower = 0.01, upper = 0.9, method = "Brent", ## "L-BFGS-B",
                    dat = dataGridComp[3, ], refVal = refpointsComp[1 , "LmaxPRP"])

SPRcorr1.5 <- optim(par = 0.3, fn = objFun, lower = 0.01, upper = 0.9, method = "Brent", ## "L-BFGS-B",
                    dat = dataGridComp[4, ], refVal = refpointsComp[2 , "LmaxPRP"])


print(SPRcorr0.5$par)
print(SPRcorr1.5$par)


## ####################################################################################################
## Effect on Lmax5% value of taking CVLinf into account (0 versus 0.1):
multLm <- 0.5
multLc <- 1                             # Lm=Lc=0.5Linf

message(ifelse(multLc == 1,
               paste0("Lm=Lc=", multLm, "Linf"),
               paste0("Lm=", multLm, "Linf_Lc=",
                      multLc, "Lm")))

## Space to explore:
MK <- c(0.5, 1.5) ##seq(0.2, 2.5, 0.1)
CVLinf <- c(0.000001, 0.1)## seq(0.03, .3, by = 0.03))
Linf <- 100
Mpow <- 0
Lm50 <- Linf * multLm
Lm95 <- Linf * multLm + 1
Ls50 <- Lm50 * multLc
Ls95 <- Lm50 * multLc + 1
BinWidth <- 1

## possible combinations of factors:
dataGrid <- expand.grid(MK = MK, Linf = Linf, Mpow = Mpow,
                        Lm50 = Lm50, Lm95 = Lm95, Ls50 = Ls50, Ls95 = Ls95,
                        BinWidth = BinWidth, CVLinf = CVLinf, P = 0.05)

## Calculate Lmax5% for combinations of factors:
pb <- ProgressBar(max=50, stepLength = (50/nrow(dataGrid)) + sqrt(.Machine$double.eps), ticks = 10)
reset(pb)
refpoints <- by(data = dataGrid,
                INDICES = as.list(dataGrid),
                FUN = function(x)
                {
                    increase(pb)
                    suppressWarnings(res <- do.call(what = LmaxPrefpointLBSPR, args = as.list(x)))
                    ## Appending results as supplementary columns:
                    cbind(x, as.matrix(t(res$res)))
                }, simplify = FALSE)

refpoints <- do.call(rbind, refpoints)


## Ordered results:
print(refpoints[do.call(order, as.list(refpoints[ , c("MK", "CVLinf")])) ,
                c("MK", "CVLinf", "LmaxPRP")])

## Increase (%) in Lmax5% between CVLinf = 1e-6 (~0) and CVLinf = 0.1, at M/k = 0.5:
message("## Increase (%) in Lmax5% between CVLinf = 1e-6 (~0) and CVLinf = 0.1, at M/k = 0.5:")
print(round(100 * diff(refpoints$LmaxPRP[c(1, 3)]) / refpoints$LmaxPRP[1], 1))

## Increase (%) in Lmax5% between CVLinf = 1e-6 (~0) and CVLinf = 0.1, at M/k = 1.5:
message("## Increase (%) in Lmax5% between CVLinf = 1e-6 (~0) and CVLinf = 0.1, at M/k = 1.5:")
print(round(100 * diff(refpoints$LmaxPRP[c(2, 4)]) / refpoints$LmaxPRP[2], 1))



### Local Variables:
### ispell-local-dictionary: "english"
### fill-column: 120
### End:
