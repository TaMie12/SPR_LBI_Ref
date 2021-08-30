
#########################################################################################################################
### File: calculate_SPR_LBI_ref.R
###
### Created: 16/08/2021
### Author: Tanja Miethe
###
### R 3.6.2
#########################################################################################################################
### Description:
###   Functions to calculate SPR-based reference points for LBIs: Lmax5% and Lmean based on Miethe et al. (2019) using
###   either determinsitic equations as given in paper (with CVLinf = 0, knife-edge selectivity and maturity ogive) or
###   using the simulation tool in R-library LBSPR with setting CVLinf and define selectivity and maturity ogives.
###
###  Miethe, T., Reecht, Y., and Dobby, H. 2019. Reference points for length-based indicator Lmax5% to support
###  assessment of data-limited stocks and fisheries.  ICES Journal of Marine Science, 76:7, 2125-2139.
#########################################################################################################################


## determinsitic reference points (RP) with CVLinf = 0, knife-edged maturity/selectivity, const M
source(file = "SPR_LBI_Ref_fun1.R")

## fill in basic life history parameters
## asymptotic length, vB growth, cm
Linf <- 1000

## growth parameter vB growth
k <- 0.14

## natural mortality
M <- 0.2

## b parameter from length weight relationship W = aL^b
b <- 3

## length at 50% maturity (knife-edged ogive), cm
Lmat <- 700

## length at 50% selectivity (length at first capture, knife-edged), cm
Lc <- 400

## Spawning potential ratio (SPR) reference point 40%
SPRref <- 0.4

## calculate reference points,
## returns Lmax5% RP and unfished value (UF, F = 0), Lmean RP and unfished value estimated using LBSPR,
## unfished values (UF), F/M and F at 40%SPR, M/k, SPR used:
res1 <- Lmax5_Lmean_LBSPR1(Linf = Linf, M = M, k = k, Wbeta = b,
                           Lm50 = Lmat, Ls50 = Lc, SPRref = SPRref)

res1


######################################################################################################################

## reference point using LBSPR simulation tool

library(LBSPR)

source(file = "SPR_LBI_Ref_fun2.R")

## basic life history parameters as above
## CV for Linf
CVLinf <- 0.000001

## selectivity ogive, length 5%, length 95%
Ls5 <- Lc - 0.01
Ls95 <- Lc + 0.01

## maturity ogive, length 95%
Lm95 <- Lmat+0.01

res2 <- LmaxP_Lmean_ref_LBSPR2(M, k, MK = M / k,
                               Linf = Linf, Lm50 = Lmat, Lm95 = Lm95,
                               Ls5 = Ls5, Ls50 = Lc, Ls95 = Ls95,
                               P = 0.05,
                               SPRref = SPRref,
                               FMref = NULL,
                               Mpow = 0, Wbeta = b, Walpha = 1e-04,
                               FecB = b, R0 = 10000,
                               CVLinf = CVLinf, Steepness = 0.7,
                               MLL = numeric(0), sdLegal = numeric(0),
                               fDisc = numeric(0),
                               BinWidth = numeric(0),
                               Control = list(modtype = "GTG", ngtg = 101))
res2$res

