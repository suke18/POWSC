########################################################
## Power function, compute the stratified power-realted quantities
##
## Return values:
## TD: number of True Discoveries in each stratum
## FD: number of False Discoveries in each stratum
## power: within strata, proportion of TD out of total DE
## alpha.nomial: cutoff of alpha on raw p-values
## alpha: empirical p-value in each stratum
## alpha.marginal: overall empirical p-value
########################################################
PowerEst = function(fdr, alpha, Zg, Zg2, xgr){
    ## p is input nominal p-value or FDR.
    ## alpha is cutoff of p
    ## !Zg is indicator for TN,  Zg2 is indicators for TP,
    ## xgr is grouping in covariate

    ix.D = fdr <= alpha
    N = sum(ix.D) ## this is the total discovery
    N.stratified = tapply(ix.D, xgr, sum)

    ##  TD (called DE genes)
    id.TP = Zg2==1
    TD = tapply(fdr[id.TP] <= alpha, xgr[id.TP], sum)
    TD[is.na(TD)] = 0

    ##  FD
    id.TN = Zg==0
    FD = tapply(fdr[id.TN] <= alpha, xgr[id.TN], sum)
    FD[is.na(FD)] = 0

    ## type I error
    alpha = as.vector(FD/table(xgr[id.TN]))
    alpha.marginal = sum(FD)/sum(id.TN)

    ## power & precision
    power=as.vector(TD/table(xgr[id.TP]))
    power.marginal=sum(TD,na.rm=TRUE)/sum(id.TP)

    ## FDR
    FDR = FD / N.stratified
    FDR.marginal = sum(FD, na.rm=TRUE) / N

    ## conditional truths (true DE genes)
    CT = table(xgr[id.TP])

    list(CD=TD,FD=FD,TD=CT,alpha.nominal=alpha,
         alpha=alpha, alpha.marginal=alpha.marginal,
         power=power, power.marginal=power.marginal,
         FDR=FDR, FDR.marginal=FDR.marginal)
}

#' Run DE analysis by using MAST. Here we output two result tables corresponding to two forms of DE genes.
#' These parameters include four gene-wise parameters and two cell-wise parameters.
#'
#' @param DErslt is from the DE analysis by MAST
#' @param simData is the corresponding simulated scRNA-seq dataset (SingCellExperiment)
#' @param alpha is the cutoff for the fdr which can be modified
#' @param delta or the lfc is the cutoff (=0.5) used to determined the high DE genes for Form II
#' @param strata can be modified by the user. By default, it is (0, 10], (10, 20], (20, 40], (40, 80], (80, Inf]
#' @return a list of metrics for power analysis such as: stratified targeted power and marginal power.
#' @examples
#' data("es_mef_sce")
#' sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
#' estParas = Est2Phase(sce)
#' simData = Simulate2SCE(n=100, estParas1 = estParas, estParas2 = estParas)
#' DErslt = runDE(simData$sce)
#' Disc_cont = Power_Cont(DErslt, simData)
#' @export Power_Cont
## Continous case corresponding to the Phase II DE, delta means lfc
Power_Cont = function(DErslt, simData, alpha = 0.1, delta = 0.5, strata = c(0,10,2^(1:4)*10,Inf)){
    fdrvec = DErslt$cont$fdr
    lfc = simData$lfc
    ngenes = nrow(simData$sce)
    DEid = simData$ix.DE2
    Zg = Zg2 = rep(0, ngenes)
    Zg[DEid] = 1
    ix = which(abs(lfc) > delta)
    Zg2[DEid[ix]] = 1
    sce = simData$sce
    Y = round(assays(sce)[[1]])
    # sizeF = colSums(Y)
    # sizeF = sizeF/median(sizeF)
    # X.bar = rowMeans(sweep(Y,2,sizeF,FUN="/"))
    X.bar = rowMeans(Y)
    ix.keep = which(X.bar>0)
    xgr = cut(X.bar[ix.keep], strata)

    # lev = levels(xgr)
    # ix.keep = ix.keep[!(xgr %in% lev[strata.filtered])]
    # ## recut
    # xgr = cut(X.bar[ix.keep], strata[-strata.filtered])

    # Interested genes
    Zg = Zg[ix.keep]; Zg2 = Zg2[ix.keep]
    fdrvec = fdrvec[ix.keep]

    pow = PowerEst(fdrvec, alpha, Zg, Zg2, xgr=xgr)
    return(pow)
}



#' Run DE analysis by using MAST. Here we output two result tables corresponding to two forms of DE genes.
#' These parameters include four gene-wise parameters and two cell-wise parameters.
#'
#' @param DErslt is from the DE analysis by MAST
#' @param simData is the corresponding simulated scRNA-seq dataset (SingCellExperiment)
#' @param alpha is the cutoff for the fdr which can be modified
#' @param delta or the zero ratio change is the cutoff (=0.1) used to determined the high DE genes for Form II
#' @param strata can be modified by the user. By default, it is (0, 0.2], (0.2, 0.4], (0.4, 0.6], (0.6, 0.8], (0.8, 1]
#' @return a list of metrics for power analysis such as: stratified targeted power and marginal power.
#' @examples
#' estPower1 = Power_Disc(de, simData = simData)
#' @export Power_Disc
## Discreate case corresponding to the Phase I DE, delta means pi.df
Power_Disc = function(DErslt, simData, alpha = 0.1, delta = 0.1, strata = seq(0, 1, by = 0.2)){
    fdrvec = DErslt$disc$fdr
    pi.df = simData$pi.df
    ngenes = nrow(simData$sce)
    DEid = simData$ix.DE1
    Zg = Zg2 = rep(0, ngenes)
    Zg[DEid] = 1
    ix = which(abs(pi.df) > delta)
    Zg2[DEid[ix]] = 1
    sce = simData$sce
    ntotal = ncol(sce)
    Y = round(assays(sce)[[1]])
    rate0 = rowMeans(Y==0)
    ix.keep = intersect(which(rate0 < 0.99),  which(rate0 > 0.01)) # too small none 0 rate cannot detect
    xgr = cut(rate0[ix.keep], strata)

    # lev = levels(xgr)
    # ix.keep = ix.keep[!(xgr %in% lev[strata.filtered])]
    # ## recut
    # xgr = cut(X.bar[ix.keep], strata[-strata.filtered])

    # Interested genes
    Zg = Zg[ix.keep]; Zg2 = Zg2[ix.keep]
    fdrvec = fdrvec[ix.keep]

    pow = PowerEst(fdrvec, alpha, Zg, Zg2, xgr=xgr)
    return(pow)
}






