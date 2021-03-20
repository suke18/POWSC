#' Estimate characterized parameters for a given scRNA-seq data (SingleCellExperiment object or a count matrix).
#'
#' These parameters include four gene-wise parameters and two cell-wise parameters.
#'
#' @param sce SingleCellExperiment object with assays(sce)[[1]] is the count matrix or input directly
#' @param low.prob lower bound probability for phase I
#' @return a list of needed estimated parameters
#' @examples
#' data("es_mef_sce")
#' sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
#' set.seed(123)
#' rix = sample(1:nrow(sce), 500)
#' sce = sce[rix, ]
#' estParas = Est2Phase(sce)
#' @export Est2Phase
Est2Phase = function(sce, low.prob=0.99){
    # check the input whether it is an SingleCellExperiment Object or a matrix
    if (is(sce, "SingleCellExperiment")){
        Y = round(assays(sce)[[1]])
    }else{Y = round(sce)}
    ## ## initial estimate of prob(X\in bg)
    Cell0=colMeans(Y==0) # each cell has this percentage 0
    par1= apply(Y,2,function(yy) {
        yyy=yy[yy<=15]
        ### consider (0s, 3, 4)
        if (length(table(yyy)) <= 2){
            # print(which(colSums(Y==yy) == nrow(Y)))
            yyy = yy[yy<=20]
        }
        RobustPoi0(yyy)
        })
    pi0.hat=Cell0/(par1[1,]+(1-par1[1,])*dpois(0,par1[2,]))
    if (any((pi0.hat > 1))) {warning("Zero proportion is greater than estimation.")}
    pi0.hat <- pmin(pi0.hat, 1)
    prob0=pi0.hat*par1[1,]+ pi0.hat*(1-par1[1,])*dpois(0,par1[2,]) ## ZIP prob at 0
    ## First round
    ## get the 1-low.prob quantile of ZIP
    x0=qpois(pmax(1-(1-low.prob)/(1-par1[1,]),0),par1[2,])
    Z= sweep(Y,2,x0)>0 # indicate if a gene is > bg
    L=colSums(Y*Z)/1e6 # so far it is like simple total..

    mu.g1=log2(rowSums(Z*Y)/rowSums(sweep(Z,2,L,FUN="*")))
    mu.g1[is.na(mu.g1)]=0 ## if allZ is 0, it gets NA,
    ### but we should shrink mu.g1 as well since some mu.g1 is estimated by only a few observations
    ## leave it here for now.
    n.g1=rowSums(Z)
    y1=log2(sweep(Y,2,L,FUN="/")+1) #like TPM**
    s.g1=sqrt(rowSums(Z*sweep(y1,1,mu.g1)^2)/(n.g1-1)) ## TPM type of SD
    mu.g2 = shrink.mu(mu.g1,s.g1,n.g1)
    ## get sd.g
    res.g1=log2(sweep(Y,2,L,FUN="/")+1)-mu.g1
    ## mad of those res.g1 that are associated with Z==1
    tmp=array(0,dim=c(dim(res.g1),2))
    tmp[,,1]=res.g1;tmp[,,2]=Z
    sd.g1=apply(tmp,1,function(xx) my.mad(xx[xx[,2]==1,1]))
    sd.g1[is.na(sd.g1)]=0## if all bg, there's no info about fg sd
    ## add a shrinkage for sd.g1
    sd.prior=squeezeVar(sd.g1^2,n.g1-1)
    sd.g2=sqrt(sd.prior$var.post)
    #####  gene specific bg. Z_gi
    den.fg = den.bg = NA*Y
    for(i in 1:ncol(Y)){
        den.bg[,i]=dZinf.pois(Y[,i], par1[1,i], par1[2,i])
        den.fg[,i]=dLNP2(x=Y[,i], mu=mu.g1, sigma=sd.g2, l=L[i])
    }
    Z.fg=sweep(den.fg,2,1-pi0.hat,FUN="*")
    Z.bg=sweep(den.bg,2,pi0.hat,FUN="*")
    post.Z=Z.fg/(Z.fg+Z.bg)
    post.Z[is.na(post.Z)] <- 1

    ### if I shrink mu.g
    den.fg2 = NA*Y
    for (i in 1:ncol(Y)){
        den.fg2[,i]= dLNP2(x=Y[,i], mu=mu.g2, sigma=sd.g2, l=L[i])
    }
    Z.fg2=sweep(den.fg2,2,1-pi0.hat,FUN="*")
    post.Z2=Z.fg2/(Z.fg2+Z.bg)
    post.Z2[is.na(post.Z2)] <- 1

    pi.g = rowMeans(post.Z2)
    est_rslt = list(exprs = Y, pi.g = pi.g, p0 = par1[1,], lambda = par1[2,],
                    mu = mu.g2, sd = sd.g2, sf = L)
    return(est_rslt)
}



