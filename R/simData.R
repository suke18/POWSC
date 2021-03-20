rzip = function(n, p0, lambda) {
    y = rep(0, n)
    ix = rbinom(n, size=1, prob=p0)
    y[!ix]  = rpois(sum(!ix), lambda)
    y
}

### generate LNP random variable
## The real data is not exactly LNP, the right tail is much shorter
## Using LNP to simulate there are too many large numbers.
## How to shave the right tail??
rLNP = function(n, mu, sigma, sf) {
    theta = 2^rnorm(n, mean=mu, sd=sigma)
    ## should shave the right tail of theta a bit to avoid excessive large number??

    y = suppressWarnings(rpois(n, theta*sf))
    if (sum(is.na(y)) > 0){
        m = max(y[!is.na(y)])
        y[is.na(y)] = m
    }
    return(y)
}

## main function for generating the count matrix
GenerateCountMatrix = function(pi.g, p0, lambda, mu, sigma, sf){

    stopifnot(identical(length(p0), length(lambda), length(sf)),
              identical(length(pi.g), length(mu), length(sigma)))

    N = length(p0)
    G = length(mu)

    ## simulate Z, indicator for being in foreground.
    Z = matrix(rbinom(N*G, size=1, prob=pi.g), ncol=N)

    ## loop over cells to generate counts
    Y = matrix(0, nrow=G, ncol=N)
    for(i in 1:N) {
        ix <- Z[,i] == 0
        Y[ix, i]  = rzip(sum(ix), p0[i], lambda[i])
        Y[!ix, i]  = rLNP(sum(!ix), mu[!ix], sigma[!ix], sf[i])
    }

    return(Y)
}






##################################################################
####### Function 1 considering two (pair-wised) conditions #######
##################################################################

#' Simulate the data for two-group comparison; e.g., treatment v.s. control
#' It simulates the DE changes in two forms corresponding two types of DE genes
#'
#' @param n the number of total cells for two groups
#' @param perDE percentage of DE genes
#' @param estParas1 the set of parameters corresponding to cell type I
#' @param estParas2 the set of parameters corresponding to cell type II
#' @return a list of metrics recording the changes in the generated data: such as the DE gene indices for Form I and II DE genes, and simulated expression data in singlecellexperiment format.
#' @examples
#' data("es_mef_sce")
#' sce = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
#' set.seed(123)
#' rix = sample(1:nrow(sce), 500)
#' sce = sce[rix, ]
#' estParas = Est2Phase(sce)
#' simData = Simulate2SCE(n=100, estParas1 = estParas, estParas2 = estParas)
#' @export Simulate2SCE
Simulate2SCE = function(n = 100, perDE = 0.05, estParas1, estParas2) {
    # equally divide the sample size
    n1 = n2 = round(n / 2)
    # Parameters from estimations
    pi.g1 = estParas1$pi.g; pi.g2 = estParas2$pi.g
    p01 = estParas1$p0; p02 = estParas2$p0
    lambda1 = estParas1$lambda; lambda2 = estParas2$lambda
    mu1 = estParas1$mu; mu2 = estParas2$mu
    sigma1 = estParas1$sd; sigma2 = estParas2$sd
    sf1 = estParas1$sf; sf2 = estParas2$sf

    if (length(pi.g1) == length(pi.g2)){
        ngene = length(pi.g1)
    }else{
        stop("Please Input Proper Parameter Estimation Objects")
    }
    ## simulate index for DE genes.
    ## I'm picking from genes with large means.
    ## If DE genes are for genes with very small mean, they won't be detected.
    n0 = max(sum(mu1 > 3), sum(mu2 > 3))
    nDE1 = nDE2 = ngene * perDE
    ix1.highGenes = order(mu1, decreasing=TRUE)[1:n0]
    ix2.highGenes = order(mu2, decreasing=TRUE)[1:n0]

    ix.DE1 = sample(union(ix1.highGenes, ix2.highGenes), nDE1)
    ix.DE2 = sample(union(ix1.highGenes, ix2.highGenes), nDE2)

    ix.DEs = union(ix.DE1, ix.DE2)
    ## general parameters in both groups
    sf1 = sample(sf1, n1, replace=TRUE)
    sf2 = sample(sf2, n2, replace=TRUE)
    p0.1 = sample(p01, n1, replace=TRUE )
    p0.2 = sample(p02, n2, replace=TRUE)
    lambda1 = sample(lambda1, n1, replace=TRUE)
    lambda2= sample(lambda2, n2, replace=TRUE)

    ## For Phase I DEGs zero ratio
    tmp = pi.g2[ix.DE1]
    tmp[tmp < 0.5] = tmp[tmp < 0.5] + runif(sum(tmp<0.5), 0.1, 0.3)
    tmp[tmp >= 0.5] = tmp[tmp >= 0.5] - runif(sum(tmp >=0.5), 0.1, 0.3)
    pi.g2[ix.DE1] = tmp
    pi.df = tmp - pi.g1[ix.DE1]

    ## For Phase II DEGs lfc
    tmp = c(rnorm(1000, mean=-1, sd=1), rnorm(1000, mean=1, sd=1))
    mu.diff = sample(tmp, length(ix.DE2))
    mu2[ix.DE2] = mu1[ix.DE2] + mu.diff
    lfc  = mu.diff

    ## generate counts in two groups
    y1 = GenerateCountMatrix(pi.g1, p0.1, lambda1, mu1, sigma1, sf1)
    y2 = GenerateCountMatrix(pi.g2, p0.2, lambda2, mu2, sigma2, sf2)
    y = cbind(y1, y2)
    rownames(y) = paste0("g", 1:nrow(y))
    celltypes = rep(paste0("celltype", c(1,2)), c(n1, n2))
    sce = SingleCellExperiment(
        assays = list(counts = y),
        rowData = data.frame(geneNames = rownames(y), stringsAsFactors = FALSE),
        colData = data.frame(cellTypes = celltypes, stringsAsFactors = FALSE)
    )
    DEGs = paste0("g", union(ix.DE1, ix.DE2))
    list(ix.DE1=ix.DE1, ix.DE2=ix.DE2, ix.DEs = ix.DEs, DEGs = DEGs, sce = sce,
         pi.g1=pi.g1, pi.g2=pi.g2, mu1=mu1, mu2=mu2, lfc = lfc, pi.df=pi.df, ngenes = ngene)
}






##################################################################
##### Function 2 considering a mixture cell type conditions ######
##################################################################

#' Simulate the data for multiple-group comparisons; e.g., different cell types in blood
#' It simulates the DE changes in two forms corresponding two types of DE genes
#'
#' @param n the number of total cells for multiple groups; e.g., 1000, 2000, and etc.
#' @param estParas_set a set of parameters corresponding to different cell types.
#' @param multiProb a vector of probilities correponding to each cell type. It is not necessary to sum up to 1 because POWSC will normalize internally.
#' @param delta1 the minimum of expression change used to determine the Form I DE.
#' @param delta2 the minimum of log fold change used to determine the Form II DE.
#' @return a list of simulated datasets. Each dataset corresponds to a pair-wise comparison including a series of metrics such as the DE gene indices for Form I and II DE genes, and simulated expression data in singlecellexperiment format.
#' @examples
#' data("es_mef_sce")
#' set.seed(123)
#' rix = sample(1:nrow(es_mef_sce), 500)
#' es_mef_sce = es_mef_sce[rix, ]
#' sce1 = es_mef_sce[, colData(es_mef_sce)$cellTypes == "fibro"]
#' estParas1 = Est2Phase(sce1)
#' sce2 = es_mef_sce[, colData(es_mef_sce)$cellTypes == "stemCell"]
#' estParas2 = Est2Phase(sce2)
#' estParas_set = list(celltype1 = estParas1, celltype2 = estParas1, celltype3 =estParas2)
#' multiProb = c(0.2, 0.3, 0.5)
#' simData = SimulateMultiSCEs(n=200, estParas_set = estParas_set, multiProb = multiProb)
#' @export SimulateMultiSCEs
SimulateMultiSCEs = function(n = 1000, estParas_set, multiProb, delta1 = 0.1, delta2 = 0.5) {
    ## Initialize the simAll object
    celltypeNames = names(estParas_set)
    ncelltype = length(multiProb)
    SimulateMultiSCEs = NULL
    multiSize = as.vector(rmultinom(1, n, multiProb))
    geneNames = rownames(estParas_set[[1]]$exprs)

    ## pairwisely simulate the data
    for (i in 1:ncelltype){
        for (j in i:ncelltype){
            if (j > i){
                comID = c(i, j)
                compName = paste0(celltypeNames[i], "_vs_", celltypeNames[j])
                ## two set of parameters
                Paras_set1 = estParas_set[[i]]; Paras_set2 = estParas_set[[j]]
                pi.g1 = Paras_set1$pi.g; pi.g2 = Paras_set2$pi.g
                p01 = Paras_set1$p0; p02 = Paras_set2$p0
                lambda1 = Paras_set1$lambda; lambda2 = Paras_set2$lambda
                mu1 = Paras_set1$mu; mu2 = Paras_set2$mu
                sigma1 = Paras_set1$sd; sigma2 = Paras_set2$sd
                sf1 = Paras_set1$sf; sf2 = Paras_set2$sf
                n1 = multiSize[i]; n2 = multiSize[j]
                ## simulation of cell-wise parameters
                sf1 = sample(sf1, n1, replace=TRUE)
                p01 = sample(p01, n1, replace=TRUE )
                lambda1 = sample(lambda1, n1, replace=TRUE)

                sf2 = sample(sf2, n2, replace=TRUE)
                p02 = sample(p02, n2, replace=TRUE )
                lambda2 = sample(lambda2, n2, replace=TRUE)

                Y1 = GenerateCountMatrix(pi.g1, p01, lambda1, mu1, sigma1, sf1)
                Y2 = GenerateCountMatrix(pi.g2, p02, lambda2, mu2, sigma2, sf2)
                Y = cbind(Y1, Y2)

                ## Record the DE index
                tmp0ratio = pi.g2 - pi.g1
                tmplfc = mu2 - mu1
                ix.DE1 = which(abs(tmp0ratio) > delta1)
                ix.DE2 = which(abs(tmplfc) > delta2)
                celltypes = rep(celltypeNames[comID], multiSize[comID])

                ## Save to sce
                sce = SingleCellExperiment(
                    assays = list(counts = Y),
                    rowData = data.frame(geneNames = geneNames, stringsAsFactors = FALSE),
                    colData = data.frame(cellTypes = celltypes, stringsAsFactors = FALSE)
                )
                SimulateMultiSCEs[[compName]] = list(sce = sce, ix.DE1 = ix.DE1, ix.DE2 = ix.DE2,
                                                     pi.df = tmp0ratio, lfc = tmplfc)
            }
            else{
                next
            }
        }
    }
    return(SimulateMultiSCEs)
}


# SimulateMultiSCEs = function(n = 1000, perDE = 0.05, estParas_set, multiProb) {
#     ## Initialize the simAll object
#     celltypeNames = names(estParas_set)
#     ncelltype = length(multiProb)
#     simAll = vector(mode = "list", length = ncelltype)
#     ## Baseline Parameters
#     baseParas = estParas_set[[1]]
#     pi.g = baseParas$pi.g
#     p0 = baseParas$p0
#     lambda = baseParas$lambda
#     mu = baseParas$mu
#     sigma = baseParas$sd
#     sf = baseParas$sf
#     ngene = length(pi.g)
#     n0 = round(ngene*perDE)
#     # Simulate index for DE genes. Picking from genes with large means otherwise
#     # genes with very small mean won't be detected.
#     nhighGenes = sum(mu > 2)
#     ix.highGenes = order(mu, decreasing=TRUE)[1:nhighGenes]
#
#     ## Simulate the baseline cell type
#     multiSize = as.vector(rmultinom(1, n, multiProb))
#     baseSize = multiSize[1]
#     base.sf = sample(sf, baseSize, replace=TRUE)
#     base.p0 = sample(p0, baseSize, replace=TRUE )
#     base.lambda = sample(lambda, baseSize, replace=TRUE)
#     base.y = GenerateCountMatrix(pi.g, base.p0, base.lambda, mu, sigma, base.sf)
#     base.rslt = list(Y = base.y, ix.DE1 = NULL, ix.DE2 = NULL, pi.g = pi.g, mu = mu)
#     simAll[[1]] = base.rslt
#     ## Simulate number of the cell types from multinomial distribution
#     for (i in 2 : ncelltype){
#         tmpParas = estParas_set[[i]]
#         tmpSize = multiSize[i]
#         tmpLfc = tmpParas$mu - mu
#         tmp0ratio = tmpParas$pi.g - pi.g
#         ix.DE1 = order(abs(tmp0ratio), decreasing = T)[1:5000]
#         ix.DE1 = intersect(ix.DE1, ix.highGenes)[1:n0]
#         ix.DE2 = order(abs(tmpLfc), decreasing = T)[1:5000]
#         ix.DE2 = intersect(ix.DE2, ix.highGenes)[1:n0]
#         tmp.sf = sample(tmpParas$sf, tmpSize, replace = T)
#         tmp.p0 = sample(tmpParas$p0, tmpSize, replace = T)
#         tmp.lambda = sample(tmpParas$lambda, tmpSize, replace = T)
#         ## Perturbation of pi.g/zero ratio in phase I; Simulation of Form I DE
#         tmp.pi.g = pi.g
#         tmp.pi.g[ix.DE1] = tmpParas$pi.g[ix.DE1]
#         pi.df = tmp.pi.g - pi.g
#         ## Perturbation of mu/lfc in phase II; Simulation of Form II DE
#         tmp.mu = mu
#         tmp.mu[ix.DE2] = tmpParas$mu[ix.DE2]
#         lfc = tmpParas$mu - tmp.mu
#         tmp.sigma = tmpParas$sd
#         ## Generate count matrix for this cell type
#         tmp.y = GenerateCountMatrix(tmp.pi.g, tmp.p0, tmp.lambda, tmp.mu, sigma, tmp.sf)
#         tmp.rslt = list(Y = tmp.y, ix.DE1 = ix.DE1, ix.DE2 = ix.DE2, pi.g = tmp.pi.g, mu = tmp.mu)
#         simAll[[i]] = tmp.rslt
#     }
#
#     SimulateMultiSCEs = NULL
#     for (i in 1:ncelltype) {
#         for (j in i: ncelltype) {
#             if (j  > i){
#                 comID = c(i, j)
#                 comName = paste0(celltypeNames[i], "_vs_", celltypeNames[j])
#                 Y = do.call(cbind, lapply(simAll[comID], function(x) x$Y))
#                 ix.DE1 = do.call(union, lapply(simAll[comID], function(x) x$ix.DE1))
#                 ix.DE2 = do.call(union, lapply(simAll[comID], function(x) x$ix.DE2))
#                 pi.g = do.call(cbind, lapply(simAll[comID], function(x) x$pi.g))
#                 pi.df = (pi.g[,1] - pi.g[,2])[ix.DE1]
#                 mu = do.call(cbind, lapply(simAll[comID], function(x) x$mu))
#                 lfc = (mu[,1] - mu[,2])[ix.DE2]
#                 rownames(Y) = paste0("g", 1:nrow(Y))
#                 celltypes = rep(paste0("celltype", comID), multiSize[comID])
#                 sce = SingleCellExperiment(
#                     assays = list(counts = Y),
#                     rowData = data.frame(geneNames = rownames(Y), stringsAsFactors = F),
#                     colData = data.frame(cellTypes = celltypes, stringsAsFactors = F)
#                 )
#                 SimulateMultiSCEs[[comName]] = list(sce = sce, ix.DE1 = ix.DE1, ix.DE2 = ix.DE2,
#                                                     pi.df = pi.df, lfc = lfc)
#             }else{
#                 next
#             }
#         }
#     }
#     return(SimulateMultiSCEs)
# }

