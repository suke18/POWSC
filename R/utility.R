## parameter estimation
RobustPoi0 = function(x){
    tt=table(x)
    n0=sum(x==0)
    if (n0 == length(x)){
        c("p0"=1,"mu"=0)
    }else{
        if (names(tt)[1]=="0") {
            xx=as.numeric(names(tt))[-1]
            tt=as.vector(tt)[-1]
        }else{
            xx=as.numeric(names(tt))
            tt=as.vector(tt)
        }
        tt=log(tt)+log(gamma(1+xx))
        ##fit the regression without N0
        beta=lm(tt~xx,weight=1/exp(xx))$coef
        mu.hat=exp(beta[2])
        p0=(n0-exp(beta[1]))/(exp(beta[1]+mu.hat)+n0-exp(beta[1]))
        if (any(p0<0)) {warning("Proportion of zero inflation is negative")}
        p0 <- pmax(p0, 0)
        c("p0"=p0,"mu"=mu.hat)
    }
}


## Shrink the mu for sample sizeas are so small
shrink.mu=function(y,s,n){
    mu.g=rep(NA,length(y))
    k=which(n>1)
    if (length(k)<length(n)) {fill=TRUE} else {fill=FALSE}
    s=s[k];y=y[k];n=n[k]
    
    mu0=weighted.mean(y,w=n)
    
    s2.total=sum(s^2*(n-1))+sum(n*(y-mu0)^2)
    s2.total=s2.total/sum(n)
    
    s2.1=sum(s^2*(n-1))/sum(n)
    s2.0=s2.total-s2.1
    ### shrink mu
    mu.sub=  (y*n/s2.1+mu0/s2.0)/(n/s2.1+1/s2.0)
    mu.g[k]=mu.sub
    if (fill) mu.g[-k]=mu0
    
    mu.g
}

my.mad=function (x, center = median(x), constant = 1.4826, na.rm = FALSE){
    if (na.rm)
        x <- x[!is.na(x)]
    res=x-center
    constant * median(res[res>0])
}

dZinf.pois=function(x, p0, mu){
    (x==0)*(p0+(1-p0)*dpois(x,mu))+(x>0)*(1-p0)*dpois(x,mu)
}

dLNP2 <- function(x, mu, sigma, l=1){
    x.min <- pmax(0, x-0.5)
    pnorm(log2((x+0.5)/l), mu, sigma) - pnorm(log2(x.min/l), mu, sigma)
}




### for the use of SC2P
eset2Phase <- function(eset, low.prob=0.99){  ## takes eSet as input
    Y <- round(exprs(eset))
    #################################################
    ## ## initial estimate of prob(X\in bg)
    ##################################################
    Cell0=colMeans(Y==0) # each cell has this percentage 0
    par1=apply(Y,2,function(yy) {
        yy=yy[yy<=15]
        RobustPoi0(yy)}
    )
    pi0.hat=Cell0/(par1[1,]+(1-par1[1,])*dpois(0,par1[2,]))
    if (any((pi0.hat > 1))) {warning("Zero proportion is greater than estimation.")}
    pi0.hat <- pmin(pi0.hat, 1)
    prob0=pi0.hat*par1[1,]+ pi0.hat*(1-par1[1,])*dpois(0,par1[2,]) ## ZIP prob at 0
    ############################################
    ## First round
    ###########################################
    ## get the 1-low.prob quantile of ZIP
    x0=qpois(pmax(1-(1-low.prob)/(1-par1[1,]),0),par1[2,])
    Z= sweep(Y,2,x0)>0 # indicate if a gene is > bg
    L=colSums(Y*Z)/1e6 # so far it is like simple total..

    mu.g1=log2(rowSums(Z*Y)/rowSums(sweep(Z,2,L,FUN="*")))
    mu.g1[is.na(mu.g1)]=0 ## if allZ is 0, it gets NA,
    ### but we should shrink mu.g1 as well since some mu.g1 is estimated by only a few observations
    ## leave it here for now.
    n.g1=rowSums(Z)
    y1=log2(sweep(Y,2,L,FUN="/")+1) #like CPM**
    s.g1=sqrt(rowSums(Z*sweep(y1,1,mu.g1)^2)/(n.g1-1)) ## CPM type of SD
    mu.g2 = shrink.mu(mu.g1,s.g1,n.g1)
    ###############################################
    ## get sd.g
    ############################################
    res.g1=log2(sweep(Y,2,L,FUN="/")+1)-mu.g1
    ## mad of those res.g1 that are associated with Z==1
    tmp=array(0,dim=c(dim(res.g1),2))
    tmp[,,1]=res.g1;tmp[,,2]=Z
    sd.g1=apply(tmp,1,function(xx) my.mad(xx[xx[,2]==1,1]))
    sd.g1[is.na(sd.g1)]=0## if all bg, there's no info about fg sd
    ## add a shrinkage for sd.g1
    sd.prior=squeezeVar(sd.g1^2,n.g1-1)
    sd.g2=sqrt(sd.prior$var.post)
    ####################################### ########
    #####  gene specific bg. Z_gi
    #######################
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
    ##################################################
    ## compute offsets
    ##################################################
    Offset = Y*0
    Ylim=range(log2(1+Y)-mu.g1);Xlim=range(mu.g1)

    for(i in 1:ncol(Y)){
        tmp.y=log2(1+Y[,i])-mu.g2
        subset= post.Z2[,i] > .99
        tmp.Z2 = post.Z2[, i]

        ## lm1 <- loess(tmp.y~mu.g1,
        ##              weights=post.Z2[,i]*mu.g2,subset=subset,degree=1,span=.3)
        if (sum(subset) < 2)
            next
        else
            lm1 <- tryCatch(expr = loess(tmp.y ~ mu.g1, weights = tmp.Z2 * mu.g2,
                                         subset = subset, degree = 1, span = span),
                            error = function(e) loess(tmp.y ~ mu.g1, weights = tmp.Z2 * mu.g2,
                                                      subset = subset, degree = 1, span = 0.8))

        Offset[subset,i]=lm1$fitted
        ## par(mfrow=c(1,2))
        ## plot(mu.g1, log2(1+Y[,i])-mu.g1, pch=16,cex=.6,ylab="",
        ##      col=rgb(1-post.Z2[,i],0,post.Z2[,i],alpha=rowMeans(post.Z2))
        ##     ,ylim=Ylim,xlim=Xlim,main=i)
        ## points(lm1$x,lm1$fitted,col=5)

        ## plot(mu.g1, log2(1+Y[,i])-Offset[,i]-mu.g1, pch=16,cex=.6,ylab="",
        ##      col=rgb(1-post.Z2[,i],0,post.Z2[,i],alpha=rowMeans(post.Z2))
        ##     ,ylim=Ylim,xlim=Xlim,main=i)
        ## tmp.y2 <- tmp.y-Offset[,i]
        ## lm2 <- loess(tmp.y2 ~ mu.g1,
        ##              weights=post.Z2[,i]*mu.g2, subset=subset,degree=1,span=.3)
        ## points(lm2$x,lm2$fitted,col=5)
    }
    ##################################################
    ## assemble the estimators into sc2pSet object
    ##################################################
    ## add mu and sd to feature data
    fdata <- fData(eset)
    fdata2 <- as.data.frame(cbind(fdata, mu.g2, sd.g2))
    colnames(fdata2) <- c(colnames(fdata), "mean", "sd")
    fvar <- rbind(fvarMetadata(eset), "mean"="shrinkage estimated foreground mean",
                  "sd"="shrinkage estimated foreground standard deviation")
    featureData <- new("AnnotatedDataFrame", data=fdata2,
                       varMetadata=fvar)
    ## add lambda and p0 to phenoData
    pdata <- pData(eset)
    pdata2 <- as.data.frame(cbind(pdata, par1[1,], par1[2,], L))
    colnames(pdata2) <- c(colnames(pdata), "p0", "lambda", "L")
    pvar <-rbind(varMetadata(eset), "p0"="proportion of zero inflation",
                 "lambda"="mean of background poisson",
                 "L"="foreground library size")
    phenoData <- new("AnnotatedDataFrame", data=pdata2, varMetadata=pvar)

    out <- list(exprs=Y, Z=post.Z2, Offset=Offset,
               phenoData=phenoData,
               featureData=featureData,
               experimentData=experimentData(eset),
               annotation=annotation(eset))
    out
}



##################################################
twoPhaseDE0 <- function(Y, Z, X, Offset, test.which,
                        low.prob=.99){
    vars <- colnames(X);  vars0 <- vars[-test.which]
    group <- X[, test.which]
    group = as.factor(group)
    if (!is.factor(group)){ stop("The variable to be tested must be a factor") }
    group <- droplevels(group)
    X[, test.which] <- group ## just putting it back for Phase II
    Ng <- length(levels(group))
    parse1 <- parse(text= paste0("glm(yyy~", paste(vars, collapse="+"),
                                 ",data=X, family=binomial)"))
    contrast <- paste0(vars[test.which], levels(group)[Ng])

    ################################################
    ## phase 1: change of on rate
    ################################################
    Z1=( Z > low.prob)^2
    ## avgZ=tapply(1:ncol(Z), group, function(ind){
    ##     rowMeans(Z[,ind]) })
    ## avgZ=matrix(unlist(avgZ),ncol=Ng)
    n.on=rowSums(Z1)
    ind=which(n.on > 0 & n.on < ncol(Y))
    if (length(vars)==1){ ## single binary variable
        DE.z <- matrix(NA, nrow=nrow(Z), ncol=4)
        rownames(DE.z) <- rownames(Y)
        DE.z[ind, ]= t(apply(Z1[ind,], 1,function(yyy){
            fit <- eval(parse1)
            ss <- summary(fit)
            coef <- ss$coef[contrast, "Estimate"]
            pval <- pchisq(ss$null.deviance - ss$deviance,
                           df=ss$df.null - ss$df.residual,
                           lower.tail=FALSE)
            c(mean(yyy[group==levels(group)[1]]),
              mean(yyy[group==levels(group)[2]]),
              coef, pval)
        }))
        colnames(DE.z) <- c("p1", "p2", "Ph1.coef", "Ph1.pval")
    } else {
        DE.z <- matrix(NA, nrow=nrow(Z), ncol=2)
        rownames(DE.z) <- rownames(Y)
        parse0 <- parse(text= paste0("glm(yyy~", paste(vars0, collapse="+"),
                                     ",data=X, family=binomial)"))
        DE.z[ind, ]= t(apply(Z1[ind,], 1,function(yyy){
            fit1 <- eval(parse1); fit0 <- eval(parse0)
            ss1 <- summary(fit1); ss0 <- summary(fit0)
            coef <- ss1$coef[contrast, "Estimate"]
            pval <- pchisq(ss0$deviance - ss1$deviance,
                           df=ss0$df.residual - ss1$df.residual,
                           lower.tail=FALSE)
            c(coef, pval)
        }))
        colnames(DE.z) <- c("Ph1.coef", "Ph1.pval")
    }

    ##################################################
    ## phase 2: conditional FC
    ################################################
    W <- log2(Y+1) - Offset; W[!Z1] <- NA
    modelX <- eval(parse(text=paste0("model.matrix(~", paste(vars, collapse="+"),
                                     ", data=X)")))
    fit <- lmFit(W, modelX)
    fit <- eBayes(fit)
    mu1 <- fit$coef[,"(Intercept)"]
    mu2 <- rowSums(fit$coef)
    coef <- fit$coef[, contrast]
    stdev <- fit$stdev.unscaled[, contrast]
    delta <- qt(.975, fit$df.residual + fit$df.prior)*stdev*sqrt(fit$s2.post)
    ci.lo <- coef - delta; ci.hi <- coef + delta
    DE.y <- cbind(mu1, mu2, coef, ci.lo, ci.hi, fit$p.value[, contrast])
    colnames(DE.y) <- c("m1", "m2", "Ph2.coef",
                        "Ph2.ci.lo", "Ph2.ci.hi", "Ph2.pval")

    ## ################################################
    ## marginal logFC change
    sel1 <- as.numeric(group)==1
    sel2 <- as.numeric(group)==2
    logFC <- apply(log2(Y+1), 1, function(y){
        mean(y[sel2], na.rm=TRUE) - mean(y[sel1], na.rm=TRUE)
    })
    as.data.frame(cbind(DE.z, DE.y, logFC))
}


