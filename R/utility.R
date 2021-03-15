## parameter estimation
RobustPoi0 <- function(x){
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

