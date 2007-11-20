edrcv <- function(x,y,m=2,rho0=1,h0=NULL,ch=exp(.5/max(4,(dim(x)[2]))),
    crhomin=1,cm=4,method="Penalized",basis="Quadratic",cw=NULL,graph=FALSE,
    show=1,trace=FALSE,seed=1,cvsize=1,m0=min(m,2),hsm=NULL){
#
#   first prepare random grouping for Cross-validation 
#
args <- match.call()      
n <- length(y)
ind <- groups(n,cvsize,seed)

#
#  Now the CV-loop
#
require(sm)
if(is.null(hsm)) lh <- 1 else lh <- length(hsm)
cvres <- matrix(0,n,lh)
cat("Start of CV loop \n Progress:")
ngroups <- length(unique(ind))
for(i in 1:ngroups){
xtest <- x[ind==i,]
ytest <- y[ind==i]
edrestimate <- edr(x[ind!=i,],y[ind!=i],m=m,rho0=rho0,h0=h0,ch=ch,
    crhomin=crhomin,cm=cm,method=method,basis=basis,cw=cw,
    show=show,trace=trace)
for(j in 1:lh){
cvfits <- fit.edr(edrestimate,xtest,m=m0,hsm[j])
cvres[ind==i,j] <- ytest - cvfits$fhat
}
cat(" ",signif(i/ngroups*100,2),"%")
}
cat("\n")
# determine the optimal hsm
mae <- apply(abs(cvres),2,mean)
mse <- apply(cvres^2,2,mean)
if(!is.null(hsm)) hsmopt <- hsm[mae==min(mae)] else hsmopt <- NULL
cat("Estimate EDR using the full dataset")
#
#  Now produce all estimates using the full dataset
#
edrestimate <- edr(x,y,m=m,rho0=rho0,h0=h0,ch=ch,
    crhomin=crhomin,cm=cm,method=method,basis=basis,cw=cw,graph=graph,
    show=show,trace=trace)
edrestimate$call <- args
edrestimate$fhat <- fit.edr(edrestimate,x,m=m0,hsmopt)
edrestimate$cvres <- cvres
edrestimate$cvmseofh <- mse
edrestimate$cvmaeofh <- mae
edrestimate$cvmse <- if(lh>1) mse[hsm==hsmopt] else mse
edrestimate$cvmae <- if(lh>1) mae[hsm==hsmopt] else mae
edrestimate$hsm <- hsm
edrestimate$hsmopt <- hsmopt
class(edrestimate) <- "edr"
edrestimate
}

groups <- function(n,size,seed=1){
set.seed(seed)
ox <- sample(1:n,n)
ngroups <- n%/%size
ind <- numeric(n)
k<-1
while(k <= n){
for( i in 1:ngroups ){
ind[ox[k]] <- i
k <- k+1
}
}
ind
}

fit.edr <- function(x,xtest,m=1,h){
    if(!require(sm)) return("Please install package sm")
    xx<-x$x
    yy<-x$y
    n <- length(yy)
    bhat<-x$bhat
    fhat <- x$fhat
    if(m==1) ngrid<-length(yy) else ngrid<-max(100,n^(2/3))
    z <- sm::sm.regression(xx%*%t(edr.R(bhat,m)),yy,h=if(is.null(h)) sm::h.select(xx%*%t(edr.R(bhat,m)),yy,method="cv") else rep(h,m),display="none",eval.points=xtest%*%t(edr.R(bhat,m)),eval.grid=FALSE)
list(x=xtest,fhat=z$estimate)
}
