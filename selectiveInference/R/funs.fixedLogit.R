
fixedLogitLassoInf=function(x,y,beta,lambda,alpha=.1, type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                            gridrange=c(-100,100), bits=NULL, verbose=FALSE,this.call=NULL){
  
  
  type = match.arg(type)
  checkargs.xy(x,y)
  if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
  if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda") 
  checkargs.misc(beta=beta,lambda=lambda,alpha=alpha,
                 gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
  if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
    warning("Package Rmpfr is not installed, reverting to standard precision")
    bits = NULL
  }
  
  
  n=length(y)
  p=ncol(x)
  # I assume that intcpt was used
  if(length(beta)!=p+1) stop("Since family='binomial', beta must be of length ncol(x)+1, that is, it should include an intercept")
  vars = which(abs(beta[-1]) > tol.beta / sqrt(colSums(x^2)))
  nvar=length(vars)
  pv=vlo=vup=sd=rep(NA, nvar)
  ci=tailarea=matrix(NA,nvar,2)
  
  #do we need to worry about standardization?
  
  #  obj = standardize(x,y,TRUE,FALSE)
  #  x = obj$x
  #  y = obj$y
  
  # m=beta[-1]!=0  #active set
  
  bhat=c(beta[1],beta[-1][vars]) # intcpt plus active vars
  s2=sign(bhat)
  lam2m=diag(c(0,rep(lambda,nvar)))
  
  
  xm=cbind(1,x[,vars])
  xnotm=x[,-vars]
  
  etahat = xm %*% bhat
  prhat = as.vector(exp(etahat) / (1 + exp(etahat)))
  ww=prhat*(1-prhat)
  w=diag(ww)
  
  #check KKT
  z=etahat+(y-prhat)/ww
  # g=  t(x)%*%w%*%(z-etahat)/lambda # negative gradient scaled by lambda
  g=scale(t(x),FALSE,1/ww)%*%(z-etahat)/lambda # negative gradient scaled by lambda
  if (any(abs(g) > 1+tol.kkt) )
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances)"))
  
  if(length(vars)==0){
    cat("Empty model",fill=T)
    return()
  }
  if (any(sign(g[vars]) != sign(beta[-1][vars])))
    warning(paste("Solution beta does not satisfy the KKT conditions",
                  "(to within specified tolerances). You might try rerunning",
                  "glmnet with a lower setting of the",
                  "'thresh' parameter, for a more accurate convergence."))
  # warnings(paste(g[vars],beta[-1][vars]))
  
  #constraints for active variables             
  # MM=solve(t(xxm)%*%w%*%xxm)
  MM=solve(scale(t(xm),F,1/ww)%*%xm)
  gm = c(0,-g[vars]*lambda) # gradient at LASSO solution, first entry is 0 because intercept is unpenalized
  # at exact LASSO solution it should be s2[-1]
  dbeta = MM %*% gm
  
  # bbar=(bhat+lam2m%*%MM%*%s2)  # JT: this is wrong, shouldn't use sign of intercept anywhere...
  bbar = bhat - dbeta
  
  A1= matrix(-(mydiag(s2))[-1,],nrow=length(s2)-1)
  b1= (s2 * dbeta)[-1]
  V = diag(length(bbar))[-1,]
  null_value = rep(0,nvar)
  
  if (type=='full') {
    Xordered = cbind(xm,xnotm)
    
    hsigma <- 1/n*(t(Xordered)%*%w%*%Xordered)
    
    M <- matrix(InverseLinfty(hsigma,n,dim(xm)[2],verbose=F,max.try=10),ncol=p+1)[-1,] # remove intercept row
    # I <- matrix(diag(dim(xm)[2]),nrow=dim(xm)[2])
    I <- matrix(diag(dim(xm)[2])[-1,],nrow=dim(xm)[2]-1)
    if (is.null(dim(M))) Msubset <- M[-c(1,vars+1)]
    else Msubset <- M[,-c(1,vars+1)]
    # Msubset = matrix(Msubset,nrow=dim(xm)[2])
    # Msubset = matrix(Msubset,nrow=dim(xm)[2]-1)
    # print("*********")
    # print(Msubset)
    # print(I)
    # print("*********")
    V <- cbind(I,Msubset/sqrt(n)) # matrix(cbind(I,Msubset/sqrt(n)),nrow=dim(xm)[2]-1)
    
    c <- matrix(c(gm,t(xnotm)%*%w%*%xm%*%(-dbeta)),ncol=1)
    d <- -dbeta[-1] # remove intercept
    
    # d <- matrix(-dbeta,ncol=1) # remove intercept
    # d <- matrix(-dbeta*lambda,ncol=1) # remove intercept
    # print(dim(xnotm))
    # print(dim(w))
    # print(dim(xm))
    # c <- matrix(c(gm,t(xnotm)%*%w%*%xm%*%d),ncol=1)
    # c <- matrix(c(gm[-1],t(xnotm)%*%w%*%xm%*%d),ncol=1)
    # d <- d[-1,]
    
    # print(dim(c))
    # print(dim(M))
    # print(length(d))
    
    null_value = -(M%*%c/sqrt(n) - d)
    # null_value = -(M[,-1]%*%c/sqrt(n) - d)
    
    # A0 = (t(xnotm)%*%w%*%xm)/lambda
    # A0 = cbind(A0,matrix(0,nrow(A0),ncol(xnotm)))
    # A0 = rbind(A0,-A0)
    # 
    # b0 = matrix(t(xnotm)%*%w%*%(z/lambda+xm%*%MM%*%gm),ncol=1)
    # 
    # 
    # print("------")
    # print(dim(A0))
    # print(dim(b0))
    # print(length(b1))
    
    
    # b1 = rbind(1+b0,1-b0,matrix(b1,ncol=1))
    
    A0 = matrix(0,ncol(xnotm),ncol(A1))
    A0 = cbind(A0,diag(nrow(A0)))
    fill = matrix(0,nrow(A1),ncol(xnotm))
    A1 = cbind(A1,fill)
    A1 = rbind(A1,A0,-A0)
    
    b1 = matrix(c(b1,rep(lambda,2*nrow(A0))),ncol=1)
    
    # full covariance
    MMbr = (t(xnotm)%*%w%*%xnotm - t(xnotm)%*%w%*%xm%*%MM%*%t(xm)%*%w%*%xnotm) # matrix(0,ncol(xnotm),ncol(xnotm))
    MM = cbind(MM,matrix(0,nrow(MM),ncol(MMbr)))
    MMbr = cbind(matrix(0,nrow(MMbr),nrow(MM)),MMbr)
    # print(dim(MM))
    # print(dim(MMbr))
    MM = rbind(MM,MMbr)
    
    gnotm = g[-vars]*lambda
    bbar = matrix(c(bbar,gnotm),ncol=1)
    
    # print(dim(A1))
    # print(dim(b1))
    # print(dim(bbar))
    # print(dim(V))
  }
  
  
  if (is.null(dim(V))) V=matrix(V,nrow=1)
  
  tol.poly = 0.01 
  if (max((A1 %*% bbar) - b1) > tol.poly)
    stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
               "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
               "and check whether the specified value of lambda is too small",
               "(beyond the grid of values visited by glmnet).",
               "You might also try rerunning glmnet with a lower setting of the",
               "'thresh' parameter, for a more accurate convergence."))
  
  
  sign=numeric(nvar)
  coef0=numeric(nvar)
  
  
  for(j in 1:nvar){
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
    
    if (is.null(dim(V))) vj = V
    else vj = matrix(V[j,],nrow=1)
    coef0[j] = vj%*%bbar # sum(vj * bbar)
    sign[j] = sign(coef0[j])
    vj = vj * sign[j]
    # print(dim(MM))
    # print(dim(vj))
    # print(nvar)
    
    # compute p-values
    limits.info = TG.limits(bbar, A1, b1, vj, Sigma=MM)
    # if(is.null(limits.info)) return(list(pv=NULL,MM=MM,eta=vj))
    a = TG.pvalue.base(limits.info, null_value=null_value[j], bits=bits)
    pv[j] = a$pv
    vlo[j] = a$vlo # * mj # Unstandardize (mult by norm of vj)
    vup[j] = a$vup # * mj # Unstandardize (mult by norm of vj)
    sd[j] = a$sd # * mj # Unstandardize (mult by norm of vj)
    
    a = TG.interval.base(limits.info, 
                         alpha=alpha,
                         gridrange=gridrange,
                         flip=(sign[j]==-1),
                         bits=bits)
    ci[j,] = (a$int-null_value[j]) # * mj # Unstandardize (mult by norm of vj)
    tailarea[j,] = a$tailarea
  }
  
  se0 = sqrt(diag(V%*%MM%*%t(V)))
  zscore0 = (coef0+null_value)/se0
  
  out = list(type=type,lambda=lambda,pv=pv,ci=ci,
             tailarea=tailarea,vlo=vlo,vup=vup,sd=sd,
             vars=vars,alpha=alpha,coef0=coef0,zscore0=zscore0,
             call=this.call,
             info.matrix=MM) # info.matrix is output just for debugging purposes at the moment
  class(out) = "fixedLogitLassoInf"
  return(out)
  
}



print.fixedLogitLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)
  
  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))
  
  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
              round(x$coef0,3),
              round(x$zscore0,3),
              round(x$pv,3),round(x$ci,3))
  colnames(tab) = c("Var", "Coef", "Z-score", "P-value", "LowConfPt", "UpConfPt")
  if (tailarea) {
    tab = cbind(tab,round(x$tailarea,3))
    colnames(tab)[(ncol(tab)-1):ncol(tab)] = c("LowTailArea","UpTailArea")
  }
  rownames(tab) = rep("",nrow(tab))
  print(tab)
  
  cat(sprintf("\nNote: coefficients shown are %s regression coefficients\n",
              ifelse(x$type=="partial","partial","full")))
  invisible()
}



