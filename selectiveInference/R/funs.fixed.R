# Lasso inference function (for fixed lambda). Note: here we are providing inference
# for the solution of
# min 1/2 || y - \beta_0 - X \beta ||_2^2 + \lambda || \beta ||_1

fixedLassoInf <- function(x, y, beta, lambda, family=c("gaussian","binomial","cox"),intercept=TRUE, status=NULL,
                          sigma=NULL, alpha=0.1,
                          type=c("partial","full"), tol.beta=1e-5, tol.kkt=0.1,
                          gridrange=c(-100,100), bits=NULL, verbose=FALSE) {
  
  family = match.arg(family)
  this.call = match.call()
  type = match.arg(type)
  
  if(family=="binomial")  {
    # if(type!="partial") stop("Only type= partial allowed with binomial family")
    out=fixedLogitLassoInf(x,y,beta,lambda,alpha=alpha, type=type, tol.beta=tol.beta, tol.kkt=tol.kkt,
                           gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }
  else if(family=="cox")  {
    if(type!="partial") stop("Only type= partial allowed with Cox family")
    out=fixedCoxLassoInf(x,y,status,beta,lambda,alpha=alpha, type="partial",tol.beta=tol.beta,
                         tol.kkt=tol.kkt, gridrange=gridrange, bits=bits, verbose=verbose,this.call=this.call)
    return(out)
  }
  
  else{
    
    
    
    checkargs.xy(x,y)
    if (missing(beta) || is.null(beta)) stop("Must supply the solution beta")
    if (missing(lambda) || is.null(lambda)) stop("Must supply the tuning parameter value lambda")
    checkargs.misc(beta=beta,lambda=lambda,sigma=sigma,alpha=alpha,
                   gridrange=gridrange,tol.beta=tol.beta,tol.kkt=tol.kkt)
    if (!is.null(bits) && !requireNamespace("Rmpfr",quietly=TRUE)) {
      warning("Package Rmpfr is not installed, reverting to standard precision")
      bits = NULL
    }
    
    # Estimate sigma
    if (is.null(sigma)) {
      if (n >= 2*p) {
        oo = intercept
        sigma = sqrt(sum(lsfit(x,y,intercept=oo)$res^2)/(n-p-oo))
      }
      else {
        sigma = sd(y)
        warning(paste(sprintf("p > n/2, and sd(y) = %0.3f used as an estimate of sigma;",sigma),
                      "you may want to use the estimateSigma function"))
      }
    }
    
    n = nrow(x)
    p = ncol(x)
    beta = as.numeric(beta)
    if (intercept == F & length(beta) != p) stop("Since intercept=FALSE, beta must have length equal to ncol(x)")
    if (intercept == T & length(beta) != p+1) stop("Since intercept=TRUE, beta must have length equal to ncol(x)+1")
    
    
    
    if (intercept == T) {
      bbeta = beta[-1]
      # m=beta[-1]!=0  #active set
      vars = which(abs(bbeta) > tol.beta * sqrt(n / colSums(x^2)))
      xm=cbind(1,x[,vars])
      bhat=c(beta[1],bbeta[vars]) # intcpt plus active vars
    } else {
      bbeta = beta
      # m=beta!=0  #active set
      vars = which(abs(bbeta) > tol.beta * sqrt(n / colSums(x^2)))
      bhat=bbeta[vars] # active vars
      xm=x[,vars]
    }
    xnotm=x[,-vars]
    
    # vars = which(abs(bbeta) > tol.beta / sqrt(colSums(x^2)))
    nvar = length(vars)
    if(nvar==0){
      cat("Empty model",fill=T)
      return()
    }
    
    pv=vlo=vup=sd=rep(NA, nvar)
    vmat = matrix(0,nvar,n)
    ci=tailarea=matrix(NA,nvar,2)
    
    # If glmnet was run with an intercept term, center x and y
    if (intercept==TRUE) {
      obj = standardize(x,y,TRUE,FALSE)
      x = obj$x
      y = obj$y
    }
    
    s2=sign(bhat)
    
    #check KKT
    g = t(x)%*%(y-xm%*%bhat)/lambda # negative gradient scaled by lambda
    # print(g[which(abs(g)>1)])
    if (any(abs(g) > 1+tol.kkt) )
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances)"))
    
    if (any(sign(g[vars]) != sign(bbeta[vars])))
      warning(paste("Solution beta does not satisfy the KKT conditions",
                    "(to within specified tolerances). You might try rerunning",
                    "glmnet with a lower setting of the",
                    "'thresh' parameter, for a more accurate convergence."))
    
    MM = pinv(crossprod(xm))/sigma^2
    # gradient at LASSO solution, first entry is 0 because intercept is unpenalized
    # at exact LASSO solution it should be s2[-1]
    if (intercept == T) gm = c(0,-g[vars]*lambda)
    else gm = -g[vars]*lambda
    
    dbeta = MM%*%gm
    
    bbar = bhat - dbeta
    
    if (intercept == T) {
      A1 = -(mydiag(s2))[-1,]
      b1 = (s2*dbeta)[-1]
      V = diag(length(bbar))[-1,]
    } else {
      A1 = -(mydiag(s2))
      b1 = (s2*dbeta)
      V = diag(length(bbar))
    }
    
    null_value = rep(0,nvar)
    
    # if (p > n) {
      if (type=="full") {
        
        
        # Reorder so that active set is first
        Xordered = cbind(xm,xnotm)
        
        hsigma <- 1/n*(t(Xordered)%*%Xordered)

        M <- InverseLinfty(hsigma,n,dim(xm)[2],verbose=F)[-1,] # remove intercept row
        I <- diag(dim(xm)[2])[-1,]
        if (intercept) Msubset <- M[,-c(1,vars+1)]
        else Msubset <- M[,-vars]
        B <- cbind(I,Msubset/sqrt(n))
        V = B
        
        c <- matrix(c(gm,t(xnotm)%*%xm%*%(-dbeta)),ncol=1)
        d <- -dbeta[-1] # remove intercept
        
        null_value = -(M%*%c/sqrt(n) - d)
        
        A0 = matrix(0,ncol(xnotm),ncol(A1))
        A0 = cbind(A0,diag(nrow(A0)))
        fill = matrix(0,nrow(A1),nrow(A0))
        A1 = cbind(A1,fill)
        A1 = rbind(A1,A0,-A0)
        
        b1 = matrix(c(b1,rep(lambda,2*nrow(A0))),ncol=1)
        
        # full covariance
        MMbr = (crossprod(xnotm) - t(xnotm)%*%xm%*%pinv(crossprod(xm))%*%t(xm)%*%xnotm)*sigma^2
        MM = cbind(MM,matrix(0,nrow(MM),ncol(MMbr)))
        MMbr = cbind(matrix(0,nrow(MMbr),nrow(MM)),MMbr)
        MM = rbind(MM,MMbr)

        gnotm = g[-vars]*lambda
        bbar = matrix(c(bbar,gnotm),ncol=1)
      }
    # }
    
    # Check polyhedral region
    tol.poly = 0.01
    if (max((A1 %*% bbar) - b1) > tol.poly)
      stop(paste("Polyhedral constraints not satisfied; you must recompute beta",
                 "more accurately. With glmnet, make sure to use exact=TRUE in coef(),",
                 "and check whether the specified value of lambda is too small",
                 "(beyond the grid of values visited by glmnet).",
                 "You might also try rerunning glmnet with a lower setting of the",
                 "'thresh' parameter, for a more accurate convergence."))
    
    sign = numeric(nvar)
    coef0 = numeric(nvar)
    
    for (j in 1:nvar) {
      if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))
      
      vj = V[j,]
      # mj = sqrt(sum(vj^2))
      # vj = vj/mj
      coef0[j] = sum(vj * bbar)
      sign[j] = sign(coef0[j])
      vj = vj * sign[j]
      
      limits.info = TG.limits(bbar, A1, b1, vj, Sigma=MM)
      if(is.null(limits.info)) return(list(pv=NULL,MM=MM,eta=vj))
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
    
    out = list(type=type,lambda=lambda,pv=pv,ci=ci,
               tailarea=tailarea,vlo=vlo,vup=vup,y=y,
               vars=vars,sign=sign,sigma=sigma,alpha=alpha,
               sd=sd,
               coef0=coef0,
               call=this.call)
    class(out) = "fixedLassoInf"
    return(out)
  }
}

#############################


fixedLasso.poly=
  function(x, y, beta, lambda, a, inactive = FALSE) {
    xa = x[,a,drop=F]
    xac = x[,!a,drop=F]
    xai = pinv(crossprod(xa))
    xap = xai %*% t(xa)
    za = sign(beta[a])
    if (length(za)>1) dz = diag(za)
    if (length(za)==1) dz = matrix(za,1,1)
    
    if (inactive) {
      P = diag(1,nrow(xa)) - xa %*% xap
      
      G = -rbind(
        1/lambda * t(xac) %*% P,
        -1/lambda * t(xac) %*% P,
        -dz %*% xap
      )
      lambda2=lambda
      if(length(lambda)>1) lambda2=lambda[a]
      u = -c(
        1 - t(xac) %*% t(xap) %*% za,
        1 + t(xac) %*% t(xap) %*% za,
        -lambda2 * dz %*% xai %*% za)
    } else {
      G = -rbind(
        #   1/lambda * t(xac) %*% P,
        # -1/lambda * t(xac) %*% P,
        -dz %*% xap
      )
      lambda2=lambda
      if(length(lambda)>1) lambda2=lambda[a]
      u = -c(
        #   1 - t(xac) %*% t(xap) %*% za,
        #   1 + t(xac) %*% t(xap) %*% za,
        -lambda2 * dz %*% xai %*% za)
    }
    
    return(list(G=G,u=u))
  }

##############################

### Functions borrowed and slightly modified from lasso_inference.R

## Approximates inverse covariance matrix theta
InverseLinfty <- function(sigma, n, e, resol=1.2, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE, max.try=10) {
    isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }
  
  p <- nrow(sigma);
  M <- matrix(0, e, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:e) {
    if ((i %% xp)==0){
      xperc = xperc+10;
      if (verbose) {
        print(paste(xperc,"% done",sep="")); }
    }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;

    output = NULL

    while ((mu.stop != 1) && (try.no<max.try) ){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, soln_result=output) # uses a warm start
      beta <- output$soln
      iter <- output$iter
      if (isgiven==1) {
        mu.stop <- 1
      }
      else{
        if (try.no==1){
          if (iter == (maxiter+1)){
            incr <- 1;
            mu <- mu*resol;
          } else {
            incr <- 0;
            mu <- mu/resol;
          }
        }
        if (try.no > 1){
          if ((incr == 1)&&(iter == (maxiter+1))){
            mu <- mu*resol;
          }
          if ((incr == 1)&&(iter < (maxiter+1))){
            mu.stop <- 1;
          }
          if ((incr == 0)&&(iter < (maxiter+1))){
            mu <- mu/resol;
          }
          if ((incr == 0)&&(iter == (maxiter+1))){
            mu <- mu*resol;
            beta <- last.beta;
            mu.stop <- 1;
          }
        }
      }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

InverseLinftyOneRow <- function (Sigma, i, mu, maxiter=50, soln_result=NULL, kkt_tol=1.e-6, objective_tol=1.e-6,
		                 use_QP=TRUE) {

  # If soln_result is not Null, it is used as a warm start.
  # It should be a list
  # with entries "soln", "gradient", "ever_active", "nactive"

  p = nrow(Sigma)

  if (is.null(soln_result)) {
     soln = rep(0, p)
     ever_active = rep(0, p)
     ever_active[1] = i      # 1-based
     ever_active = as.integer(ever_active)
     nactive = as.integer(1)
     if (use_QP) {
          linear_func = rep(0, p)
	  linear_func[i] = -1
	  linear_func = as.numeric(linear_func)
	  gradient = 1. * linear_func 
     } else {
          gradient = rep(0, p)
     }
  }
  else {
     soln = soln_result$soln
     gradient = soln_result$gradient  
     ever_active = as.integer(soln_result$ever_active)
     nactive = as.integer(soln_result$nactive)
     if (use_QP) { 
         linear_func = soln_result$linear_func
     }
  }

  if (use_QP) {
      result = solve_QP(Sigma, mu, maxiter, soln, linear_func, gradient, ever_active, nactive, kkt_tol, objective_tol) 
  } else {
      result = find_one_row_debiasingM(Sigma, i, mu, maxiter, soln, gradient, ever_active, nactive, kkt_tol, objective_tol) # C function uses 1-based indexing for active set
  }

  # Check feasibility

  if (!result$kkt_check) {
     warning("Solution for row of M does not seem to be feasible")
  } 

  return(result)

}

##############################

print.fixedLassoInf <- function(x, tailarea=TRUE, ...) {
  cat("\nCall:\n")
  dput(x$call)

  cat(sprintf("\nStandard deviation of noise (specified or estimated) sigma = %0.3f\n",
              x$sigma))

  cat(sprintf("\nTesting results at lambda = %0.3f, with alpha = %0.3f\n",x$lambda,x$alpha))
  cat("",fill=T)
  tab = cbind(x$vars,
    round(x$coef0,3),
    round(x$coef0 / x$sd,3),
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

#estimateLambda <- function(x, sigma, nsamp=1000){
#  checkargs.xy(x,rep(0,nrow(x)))
#  if(nsamp < 10) stop("More Monte Carlo samples required for estimation")
#  if (length(sigma)!=1) stop("sigma should be a number > 0")
 # if (sigma<=0) stop("sigma should be a number > 0")

 # n = nrow(x)
 # eps = sigma*matrix(rnorm(nsamp*n),n,nsamp)
 # lambda = 2*mean(apply(t(x)%*%eps,2,max))
 # return(lambda)
#}

