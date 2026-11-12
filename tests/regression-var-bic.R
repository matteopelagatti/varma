library(varma)
internal <- function(name) get(name,asNamespace("varma"))
check <- function(name,expr) { force(expr); cat("PASS:",name,"\n") }
close <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),"error"))

check("BIC uses common observations and agrees with independent OLS",{
  set.seed(781);Y <- sim_varma(varma_paper_designs()$model2,180)
  P <- 8L; m <- ncol(Y); N <- nrow(Y)-P
  for(icpt in c(FALSE,TRUE)) {
    a <- fit_var_bic(Y,p_max=P,intercept=icpt)
    bic <- vapply(seq_len(P),function(p) {
      # Independent construction using embed, not the package lag builder.
      Z <- embed(Y,P+1L)
      X <- Z[,m+seq_len(m*p),drop=FALSE]
      if(icpt) X<-cbind(X,1)
      E <- qr.resid(qr(X),Z[,seq_len(m),drop=FALSE])
      as.numeric(determinant(crossprod(E)/N,logarithm=TRUE)$modulus)+
        log(N)/N*(m*m*p+m*as.integer(icpt))
    },numeric(1))
    close(a$lag_selection$bic,bic)
    stopifnot(a$selected_order==which.min(bic),all(a$lag_selection$nobs==N),
      a$nobs==nrow(Y)-a$selected_order,a$lag_selection_nobs==N)
    p <- a$selected_order;Z <- embed(Y,p+1L);X <- Z[,-seq_len(m),drop=FALSE]
    if(icpt) X<-cbind(1,X)
    coef <- qr.coef(qr(X),Z[,seq_len(m),drop=FALSE])
    dynamic <- if(icpt) coef[-1,,drop=FALSE] else coef
    close(a$ar,array(t(dynamic),c(m,m,p)))
    E <- qr.resid(qr(X),Z[,seq_len(m),drop=FALSE])
    close(a$cov,crossprod(E)/nrow(E))
    # Independent recursive forecasts from the estimated coefficient matrices.
    history <- Y; pred <- matrix(0,4,m)
    for(h in 1:4) {
      value <- if(icpt) coef[1,] else numeric(m)
      for(j in seq_len(p)) value <- value+a$ar[,,j]%*%history[nrow(history)-j+1L,]
      pred[h,] <- value; history <- rbind(history,as.numeric(value))
    }
    close(predict(a,4,cov=FALSE),pred)
  }
})

check("automatic bounds grow, are capped safely, and remain configurable",{
  set.seed(88)
  bounds <- vapply(c(100,200,400),function(n)
    fit_var_bic(matrix(rnorm(n*3),n,3),intercept=FALSE)$lag_selection_p_max,integer(1))
  stopifnot(identical(bounds,c(10L,13L,15L)))
  Y <- matrix(rnorm(100),20,5)
  a <- fit_var_bic(Y,intercept=TRUE)
  stopifnot(a$lag_selection_p_max==2L,a$lag_selection_capped)
  fails(fit_var_bic(Y,p_max=3,intercept=TRUE))
  b <- fit_var_bic(matrix(rnorm(200),100,2),lag_power=1)
  stopifnot(b$lag_selection_p_max==ceiling(log(100)))
  fails(fit_var_bic(Y,p_min=4,p_max=2));fails(fit_var_bic(Y,p_max=1.5))
  fails(fit_var_bic(Y,lag_power=NA));fails(fit_var_bic(Y,intercept=1))
})

check("white noise, scalar series, singular data and unstable OLS are handled",{
  set.seed(119);Y <- matrix(rnorm(300),150,2)
  for(icpt in c(FALSE,TRUE)) {
    w <- fit_var_bic(Y,p_min=0,p_max=0,intercept=icpt)
    stopifnot(w$selected_order==0,all(dim(w)==c(2,0,0)))
    target <- if(icpt) colMeans(Y) else numeric(2)
    close(predict(w,3,cov=FALSE),matrix(rep(target,each=3),3,2))
    close(w$cov,crossprod(sweep(Y,2,target,"-"))/nrow(Y))
  }
  a <- fit_var_bic(Y[,1],p_max=4)
  stopifnot(all(is.finite(irf(a,5))),all(is.finite(predict(a,3,cov=FALSE))))
  fails(fit_var_bic(cbind(Y[,1],Y[,1]),p_max=3))
  trend <- matrix(1.05^(1:100),ncol=1)
  a <- fit_var_bic(trend,p_min=1,p_max=1,intercept=FALSE)
  stopifnot(!is_stable(a)[1]);close(a$ar,1.05)
})

check("bootstrap repeats VAR lag selection and summaries disclose boundary selection",{
  set.seed(61);Y <- sim_varma(varma_paper_designs()$model1,200)
  a <- fit_var_bic(Y,p_max=5,intercept=FALSE)
  b <- bootstrap_varma(a,nsim=1,n=150,parallel=FALSE)[[1]]
  stopifnot(b$estimation_method=="fit_var_bic",b$lag_selection_nobs==145L,
    nrow(b$lag_selection)==5L)
  d <- data.frame(selected_p=c(1L,5L,NA_integer_),var_p_max=c(5L,5L,NA_integer_),
    var_bound_capped=c(FALSE,FALSE,NA))
  s <- internal(".varma_bic_selection_summary")(d,"model1",100)
  stopifnot(s$summary$attempted==3,s$summary$returned==2,
    s$summary$mean_p==3,s$summary$at_max_fraction==.5,
    sum(s$frequencies$count)==3)
})
cat("All VAR BIC checks passed.\n")
