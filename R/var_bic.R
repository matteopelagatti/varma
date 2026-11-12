#' Fit a VAR with lag order selected by BIC
#'
#' Selects a pure VAR by Gaussian BIC on a common sample and refits the
#' selected order by equation-by-equation ordinary least squares. No VARMA
#' population orders, parameters or future observations enter the selection.
#' @param Y finite numeric n by m data matrix; rows are observations.
#' @param p_max maximum candidate lag, or NULL for ceiling(log(n)^lag_power),
#' capped to leave more than m residual degrees of freedom per equation.
#' An explicitly supplied infeasible bound raises an error.
#' @param intercept logical; include the same intercept specification in every
#' candidate. The simulation default is FALSE, matching the other estimators.
#' @param p_min smallest candidate lag (default 1); use 0 to include white noise.
#' @param lag_power positive exponent for the automatic natural-logarithmic
#' bound (default 1.5). Unused for an explicitly supplied p_max.
#' @details Let P be the maximum candidate lag and N=n-P. Every candidate is
#' fitted on observations P+1 through n. The selection criterion is
#' log(det(S_p)) + log(N)/N * (m*m*p + m*intercept), where S_p is the
#' residual cross-product divided by N. Covariance parameter counts are
#' constant across candidates and are omitted. Exact ties favor the smaller
#' order. Rank-deficient regressions and singular covariance estimates are
#' recorded as inadmissible candidates, never assigned an artificially good BIC.
#'
#' After selection, the model is refitted on observations p+1 through n.
#' This matches the common-sample SC criterion and subsequent OLS refit used
#' in conventional VAR implementations. The default search begins at p=1.
#' Stability is reported, not imposed by rejecting or shrinking OLS estimates.
#'
#' The automatic bound uses a log-power envelope motivated by the discussion
#' of VAR approximations to VARMA models in Kuersteiner (2005, Section 2).
#' The exponent 1.5 is a configurable implementation choice, not a uniquely
#' prescribed constant or a theorem guaranteeing BIC post-selection inference.
#' In particular, a growing upper bound alone does not guarantee that the
#' selected lag removes approximation bias at the rate required for inference.
#' @return A varma object with q=0, selected_order, lag_selection (a table of
#' candidate BICs and diagnostics), lag_selection_nobs, lag_selection_p_max,
#' lag_selection_capped, and refit_args. The reported model loglik uses the
#' refit sample; it need not reproduce the common-sample selection BIC.
#' Bootstrap refits repeat BIC selection with the same search rule.
#' @references Kuersteiner, G. M. (2005). Automatic inference for infinite order
#' vector autoregressions. Econometric Theory 21(1), 85--115.
#' \url{https://www.econ.umd.edu/sites/www.econ.umd.edu/files/pubs/automaticvarinf10.pdf}
#' @export
#' @examples
#' set.seed(24)
#' y <- matrix(rnorm(400),200,2)
#' fit <- fit_var_bic(y)
#' fit$selected_order
#' fit$lag_selection
fit_var_bic <- function(Y,p_max=NULL,intercept=TRUE,p_min=1L,lag_power=1.5) {
  Y <- as.matrix(Y)
  if(!is.numeric(Y) || is.complex(Y) || !length(Y) || any(!is.finite(Y)))
    stop("Y must be a nonempty finite real numeric matrix")
  if(length(intercept)!=1L || is.na(intercept) || !is.logical(intercept))
    stop("intercept must be TRUE or FALSE")
  p_min <- .varma_integer(p_min,"p_min",0L)
  if(!is.numeric(lag_power) || length(lag_power)!=1L || !is.finite(lag_power) || lag_power<=0)
    stop("lag_power must be positive and finite")
  n <- nrow(Y); m <- ncol(Y); d <- as.integer(intercept)
  # N-(m*P+d)>m guarantees room for a full-rank innovation covariance.
  cap <- floor((n-d-m-1L)/(m+1L))
  if(cap<p_min) stop("Too few observations for p_min and a nonsingular residual covariance")
  automatic <- is.null(p_max)
  proposed <- if(automatic) ceiling(log(n)^lag_power) else .varma_integer(p_max,"p_max",0L)
  P <- as.integer(if(automatic) min(proposed,cap) else proposed)
  if(P<p_min || P>cap) stop("p_max must be >= p_min and leave more than m residual degrees of freedom")
  rows <- seq.int(P+1L,n); N <- length(rows)
  lagged <- .varma_lags(Y,P)
  candidates <- seq.int(p_min,P)
  selection <- data.frame(p=candidates,nobs=N,ncoef=m*m*candidates+m*d,
    rank=NA_integer_,logdet=NA_real_,bic=Inf,admissible=FALSE,error="")
  for(j in seq_along(candidates)) {
    p <- candidates[j]
    X <- lagged[rows,seq_len(m*p),drop=FALSE]
    if(intercept) X <- cbind(1,X)
    if(ncol(X)) {
      reg <- .varma_regression(X,Y[rows,,drop=FALSE])
      selection$rank[j] <- reg$rank
      if(reg$rank<ncol(X)) { selection$error[j]<-"Rank-deficient regression"; next }
      E <- reg$residuals
    } else { E <- Y[rows,,drop=FALSE]; selection$rank[j] <- 0L }
    S <- crossprod(E)/N
    ch <- tryCatch(chol(S),error=function(e) NULL)
    if(is.null(ch) || any(!is.finite(ch)) || any(diag(ch)<=0)) {
      selection$error[j] <- "Singular or nonfinite residual covariance"; next
    }
    selection$logdet[j] <- 2*sum(log(diag(ch)))
    selection$bic[j] <- selection$logdet[j]+log(N)/N*selection$ncoef[j]
    selection$admissible[j] <- is.finite(selection$bic[j])
  }
  if(!any(selection$admissible)) stop("No admissible VAR candidate: check data rank and lag bounds")
  p <- selection$p[which.min(selection$bic)]
  if(p>0L) {
    # q=0 in HR uses ordinary least squares, with no long-VAR pilot or iteration.
    fit <- fit_varma_ihr(Y,p=p,q=0L,intercept=intercept)
  } else {
    mu <- if(intercept) colMeans(Y) else numeric(m)
    fit <- varma(mean=if(intercept) mu else NULL,cov=diag(m))
    fit <- .varma_finish(fit,Y,"fit_var_bic",m*d,start=0L,
                         rank=d,convergence=0L,iterations=0L)
  }
  if(!is.finite(fit$loglik)) stop("Selected VAR refit has an invalid covariance")
  fit$estimation_method <- "fit_var_bic"
  fit$selected_order <- p
  fit$lag_selection <- selection
  fit$lag_selection_nobs <- N; fit$lag_selection_p_max <- P
  fit$lag_selection_capped <- automatic && proposed>P
  fit$refit_args <- list(p_max=p_max,intercept=intercept,p_min=p_min,lag_power=lag_power)
  fit
}

# Shared summaries used by the package and the standalone R report collector.
.varma_bic_selection_summary <- function(d,design,n) {
  p <- d$selected_p; bound <- d$var_p_max
  ok <- is.finite(p); returned <- sum(ok); attempted <- nrow(d)
  hits <- sum(ok & p==bound,na.rm=TRUE)
  limits <- unique(bound[is.finite(bound)])
  summary <- data.frame(design=design,n=n,method="VAR_BIC",attempted=attempted,
    returned=returned,mean_p=if(returned) mean(p[ok]) else NA_real_,
    median_p=if(returned) stats::median(p[ok]) else NA_real_,
    p_max=if(length(limits)==1L) limits else NA_integer_,
    at_max=hits,at_max_fraction=if(returned) hits/returned else NA_real_,
    capped=sum(d$var_bound_capped %in% TRUE))
  labels <- ifelse(ok,as.character(p),"unavailable_after_failed_fit")
  counts <- table(labels)
  frequencies <- data.frame(design=design,n=n,method="VAR_BIC",selected_p=names(counts),
    count=as.integer(counts),attempted=attempted,fraction=as.integer(counts)/attempted)
  list(summary=summary,frequencies=frequencies)
}
