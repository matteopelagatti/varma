.varma_paper_rng <- function() list(kind=RNGkind(),
  seed=if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) get(".Random.seed",.GlobalEnv) else NULL)
.varma_paper_restore_rng <- function(state) {
  do.call(RNGkind,as.list(state$kind))
  if(is.null(state$seed)) {
    if(exists(".Random.seed",.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)
  } else assign(".Random.seed",state$seed,.GlobalEnv)
}

#' Reproduce the four population models used in the paper
#'
#' Recreates the full-precision, seeded designs in inst/scripts/sim2.R without
#' changing the caller's random-number state. Means are zero. Model 1 has
#' nilpotent coefficients, Model 3 a common left factor, and Models 2 and 4
#' are identified at their supplied orders.
#' @return A named list of four varma objects.
#' @export
varma_paper_designs <- function() {
  saved <- .varma_paper_rng(); on.exit(.varma_paper_restore_rng(saved))
  RNGkind("Mersenne-Twister","Inversion","Rejection")
  set.seed(1); a <- rvarma(2,1,1)
  a$ar[c(1,2,4)] <- a$ma[c(1,2,4)] <- 0
  a$cov[1,2] <- a$cov[2,1] <- .3
  set.seed(2); b <- rvarma(2,1,1)
  b$cov[1,2] <- b$cov[2,1] <- .3
  set.seed(3); left <- rvarma(2,1,0,dist=stats::runif)
  reduced <- rvarma(2,1,1,dist=stats::runif); c <- rvarma(2,2,2)
  c$ar[,,1] <- left$ar[,,1]+reduced$ar[,,1]
  c$ar[,,2] <- -left$ar[,,1]%*%reduced$ar[,,1]
  c$ma[,,1] <- -(left$ar[,,1]+reduced$ma[,,1])
  c$ma[,,2] <- left$ar[,,1]%*%reduced$ma[,,1]
  c$cov[1,2] <- c$cov[2,1] <- .3
  set.seed(4); d <- rvarma(3,2,2)
  d$cov[1,2] <- d$cov[2,1] <- .3
  d$cov[1,3] <- d$cov[3,1] <- .2
  d$cov[2,3] <- d$cov[3,2] <- -.1
  list(model1=a,model2=b,model3=c,model4=d)
}

.varma_paper_method_names <- c("ECH_oracle","ECH_selected","IHR","HR1","NET","RBY","DPF","MLE")

.varma_paper_controls <- function(methods,estimator_args,mle_backend) {
  defaults <- list(ECH_oracle=list(control=list(maxit=500L,reltol=1e-8)),
    ECH_selected=list(plag=5L,crit=.05,control=list(maxit=500L,reltol=1e-8)),
    IHR=list(),HR1=list(maxit=1L),NET=list(parallel=FALSE),RBY=list(),DPF=list(),
    MLE=list(control=list(maxit=500L,reltol=1e-8)))
  if(!is.list(estimator_args) || (length(estimator_args) &&
      (is.null(names(estimator_args)) || anyDuplicated(names(estimator_args)) ||
       any(!names(estimator_args)%in%.varma_paper_method_names))))
    stop("estimator_args must be a named list indexed by method")
  for(nm in names(estimator_args)) {
    x <- estimator_args[[nm]]
    if(!is.list(x) || (length(x) && (is.null(names(x)) || anyDuplicated(names(x)) ||
        any(!nzchar(names(x)))))) stop("Each estimator_args entry must be a named list")
    if(any(names(x)%in%c("Y","p","q","intercept","kidx","ret","parallel","n_cores")))
      stop("estimator_args cannot override data, orders, indices, intercept, return type, or parallelism")
    if(nm=="HR1" && "maxit"%in%names(x) && !identical(as.integer(x$maxit),1L))
      stop("HR1 must retain maxit=1")
    fn <- .varma_paper_fit_function(nm,mle_backend)
    if(any(!names(x)%in%names(formals(fn)))) stop("Unknown estimator argument for ",nm)
    defaults[[nm]] <- utils::modifyList(defaults[[nm]],x,keep.null=TRUE)
  }
  defaults[methods]
}

.varma_paper_fit_function <- function(method,backend) switch(method,
  ECH_oracle=fit_varma_echelon,ECH_selected=fit_varma_echelon,
  IHR=fit_varma_ihr,HR1=fit_varma_ihr,NET=fit_varma_net,RBY=fit_varma_rby,
  DPF=fit_varma_dpf,MLE=switch(backend,skf=fit_varma_skf,kfas=fit_varma_kfas,fkf=fit_varma_fkf))

# Exact Gaussian forecast conditional on the same observed training sample.
# It never receives held-out observations or the simulated innovations.
.varma_paper_oracle <- function(model,Y,h,ss=NULL,P=NULL) {
  m <- ncol(Y); n <- nrow(Y)
  if(is.null(ss)) ss <- varma_to_ss(model)
  d <- nrow(ss$T)
  if(is.null(P)) P <- solve_discrete_lyapunov(ss$T,ss$R%*%ss$Q%*%t(ss$R))
  aa <- matrix(0,d,n+1L); PP <- array(0,c(d,d,n+1L)); PP[,,1L] <- P
  errors <- matrix(0,m,n); FF <- array(0,c(m,m,n))
  mu <- .varma_process_mean(model)
  kalman(Matrix::Matrix(ss$T,sparse=TRUE),Matrix::Matrix(ss$R,sparse=TRUE),ss$Q,
         aa,PP,t(sweep(Y,2,mu,"-")),errors,FF)
  out <- matrix(0,h,m); state <- aa[,n+1L]
  for(j in seq_len(h)) { out[j,] <- as.numeric(ss$Z%*%state)+mu; state <- ss$T%*%state }
  out
}

.varma_paper_capture <- function(expr) {
  warnings <- character(); started <- proc.time()[["elapsed"]]
  value <- tryCatch(withCallingHandlers(force(expr),warning=function(w) {
    warnings <<- unique(c(warnings,conditionMessage(w))); invokeRestart("muffleWarning")
  }),error=identity)
  list(value=value,warnings=paste(warnings,collapse=" | "),seconds=proc.time()[["elapsed"]]-started)
}

.varma_paper_replication <- function(seed,replication,cell,config) {
  assign(".Random.seed",seed,.GlobalEnv)
  model <- config$designs[[cell$design]]; n <- cell$n; h <- config$h
  orders <- as.integer(dim(model)); m <- orders[1L]; p <- orders[2L]; q <- orders[3L]
  t0 <- proc.time()[["elapsed"]]
  sample <- sim_varma(model,n+h,burn_in=config$burn_in)
  Y <- sample[seq_len(n),,drop=FALSE]; actual <- sample[n+seq_len(h),,drop=FALSE]
  simulation_seconds <- proc.time()[["elapsed"]]-t0
  oracle <- .varma_paper_capture(.varma_paper_oracle(model,Y,h,cell$ss,cell$P))
  if(inherits(oracle$value,"error") || any(!is.finite(oracle$value)))
    stop("True-model forecast failed in ",cell$design,", replication ",replication)
  truth <- cell$irf
  results <- list()
  for(method in c("DGP",config$methods)) {
    if(method=="DGP") {
      fitted <- list(value=model,warnings="",seconds=oracle$seconds)
    } else {
      args <- list(Y=Y,intercept=config$intercept)
      if(method%in%c("ECH_oracle","ECH_selected")) {
        args$kidx <- if(method=="ECH_oracle") cell$indices else NULL
      } else {
        args$p <- if(method=="DPF") p+(m-1L)*q else p
        args$q <- if(method=="DPF") m*q else q
      }
      args <- c(args,config$estimator_args[[method]])
      fitted <- .varma_paper_capture(do.call(.varma_paper_fit_function(method,config$mle_backend),args))
    }
    fit <- fitted$value; good <- !inherits(fit,"error") && inherits(fit,"varma")
    eirf <- array(NA_real_,c(m,m,config$irf_lags))
    efc <- gap <- matrix(NA_real_,h,m)
    status <- list(method=method,fit_ok=good,irf_ok=FALSE,forecast_ok=FALSE,
      fit_seconds=fitted$seconds,evaluation_seconds=0,convergence=NA_integer_,
      convergence_type=if(method%in%c("DGP","HR1","DPF")) "not_applicable" else
        if(method=="NET") "IHR_pilot" else "estimator",
      ar_stable=NA,ma_invertible=NA,qualified=FALSE,indices="",index_correct=NA,
      warnings=fitted$warnings,error=if(inherits(fit,"error")) conditionMessage(fit) else
        if(!good) "Estimator did not return a varma object" else "")
    if(good) {
      t1 <- proc.time()[["elapsed"]]
      stable <- tryCatch(is_stable(fit),error=function(e) c(NA,NA))
      status$ar_stable <- as.logical(stable[1L]); status$ma_invertible <- as.logical(stable[2L])
      if(length(fit$convergence)) status$convergence <- as.integer(fit$convergence)[1L]
      status$qualified <- isTRUE(status$ar_stable) && isTRUE(status$ma_invertible) &&
        (status$convergence_type=="not_applicable" || isTRUE(status$convergence==0L))
      if(length(fit$kronecker_indices)) {
        status$indices <- paste(fit$kronecker_indices,collapse=",")
        status$index_correct <- identical(as.integer(fit$kronecker_indices),as.integer(cell$indices))
      }
      response <- .varma_paper_capture(if(method=="DGP") truth else irf(fit,config$irf_lags,orth="none"))
      if(!inherits(response$value,"error")) {
        eirf <- array((response$value-truth)[,,-1L,drop=FALSE],c(m,m,config$irf_lags))
        status$irf_ok <- all(is.finite(eirf^2))
      }
      forecast <- .varma_paper_capture(if(method=="DGP") oracle$value else stats::predict(fit,n.ahead=h,cov=FALSE))
      if(!inherits(forecast$value,"error") && identical(dim(forecast$value),c(h,m))) {
        efc <- actual-forecast$value; gap <- forecast$value-oracle$value
        status$forecast_ok <- all(is.finite(efc^2)) && all(is.finite(gap^2))
      }
      if(!status$irf_ok) status$error <- paste(status$error,"IRF evaluation failed or produced nonfinite loss")
      if(!status$forecast_ok) status$error <- paste(status$error,"Forecast evaluation failed or produced nonfinite loss")
      status$warnings <- paste(Filter(nzchar,c(status$warnings,response$warnings,forecast$warnings)),collapse=" | ")
      status$evaluation_seconds <- proc.time()[["elapsed"]]-t1
    }
    results[[method]] <- list(status=status,irf_error=eirf,forecast_error=efc,forecast_gap=gap,
      fit=if(config$keep_fits && good && method!="DGP") fit else NULL)
  }
  list(replication=replication,seed=seed,simulation_seconds=simulation_seconds,
       elapsed_seconds=proc.time()[["elapsed"]]-t0,actual=actual,oracle=oracle$value,results=results)
}

.varma_paper_chunk <- function(job,config,config_id) {
  result <- lapply(seq_along(job$replications),function(i)
    .varma_paper_replication(job$seeds[[i]],job$replications[i],job$cell,config))
  payload <- list(schema=1L,config_id=config_id,cell=job$cell[c("id","design","n","indices")],
                  replications=job$replications,results=result)
  temporary <- paste0(job$file,".partial-",Sys.getpid())
  saveRDS(payload,temporary,compress=TRUE)
  if(file.exists(job$file) || !file.rename(temporary,job$file))
    stop("Could not publish checkpoint: ",job$file)
  invisible(job$file)
}

.varma_paper_hash <- function(x) {
  path <- tempfile("varma-config-"); on.exit(unlink(path))
  saveRDS(x,path,compress=FALSE,version=2)
  unname(tools::md5sum(path))
}

#' Run the paper's simulation experiment with resumable checkpoints
#'
#' Repeats the four original Gaussian designs, retaining the original estimation
#' competitors and adding echelon estimation with data-selected indices. All
#' methods see the same training sample in each replication. The DGP benchmark
#' forecasts from known parameters conditional on that same observed history.
#' @param output_dir directory for a manifest, checkpoints, and CSV summaries.
#' @param nsim number of replications per design and sample size (default 1000).
#' @param sample_sizes training sample sizes (default c(100,200,400)).
#' @param designs named list of stable, invertible varma population models.
#' @param h maximum forecast horizon (default 5).
#' @param irf_lags number of positive response lags (default 5); lag zero is excluded.
#' @param methods estimation methods; defaults to all eight paper competitors.
#' DGP forecasts are always included separately.
#' @param seed integer seed. Independent L'Ecuyer streams are assigned to cells
#' and substreams to replications, preserving results across workers and resumes.
#' @param workers number of parallel R processes; 1 runs sequentially.
#' @param chunk_size replications per checkpoint (default 10).
#' @param burn_in simulated observations discarded before training (default 100,
#' as in the original script).
#' @param intercept logical; estimate intercepts (default FALSE, as in the paper).
#' @param mle_backend stationary Gaussian likelihood backend; default "skf" is
#' the compiled Kalman objective. Alternatives are "fkf" and "kfas".
#' @param estimator_args named list of method-specific argument lists. MLE and
#' both echelon estimators default to control=list(maxit=500,reltol=1e-8).
#' Other methods retain their package defaults; HR1 uses one iteration.
#' @param resume reuse compatible completed checkpoints (default TRUE).
#' @param keep_fits retain fitted objects (default FALSE). Raw IRF errors,
#' forecast errors, forecasts relative to the DGP, and diagnostics are always saved.
#' @param verbose report progress after each batch of checkpoints.
#' @details The default experiment has 12 cells and 96,000 estimated models,
#' plus 12,000 true-parameter forecast calculations. DPF uses the original
#' expanded orders p+(m-1)q and mq. Oracle echelon receives population Kronecker
#' indices; feasible echelon selects them using training data only. Neither is
#' allowed access to held-out outcomes.
#'
#' Errors and warnings are recorded rather than silently discarding failures.
#' Primary accuracy summaries use all available finite results, without trimming
#' large errors. Additional summaries require stability, invertibility and
#' reported convergence (where applicable). HR1 is deliberately one-step, so
#' its IHR stopping code is not treated as an optimization failure. NET's
#' convergence flag refers to its IHR pilot. No automatic retries or fallback
#' estimators change the named procedures.
#'
#' A checkpoint is published only after its chunk finishes. Resume validates
#' the simulation configuration; workers may change, but other settings must
#' match. Run only one driver per output directory. At most an unfinished
#' chunk per worker needs repeating after interruption. The caller's RNG state
#' and future plan are restored. Use an installed version of varma on every
#' worker. Avoid nested or multithreaded BLAS parallelism when using many workers.
#' @return A list with output_dir, config, execution_seconds, and summary tables.
#' @seealso \code{\link{summarize_varma_paper_simulation}},
#' \code{\link{benchmark_varma_paper_simulation}}
#' @export
#' @examples
#' \dontrun{
#' run_varma_paper_simulation("paper-results", workers=4)
#' # Resume after interruption, possibly with a different number of workers:
#' run_varma_paper_simulation("paper-results", workers=2)
#' }
run_varma_paper_simulation <- function(output_dir,nsim=1000L,sample_sizes=c(100L,200L,400L),
  designs=varma_paper_designs(),h=5L,irf_lags=5L,methods=.varma_paper_method_names,
  seed=20261110L,workers=1L,chunk_size=10L,burn_in=100L,intercept=FALSE,
  mle_backend=c("skf","fkf","kfas"),estimator_args=list(),resume=TRUE,keep_fits=FALSE,verbose=TRUE) {
  saved <- .varma_paper_rng(); on.exit(.varma_paper_restore_rng(saved),add=TRUE)
  nsim <- .varma_integer(nsim,"nsim"); h <- .varma_integer(h,"h")
  irf_lags <- .varma_integer(irf_lags,"irf_lags"); workers <- .varma_integer(workers,"workers")
  chunk_size <- .varma_integer(chunk_size,"chunk_size"); burn_in <- .varma_integer(burn_in,"burn_in",0L)
  seed <- .varma_integer(seed,"seed",0L)
  if(!is.numeric(sample_sizes) || !length(sample_sizes) || any(!is.finite(sample_sizes)) ||
      any(sample_sizes<2) || any(sample_sizes!=as.integer(sample_sizes)) || anyDuplicated(sample_sizes))
    stop("sample_sizes must contain distinct integers >= 2")
  sample_sizes <- as.integer(sample_sizes)
  if(!is.character(methods) || !length(methods) || anyNA(methods) || anyDuplicated(methods) ||
      any(!methods%in%.varma_paper_method_names)) stop("Unknown or repeated estimation methods")
  for(nm in c("intercept","resume","keep_fits","verbose")) {
    value <- get(nm)
    if(!is.logical(value) || length(value)!=1L || is.na(value)) stop(nm," must be TRUE or FALSE")
  }
  if(!is.character(output_dir) || length(output_dir)!=1L || is.na(output_dir) || !nzchar(output_dir))
    stop("output_dir must be a directory path")
  if(!is.list(designs) || !length(designs) || is.null(names(designs)) ||
      anyNA(names(designs)) || any(!nzchar(names(designs))) || anyDuplicated(names(designs)))
    stop("designs must be a named list with distinct nonempty names")
  for(model in designs) {
    if(!inherits(model,"varma") || any(!is_stable(model))) stop("Every design must be a stable and invertible varma model")
    if(sum(dim(model)[-1L])==0L) stop("Supply dynamic designs; white-noise fitted echelon models remain supported")
    if(is.null(model$cov) || any(!is.finite(model$cov)) ||
       !isTRUE(all.equal(model$cov,t(model$cov)))) stop("Every design needs a symmetric innovation covariance")
    chol(model$cov)
  }
  mle_backend <- match.arg(mle_backend)
  controls <- .varma_paper_controls(methods,estimator_args,mle_backend)
  config <- list(schema=1L,designs=designs,nsim=nsim,sample_sizes=sample_sizes,h=h,
    irf_lags=irf_lags,methods=methods,seed=seed,chunk_size=chunk_size,burn_in=burn_in,
    intercept=intercept,mle_backend=mle_backend,estimator_args=controls,keep_fits=keep_fits,
    package_version=as.character(utils::packageVersion("varma")))
  config_id <- .varma_paper_hash(config)
  if(!dir.exists(output_dir) && !dir.create(output_dir,recursive=TRUE)) stop("Could not create output_dir")
  output_dir <- normalizePath(output_dir,winslash="/",mustWork=TRUE)
  manifest_file <- file.path(output_dir,"manifest.rds")
  if(file.exists(manifest_file)) {
    if(!resume) stop("An experiment already exists here; use resume=TRUE or a new directory")
    old <- readRDS(manifest_file)
    if(!identical(old$config_id,config_id)) stop("Checkpoint configuration differs; use the original settings or a new directory")
  } else {
    if(length(list.files(output_dir,pattern="^chunk-.*[.]rds$"))) stop("Checkpoints exist without a manifest")
    temporary <- paste0(manifest_file,".partial-",Sys.getpid())
    saveRDS(list(config=config,config_id=config_id,created=Sys.time(),session=utils::sessionInfo(),
                 system=Sys.info(),logical_cores=parallel::detectCores()),temporary)
    if(!file.rename(temporary,manifest_file)) stop("Could not publish the manifest")
  }
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed); cell_seed <- get(".Random.seed",.GlobalEnv)
  jobs <- list(); cell_id <- 0L
  for(n in sample_sizes) for(design in names(designs)) {
    cell_id <- cell_id+1L; model <- designs[[design]]
    ss <- varma_to_ss(model)
    cell <- list(id=cell_id,design=design,n=n,indices=kronecker_indices(model),ss=ss,
      P=solve_discrete_lyapunov(ss$T,ss$R%*%ss$Q%*%t(ss$R)),irf=irf(model,irf_lags,orth="none"))
    seeds <- vector("list",nsim); next_seed <- cell_seed
    for(b in seq_len(nsim)) { seeds[[b]] <- next_seed; next_seed <- parallel::nextRNGSubStream(next_seed) }
    cell_seed <- parallel::nextRNGStream(cell_seed)
    for(first in seq.int(1L,nsim,by=chunk_size)) {
      ids <- seq.int(first,min(nsim,first+chunk_size-1L))
      file <- file.path(output_dir,sprintf("chunk-%03d-%06d-%06d.rds",cell_id,min(ids),max(ids)))
      if(file.exists(file)) {
        existing <- tryCatch(readRDS(file),error=function(e) stop("Unreadable checkpoint: ",file))
        if(!identical(existing$config_id,config_id) || !identical(existing$replications,ids) ||
            !identical(existing$cell$design,design) || existing$cell$n!=n || length(existing$results)!=length(ids))
          stop("Invalid checkpoint: ",file)
      } else jobs[[length(jobs)+1L]] <- list(cell=cell,replications=ids,seeds=seeds[ids],file=file)
    }
  }
  started <- proc.time()[["elapsed"]]; total <- length(jobs)
  if(workers>1L && total) {
    old_plan <- future::plan(); on.exit(future::plan(old_plan),add=TRUE,after=FALSE)
    future::plan(future::multisession,workers=workers)
  }
  if(verbose) message(total," checkpoint chunks pending; ",nsim*length(designs)*length(sample_sizes),
                     " replications in the complete experiment")
  if(total) for(first in seq.int(1L,total,by=max(1L,workers*2L))) {
    indices <- seq.int(first,min(total,first+workers*2L-1L))
    batch <- jobs[indices]
    if(workers>1L) future.apply::future_lapply(batch,.varma_paper_chunk,config=config,config_id=config_id,
      future.seed=lapply(batch,function(x) x$seeds[[1L]]),future.packages="varma",future.scheduling=Inf)
    else lapply(batch,.varma_paper_chunk,config=config,config_id=config_id)
    if(verbose) message(sprintf("Completed %d/%d pending chunks; elapsed %.1f minutes",
                                max(indices),total,(proc.time()[["elapsed"]]-started)/60))
  }
  elapsed <- proc.time()[["elapsed"]]-started
  execution <- list(started=Sys.time()-elapsed,finished=Sys.time(),seconds=elapsed,
    workers=workers,new_chunks=total,session=utils::sessionInfo())
  saveRDS(execution,file.path(output_dir,paste0("execution-",format(Sys.time(),"%Y%m%d-%H%M%S"),"-",Sys.getpid(),".rds")))
  summary <- summarize_varma_paper_simulation(output_dir)
  invisible(list(output_dir=output_dir,config=config,execution_seconds=elapsed,summary=summary))
}

.varma_paper_loss_summary <- function(loss,reference=NULL,trim=.01) {
  good <- is.finite(loss); x <- loss[good]; n <- length(x)
  mse <- if(n) mean(x) else NA_real_; rmse <- sqrt(mse)
  se <- if(n>1L) stats::sd(x)/sqrt(n) else NA_real_
  ans <- list(n_available=n,mse=mse,rmse=rmse,mse_mcse=se,
    rmse_mcse=if(isTRUE(rmse==0)) 0 else se/(2*rmse),
    trimmed_rmse=if(n) sqrt(mean(x,trim=trim)) else NA_real_,
    n_paired=0L,paired_mse=NA_real_,reference_mse=NA_real_,relative_mse=NA_real_,relative_mse_mcse=NA_real_)
  if(!is.null(reference)) {
    paired <- good & is.finite(reference); np <- sum(paired)
    ans$n_paired <- np
    if(np) {
      a <- loss[paired]; b <- reference[paired]; mb <- mean(b)
      ans$paired_mse <- mean(a); ans$reference_mse <- mb
      if(mb>0) {
        ratio <- mean(a)/mb; ans$relative_mse <- ratio
        if(np>1) ans$relative_mse_mcse <- stats::sd(a-ratio*b)/(sqrt(np)*mb)
      }
    }
  }
  ans
}

#' Summarize complete or partially completed paper simulations
#'
#' Reads checkpoints one design/sample-size cell at a time. Primary summaries
#' use all available finite losses without trimming. A separate qualified subset
#' requires stable, invertible estimates and reported convergence when applicable.
#' @param output_dir directory produced by run_varma_paper_simulation.
#' @param trim symmetric trimming proportion for an additional, explicitly
#' labeled trimmed RMSE column (default 0.01); primary RMSE is always untrimmed.
#' @param write_csv write summary CSV files into output_dir (default TRUE).
#' @details IRF MSE averages squared reduced-form coefficient errors over entries
#' and positive lags. Horizon zero in the summary means pooling all positive
#' lags. Forecast MSE averages squared held-out forecast errors. Variable zero
#' means pooling variables within a replication, not treating them as independent
#' Monte Carlo draws. Forecast-gap MSE measures distance from the true-model
#' conditional forecast, isolating forecast estimation error.
#'
#' Monte Carlo standard errors use independent replication-level losses; RMSE
#' errors use the delta method. Relative MSE compares methods on common usable
#' replications: IRFs against oracle echelon, forecasts against DGP. The paired
#' reference denominator and sample count are reported explicitly. Quantities
#' unavailable because of failed fits remain missing and have separate failure
#' counts. These conditional summaries are not unconditional risks of a procedure
#' that sometimes fails. No observations are removed based on error magnitude.
#' @return Lists named irf, forecast, forecast_gap, failures, timing, indices,
#' and diagnostics, also written as CSV files when requested.
#' @export
summarize_varma_paper_simulation <- function(output_dir,trim=.01,write_csv=TRUE) {
  if(length(trim)!=1L || !is.finite(trim) || trim<0 || trim>=.5) stop("trim must be in [0,0.5)")
  manifest <- readRDS(file.path(output_dir,"manifest.rds")); config <- manifest$config
  methods <- c("DGP",config$methods); K <- length(methods)
  out <- stats::setNames(vector("list",7L),c("irf","forecast","forecast_gap","failures","timing","indices","diagnostics"))
  cell_id <- 0L
  for(n in config$sample_sizes) for(design in names(config$designs)) {
    cell_id <- cell_id+1L; m <- as.integer(dim(config$designs[[design]])[1L])
    files <- list.files(output_dir,pattern=sprintf("^chunk-%03d-[0-9]+-[0-9]+[.]rds$",cell_id),full.names=TRUE)
    if(!length(files)) next
    irf_loss <- array(NA_real_,c(config$nsim,K,config$irf_lags+1L))
    forecast_loss <- gap_loss <- array(NA_real_,c(config$nsim,K,config$h,m+1L))
    qualified <- matrix(FALSE,config$nsim,K); present <- rep(FALSE,config$nsim)
    diagnostics <- list()
    for(file in files) {
      chunk <- readRDS(file)
      if(!identical(chunk$config_id,manifest$config_id)) stop("Incompatible checkpoint: ",file)
      for(rep in chunk$results) {
        b <- rep$replication
        if(b<1L || b>config$nsim || present[b]) stop("Duplicate or invalid replication in ",file)
        present[b] <- TRUE
        for(k in seq_along(methods)) {
          result <- rep$results[[methods[k]]]; status <- result$status
          qualified[b,k] <- isTRUE(status$qualified)
          if(isTRUE(status$irf_ok)) {
            squared <- result$irf_error^2
            irf_loss[b,k,1L] <- mean(squared)
            irf_loss[b,k,-1L] <- apply(squared,3L,mean)
          }
          if(isTRUE(status$forecast_ok)) {
            squared <- result$forecast_error^2; gaps <- result$forecast_gap^2
            forecast_loss[b,k,,1L] <- rowMeans(squared); forecast_loss[b,k,,-1L] <- squared
            gap_loss[b,k,,1L] <- rowMeans(gaps); gap_loss[b,k,,-1L] <- gaps
          }
          diagnostics[[length(diagnostics)+1L]] <- as.data.frame(c(list(design=design,n=n,
            replication=b),status),stringsAsFactors=FALSE)
        }
      }
    }
    diag <- do.call(rbind,diagnostics); attempted <- sum(present)
    out$diagnostics[[length(out$diagnostics)+1L]] <- diag
    for(k in seq_along(methods)) {
      method <- methods[k]; d <- diag[diag$method==method,,drop=FALSE]
      applicable <- d$convergence_type!="not_applicable" & d$fit_ok
      out$failures[[length(out$failures)+1L]] <- data.frame(design=design,n=n,method=method,
        planned=config$nsim,attempted=attempted,fit_failures=sum(!d$fit_ok),
        irf_failures=sum(!d$irf_ok),forecast_failures=sum(!d$forecast_ok),
        nonconverged=sum(applicable & !is.na(d$convergence) & d$convergence!=0L),
        convergence_unknown=sum(applicable & is.na(d$convergence)),
        unstable_ar=sum(d$fit_ok & !is.na(d$ar_stable) & !d$ar_stable),
        noninvertible_ma=sum(d$fit_ok & !is.na(d$ma_invertible) & !d$ma_invertible),
        qualified=sum(d$qualified),with_warnings=sum(nzchar(d$warnings)))
      out$timing[[length(out$timing)+1L]] <- data.frame(design=design,n=n,method=method,
        attempted=attempted,mean_fit_seconds=mean(d$fit_seconds),median_fit_seconds=stats::median(d$fit_seconds),
        p90_fit_seconds=as.numeric(stats::quantile(d$fit_seconds,.9)),total_fit_seconds=sum(d$fit_seconds),
        mean_evaluation_seconds=mean(d$evaluation_seconds))
      if(method%in%c("ECH_oracle","ECH_selected")) {
        labels <- ifelse(nzchar(d$indices),d$indices,"unavailable_after_failed_fit")
        counts <- table(labels)
        out$indices[[length(out$indices)+1L]] <- data.frame(design=design,n=n,method=method,
          indices=names(counts),count=as.integer(counts),attempted=attempted,
          fraction=as.integer(counts)/attempted)
      }
      for(subset in c("available","qualified")) {
        mask <- if(subset=="qualified") qualified[,k] else present
        for(lag in 0:config$irf_lags) {
          loss <- irf_loss[,k,lag+1L]; loss[!mask] <- NA_real_
          ref <- NULL; r <- match("ECH_oracle",methods)
          if(!is.na(r)) { ref <- irf_loss[,r,lag+1L]; if(subset=="qualified") ref[!qualified[,r]] <- NA_real_ }
          out$irf[[length(out$irf)+1L]] <- as.data.frame(c(list(design=design,n=n,method=method,
            subset=subset,horizon=lag,attempted=attempted,reference="ECH_oracle"),
            .varma_paper_loss_summary(loss,ref,trim)),stringsAsFactors=FALSE)
        }
        for(horizon in seq_len(config$h)) for(variable in 0:m) for(target in c("forecast","forecast_gap")) {
          loss <- if(target=="forecast") forecast_loss[,k,horizon,variable+1L] else gap_loss[,k,horizon,variable+1L]
          loss[!mask] <- NA_real_
          ref <- if(target=="forecast") forecast_loss[,1L,horizon,variable+1L] else NULL
          out[[target]][[length(out[[target]])+1L]] <- as.data.frame(c(list(design=design,n=n,
            method=method,subset=subset,horizon=horizon,variable=variable,attempted=attempted,
            reference=if(target=="forecast") "DGP" else NA_character_),
            .varma_paper_loss_summary(loss,ref,trim)),stringsAsFactors=FALSE)
        }
      }
    }
  }
  out <- lapply(out,function(x) if(length(x)) do.call(rbind,x) else data.frame())
  if(write_csv) for(nm in names(out)) utils::write.csv(out[[nm]],file.path(output_dir,paste0(nm,"-summary.csv")),row.names=FALSE)
  out
}

#' Benchmark and project the duration of the paper simulation
#'
#' Runs a small experiment with the same estimation controls as the full run.
#' @param output_dir a new directory for pilot checkpoints and timing projections.
#' @param pilot_nsim replications per cell for timing (default 5; use at least 20
#' for a more stable estimate when practical).
#' @param full_nsim replications per cell in the intended full run (default 1000).
#' @param workers actual number of worker processes to benchmark (default 1).
#' @param ... further arguments to run_varma_paper_simulation, excluding nsim.
#' @details Sequential projections sum mean elapsed replication times by cell,
#' including failures and all procedures. Their Monte Carlo standard error
#' reflects pilot variability, not changes in hardware, thermal throttling, or
#' contention. Parallel projections extrapolate the measured worker throughput;
#' they do not assume perfect scaling by the logical core count. Small pilots
#' can miss rare slow optimizations. Full estimator iteration budgets are used.
#' @return A list with cell_times, task_seconds, serial_seconds,
#' serial_mcse_seconds, throughput_seconds, and pilot_wall_seconds. Serial
#' projections are NA for parallel pilots because per-task elapsed times under
#' CPU contention do not estimate single-worker performance. CSV/RDS projections
#' are saved.
#' @export
benchmark_varma_paper_simulation <- function(output_dir,pilot_nsim=5L,full_nsim=1000L,workers=1L,...) {
  full_nsim <- .varma_integer(full_nsim,"full_nsim")
  pilot <- run_varma_paper_simulation(output_dir,nsim=pilot_nsim,workers=workers,resume=FALSE,...)
  files <- list.files(output_dir,pattern="^chunk-.*[.]rds$",full.names=TRUE)
  times <- do.call(rbind,lapply(files,function(file) {
    x <- readRDS(file)
    data.frame(design=x$cell$design,n=x$cell$n,
      seconds=vapply(x$results,function(r) r$elapsed_seconds,numeric(1)))
  }))
  cells <- split(times,interaction(times$design,times$n,drop=TRUE))
  cell_times <- do.call(rbind,lapply(cells,function(x) data.frame(design=x$design[1],n=x$n[1],
    pilot_nsim=nrow(x),mean_seconds=mean(x$seconds),median_seconds=stats::median(x$seconds),
    max_seconds=max(x$seconds),sd_seconds=stats::sd(x$seconds),
    projected_seconds=mean(x$seconds)*full_nsim)))
  rownames(cell_times) <- NULL
  out <- list(cell_times=cell_times,task_seconds=sum(cell_times$projected_seconds),
    serial_seconds=if(workers==1L) sum(cell_times$projected_seconds) else NA_real_,
    serial_mcse_seconds=if(workers==1L) full_nsim*sqrt(sum(cell_times$sd_seconds^2/cell_times$pilot_nsim)) else NA_real_,
    throughput_seconds=pilot$execution_seconds*full_nsim/pilot_nsim,
    pilot_wall_seconds=pilot$execution_seconds,workers=workers,full_nsim=full_nsim)
  utils::write.csv(cell_times,file.path(output_dir,"runtime-projection.csv"),row.names=FALSE)
  saveRDS(out,file.path(output_dir,"runtime-projection.rds"))
  out
}
