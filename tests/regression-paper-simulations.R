library(varma)
internal <- function(name) get(name,asNamespace("varma"))
check <- function(name,expr) { force(expr); cat("PASS:",name,"\n") }
close <- function(a,b,tol=1e-9) stopifnot(isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=tol)))
fails <- function(expr) stopifnot(inherits(tryCatch({force(expr);NULL},error=identity),"error"))
root <- tempfile("varma-paper-tests-"); dir.create(root)

check("paper models reproduce the recorded designs and preserve RNG state",{
  set.seed(916);seed_before<-.Random.seed;kind_before<-RNGkind()
  models<-varma_paper_designs()
  stopifnot(identical(seed_before,.Random.seed),identical(kind_before,RNGkind()))
  close(models$model1$ar[1,2,1],.145706726703793,1e-13)
  close(models$model1$ma[1,2,1],.730055906175879,1e-13)
  close(models$model2$ar[,,1],matrix(c(-.630235480144620,.404748071916401,.146652669645846,-.663896159268916),2),1e-13)
  expected<-list(c(1L,0L),c(1L,1L),c(1L,1L),c(2L,2L,2L))
  for(i in seq_along(models)) stopifnot(identical(kronecker_indices(models[[i]]),expected[[i]]),all(is_stable(models[[i]])))
  stopifnot(identical(vapply(models,is_identified,logical(1)),c(model1=FALSE,model2=TRUE,model3=FALSE,model4=TRUE)))
})

check("true-model forecasts equal independent Gaussian conditioning",{
  f<-models$model2;n<-8L;h<-3L;m<-2L
  set.seed(145);Y<-sim_varma(f,n)
  G<-autocov_ss(f,n+h)
  V<-matrix(0,n*m,n*m);C<-matrix(0,h*m,n*m)
  for(i in seq_len(n)) for(j in seq_len(n))
    V[(i-1L)*m+seq_len(m),(j-1L)*m+seq_len(m)]<-if(i>=j) G[,,i-j+1L] else t(G[,,j-i+1L])
  for(i in seq_len(h)) for(j in seq_len(n))
    C[(i-1L)*m+seq_len(m),(j-1L)*m+seq_len(m)]<-G[,,n+i-j+1L]
  expected<-matrix(C%*%solve(V,as.numeric(t(Y))),h,m,byrow=TRUE)
  close(internal(".varma_paper_oracle")(f,Y,h),expected,1e-8)
})

check("loss summaries do not trim primary risk and use paired references",{
  x<-c(1,4,1e6,NA_real_);ref<-c(2,NA,5,8)
  s<-internal(".varma_paper_loss_summary")(x,ref)
  close(s$mse,mean(x,na.rm=TRUE));close(s$rmse,sqrt(mean(x,na.rm=TRUE)))
  close(s$mse_mcse,sd(x,na.rm=TRUE)/sqrt(3));stopifnot(s$n_available==3L,s$n_paired==2L)
  close(s$relative_mse,mean(c(1,1e6))/mean(c(2,5)))
})

args<-list(nsim=3L,sample_sizes=100L,designs=models[1],h=2L,irf_lags=2L,
  methods=c("HR1","ECH_oracle","ECH_selected","NET","MLE"),seed=772L,chunk_size=1L,
  estimator_args=list(MLE=list(control=list(maxit=25))),keep_fits=TRUE,verbose=FALSE)
dir1<-file.path(root,"sequential")
check("simulation records shared data, raw losses, convergence, and valid summaries",{
  set.seed(37);saved_seed<-.Random.seed
  result<-do.call(run_varma_paper_simulation,c(list(output_dir=dir1),args))
  stopifnot(identical(saved_seed,.Random.seed))
  paths<-list.files(dir1,pattern="^chunk-.*rds$",full.names=TRUE)
  records<-lapply(paths,readRDS)
  for(chunk in records) for(rep in chunk$results) {
    fits<-lapply(rep$results[-1],function(x) x$fit)
    stopifnot(all(vapply(fits,function(x) nrow(x$y)==100L,logical(1))))
    for(f in fits[-1]) close(f$y,fits[[1]]$y)
    stopifnot(all(rep$results$DGP$irf_error==0),all(rep$results$DGP$forecast_gap==0))
    close(rep$results$DGP$forecast_error,rep$actual-rep$oracle)
  }
  loss<-vapply(records,function(x) mean(x$results[[1]]$results$HR1$irf_error^2),numeric(1))
  row<-subset(result$summary$irf,method=="HR1" & subset=="available" & horizon==0)
  close(row$mse,mean(loss));close(row$mse_mcse,sd(loss)/sqrt(3))
  loss<-vapply(records,function(x) mean(x$results[[1]]$results$HR1$forecast_error[1,]^2),numeric(1))
  row<-subset(result$summary$forecast,method=="HR1" & subset=="available" & horizon==1 & variable==0)
  close(row$mse,mean(loss));close(row$mse_mcse,sd(loss)/sqrt(3))
})

check("completed and partial resumes preserve numerical results",{
  times<-file.info(paths)$mtime
  resumed<-do.call(run_varma_paper_simulation,c(list(output_dir=dir1),args))
  stopifnot(identical(times,file.info(paths)$mtime))
  close(resumed$summary$irf$mse,result$summary$irf$mse)
  dir2<-file.path(root,"resumed");dir.create(dir2)
  file.copy(file.path(dir1,"manifest.rds"),dir2)
  file.copy(paths[1],dir2)
  partial<-do.call(run_varma_paper_simulation,c(list(output_dir=dir2),args))
  close(partial$summary$irf$mse,result$summary$irf$mse)
  close(partial$summary$forecast$mse,result$summary$forecast$mse)
  changed<-args;changed$nsim<-4L
  fails(do.call(run_varma_paper_simulation,c(list(output_dir=dir1),changed)))
  fails(do.call(run_varma_paper_simulation,c(list(output_dir=dir1,resume=FALSE),args)))
})

check("failed estimators remain in attempted counts with explicit diagnostics",{
  bad<-run_varma_paper_simulation(file.path(root,"failures"),nsim=2L,sample_sizes=100L,
    designs=models[1],methods="RBY",estimator_args=list(RBY=list(r=1000L)),verbose=FALSE)
  row<-subset(bad$summary$failures,method=="RBY")
  stopifnot(row$attempted==2L,row$fit_failures==2L,row$forecast_failures==2L)
  stopifnot(all(is.na(subset(bad$summary$irf,method=="RBY")$mse)))
  stopifnot(all(nzchar(subset(bad$summary$diagnostics,method=="RBY")$error)))
})

if(identical(Sys.getenv("VARMA_TEST_PARALLEL"),"1")) check("parallel workers preserve the same samples and results",{
  old_plan<-future::plan();set.seed(319);saved_seed<-.Random.seed
  parallel_result<-do.call(run_varma_paper_simulation,c(list(output_dir=file.path(root,"parallel"),workers=2L),args))
  stopifnot(identical(saved_seed,.Random.seed),identical(class(future::plan()),class(old_plan)))
  close(parallel_result$summary$irf$mse,result$summary$irf$mse,1e-8)
  close(parallel_result$summary$forecast$mse,result$summary$forecast$mse,1e-8)
})
cat("All paper-simulation checks passed.\n")
