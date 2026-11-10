# Developer-only comparison with MTS 1.2.1; not a package dependency or test.
# Run from the repository root after installing MTS in a development library.
Rcpp::sourceCpp("src/rcpp_functions.cpp",showOutput=FALSE,
  cacheDir=file.path(Sys.getenv("TEMP",tempdir()),"varma-regression-source-cache"))
for(file in c("estimation_helpers.R","functions.R","echelon.R")) source(file.path("R",file))
stopifnot(requireNamespace("MTS",quietly=TRUE))
simulate <- function(seed,m,n,ar,ma=rep(.15,m)) {
  set.seed(seed)
  f<-varma(ar=array(diag(ar,m),c(m,m,1)),ma=array(diag(ma,m),c(m,m,1)),cov=diag(m))
  sim_varma(f,n)
}
set.seed(17);white2<-matrix(rnorm(1000),500,2)
set.seed(19);white3<-matrix(rnorm(1200),400,3)
cases<-list(white2=white2,white3=white3,
  dynamic2=simulate(332,2,1000,c(.65,.5),c(.2,.1)),
  dynamic3=simulate(672,3,700,c(.7,.6,.5)),
  dynamic5=simulate(476,5,1000,seq(.5,.7,length.out=5)),
  zero_first=simulate(845,2,800,c(0,.8),c(0,.2)))
comparison<-list();fixtures<-list()
for(nm in names(cases)) {
  y<-cases[[nm]]
  cutoff<-if(grepl("^white",nm)) .01 else .05
  invisible(capture.output(reference<-MTS::Kronid(scale(y),crit=cutoff)))
  current<-select_kronecker_indices(y,crit=cutoff)
  difference<-max(abs(current$tests-reference$tests))
  cat(nm, "indices", current$index, "maximum test difference", difference, "\n")
  stopifnot(identical(as.integer(reference$index),as.integer(current$index)),
    isTRUE(all.equal(as.numeric(current$tests),as.numeric(reference$tests),tolerance=1e-10)))
  fixtures[[nm]]<-list(Y=y,plag=5L,crit=cutoff,index=reference$index,tests=reference$tests)
  comparison[[nm]]<-data.frame(case=nm,n=nrow(y),m=ncol(y),
    indices=paste(current$index,collapse=","),max_test_difference=difference)
}
patterns<-lapply(list(c(0,0),c(2,0),c(0,2),c(2,1,0),c(0,1,2),2L),function(k)
  list(indices=k,spec=MTS::Kronspec(k,output=FALSE)))
dir.create("tests/fixtures",showWarnings=FALSE)
saveRDS(list(source="MTS 1.2.1, Kronid and Kronspec; numerical outputs only",
  url="https://cran.r-project.org/package=MTS",cases=fixtures,patterns=patterns),
  "tests/fixtures/kronid-reference.rds",compress="xz")
write.csv(do.call(rbind,comparison),"audit/2026-11-10/kronid/agreement.csv",row.names=FALSE)
print(do.call(rbind,comparison))

cat("\nKnown MTS edge cases (intentional differences):\n")
zero_last<-cases$zero_first[,2:1]
ours<-select_kronecker_indices(zero_last)
invisible(capture.output(old<-MTS::Kronid(scale(zero_last))))
print(ours$details)
cat("MTS degrees of freedom:",old$tests[2,],"\n")
cat("Internal degrees of freedom:",ours$tests[2,],"\n")
uni<-cases$dynamic2[,1,drop=FALSE]
print(select_kronecker_indices(uni)$index)
old_uni<-tryCatch({invisible(capture.output(MTS::Kronid(scale(uni))));"OK"},error=conditionMessage)
cat("MTS univariate result:",old_uni,"\n")
capture.output(sessionInfo(),file="audit/2026-11-10/kronid/session.txt")
