library(varma)
set.seed(332)
model <- varma(ar=array(diag(c(.65,.5)),c(2,2,1)),
               ma=array(diag(c(.2,.1)),c(2,2,1)),cov=diag(2))
Y <- sim_varma(model,n=1000)

# Indices refer to the existing order of the columns of Y.
known <- fit_varma_echelon(Y,kidx=c(1,1))
known$convergence
known$echelon$Phi0
predict(known,n.ahead=5)
irf(known,maxlag=10)

# Automatic selection runs internally, without MTS.
selected <- fit_varma_echelon(Y,plag=5,crit=.05)
selected$kronecker_indices
selected$selection$details

# Selection can also be run independently of coefficient estimation.
select_kronecker_indices(Y,plag=5,crit=.05)

# Faster stationary covariance and autocovariances, with no Kronecker system.
autocov_ss(known,maxlag=20)
ss <- varma_to_ss(known)
P <- solve_discrete_lyapunov(ss$T,ss$R %*% ss$Q %*% t(ss$R))
max(abs(P-ss$T %*% P %*% t(ss$T)-ss$R %*% ss$Q %*% t(ss$R)))

# Bootstrap refits retain the supplied/selected indices unless kidx=NULL is passed.
# bootstrap_varma(known,nsim=100)
# bootstrap_varma(known,nsim=100,kidx=NULL)
