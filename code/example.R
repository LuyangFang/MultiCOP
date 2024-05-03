rm(list=ls())
setwd("~/Desktop/microbiome/code/github")

source("utils.R")
source("main.R")
library(parallel)
# library(ggplot2)


##--
# function to generate synthetic data
sce2 <- function(n, p, q, rho, sigma_epsilon){
  # Generate X
  Sigma_X <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma_X[i, j] <- rho^(abs(i - j))
    }
  }
  X <- mvtnorm::rmvnorm(n = n, mean = rep(0, p), sigma = Sigma_X)
  
  epsilon = sigma_epsilon * rnorm(n*q, 0, 1)
  epsilon = matrix(epsilon, nrow=n, ncol=q)
  
  # Calculate Y based on the linear model
  bias = ( X[,1] + X[,2] ) / ( 0.5 + (1.5+X[,1])^2 )
  Y <- epsilon
  for (i in 1:2){
    Y[,i] = Y[,i] + bias
  }
  
  return(list(X = X, Y = Y))
}




##--- simulation scenario 2, linear setting, rho=0.5, sigma=3, n=100 ----
n = 100
p = 5
q = 6
rho = 0.5
sigma_epsilon = 3

my.range = 100;
m = 100
alpha.in.list = c(0.90, 0.95, 0.99)
alpha.out.list = c(0.85, 0.90, 0.95) 
rep = 30

##-- generate data:
# dat_list_s2 = list()
# for (rr in 1:rep){
#   sce2_dat = sce2(n, p, q, rho, sigma_epsilon)
#   dat_list_s2[[rr]] = sce2_dat
# }
# save(dat_list_s2, file=paste0("example_data.rdata"))
load(paste0("example_data.rdata"))



##--- MultiCOP:
rep = 30
results_list = list()
for (rr in 1:rep){
  print(paste("-------repetition",rr))
  # generate data:
  dat = dat_list_s2[[rr]]
  X = dat$X; Y = dat$Y
  # get results:
  result = get_result(X, Y, m, alpha.in.list, alpha.out.list, seed=rr, k0=5)
  results_list[[rr]] = result
}

save(results_list, file="example_results.rdata")
load(paste0("example_results.rdata"))



## get score:
true_X = rep(0,p); true_X[1:2] = 1
true_Y = rep(0,q); true_Y[1:2] = 1

FR_X = matrix(NA, nrow=rep, ncol=2); colnames(FR_X) = c("FPR_X", "FNR_X"); rownames(FR_X) = paste("rep", 1:rep)
FR_Y = matrix(NA, nrow=rep, ncol=2); colnames(FR_Y) = c("FPR_Y", "FNR_Y"); rownames(FR_Y) = paste("rep", 1:rep)
for (rr in 1:rep){
  X_sub = results_list[[rr]]$X_sub
  Y_sub = results_list[[rr]]$Y_sub
  pred_X = rep(0,p); pred_X[as.integer(X_sub)] = 1
  pred_Y = rep(0,q); pred_Y[as.integer(Y_sub)] = 1
  FR_X[rr,] = get_rate(true_X, pred_X)
  FR_Y[rr,] = get_rate(true_Y, pred_Y)
}

apply(FR_X, 2, mean) # 0     0 
apply(FR_X, 2, sd) # 0     0 
apply(FR_Y, 2, mean) # 0.05  0.10 
apply(FR_Y, 2, sd) # 0.1017095 0.2034191 






