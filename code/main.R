##==============
get_result = function(X, Y, m, alpha.in.list, alpha.out.list, seed=2023, k0=10){
  cl <- makeCluster(detectCores() - 3)
  p = ncol(X); q = ncol(Y); k0 = k0
  #--- select from X:
  time0 <- proc.time()
  scalar.Y = scalar(Y, m=m)
  
  result.X <- mclapply(1:m, function(i) {
    suppressWarnings(step.multicop.x(i, X, scalar.Y, alpha.in.list, alpha.out.list, min(k0,p)))
  }, mc.cores = detectCores() - 3)
  
  X_sub = select.idx(result.X)
  cat("The selected subset in X: ", X_sub, '\n')
  
  duration <- proc.time() - time0
  elapsed_minutes_X <- duration['elapsed'] / 60
  print(paste("Time used for X in minutes:", elapsed_minutes_X)) # 1.68 min
  
  #--- select from Y:
  time0 <- proc.time()
  
  if(length(X_sub)==1){
    X.selected = as.matrix(X[,as.numeric(X_sub)], ncol=1)
  }else{
    X.selected = X[,as.numeric(X_sub)]
  }
  scalar.X = scalar(X.selected, m)
  result.Y <- mclapply(1:m, function(i) {
    suppressWarnings(step.multicop.y(i, Y, scalar.X, alpha.in.list, alpha.out.list, min(k0,q)))
  }, mc.cores = detectCores() - 3)
  
  stopCluster(cl)
  
  Y_sub = select.idx(result.Y)
  cat("The selected subset in Y: ", Y_sub, '\n')
  
  duration <- proc.time() - time0
  elapsed_minutes_Y <- duration['elapsed'] / 60
  print(paste("Time used for Y in minutes:", elapsed_minutes_Y)) 
  
  result = list(result_X=result.X, result_Y=result.Y, 
                time_X=elapsed_minutes_X, time_Y=elapsed_minutes_Y, 
                X_sub=X_sub, Y_sub=Y_sub)
  return(result)
}




get_rate = function(true, pred){
  conf_matrix <- table(true, pred)
  if (nrow(conf_matrix)==1){
    fpr = 0
    fnr = conf_matrix[1, 1] / (conf_matrix[1, 1] + conf_matrix[1, 2])
  }else{
    fpr <- conf_matrix[1, 2] / (conf_matrix[1, 1] + conf_matrix[1, 2])
    fnr <- conf_matrix[2, 1] / (conf_matrix[2, 1] + conf_matrix[2, 2])
  }
  return(c(fpr, fnr))
}
