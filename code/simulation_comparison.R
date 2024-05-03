library(tidyverse)
library(PMA) # CCA
library(plsVarSel) # PLS
library(pls)
library(glmnet) # multiLasso
library(rrpack) # RRR

#### Calculate FPR and FNR
get_rate <- function(true, pred)
{
  if(sum(pred) == length(pred))
  {
    fpr <- 1
    fnr <- (length(true) - sum(true))/length(true)
    return(c(fpr, fnr))
  }
  else if(sum(pred) == 0)
  {
    fpr <- sum(true)/length(true)
    fnr <- 1
    return(c(fpr, fnr))
  }
  else
  {
    conf_matrix <- table(true, pred)
    if(sum(true) != length(true))
    {
      fpr <- conf_matrix[1, 2] / (conf_matrix[1, 1] + conf_matrix[1, 2])
      fnr <- conf_matrix[2, 1] / (conf_matrix[2, 1] + conf_matrix[2, 2])
      return(c(fpr, fnr))  
    }
    else
    {
      fpr <- 0
      fnr <- sum(pred)/length(pred)
      return(c(fpr, fnr))
    }
  }
}


#### Main fuction
get_result <- function(data, result, true_x, true_y)
{
  result_fpr_x <- matrix(NA, nrow = length(result), ncol = 5)
  result_fpr_y <- matrix(NA, nrow = length(result), ncol = 5)
  result_fnr_x <- matrix(NA, nrow = length(result), ncol = 5)
  result_fnr_y <- matrix(NA, nrow = length(result), ncol = 5)
  colnames(result_fpr_x) = colnames(result_fpr_y) = colnames(result_fnr_x) = colnames(result_fnr_y) <- 
    c("MultiCOP", "CCA", "PLS", "MultiLasso", "RRR")
  
  for(i in 1:length(result))
  {
    #### MultiCOP
    result_rep <- result[[i]]
    index_x <- as.numeric(result_rep$X_sub)
    pred_x <- rep(0, length(true_x))
    pred_x[index_x] <- 1
    index_x_count <- sum(pred_x)
    rate <- get_rate(true_x, pred_x)
    result_fpr_x[i, 1] <- rate[1]
    result_fnr_x[i, 1] <- rate[2]

    index_y <- as.numeric(result_rep$Y_sub)
    pred_y <- rep(0, length(true_y))
    pred_y[index_y] <- 1
    index_y_count <- sum(pred_y)
    rate <- get_rate(true_y, pred_y)
    result_fpr_y[i, 1] <- rate[1]
    result_fnr_y[i, 1] <- rate[2]
    
    #### Load data
    data_rep <- data[[i]]
    X <- data_rep$X
    Y <- data_rep$Y
    
    #### CCA
    result_cca <- CCA(x = X, z = Y, typex = "standard", typez = "standard")
    index_x <- which(result_cca$u != 0)
    pred_x <- rep(0, length(true_x))
    pred_x[index_x] <- 1
    rate <- get_rate(true_x, pred_x)
    result_fpr_x[i, 2] <- rate[1]
    result_fnr_x[i, 2] <- rate[2]
    
    index_y <- which(result_cca$v != 0)
    pred_y <- rep(0, length(true_y))
    pred_y[index_y] <- 1
    rate <- get_rate(true_y, pred_y)
    result_fpr_y[i, 2] <- rate[1]
    result_fnr_y[i, 2] <- rate[2]
    
    #### PLS
    pls_model <- plsr(Y ~ X, ncomp = min(ncol(X), nrow(X)))
    # validationplot(pls_model, val.type = "MSEP")
    vip_scores <- VIP(pls_model, opt.comp = min(ncol(X), nrow(X)))
    index_x <- which(vip_scores > 1)
    pred_x <- rep(0, length(true_x))
    pred_x[index_x] <- 1
    rate <- get_rate(true_x, pred_x)
    result_fpr_x[i, 3] <- rate[1]
    result_fnr_x[i, 3] <- rate[2]
    
    pls_model <- plsr(X ~ Y, ncomp = min(ncol(Y), nrow(Y)))
    # validationplot(pls_model, val.type = "MSEP")
    vip_scores <- VIP(pls_model, opt.comp = min(ncol(Y), nrow(Y)))
    index_y <- which(vip_scores > 1)
    pred_y <- rep(0, length(true_y))
    pred_y[index_y] <- 1
    rate <- get_rate(true_y, pred_y)
    result_fpr_y[i, 3] <- rate[1]
    result_fnr_y[i, 3] <- rate[2]
    
    #### multivariate Lasso
    cv_lasso <- cv.glmnet(X, Y, family = "mgaussian", alpha = 1, type.measure = "mse")
    lambda_opt <- cv_lasso$lambda.min
    lasso_model_opt <- glmnet(X, Y, family = "mgaussian", alpha = 1, lambda = lambda_opt)
    coefficients <- coef(lasso_model_opt, s = lambda_opt)  # Extract coefficients at the optimal lambda
    coeff_matrices <- lapply(coefficients, as.matrix)
    combined_matrix <- do.call(cbind, coeff_matrices)
    combined_matrix <- combined_matrix[-1, ] # remove intercept
    pred_x <- rep(0, length(true_x))
    index_x <- which(rowSums(abs(combined_matrix)) > 0)
    pred_x[index_x] <- 1
    rate <- get_rate(true_x, pred_x)
    result_fpr_x[i, 4] <- rate[1]
    result_fnr_x[i, 4] <- rate[2]
    
    cv_lasso <- cv.glmnet(Y, X, family = "mgaussian", alpha = 1, type.measure = "mse")
    lambda_opt <- cv_lasso$lambda.min
    lasso_model_opt <- glmnet(Y, X, family = "mgaussian", alpha = 1, lambda = lambda_opt)
    coefficients <- coef(lasso_model_opt, s = lambda_opt)  # Extract coefficients at the optimal lambda
    coeff_matrices <- lapply(coefficients, as.matrix)
    combined_matrix <- do.call(cbind, coeff_matrices)
    combined_matrix <- combined_matrix[-1, ] # remove intercept
    pred_y <- rep(0, length(true_y))
    index_y <- which(rowSums(abs(combined_matrix)) > 0)
    pred_y[index_y] <- 1
    rate <- get_rate(true_y, pred_y)
    result_fpr_y[i, 4] <- rate[1]
    result_fnr_y[i, 4] <- rate[2]
    
    #### reduced rank regression
    rrr_model <- rrr(Y, X)
    coefficient <- rrr_model$coef
    sum_order <- order(rowSums(abs(coefficient)), decreasing = T)
    pred_x <- rep(0, length(true_x))
    index_x <- sum_order[1:index_x_count]
    pred_x[index_x] <- 1
    rate <- get_rate(true_x, pred_x)
    result_fpr_x[i, 5] <- rate[1]
    result_fnr_x[i, 5] <- rate[2]

    rrr_model <- rrr(X, Y)
    coefficient <- rrr_model$coef
    sum_order <- order(rowSums(abs(coefficient)), decreasing = T)
    pred_y <- rep(0, length(true_y))
    index_y <- sum_order[1:index_y_count]
    pred_y[index_y] <- 1
    rate <- get_rate(true_y, pred_y)
    result_fpr_y[i, 5] <- rate[1]
    result_fnr_y[i, 5] <- rate[2]
  }
  
  final_result <- matrix(NA, nrow = 8, ncol = 5)
  colnames(final_result) <- colnames(result_fnr_x)
  rownames(final_result) <- c("fpr_x", "fpr_x_sd", "fnr_x", "fnr_x_sd",
                              "fpr_y", "fpr_y_sd", "fnr_y", "fnr_y_sd")
  final_result[1, ] <- round(colMeans(result_fpr_x), 2)
  final_result[2, ] <- round(apply(result_fpr_x, MARGIN = 2, FUN = sd), 3)
  final_result[3, ] <- round(colMeans(result_fnr_x), 2)
  final_result[4, ] <- round(apply(result_fnr_x, MARGIN = 2, FUN = sd), 3)
  final_result[5, ] <- round(colMeans(result_fpr_y), 2)
  final_result[6, ] <- round(apply(result_fpr_y, MARGIN = 2, FUN = sd), 3)
  final_result[7, ] <- round(colMeans(result_fnr_y), 2)
  final_result[8, ] <- round(apply(result_fnr_y, MARGIN = 2, FUN = sd), 3)
  
  return(final_result)
}


work_path <- "path/to/the/project"
index <- 1

file_path <- paste0(work_path, "data_sim/")
file_path_list <- list.files(file_path)
load(paste0(file_path, file_path_list[index]))

file_path <- paste0(work_path, "result_sim/")
file_path_list <- list.files(file_path)
load(paste0(file_path, file_path_list[index]))


data <- dat_list_s2
result <- results_list
#### Change the actual selected features according to the simulation setting
true_x <- c(rep(1, 3), rep(0, 2))
true_y <- c(rep(1, 3), rep(0, 3))
result <- get_result(data, result, true_x, true_y)
print(result)


