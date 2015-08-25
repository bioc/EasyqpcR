totData <- function (data, r, base, nSpl, nCTL, na.rm = na.rm) 
{
  
  if (nCTL >= 2) {
    data1 <- (colProds(data[1:nCTL, ], na.rm = na.rm))^(1/length(1:nCTL))
    data1 <- as.data.frame(data1)
  }
  else {
    data1 <- (data[1:nCTL, ])
    data1 <- as.data.frame(data1)
  }
  if (nCTL >= 2) {
    data1 <- data1
  }
  else {
    data1 <- as.numeric(data1)
  }
  data1 <- as.data.frame(data1)
  data2 <- (data1[rep(1:(ncol(data)), each = nrow(data)), ])
  data3 <- as.data.frame(matrix(data2, ncol = length(1:ncol(data)), 
                                byrow = FALSE))
  rownames(data3) <- rownames(data)
  colnames(data3) <- colnames(data)
  data4 <- data/data3
  
  colnames_removing_prefix <- function(df, prefix) {
    names <- colnames(df)
    indices <- (substr(names, 1, nchar(prefix)) == prefix)
    names[indices] <- substr(names[indices], nchar(prefix) + 
                               1, nchar(names[indices]))
  }
  if (nchar(colnames(data[1])) <= 28) {
    colnames(data) <- colnames(data)
  }
  else {
    colnames(data) <- colnames_removing_prefix(data, "NRQs.normalized.to.control.")
  }
  
  rownames(data4) <- rownames(data)
  colnames(data4) <- colnames(data)   
  
  x <- data4
  x2 <- x[order(rownames(x)), ]
  
  if (base == 2) {
    x <- log2(x)
  }
  else {
    x <- log10(x)
  }
  
  exp_avg <- aggregate(x, by = list(rep(1:r, each = (nrow(x)/r))), 
                       mean, na.rm = na.rm)[1:r, -1]
  
  meancentered <- x - exp_avg[rep(1:nrow(exp_avg), 
                                  each = (nrow(x)/r)), ]
  
  exp_sd <- aggregate(x, by = list(rep(1:r, each = (nrow(x)/r))), 
                      sd, na.rm = na.rm)[, -1]
  
  autoscaled_data <- meancentered/exp_sd[rep(1:nrow(exp_sd), 
                                             each = (nrow(x)/r)), ]
  
  mean_exp_sd <- as.data.frame(t(colMeans(exp_sd)))
  mean_exp_sd_rep <- mean_exp_sd[rep(1:nrow(mean_exp_sd), 
                                     each = nSpl), ]
  autoscaled_data_multiplied_by_mean_exp_sd <- autoscaled_data*mean_exp_sd_rep
  
  testedit <- NULL
  totmean <- NULL
  totsd <- NULL
  sem <- function(x) sqrt(var(x)/length(x))
  
  for(i in 1:(nSpl/r)){
    testedit[[i]] <- seq(i, nSpl, by = nSpl/r)
    totmean[[i]] <- colMeans(autoscaled_data_multiplied_by_mean_exp_sd[c(testedit[[i]]), ])
    totsd[[i]] <- apply(autoscaled_data_multiplied_by_mean_exp_sd[c(testedit[[i]]), ], 2, sem)
  }  
  
  autoscaled_totmean_final <- data.frame(matrix(unlist(totmean), 
                                                nrow = nSpl/r, byrow = TRUE), 
                                         stringsAsFactors = FALSE)
  rownames(autoscaled_totmean_final) <- rownames(autoscaled_data_multiplied_by_mean_exp_sd[1:(nSpl/r), ]) 
  colnames(autoscaled_totmean_final) <- colnames(autoscaled_data_multiplied_by_mean_exp_sd)                                                      
  
  autoscaled_totsd_final <- data.frame(matrix(unlist(totsd), nrow = nSpl/r, 
                                              byrow = TRUE), stringsAsFactors = FALSE)
  rownames(autoscaled_totsd_final) <- rownames(autoscaled_data_multiplied_by_mean_exp_sd[1:(nSpl/r), ]) 
  colnames(autoscaled_totsd_final) <- colnames(autoscaled_data_multiplied_by_mean_exp_sd)                                                      
  
  lower95CI <- autoscaled_totmean_final - qt((1 - 0.05/2), r-1)*autoscaled_totsd_final
  upper95CI <- autoscaled_totmean_final + qt((1 - 0.05/2), r-1)*autoscaled_totsd_final
  
  if (base == 2) {
    lin_autoscaled_data_multiplied_by_mean_exp_sd <- 2^autoscaled_data_multiplied_by_mean_exp_sd
    lin_autoscaled_totmean_final <- 2^autoscaled_totmean_final
    lin_lower95CI <- 2^lower95CI
    lin_upper95CI <- 2^upper95CI
  }
  else {
    lin_autoscaled_data_multiplied_by_mean_exp_sd <- 10^autoscaled_data_multiplied_by_mean_exp_sd
    lin_autoscaled_totmean_final <- 10^autoscaled_totmean_final
    lin_lower95CI <- 10^lower95CI
    lin_upper95CI <- 10^upper95CI
  }
  
  
  return(list('Linearized Standardized data according to Willems procedure (2008)' = lin_autoscaled_data_multiplied_by_mean_exp_sd,
              `Mean` = lin_autoscaled_totmean_final, 
              `Lower Confidence Interval 95%` = lin_lower95CI,
              `Upper Confidence Interval 95%` = lin_upper95CI))
  
}
