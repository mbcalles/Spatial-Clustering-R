localmoran_mc <- function(nsim,x,sw,sig_threshold){

  n <- length(x)#number of polylines
  mew <- mean(x) #mean of attribute over all polylines
  y <- x-mew #centralized
  p_var <- sum((x-mew)^2)*(1/n)#population variance
  psd <- sqrt(p_var)#population standard deviation
  z <- y/psd #standardized values vector (how far does value fall from the mean)
  lag_z <- lag.listw(sw,z,zero.policy = TRUE) #Lagged Local Moran's I values for neighbours
  Ii <- z*lag_z

  sim <- matrix(0, ncol = nsim, nrow = n)
  for (i in 1:nsim) {
    x_s <- sample(x) #Permute values of x
    y_s <- x_s-mew #centralized
    z_s <- y_s/psd #standardized values vector (how far does value fall from the mean)
    lag_z_s <- lag.listw(sw,z_s,zero.policy = TRUE)
    sim[, i] <- z*lag_z_s
  }
  xrank <- apply(cbind(Ii, sim), 1, function(x) rank(x)[1])
  diff <- nsim - xrank
  diff <- ifelse(diff > 0, diff, 0)
  pval <- punif((diff + 1)/(nsim + 1))

  type = factor(case_when(pval <=sig_threshold & z>0 & lag_z>0 ~ "High-High",
                          pval <=sig_threshold & z>0 & lag_z<0 ~ "High-Low",
                          pval <=sig_threshold & z<0 & lag_z>0 ~ "Low-High",
                          pval <=sig_threshold & z<0 & lag_z<0 ~ "Low-Low",
                          TRUE ~ "Insig."),
                levels = c("High-High","Low-Low","Low-High","High-Low","Insig."))#significant clusters

  return(data.frame(Ii = Ii,z=z,lag_z = lag_z,pval=pval,cluster_type = type))

}
