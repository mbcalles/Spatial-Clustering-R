localmoran_cr_mc <- function(event,exposure,sw,nsim,sig_threshold){
  ## Constant Risk
  r <- sum(event)/sum(exposure)
  rni <- r * exposure
  z <- (event - rni)/sqrt(rni)
  lag_z <- lag.listw(sw, z)
  Ii_cr <- z * lag_z
  n <- length(event)

  ti <- exposure/sum(exposure)
  sims <- matrix(0, ncol = nsim, nrow = n)
  for (i in 1:nsim) {
    # y <- rpois(length(event), lambda = rni)
    y <- sample.int(n=n,size = sum(event),prob=ti,replace = TRUE) #reshuffle all events amongst spatial units based on constant risk hypothesis, with different at risk populations (e.g. what would you expect to see given same risk per unit of exposure)
    count_y <- aggregate(y,by=list(y),FUN = length)
    y <- left_join(data.frame(Group.1=1:n),count_y,by = "Group.1") %>%
      tidyr::replace_na(list(x=0)) %>% pull(x)
    zi <- (y - rni)/sqrt(rni)
    lag_zi <- lag.listw(sw, zi)
    sims[, i] <- zi * lag_zi
  }
  xrank <- apply(cbind(Ii_cr, sims), 1, function(x) rank(x)[1])
  diff <- nsim - xrank
  diff <- ifelse(diff > 0, diff, 0)
  pval <- punif((diff + 1)/(nsim + 1))

  type = factor(case_when(pval <=sig_threshold & z>0 & lag_z>0 ~ "High-High",
                          pval <=sig_threshold & z>0 & lag_z<0 ~ "High-Low",
                          pval <=sig_threshold & z<0 & lag_z>0 ~ "Low-High",
                          pval <=sig_threshold & z<0 & lag_z<0 ~ "Low-Low",
                          TRUE ~ "Insig."),
                levels = c("High-High","Low-Low","Low-High","High-Low","Insig."))#significant clusters

  return(data.frame(Ii_cr = Ii_cr,z= z,lag_z = lag_z, pval=pval,cluster_type=type))

}
