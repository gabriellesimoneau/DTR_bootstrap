dbalpha <- function(data, psin, B1 = 500, B2 = 500)
{
  # grid alpha
  grid <- seq(0.025, 1, by = 0.025)
  maxit <- length(grid)
  
  # initialize
  coverage <- 0
  it <- 0
  while(coverage < 0.95 & it < maxit)
  {
    it <- it + 1
    alpha <- grid[it]
    est <- matrix(NA, ncol = B2 + 3, nrow = B1)
    
    for(j in 1:B1) # loop over B1 first stage bootstrap samples
    {
      # draw a n-out-of-n bootstrap sample
      index <- sample(1:n, n, replace = TRUE)
      boot1 <- as.data.frame(data[index,])
      
      # fit the model to b1-th bootstrap sample
      proba <- list(as.vector(rep(0.5,n)))
      res1 <- try(dtrreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = boot1))
      esb1 <- try(extract(res1))
      est[j,1] <- esb1[2]
      
      # estimate m for each b1 bootstrap sample
      t2 <- esb1[2] + boot1$O2 * esb1[3] + boot1$A1 * esb1[4]
      phat <- length(t2[which(abs(t2) < 0.1)])/n # subjective threshold
      est[j,2] <- phat
      
      # resampling size
      m <- n^((1 + alpha*(1-phat))/(1 + alpha)) 
      est[j,3] <- m
      
      # probability treatment with m
      proba <- list(as.vector(rep(0.5,floor(m))))
      
      for(k in 1:B2) # loop over B2 second stage bootstrap samples
      {
        # resample with replacement 
        index <- sample(1:n, floor(m), replace = TRUE)
        boot2 <- boot1[index,]
        
        # fit the model to bootstrap sample
        res2 <- try(dtrreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(boot2)))
        esb2 <- try(extract(res2))
        
        # save the (b1,b2) bootstrap estimates j in the (k+3) column
        est[j,k + 3] <- esb2[2]
      }
    }
    # for each B1 first stage bootstrap sample, calculate 0.025, 0.975 percentiles
    ul <- t(apply(est, 1, percentile))
    temp <- cbind(est[,1], est[,3], ul[,1], ul[,2])
    CI <- t(apply(temp, 1, dbCI))
    temp1 <- cbind(CI, rep(psin,B1))
    ind <- ifelse(temp1[,1] < temp1[,3] & temp1[,2] > temp1[,3], 1, 0)
    coverage <- mean(ind)
  }
  return(as.numeric(alpha))
}