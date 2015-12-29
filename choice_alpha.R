###################################################################################
#library(DTRreg)
source("dtrreg_fun.R")

expit <- function(x) exp(x)/(1+exp(x))

# gamma parameters following Chakraborty et al (2013) to control for irregularity in the generated data
g <- matrix(NA, nrow = 9, ncol = 7)
g[1,] <- c(0,0,0,0,0,0,0)
g[2,] <- c(0,0,0,0,0.01,0,0)
g[3,] <- c(0,0,-0.5,0,0.5,0,-0.5)
g[4,] <- c(0,0,-0.5,0,0.99,0,-0.98)
g[5,] <- c(0,0,-0.5,0,1,0.5,-0.5)
g[6,] <- c(0,0,-0.5,0,0.25,0.5,0.5)
g[7,] <- c(0,0,-0.25,0,0.75,0.5,0.5)
g[8,] <- c(0,0,0,0,1,0,-1)
g[9,] <- c(0,0,0,0,0.25,0,-0.24)

# delta parameters following Chakraborty et al (2013) to control for irregularity in the generated data
d <- matrix(NA, nrow = 9, ncol = 2)
d[1,] <- c(0.5,0.5)
d[2,] <- c(0.5,0.5)
d[3,] <- c(0.5,0.5)
d[4,] <- c(0.5,0.5)
d[5,] <- c(1,0)
d[6,] <- c(0.1,0.1)
d[7,] <- c(0.1,0.1)
d[8,] <- c(0,0)
d[9,] <- c(0,0)

# extract is a function to extract the treatment estimates phi_1 (effect of treatment first stage) 
#   and phi_2 (effect of treatment second stage) from a DTRreg mdoel fit
extract <-function(out)
{
  psi11 <-out["psi"][[1]][[1]][1] 
  psi21 <-out["psi"][[1]][[2]][1]
  B.o2 <- out["psi"][[1]][[2]][2]
  B.a1 <- out["psi"][[1]][[2]][3]
  return(c(psi11, psi21, B.o2, B.a1))
}

# percentile is to extract the eta-level percentiles from step 4 of double bootstrap
percentile <- function(x)
{
  phi1n <- x[1] 
  m <- floor(x[3])
  dis <- x[4:(B2+3)]-phi1n
  quan <- sqrt(m) * quantile(dis, probs = c(0.025, 0.975))
  return(quan)
}

# dbCI to construct the double bootstrap sample from step 4
dbCI <- function(x)
{
  ph1n <- x[1]
  m <- x[2]
  l <- x[3]
  u <- x[4]
  CI <- c(ph1n - u/sqrt(m), ph1n - l/sqrt(m))
  return(CI)
}
######################### m-out-of-n bootstrap : adaptive choice of alpha #############################

#### choose alpha (will be different for each scenario)

# scenario id
sc <- seq(1,9)
# number of simulated dataset
Nsimul <- 1000 
# number of first-stage bootstrap samples following Chakraborty et al (2013)
B1 <- 500
# number of second-stage bootstrap samples
B2 <- 500
# sample size
n <- 300

# model specification
blip.model <- list(~ O1, ~ O2 + A1)
treat.model <- list(A1~1, A2~1) 
tf.model <- list(~ O1, ~ O1 + A1 + O1*A1)

# allocate space
phi1n <- NA 
phi1nb1 <- NA
selected.alpha <- rep(NA, 9)

# grid alpha
al <- seq(0.025, 1, by = 0.025)
maxit <- length(al)

for(i in 1:9)
{
  coverage <- 0
  it <- 0
  while(coverage < 0.95 & it < maxit)
  {
    it <- it + 1
    alpha <- al[it]
    
    # reset estimates to NA for new alpha
    # estm[,1] = phi(b1), estm[,2] = phat, estm[,3]=m(b1)
    # remaining estm columns = phi(b1,b2)
    estm <- matrix(NA, ncol = B2 + 3, nrow = B1)
    
    # generate data for each scenario
    # treatment A1, A2: P(Aj = 1) = P(Aj = 0) = 0.5
    A1 <- rbinom(n, size = 1, prob = 0.5)
    A2 <- rbinom(n, size = 1, prob = 0.5)
    
    # treatment A1 coded as -1,1 so I don't have to adapt the delta_1 and delta_2 parameters
    A1.min <- 2*A1 - 1
    
    # covariates O1, O2: coded as -1, 1, where O2 depends on A1, O1 and (delta_1,delta_2)
    O1 <- 2*rbinom(n, size = 1, prob = 0.5) - 1
    O2 <- 2*rbinom(n, size = 1, prob = expit(d[sc[i],1]*O1 + d[sc[i],2]*A1.min)) - 1
    
    # generated outcome Y2 (Y1 set to 0), using parameters (gamma_1,...,gamma_7)
    Y2 <- g[sc[i],1] + g[sc[i],2]*O1 + g[sc[i],3]*A1 + g[sc[i],4]*O1*A1 + g[sc[i],5]*A2 + g[sc[i],6]*O2*A2 + g[sc[i],7]*A1*A2 + rnorm(n)
    
    # generated dataset
    complete <- cbind(A1, A2, O1, O2, Y2)
    
    # estimate phi_1n with all observations
    proba <- list(as.vector(rep(0.5,n)))
    res.n <- try(dtrreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(complete)))
    es <- try(extract(res.n))
    phi1n <- es[2]
    
    for(j in 1:B1) # loop over B1 first stage bootstrap samples
    {
      # draw a n-out-of-n bootstrap sample
      index <- sample(1:n, n, replace = TRUE)
      boot1 <- as.data.frame(complete[index,])
      
      # fit the model to b1-th bootstrap sample
      proba <- list(as.vector(rep(0.5,n)))
      res1 <- try(dtrreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = boot1))
      esb1 <- try(extract(res1))
      estm[j,1] <- esb1[2]
      
      # estimate m for each b1 bootstrap sample
      t2 <- esb1[2] + boot1$O2 * esb1[3] + boot1$A1 * esb1[4]
      phat <- length(t2[which(abs(t2) < 0.1)])/n # subjective threshold
      estm[j,2] <- phat
      
      # resampling size
      m <- n^((1 + alpha*(1-phat))/(1 + alpha)) 
      estm[j,3] <- m
      
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
        estm[j,k + 3] <- esb2[2]
      }
    }
    
    # for each B1 first stage bootstrap sample, calculate 0.025, 0.975 percentiles
    ul <- t(apply(estm, 1, percentile))
    temp <- cbind(estm[,1], estm[,3], ul[,1], ul[,2])
    CI <- t(apply(temp, 1, dbCI))
    temp1 <- cbind(CI, rep(phi1n,B1))
    ind <- ifelse(temp1[,1] < temp1[,3] & temp1[,2] > temp1[,3], 1, 0)
    coverage <- mean(ind)
    print(c(alpha,coverage))
  }
  # return alpha for scenario
  selected.alpha[i] <- alpha
}
write.csv(selected.alpha, file = "alpha.csv")