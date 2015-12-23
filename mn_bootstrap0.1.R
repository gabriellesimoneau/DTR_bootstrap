###################################################################################
library(DTRreg)
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
g[8,] <- c(0,0,0,0,0.25,0,0.25)
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
  psi11 <-out["psi"][[1]][[1]][2] 
  psi21 <-out["psi"][[1]][[2]][3]
  int2 <- out["psi"][[1]][[2]][1]
  B.o2 <- out["psi"][[1]][[2]][2]
  return(c(psi11, psi21, int, B.o2))
}

######################### m-out-of-n bootstrap : fixed alpha #############################

# fixed alpha 
alpha <- 0.1
# scenario id
sc <- seq(1,9)
# number of simulated dataset
Nsimul <- 1000 
# number of boostrap samples
Nboot <- 1000
# sample size
n <- 300

# model specification
blip.model <- list(~ O1, ~ O2 + A1)
proba <- list(as.vector(rep(0.5,n)))
treat.model <- list(A1~1, A2~1) 
tf.model <- list(~ O1, ~ O1 + A1 + O1*A1)

# allocate space: estm[[1]] -> Nsimul by Nboot+1 matrix, estm[[2]] -> Nsimul by Nboot+1 matrix,
#                 estm[[1]] -> save phi_1 treatment effect stage 1, estm[[2]] -> ave phi_1 treatment effect stage 2
#                 first column of estm[[1]] and estm[[2]] -> estimates using all observations (for each simulated dataset)
#                 second column of estm[[1]] and estm[[2]] -> estimate of nonregularity p hat
#                 third column of estm[[1]] and estm[[2]] -> resample size b
#                 next columns -> 1000 m-out-of n bootstrap estimates (for each simulated dataset)
estm <- vector(mode = "list", length = 2)

for(i in sc)
{
  # reset estimates to NA for new scenario
  for(k in 1:2)
  {
    estm[[k]] <- matrix(NA, nrow = Nsimul , ncol = Nboot + 3)
  }
  for(s in 1:Nsimul)
  {
    # treatment A1, A2: P(Aj = 1) = P(Aj = 0) = 0.5
    A1 <- rbinom(n, size = 1, prob = 0.5)
    A2 <- rbinom(n, size = 1, prob = 0.5)
    
    # treatment A1 coded as -1,1 so I don't have to adapt the delta_1 and delta_2 parameters
    A1.min <- 2*A1 - 1
    
    # covariates O1, O2: coded as -1, 1, where O2 depends on A1, O1 and (delta_1,delta_2)
    O1 <- 2*rbinom(n, size = 1, prob = 0.5) - 1
    O2 <- 2*rbinom(n, size = 1, prob = expit(d[sc,1]*O1 + d[sc,2]*A1.min)) - 1
    
    # generated outcome Y2 (Y1 set to 0), using parameters (gamma_1,...,gamma_7)
    Y2 <- g[sc,1] + g[sc,2]*O1 + g[sc,3]*A1 + g[sc,4]*O1*A1 + g[sc,5]*A2 + g[sc,6]*O2*A2 + g[sc,7]*A1*A2 + rnorm(n)
    
    # generated dataset
    complete <- cbind(A1, A2, O1, O2, Y2)
    
    # fit dWOLS to the generated dataset, using all n=300 observations
    res.n <- try(DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(complete)))
    es <- try(extract(res.n))
    
    # save estimates using all observations in the first column
    estm[[1]][s,1] <- es[1]
    estm[[2]][s,1] <- es[2]
    
    # estimate of nonregularity
    t2 <- es[3] + O2*es[4] + A1*es[2]
    phat <- length(t2[which(abs(t2) < 0.1)])/n # subjective threshold
    estm[[1]][s,2] <- phat
    estm[[2]][s,2] <- phat
    
    # resampling size
    m <- n^((1 + alpha*(1-phat))/(1 + alpha))
    estm[[1]][s,3] <- m
    estm[[2]][s,3] <- m
    
    # bootstrap resampling + estimate
    for(b in 1:Nboot)
    {
      # resample with replacement 
      index <- sample(1:n, m, replace = TRUE)
      boot <- complete[index,]
      
      # fit the model to bootstrap sample
      res <- try(DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(boot)))
      esb <- try(extract(res))
      
      # save bootstrap estimates i in the (i+1) column
      estm[[1]][s, b + 3] <- esb[1]
      estm[[2]][s, b + 3] <- esb[2]
    }
  }
  # linux command to save results of the simulations in a CSV file
  name1 <- paste("mn_psi1_scenario", paste(sc[i]),".csv",sep ="")
  name2 <- paste("mn_psi2_scenario", paste(sc[i]),".csv",sep ="")
  
  write.csv(estm[[1]], file = name1)
  write.csv(estm[[2]], file = name2)
}






















