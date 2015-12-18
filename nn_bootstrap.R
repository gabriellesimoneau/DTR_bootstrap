###################################################################################
library(DTRreg)
expit <- function(x) exp(x)/(1+exp(x))

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

extract <-function(out)
{
  psi11 <-out["psi"][[1]][[1]][2] # gamma
  psi21 <-out["psi"][[1]][[2]][3]
  return(c(psi11, psi21))
}

######################### regular n-out-of-n bootstrap #############################

# scenario 5
sc <- 5

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

# allocate space
est <- vector(mode = "list", length = 2)
for(i in 1:2)
{
  est[[i]] <- matrix(NA, nrow = Nsimul, ncol = 1 + Nboot)
}

for(s in 1:Nsimul)
{
  # generate data
  A1 <- rbinom(n, size = 1, prob = 0.5)
  A2 <- rbinom(n, size = 1, prob = 0.5)
  
  A1.min <- 2*A1 - 1
  
  O1 <- 2*rbinom(n, size = 1, prob = 0.5) - 1
  O2 <- 2*rbinom(n, size = 1, prob = expit(d[sc,1]*O1 + d[sc,2]*A1.min)) - 1
  
  Y2 <- g[sc,1] + g[sc,2]*O1 + g[sc,3]*A1 + g[sc,4]*O1*A1 + g[sc,5]*A2 + g[sc,6]*O2*A2 + g[sc,7]*A1*A2 + rnorm(n)

  complete <- cbind(id,A1, A2, O1, O2, Y2)
  
  res.n <- DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(complete))
  es <- extract(res.n)
  
  est[[1]][s,1] <- es[1]
  est[[2]][s,1] <- es[2]
  
  for(b in 1:Nboot)
  {
    index <- sample(1:n, n,replace = TRUE)
    boot <- complete[index,]
    
    res <- DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols", data = as.data.frame(boot))
    esb <- extract(res)
    est[[1]][s, b + 1] <- esb[1]
    est[[2]][s, b + 1] <- esb[2]
  }
}

