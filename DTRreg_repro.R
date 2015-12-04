library(DTRreg)

#### simulated example from "Introducing R package DTRreg" ####
set.seed(1)
expit <- function(x){1/(1+exp(-x))}
n <- 10000
x1 <- rnorm(n)
x2 <- rnorm(n)
a1 <- rbinom(n,1,expit(x1))
a2 <- rbinom(n,1,expit(x2))
gamma1 <- a1*(1+x1) # blip stage 1
gamma2 <- a2*(1+x2) # blip stage 2
y <- exp(x1)+exp(x2)+gamma1+gamma2+rnorm(n)

blip.mod <- list(~x1,~x2)
treat.mod <- list(a1~x1,a2~x2)
tf.mod <- list(~x1,~x2)
weight.fun <- function(w){1/w}

obj=DTRreg(y,blip.mod,treat.mod,tf.mod, method="dwols",weight=weight.fun)
summary(obj)

#### simulation result in Table 1 ####
set.seed(1)
expit <- function(x){1/(1+exp(-x))}
weight.fun <- function(w){1/w}
extract <-function(out)
{
  psi10 <-out$psi[[1]][1]
  psi11 <-out$psi[[1]][2]
  psi20 <-out$psi[[2]][1]
  psi21 <-out$psi[[2]][2]
  return(c(psi10,psi11,psi20,psi21))
}

N <- 1000 
n <- 10000
out <- matrix(NA,nrow=N,ncol=16)

for(i in 1:N)
{
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  a1 <- rbinom(n,1,expit(x1))
  a2 <- rbinom(n,1,expit(x2))
  gamma1 <- a1*(1+x1) # blip stage 1
  gamma2 <- a2*(1+x2) # blip stage 2
  y <- exp(x1)+exp(x2)+gamma1+gamma2+rnorm(n)
  
  # both correct
  blip.mod <- list(~x1,~x2)
  treat.mod <- list(a1~x1,a2~x2)
  tf.mod <- list(~exp(x1),~exp(x2))
  out1 <- DTRreg(y,blip.mod,treat.mod,tf.mod, method="dwols",weight=weight.fun)
  
  # treatment correct / treatment-free misspecified
  tf.mod <- list(~x1,~x2)
  out2 <- DTRreg(y,blip.mod,treat.mod,tf.mod, method="dwols",weight=weight.fun)
  
  # treatment-free correct / treatment misspecified
  treat.mod <- list(a1~0.5,a2~0.5) ## problem HERE 
  tf.mod <- list(~exp(x1),~exp(x2))
  out3 <- DTRreg(y,blip.mod,treat.mod,tf.mod, method="dwols",weight=weight.fun)
  
  # neither correct
  tf.mod <- list(~x1,~x2)
  out4 <- DTRreg(y,blip.mod,treat.mod,tf.mod, method="dwols",weight=weight.fun)
 
  # result simulation
  matrix[i,] <-lapply(list(out1,out2,out3,out4),extract)
}



