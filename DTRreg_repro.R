library(DTRreg)

################### simulated example 4.1 from "Introducing R package DTRreg" ###################

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

################### simulated example 4.2 from "Introducing R package DTRreg" ###################

set.seed(1)
n <- 10000
x1 <- rnorm(n)
x2 <- rnorm(n)
a1 <- rnorm(n, x1, 0.5)
a2 <- rnorm(n, x2, 0.5)
gamma1 <- a1 * (1 + x1 - a1)
gamma2 <- a2 * (1 + x2 - a2)
y <- exp(x1) + exp(x2) + gamma1 + gamma2 + rnorm(n)

blip.mod <- list(~ x1 + a1, ~ x2 + a2)
treat.mod <- list(a1 ~ x1, a2 ~ x2)
tf.mod <- list( ~ x1 ,~ x1 + a1 + x2)

obj=DTRreg(y,blip.mod,treat.mod,tf.mod,var.estim = "bootstrap", B = 200, treat.range = c(-5,5))
summary(obj)

############################ simulation result in Table 1 ################################# 
set.seed(1)
expit <- function(x){1/(1+exp(-x))}
weight.fun <- function(w){1/w}
extract <-function(out)
{
  psi10 <-out["psi"][[1]][[1]][1]
  psi11 <-out["psi"][[1]][[1]][2]
  psi20 <-out["psi"][[1]][[2]][1]
  psi21 <-out["psi"][[1]][[2]][2]
  return(c(psi10,psi11,psi20,psi21))
}

N <- 1000 
n <- 10000
result <- matrix(NA,nrow=N,ncol=16*3-2*4)

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
  out1.dwols <- DTRreg(y, blip.mod, treat.mod, tf.mod, method = "dwols", weight = weight.fun)
  out1.g <- DTRreg(y, blip.mod, treat.mod, tf.mod, method = "gest", weight = weight.fun)
  # qlearn -> treatment free-model well specified
  out1.qlearn <- DTRreg(y, blip.mod, treat.mod, tf.mod, method = "qlearn") # qlearn=unweighted
  
  # treatment correct / treatment-free misspecified
  tf.mod <- list(~x1,~x2)
  out2.dwols <- DTRreg(y, blip.mod, treat.mod, tf.mod, method = "dwols", weight = weight.fun)
  out2.g <- DTRreg(y, blip.mod, treat.mod, tf.mod, method = "gest", weight = weight.fun)
  # qlearn -> treatment free-model well specified
  out2.qlearn <- DTRreg(y, blip.mod, treat.mod, tf.mod, method = "qlearn")
  
  # treatment-free correct / treatment misspecified
  # need to use treat.mod.man to specify treatment probability of 0.5 at each stage / for each individual
  treat.mod <- list(a1~1,a2~1) # dummy treatment model 
  proba <- list(as.vector(rep(0.5,n)))
  tf.mod <- list(~exp(x1),~exp(x2))
  out3.dwols <- DTRreg(y, blip.mod, treat.mod, tf.mod, method="dwols", weight = weight.fun, treat.mod.man = rep(proba,2))
  out3.g <- DTRreg(y, blip.mod, treat.mod, tf.mod, method="gest", weight = weight.fun, treat.mod.man = rep(proba,2))
  
  # neither correct
  tf.mod <- list(~x1,~x2)
  out4.dwols <- DTRreg(y, blip.mod, treat.mod, tf.mod, method = "dwols", weight = weight.fun, treat.mod.man = rep(proba,2))
  out4.g <- DTRreg(y, blip.mod, treat.mod, tf.mod, method="gest", weight = weight.fun, treat.mod.man = rep(proba,2))

  # result simulation
  result[i,] <- unlist(lapply(list(out1.dwols,out1.g,out1.qlearn,out2.dwols,out2.g,out2.qlearn,out3.dwols,out3.g,out4.dwols,out4.g),extract))
  print(i)
}

# Rebuild TABLE 1 from the paper -- Results are close enough
gesti=cbind(result[,5:8],result[,17:20],result[,29:32],result[,37:40])
dwols=cbind(result[,1:4],result[,13:16],result[,25:28],result[,33:36])
qle=cbind(result[,9:12],result[,21:24])
methods=c("G-estimation","","","","dWOLS","","","","Q-learning","")
models=c("Both correct","Treatment correct","TF correct","Neither correct","Both correct","Treatment correct","TF correct","Neither correct","TF correct","TF incorrect")

TABLE1=as.data.frame(cbind(methods,models,matrix(round(c(apply(gesti,2,mean),apply(dwols,2,mean),apply(qle,2,mean)),3),nrow=10,byrow=TRUE)))
colnames(TABLE1)[3:6]=c("Stage 1:int","Stage 1:psi","Stage 2:int","Stage 2:psi")
TABLE1
