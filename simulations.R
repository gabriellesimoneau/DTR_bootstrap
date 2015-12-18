###################################################################################
library(DTRreg)
expit <- function(x) exp(x)/(1+exp(x))
# m-out-of-n bootstrap

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

################################### scenario 3 - nonregular ###################################
# dWOLS notation -- treatment coded 0,1
n <- 300
A1 <- rbinom(n, size = 1, prob = 0.5)
A2 <- rbinom(n, size = 1, prob = 0.5)

A1.min <- 2*A1 - 1

O1 <- 2*rbinom(n, size = 1, prob = 0.5) - 1
O2 <- 2*rbinom(n, size = 1, prob = expit(d[3,1]*O1 + d[3,2]*A1.min)) - 1

Y1 <- rep(0, n)
Y2 <- g[3,1] + g[3,2]*O1 + g[3,3]*A1 + g[3,4]*O1*A1 + g[3,5]*A2 + g[3,6]*O2*A2 + g[3,7]*A1*A2 + rnorm(n)

# fit the model
blip.model <- list(~ O1, ~ O2 + A1)
proba <- list(as.vector(rep(0.5,n)))
treat.model <- list(A1~1, A2~1) 
tf.model <- list(~ O1, ~ O1 + A1 + O1*A1)

s3 <- DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols")
summary(s3)

# calculate phi*x ~ p -- should be close to 0.5
int1 <- s3["psi"][[1]][[1]][1]
B.o1 <- s3["psi"][[1]][[1]][2]
int2 <- s3["psi"][[1]][[2]][1]
B.o2 <- s3["psi"][[1]][[2]][2]
B.a1 <- s3["psi"][[1]][[2]][3]

psi <- int2 + O2*B.o2 + A1*B.a1
psi

length(psi[which(abs(psi) < 0.1)])/n 

################################### scenario 5 - nonregular ###################################
A1 <- rbinom(n, size = 1, prob = 0.5)
A2 <- rbinom(n, size = 1, prob = 0.5)

A1.min <- 2*A1 - 1

O1 <- 2*rbinom(n, size = 1, prob = 0.5) - 1
O2 <- 2*rbinom(n, size = 1, prob = expit(d[5,1]*O1 + d[5,2]*A1.min)) - 1

Y1 <- rep(0, n)
Y2 <- g[5,1] + g[5,2]*O1 + g[5,3]*A1 + g[5,4]*O1*A1 + g[5,5]*A2 + g[5,6]*O2*A2 + g[5,7]*A1*A2 + rnorm(n)

blip.model <- list(~ O1, ~ O2 + A1)
proba <- list(as.vector(rep(0.5,n)))
treat.model <- list(A1~1, A2~1) 
tf.model <- list(~ O1, ~ O1 + A1 + O1*A1)

s5 <- DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols")
summary(s5)

# calculate phi*x ~ p -- should be close to 0.25
int1 <- s5["psi"][[1]][[1]][1]
B.o1 <- s5["psi"][[1]][[1]][2]
int2 <- s5["psi"][[1]][[2]][1]
B.o2 <- s5["psi"][[1]][[2]][2]
B.a1 <- s5["psi"][[1]][[2]][3]

psi <- int2 + O2*B.o2 + A1*B.a1
quantile(abs(psi), probs = 0.25)
psi
length(psi[which(abs(psi) < 0.1)])/n 













