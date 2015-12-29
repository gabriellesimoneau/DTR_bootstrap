###################################################################################
# to apply the m-out-of-n bootstrap, we need to have an estimate of the "amount of irregularity"
# in the data
# Irregularity occurs when the second stage treatment has a very small effect on the treatment 
# decision
#
# Following Chakraborty et al (2013), this occurs when phi*x is close to zero because then the function 
# f(x) = Indicator(phi*x >0) is not differentiable
#
# Goal here: how to estimate p = (probability that phi*x will be close to zero)
# Idea:
#       - fit the dWOLS model to get estimates of phi. Then, using these estimates, calculate
#         \hat phi*x for each observations. Then, \hat p = proportion of observations that have
#         \hat phi*x "close" to zero.        
#       - "close" to zero is subjective. Might want to vary the threshold.
# Other idea (as suggested by Wallace et al in JSS):
#       - "non-regularity occurs when optimal treatment is not unique. (...) [estimating irregularity]
#          involves identifying the proportion of subjects for whom, when all possible blip parameter
#          values within their respective confidence sets are considered, both treatment and 
#          non-treatment could be recommended."
#
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

# scenario
sc <- seq(1,9)

################################### scenario 3 - nonregular ###################################

n <- 300
i <- 3

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

# model specification
blip.model <- list(~ O1, ~ O2 + A1)
proba <- list(as.vector(rep(0.5,n)))
treat.model <- list(A1 ~ 1, A2 ~ 1) 
tf.model <- list(~ O1, ~ O1 + A1 + O1*A1)

# fit dWOLS to the generated dataset, using all n=300 observations
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

# estimate of p, the probability of generating data with gamma5 + gamma6*O2 + gamma7*A1 close to zero
#   in scenario 3, this probability should be close to 0.5
# try different threshold to quantify "close to zero"
#   the estimates of p varies a lot depending on the data
length(psi[which(abs(psi) < 0.1)])/n 
length(psi[which(abs(psi) < 0.15)])/n 

# probability of generating patient history such that g5*A2 + g6*O2*A2 + g7*A1*A2 = 0
#   this is, following the paper where A1,A2 are coded {-1,1} but this specificiation of p
#   is not relevant when A1,A2 are coded {0,1} because p will always be 0.5
gg <- g[sc[i],5]*A2 + g[sc[i],6]*O2*A2 + g[sc[i],7]*A1*A2
length(which(gg==0))/n



################################### scenario 5 - nonregular ###################################

i <- 5

# treatment A1, A2: P(Aj = 1) = P(Aj = 0) = 0.5
A1 <- rbinom(n, size = 1, prob = 0.5)
A2 <- rbinom(n, size = 1, prob = 0.5)

# treatment A1 coded as -1,1 so I don't have to adapt the delta_1 and delta_2 parameters
A1.min <- 2*A1 - 1

# covariates O1, O2: coded as -1, 1, where O2 depends on A1, O1 and (delta_1,delta_2)
O1 <- 2*rbinom(n, size = 1, prob = 0.5) - 1
O2 <- 2*rbinom(n, size = 1, prob = expit(d[3,1]*O1 + d[3,2]*A1.min)) - 1

# generated outcome Y2 (Y1 set to 0), using parameters (gamma_1,...,gamma_7)
Y1 <- rep(0, n)
Y2 <- g[5,1] + g[5,2]*O1 + g[5,3]*A1 + g[5,4]*O1*A1 + g[5,5]*A2 + g[5,6]*O2*A2 + g[5,7]*A1*A2 + rnorm(n)

# model specification
blip.model <- list(~ O1, ~ O2 + A1)
proba <- list(as.vector(rep(0.5,n)))
treat.model <- list(A1~1, A2~1) 
tf.model <- list(~ O1, ~ O1 + A1 + O1*A1)

# fit dWOLS to the generated dataset, using all n=300 observations
s5 <- DTRreg(outcome = Y2, blip.mod = blip.model, treat.mod = treat.model, tf.mod = tf.model, treat.mod.man = rep(proba,2), method = "dwols")
summary(s5)

# calculate phi*x ~ p -- should be close to 0.25
int1 <- s5["psi"][[1]][[1]][1]
B.o1 <- s5["psi"][[1]][[1]][2]
int2 <- s5["psi"][[1]][[2]][1]
B.o2 <- s5["psi"][[1]][[2]][2]
B.a1 <- s5["psi"][[1]][[2]][3]

psi <- int2 + O2*B.o2 + A1*B.a1
psi

# estimate of p, the probability of generating data with gamma5 + gamma6*O2 + gamma7*A1 close to zero
#   in scenario 5, this probability should be close to 0.25
# try different threshold to quantify "close to zero"
#   the estimates of p varies a lot depending on the data
length(psi[which(abs(psi) < 0.1)])/n 
length(psi[which(abs(psi) < 0.15)])/n 
