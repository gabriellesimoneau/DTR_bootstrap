# reproduce result tables in Chakraborty et al. (2013) with dWOLS as method of analysis

# Table 2: average bootstrap resample size
# Table 3: coverage rate for 1000 simulated dataset
# Table 4: average width for 1000 simulated dataset

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

# percentile_nn for regular bootstrap
percentile_nn <- function(x)
{
  phi1n <- x[2] 
  dis <- x[3:1002]-phi1n
  quan <- sqrt(300) * quantile(dis, probs = c(0.025, 0.975))
  return(quan)
}

# construct CI for n-out-of-n regular bootstrap
consCI_nn <- function(x)
{
  ph1n <- x[1]
  l <- x[2]
  u <- x[3]
  CI <- c(ph1n - u/sqrt(300), ph1n - l/sqrt(300))
  return(CI)
}

# percentile_mn for regular bootstrap
percentile_mn <- function(x)
{
  phi1n <- x[2] 
  m <- x[4]
  dis <- x[5:1004]-phi1n
  quan <- sqrt(m) * quantile(dis, probs = c(0.025, 0.975))
  return(quan)
}

# construct CI for n-out-of-n regular bootstrap
consCI_mn <- function(x)
{
  ph1n <- x[1]
  m <- x[2]
  l <- x[3]
  u <- x[4]
  CI <- c(ph1n - u/sqrt(m), ph1n - l/sqrt(m))
  return(CI)
}

# scenarios
sc <- seq(1,9)

######################### for psi2 #############################
pos = 5
setwd("/Users/gabriellesimoneau/Dropbox/McGill - PhD/Fall2015/MATH680/Final/code/scenario_results")

table2.2 <- matrix(NA, nrow = 4, ncol = 9)
table3.2 <- matrix(NA, nrow = 4, ncol = 9)
table4.2 <- matrix(NA, nrow = 4, ncol = 9)
rownames(table2.2) <- c("nn","mn0.05","mn0.1","mnAD")
rownames(table3.2) <- c("nn","mn0.05","mn0.1","mnAD")
rownames(table4.2) <- c("nn","mn0.05","mn0.1","mnAD")
colnames(table2.2) <- paste("sc",paste(sc),sep = "")
colnames(table3.2) <- paste("sc",paste(sc),sep = "")
colnames(table4.2) <- paste("sc",paste(sc),sep = "")

# n-out-of-n bootstrap
sc <- seq(1,9)
di <- c("nn_psi2_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.2[1,i] <- 300 # resampling size always 300
  ul <- t(apply(dat, 1, percentile_nn))
  dat2 <- cbind(dat[,2], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_nn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  print(prop.test(sum(cov),1000,p=0.95)$p.value)
  table3.2[1,i] <- mean(cov)
  table4.2[1,i] <- mean((CI[,2]-CI[,1]))
}

# m-out-of-n bootstrap alpha=0.05
di <- c("mn0.05_psi2_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.2[2,i] <- mean(dat[,4]) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat[,2], dat[,4], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  table3.2[2,i] <- mean(cov)
  table4.2[2,i] <- mean((CI[,2]-CI[,1]))
}

# m-out-of-n bootstrap alpha=0.1
di <- c("mn0.1_psi2_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.2[3,i] <- mean(dat[,4]) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat[,2], dat[,4], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  table3.2[3,i] <- mean(cov)
  table4.2[3,i] <- mean((CI[,2]-CI[,1]))
}

# m-out-of-n bootstrap adaptive alpha
di <- c("mnad_psi2_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.2[4,i] <- mean(dat[,4]) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat[,2], dat[,4], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  table3.2[4,i] <- mean(cov)
  table4.2[4,i] <- mean((CI[,2]-CI[,1]))
}

######################### for psi1 #############################
pos = 3
setwd("/Users/gabriellesimoneau/Dropbox/McGill - PhD/Fall2015/MATH680/Final/code/scenario_results")

table2.1 <- matrix(NA, nrow = 4, ncol = 9)
table3.1 <- matrix(NA, nrow = 4, ncol = 9)
table4.1 <- matrix(NA, nrow = 4, ncol = 9)
rownames(table2.1) <- c("nn","mn0.05","mn0.1","mnAD")
rownames(table3.1) <- c("nn","mn0.05","mn0.1","mnAD")
rownames(table4.1) <- c("nn","mn0.05","mn0.1","mnAD")
colnames(table2.1) <- paste("sc",paste(sc),sep = "")
colnames(table3.1) <- paste("sc",paste(sc),sep = "")
colnames(table4.1) <- paste("sc",paste(sc),sep = "")

# n-out-of-n bootstrap
sc <- seq(1,9)
di <- c("nn_psi1_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.1[1,i] <- 300 # resampling size always 300
  ul <- t(apply(dat, 1, percentile_nn))
  dat2 <- cbind(dat[,2], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_nn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  print(prop.test(sum(cov),1000,p=0.95)$p.value)
  table3.1[1,i] <- mean(cov)
  table4.1[1,i] <- mean((CI[,2]-CI[,1]))
}

# m-out-of-n bootstrap alpha=0.05
di <- c("mn0.05_psi1_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.1[2,i] <- mean(dat[,4]) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat[,2], dat[,4], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  table3.1[2,i] <- mean(cov)
  table4.1[2,i] <- mean((CI[,2]-CI[,1]))
}

# m-out-of-n bootstrap alpha=0.1
di <- c("mn0.1_psi1_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.1[3,i] <- mean(dat[,4]) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat[,2], dat[,4], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  table3.1[3,i] <- mean(cov)
  table4.1[3,i] <- mean((CI[,2]-CI[,1]))
}

# m-out-of-n bootstrap adaptive alpha
di <- c("mnad_psi1_scenario")

for(i in 1:9)
{
  name <- paste(di,paste(sc[i]),paste(".csv"),sep = "")
  dat <- read.csv(name, header = TRUE)
  table2.1[4,i] <- mean(dat[,4]) # average m
  ul <- t(apply(dat, 1, percentile_mn))
  dat2 <- cbind(dat[,2], dat[,4], ul[,1], ul[,2])
  CI <- t(apply(dat2, 1, consCI_mn))
  cov <- ifelse(CI[,1] < g[i, pos] & CI[,2] > g[i, pos], 1, 0)
  table3.1[4,i] <- mean(cov)
  table4.1[4,i] <- mean((CI[,2]-CI[,1]))
}
