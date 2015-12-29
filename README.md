## Introduction

This is my final project for MATH680: Computation Intensive Statistics. In this project, I reproduce the simulations in *Inference for Optimal Dynamic Treatment Regimes using an Adaptive m-out-of-n Bootstrap Scheme* from Chakraborty *et al* (2013) available [here](http://onlinelibrary.wiley.com/doi/10.1111/biom.12052/abstract).

## 1. Explore DTRreg package
I followed *Dynamic Treatment Regimen Estimation Via Regression-Based Techniques: Introducing R Package DTRreg* from Wallace *et al* (submitted) to explore the DTRreg package and its main function, **DTRreg**. This function is used to estimate an optimal DTR via three methods: dWOLS, G-estimation and Q-learning. I reproduced example 4.1 which demonstrated the use of dWOLS for estimating an optimal DTR in a two-stage simulated dataset. The code worked fine, but I couldn't recover the exact value in the summary output. I also reproduced the simulations shown in Table 1. These were conducted to assess the double robustness of the G-estimation and dWOLS methods, compared to Q-learning. My results were consistent with those of the paper. The R code can be found in DTRreg_repro.R

## 2. Code the simulations: first attempt
Following Table 1, I adapted the simulation parameter settings for dWOLS, where the treatment is coded as {0,1} instead of {-1,1}. The data generating process in the simulations is useful because it allows to influence the degree of nonregularity in the simulated data. The authors present 9 different scenarios classified as *regular*, *nonregular* and *near-nonregular*. In simulations.R, I coded scenarios 3 and 5. 

## 3. Reproduce Table 2, Table 3 and Table 4 from the paper
In the paper, 7 bootstrap methods are compared in each scenario. I only compared the regular n-out-of-n bootstrap method to the three m-out-of-n bootstrap introduced in the paper (fixed alpha = 0.05, fixed alpha = 0.1, adaptive choice of alpha). The four methods are coded in separate R files: nn_bootstrap.R, mn_bootstrap0.05.R, mn_bootstrap0.1.R and mn_bootstrapAD.R. The file choice_alpha.R contains the double bootstrap algorithm to select the tuning parameter alpha for each scenario. An important point here: for the first time, I used SSH to connect to the McGill server and run my codes on computers to which I have access. I downloaded [Cyberduck](https://cyberduck.io) in order to transfer files from my computer to the McGill computer. 

## 4. Results
The file table234.R reunit the results from the different scenarios/methods into tables, as in the paper. I played with the function `paste()` in R, as I never really took time to explore this useful function.