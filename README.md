## Introduction

This is my final project for MATH680: Computation Intensive Statistics. 

## 1. Explore DTRreg package
I followed *Dynamic Treatment Regimen Estimation Via Regression-Based Techniques: Introducing R Package DTRreg* from Wallace *et al* (submitted) to explore the DTRreg package and its main function, **DTRreg**. This function is used to estimate an optimal DTR via three methods: dWOLS, G-estimation and Q-learning. I reproduced example 4.1 which demonstrated the use of dWOLS for estimating an optimal DTR in a two-stage simulated dataset. The code worked fine, but I couldn't recover the exact value in the summary output. I also reproduced the simulations shown in Table 1. These were conducted to assess the double robustness of the G-estimation and dWOLS methods, compared to Q-learning. My results were consistent with those of the paper. R code can be found in DTRreg_repro.R
