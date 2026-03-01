---
title: "BIOS 731 – Homework 3"
author: "Linlin Wu"
date: "February 25, 2026"
output: pdf_document
---

# Algorithms for Logistic Regression

---

## Overview

This assignment focuses on implementing and comparing optimization algorithms for logistic regression, and deriving and applying an EM algorithm for censored exponential survival data.

The primary objectives are:

- Derive the likelihood, gradient, and Hessian for logistic regression  
- Implement Newton’s method  
- Implement a Minorize–Maximize (MM) algorithm  
- Compare Newton, MM, `glm()`, and `optim(method = "BFGS")`  
- Derive and implement an EM algorithm for censored exponential data  
- Compare EM results with an exponential AFT model  

All computations were performed in R using an R Project structure with relative file paths via the `here()` package.

---

## Project Structure

- `analysis/` – R Markdown report and knitted output  
- `source/` – Functions for Newton, MM, EM, and supporting utilities  
- `data/` – Simulated dataset  
- `results/` – Tables and figures generated from the analysis  

All results can be reproduced by knitting `analysis/hw3.Rmd`.

---

## Problem 1 – Newton’s Method

In this section:

- The logistic regression likelihood, score function, and Hessian were derived.  
- Newton’s method was implemented directly using the closed-form update.  
- Convergence was verified numerically.  
- The concavity of the log-likelihood was confirmed.  

The implementation is located in:

- `source/newton_logistic.R`  
- `source/logistic_utils.R`  

---

## Problem 2 – MM Algorithm

In this section:

- The required inequality for constructing a minorizing function was proven.  
- The arithmetic–geometric mean inequality was applied to separate parameters.  
- Coordinate-wise updating equations were derived.  
- The MM algorithm was implemented and compared with Newton’s method.  

The implementation is located in:

- `source/mm_logistic.R`  

---

## Problem 3 – Simulation Study

Simulation setup:

- Model: logit(P(Y = 1 | X)) = β₀ + β₁X  
- β₀ = 1  
- β₁ = 0.3  
- X ~ N(0,1)  
- n = 200  
- nsim = 1  

Methods compared:

- Newton’s method  
- MM algorithm  
- `glm()`  
- `optim()` with BFGS  

For each method, the following were reported:

- Parameter estimates  
- 95% confidence intervals  
- Computation time  
- Number of iterations  

Results are saved in the `results/` folder and summarized in tables and figures.

---

## Problem 4 – EM Algorithm for Censored Exponential Data

In this section:

- An EM algorithm was derived for estimating the exponential rate parameter under right censoring.  
- The E-step and M-step updates were implemented in R.  
- The algorithm was applied to the `veteran` dataset from the `survival` package.  
- Convergence was monitored via parameter change and log-likelihood increase.  
- A 95% confidence interval was constructed using observed information.  
- Results were compared with an exponential AFT model fitted using `survreg()`.  

The implementation is located in:

- `source/em_exp_censored.R`  
- `source/problem4_run.R`  

---

## Reproducibility

To reproduce the analysis:

1. Clone the repository  
2. Open `bios731_hw3_linlin.Rproj`  
3. Knit `analysis/hw3.Rmd`  

All file paths use the `here()` package to ensure reproducibility.
