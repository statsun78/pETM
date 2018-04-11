# pETM
[R package] penalized Exponential Tilt Model

## Overview

Penalized Exponential Tilt Model: Fit a penalized exponential tilt model (ETM) to identify differentially methylated loci between cases and controls. ETM is able to detect the differences in means only, in variances only or in both means and variances. A penalized exponential tilt model using combined lasso and Laplacian penalties is applied to high-dimensional DNA methylation data from case-control association studies. When CpG sites are correlated with each others within the same gene or the same genetic region, Laplacian matrix can be imposed into the penalty function to encourage grouping effects among linked CpG sites. The  selection probability of an individual CpG site is computed based on the finite number of resamplings. 

## Installation

```
library(devtools)
install_github("statsun78/pETM")
```
