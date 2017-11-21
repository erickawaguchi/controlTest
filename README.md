# controlTest
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/controlTest)](https://CRAN.R-project.org/package=controlTest)

## Introduction
controlTest is an R package for nonparametric comparisons of survival quantiles.

## Features
- ```quantileSE```: Find the standard error for a survival quantile via bootstrap.
- ```quantileControlTest```: Compare the survival quantiles of two independent samples via bootstrap.
- ```supControlTest```: (To be added in a later version)

## Examples
```
library(controlTest)
t1 <- c(1, 63, 105, 129, 182, 216, 250, 262, 301, 301,
       342, 354, 356, 358, 380, 383, 383, 338, 394, 408, 460, 489,
       499, 523, 524, 535, 562, 569, 675, 676, 748, 778, 786, 797,
       955, 968, 1000, 1245, 1271, 1420, 1551, 1694, 2363, 2754, 2950)
t2 <- c(17, 42, 44, 48, 60, 72, 74, 95, 103, 108, 122, 144, 167, 170,
       183, 185, 193, 195, 197, 208, 234, 235, 254, 307, 315, 401, 445,
       464, 484, 528, 542, 547, 577, 580, 795, 855, 1366, 1577, 2060,
       2412, 2486, 2796, 2802, 2934, 2988)
c1 <- c(rep(1, 43), 0, 0)
c2 <- c(rep(1, 39), rep(0, 6))
quantileControlTest(t1, c1, t2, c2, q = 0.5, B = 500)
```

