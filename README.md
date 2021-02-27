# scorepeak

[![Build Status](https://github.com/ShotaOchi/scorepeak/workflows/R-CMD-check/badge.svg)](https://github.com/ShotaOchi/scorepeak/actions)
[![CRAN Version](https://www.r-pkg.org/badges/version/scorepeak)](https://cran.r-project.org/package=scorepeak)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/gh/ShotaOchi/scorepeak/branch/master/graph/badge.svg)](https://codecov.io/gh/ShotaOchi/scorepeak)

## Introduction
scorepeak is an R package for peak detection in univariate time series.

scorepeak provides peak functions, which enable us to detect peaks.

## Installation
you can install scorepeak from GitHub.

Run the following R code to install scorepeak.
```r
devtools::install_github("ShotaOchi/scorepeak")
```

## Simple Demo

A simple demo of peak detection is shown below.

```r
library(scorepeak)
data("ecgca102")
lp <- detect_localmaxima(ecgca102, 13)
score <- score_type1(ecgca102, 51)
detected <- score > 0.03 & lp
plot(ecgca102, type = "l")
points(which(detected), ecgca102[detected], pch = 19, col = "red")
```

## Contribution
You're welcome to create issues for any bug report or suggestion on the [issues page](https://github.com/ShotaOchi/scorepeak/issues).

You can also fork this repository and send me a pull request for bug fixes or additional features.