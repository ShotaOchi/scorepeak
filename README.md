# scorepeak

[![Build Status](https://travis-ci.org/ShotaOchi/scorepeak.svg?branch=master)](https://travis-ci.org/ShotaOchi/scorepeak)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ShotaOchi/scorepeak?branch=master&svg=true)](https://ci.appveyor.com/project/ShotaOchi/scorepeak)
[![CRAN Version](https://www.r-pkg.org/badges/version/scorepeak)](https://cran.r-project.org/package=scorepeak)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/gh/ShotaOchi/scorepeak/branch/master/graph/badge.svg)](https://codecov.io/gh/ShotaOchi/scorepeak)

## Introduction
scorepeak is an R package for peak detection in univariate time series.

scorepeak provides peak functions, which enable us to detect peaks .

## Installation
you can install scorepeak from GitHub.

Run the following R code to install scorepeak.
```r
devtools::install_github("ShotaOchi/scorepeak")
```

## Simple Demo

```r
library(scorepeak)
data("ecgca102")
a <- detect_localmaxima(ecgca102)
plot(ecgca102, type = "l")
points(which(a), ecgca102[a], pch = 1, col = "red")
```


## Contribution
You're welcome to create issues for any bug report or suggestion on the [issues page](https://github.com/ShotaOchi/scorepeak/issues).

You can also fork this repository and send me a pull request for bug fixes or additional features.