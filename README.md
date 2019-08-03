# scorepeak

## Introduction
scorepeak is an R package for detecing peaks in univariate time series.

scorepeak provides peak functions, which enable us to detect peaks.

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