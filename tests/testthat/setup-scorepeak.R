data("ecgca102")
a <- ecgca102
b <- c(2, rep(1,4), 2:5, 4:1, rep(1,4), 3)
w <- 5
idx_lp <- which(detect_localmaxima(b, w))

expect_equal_na_allowed <- function(x, y)
{
  expect_equal(x[!is.na(x)], y[!is.na(y)])
}
