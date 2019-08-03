test_that("detect_localmaxima",
{
  correct_reflecting <- c(1, 9, 18)
  correct_periodic <- c(9, 18)
  correct_discard <- c(9)
  w_toolarge <- 2000
  w_limit <- 2 * length(a) - 1
  w_even <- 200
  expect_equal(which(detect_localmaxima(b, w, "reflecting")), correct_reflecting)
  expect_equal(which(detect_localmaxima(b, w, "periodic")), correct_periodic)
  expect_equal(which(detect_localmaxima(b, w, "discard")), correct_discard)
  expect_equal(detect_localmaxima(a, w, "reflecting"), detect_localmaxima(a, w, "r"))
  expect_equal(detect_localmaxima(a, w, "periodic"), detect_localmaxima(a, w, "p"))
  expect_equal(detect_localmaxima(a, w, "discard"), detect_localmaxima(a, w, "d"))
  expect_equal(detect_localmaxima(c(1,2)), c(FALSE, TRUE))
  expect_warning(detect_localmaxima(a, w_toolarge))
  expect_equal(detect_localmaxima(a, w_limit), detect_localmaxima(a, w_toolarge))
  expect_warning(detect_localmaxima(a, w_even))
  expect_error(detect_localmaxima(1))
  expect_error(detect_localmaxima(c(1, NaN)))
  expect_error(detect_localmaxima(a, NaN))
  expect_error(detect_localmaxima(a, rep(1, 100)))
  expect_error(detect_localmaxima(a, -3))
  expect_error(detect_localmaxima(a, boundary = NaN))
})