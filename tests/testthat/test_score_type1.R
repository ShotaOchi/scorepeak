test_that("score_type1",
{
  correct_s1_r <- c(255, 387, 439, 625)
  correct_s1_p <- c(255, 387, 439, 625)
  correct_s1_d <- c(255, 387, 439, 625)
  expect_equal(which(score_type1(a, w2) > thr & lp2), correct_s1_r)
  expect_equal(which(score_type1(a, w2, "periodic") > thr & lp2), correct_s1_p)
  expect_equal_na_allowed(which(score_type1(a, w2, "discard") > thr & lp2), correct_s1_d)
})
