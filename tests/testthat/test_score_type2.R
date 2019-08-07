test_that("score_type2",
{
  correct_s2_r <- c(255, 387, 439)
  correct_s2_p <- c(255, 387, 439)
  correct_s2_d <- c(255, 387, 439)
  expect_equal(which(score_type2(a, w2) > thr & lp2), correct_s2_r)
  expect_equal(which(score_type2(a, w2, "periodic") > thr & lp2), correct_s2_p)
  expect_equal_na_allowed(which(score_type2(a, w2, "discard") > thr & lp2), correct_s2_d)
})
