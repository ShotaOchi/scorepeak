test_that("score_type3",
{
  correct_s3_r <- c(387)
  correct_s3_p <- c(387)
  correct_s3_d <- c(387)
  expect_equal(which(score_type3(a, w2) > thr2 & lp2), correct_s3_r)
  expect_equal(which(score_type3(a, w2, "periodic") > thr2 & lp2), correct_s3_p)
  expect_equal_na_allowed(which(score_type3(a, w2, "discard") > thr2 & lp2), correct_s3_d)
})
