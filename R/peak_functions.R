#' Peak Functions for Peak Detection in Univariate Time Series
#'
#' scorepeak package provides several types of peak function. See the vignette "Introduction to scorepeak" for detail.
#'
#' @name peak_functions
#' @inheritParams detect_localmaxima
#' @return a numeric vector
#' @author Shota Ochi
#' @export
#' @examples
#' data("ecgca102")
#' plot(ecgca102, type = "l", ylim = c(-0.38, 0.53))
#' points(seq(length(ecgca102)), score_type1(ecgca102, 51), col = "red", type = "l")
#' points(seq(length(ecgca102)), score_type2(ecgca102, 51), col = "blue", type = "l")
#' points(seq(length(ecgca102)), score_type3(ecgca102, 51), col = "green", type = "l")
score_type1 <- function(data, w, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)
  
  data - (min_neighbors(data, w, "left", boundary) + min_neighbors(data, w, "right", boundary)) / 2
}

#' @rdname peak_functions
#' @export
score_type2 <- function(data, w, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)

  data - mean_neighbors(data, w, "both", boundary)
}

#' @rdname peak_functions
#' @export
score_type3 <- function(data, w, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)

  (data - pmax(mean_neighbors(data, w, "left", boundary), mean_neighbors(data, w, "right", boundary))) * sd_neighbors(data, w, "all", boundary)
}
