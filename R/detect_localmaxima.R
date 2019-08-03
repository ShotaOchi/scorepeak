#' detect local maxima in univariate time series data
#'
#' @param data a numeric vector. Length of data must be greater than 1.
#' @param w window size. w must be odd and greater than 2 and smaller than double length of data.
#' @param boundary determines how data points in the beginning and end of the time series will be treated. "reflecting", "r": reflecting boundary condition, "periodic", "p": periodic boundary condition, "discard", "d", discarding data points in the beginning and end of the time series. See the vignette "Introduction to scorepeak" for detail.
#' @return a logical vector. TRUE indicates local peak. FALSE indicates not local peak.
#' @author Shota Ochi
#' @export
#' @examples
#' data("ecgca102")
#' peaks <- detect_localmaxima(ecgca102)
#' plot(ecgca102, type = "l")
#' points(which(peaks), ecgca102[peaks], pch = 1, col = "red")
detect_localmaxima <- function(data, w = 3, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)
  
  if (boundary == "reflecting" || boundary == "r")
  {
    detect_localmaxima_reflecting_cpp(data, w)
  } else if (boundary == "periodic" || boundary == "p")
  {
    detect_localmaxima_periodic_cpp(data, w)
  } else
  {
    detect_localmaxima_discard_cpp(data, w)
  }
}
