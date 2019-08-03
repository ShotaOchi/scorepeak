assert_data <- function(data)
{
  assertNumeric(data, finite = TRUE, any.missing = FALSE, min.len = 2)
}

assert_window <- function(w)
{
  assertNumeric(w, lower = 3, finite = TRUE, any.missing = FALSE, len = 1)
  w <- as.integer(w)
  if (w %% 2 == 0)
  {
    warning(sprintf("w must be odd. w is %d. w will be treated as %d.", w, w + 1))
    w <- w + 1
  }
  w
}

assert_boundary <- function(boundary)
{
  assertChoice(boundary, c("reflecting", "periodic", "discard", "r", "p", "d"))
}

assert_side <- function(side)
{
  assertChoice(side, c("left", "right", "both", "all", "l", "r", "b", "a"))
}

assert_length_window <- function(w, data)
{
  if (w > 2 * length(data))
  {
    warning(sprintf("w must be smaller than double length of data. w is %d. w will be treated as %d.", w, 2 * length(data) - 1))
    w <- as.integer(2 * length(data) - 1)
  }
  w
}
