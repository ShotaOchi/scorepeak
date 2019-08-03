#' Building Blocks of Peak Functions
#'
#' @name building_blocks
#' @inheritParams detect_localmaxima
#' @param side determines which side of neighbors of data point will be used in calculation. "left", "l": left temporal neighbors, "right", "r": right temporal neighbors, "both", "b": left and right temporal neighbors, "all", "a": data point and its left and right temporal neighbors.
#' @return a numeric vector
#' @author Shota Ochi
#' @export
#' @examples
#' data("ecgca102")
#' max_neighbors(ecgca102, 3, "all")
#' min_neighbors(ecgca102, 3, "all")
#' mean_neighbors(ecgca102, 3, "all")
#' sd_neighbors(ecgca102, 3, "all")
max_neighbors <- function(data, w, side, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)
  assert_side(side)
  
  if (side == "left" || side == "l")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      max_neighbors_left_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      max_neighbors_left_periodic_cpp(data, w)
    } else
    {
      max_neighbors_left_discard_cpp(data, w)
    }
  } else if (side == "right" || side == "r")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      max_neighbors_right_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      max_neighbors_right_periodic_cpp(data, w)
    } else
    {
      max_neighbors_right_discard_cpp(data, w)
    }
  } else  if (side == "both" || side == "b")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      max_neighbors_both_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      max_neighbors_both_periodic_cpp(data, w)
    } else
    {
      max_neighbors_both_discard_cpp(data, w)
    }
  } else
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      max_neighbors_all_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      max_neighbors_all_periodic_cpp(data, w)
    } else
    {
      max_neighbors_all_discard_cpp(data, w)
    }
  }
}

#' @rdname building_blocks
#' @export
min_neighbors <- function(data, w, side, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)
  assert_side(side)
  
  if (side == "left" || side == "l")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      min_neighbors_left_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      min_neighbors_left_periodic_cpp(data, w)
    } else
    {
      min_neighbors_left_discard_cpp(data, w)
    }
  } else if (side == "right" || side == "r")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      min_neighbors_right_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      min_neighbors_right_periodic_cpp(data, w)
    } else
    {
      min_neighbors_right_discard_cpp(data, w)
    }
  } else  if (side == "both" || side == "b")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      min_neighbors_both_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      min_neighbors_both_periodic_cpp(data, w)
    } else
    {
      min_neighbors_both_discard_cpp(data, w)
    }
  } else
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      min_neighbors_all_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      min_neighbors_all_periodic_cpp(data, w)
    } else
    {
      min_neighbors_all_discard_cpp(data, w)
    }
  }
}

#' @rdname building_blocks
#' @export
mean_neighbors <- function(data, w, side, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)
  assert_side(side)
  
  if (side == "left" || side == "l")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      mean_neighbors_left_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      mean_neighbors_left_periodic_cpp(data, w)
    } else
    {
      mean_neighbors_left_discard_cpp(data, w)
    }
  } else if (side == "right" || side == "r")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      mean_neighbors_right_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      mean_neighbors_right_periodic_cpp(data, w)
    } else
    {
      mean_neighbors_right_discard_cpp(data, w)
    }
  } else  if (side == "both" || side == "b")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      mean_neighbors_both_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      mean_neighbors_both_periodic_cpp(data, w)
    } else
    {
      mean_neighbors_both_discard_cpp(data, w)
    }
  } else
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      mean_neighbors_all_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      mean_neighbors_all_periodic_cpp(data, w)
    } else
    {
      mean_neighbors_all_discard_cpp(data, w)
    }
  }
}

#' @rdname building_blocks
#' @export
sd_neighbors <- function(data, w, side, boundary = "reflecting")
{
  assert_data(data)
  w <- assert_window(w)
  w <- assert_length_window(w, data)
  assert_boundary(boundary)
  assert_side(side)
  
  if (side == "left" || side == "l")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      sd_neighbors_left_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      sd_neighbors_left_periodic_cpp(data, w)
    } else
    {
      sd_neighbors_left_discard_cpp(data, w)
    }
  } else if (side == "right" || side == "r")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      sd_neighbors_right_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      sd_neighbors_right_periodic_cpp(data, w)
    } else
    {
      sd_neighbors_right_discard_cpp(data, w)
    }
  } else  if (side == "both" || side == "b")
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      sd_neighbors_both_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      sd_neighbors_both_periodic_cpp(data, w)
    } else
    {
      sd_neighbors_both_discard_cpp(data, w)
    }
  } else
  {
    if (boundary == "reflecting" || boundary == "r")
    {
      sd_neighbors_all_reflecting_cpp(data, w)
    } else if (boundary == "periodic" || boundary == "p")
    {
      sd_neighbors_all_periodic_cpp(data, w)
    } else
    {
      sd_neighbors_all_discard_cpp(data, w)
    }
  }
}