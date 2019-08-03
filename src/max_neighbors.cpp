/*
 * Copyright (c) 2019, Shota Ochi <shotaochi1990@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "utils_scorepeak.h"

/* windowsize must be smaller than double length of data */
/* length of data must be greater than 1 */

// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_left_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(windowsize_half, data[windowsize_half]);
  
  //' create reflecting left window
  for (int i = windowsize_half; i > 0; --i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  
  //' store maximum of temporal left neighbors
  res[0] = max_pair.second;
  
  // store maximum of temporal left neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i - 1;
    int out_range = i - windowsize_half - 1;
    if (out_range < 0)
    {
      out_range = -out_range;
    }
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first == out_range)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_right_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(windowsize_half, data[windowsize_half]);
  
  //' create reflecting right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  
  //' store maximum of temporal right neighbors
  res[0] = max_pair.second;
  
  // store maximum of temporal right neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    if (max_range >= n)
    {
      max_range = n - max_range + n - 2; // 2 * (n - 1) - max_range
    }
    int out_range = i;
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first == out_range)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_both_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window_left;
  std::list< std::pair<int,double> > window_right;
  std::pair<int,double> max_pair_left = std::make_pair(windowsize_half, data[windowsize_half]);
  std::pair<int,double> max_pair_right = std::make_pair(windowsize_half, data[windowsize_half]);
  
  //' create reflecting left window
  for (int i = windowsize_half; i > 0; --i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window_left.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair_left.second)
    {
      max_pair_left.first = i;
      max_pair_left.second = data[i];
    }
  }
  
  //' create reflecting right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window_right.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair_right.second)
    {
      max_pair_right.first = i;
      max_pair_right.second = data[i];
    }
  }

  //' store maximum of temporal neighbots
  res[0] = max_pair_left.second; // maximum of left temporal neighbors and maximum of right temporal neighbors are same under this condition
  
  // store maximum of temporal neighbors
  for (int i = 1; i < n; ++i)
  {
    //' left window
    int max_range_left = i - 1;
    int out_range_left = i - windowsize_half - 1;
    if (out_range_left < 0)
    {
      out_range_left = -out_range_left;
    }
    window_left.pop_front();
    std::pair<int,double> tmp_pair_left = std::make_pair(max_range_left,data[max_range_left]);
    window_left.push_back(tmp_pair_left);
    if (data[max_range_left] > max_pair_left.second)
    {
      max_pair_left.first = max_range_left;
      max_pair_left.second = data[max_range_left];
    }
    if (max_pair_left.first == out_range_left)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window_left);
      max_pair_left.first = tmp_max.first;
      max_pair_left.second = tmp_max.second;
    }
    //' right window
    int max_range_right = i + windowsize_half;
    int out_range_right = i;
    if (max_range_right >= n)
    {
      max_range_right = n - max_range_right + n - 2; // 2 * (n - 1) - max_range
    }
    window_right.pop_front();
    std::pair<int,double> tmp_pair_right = std::make_pair(max_range_right,data[max_range_right]);
    window_right.push_back(tmp_pair_right);
    if (data[max_range_right] > max_pair_right.second)
    {
      max_pair_right.first = max_range_right;
      max_pair_right.second = data[max_range_right];
    }
    if (max_pair_right.first == out_range_right)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window_right);
      max_pair_right.first = tmp_max.first;
      max_pair_right.second = tmp_max.second;
    }
    res[i] = max_pair_left.second >= max_pair_right.second ? max_pair_left.second : max_pair_right.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_all_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(0, data[0]);
  
  //' create reflecting window
  for (int i = windowsize_half; i > 0; --i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  for (int i = 0; i <= windowsize_half; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  
  //' store local maximum
  res[0] = max_pair.second;
  
  //' store local maximum
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    if (max_range >= n)
    {
      max_range = n - max_range + n - 2; // 2 * (n - 1) - max_range
    }
    int out_range = i - windowsize_half - 1;
    if (out_range < 0)
    {
      out_range = -out_range;
    }
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first == out_range)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_left_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nminwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(windowsize_half, data[windowsize_half]);
  
  //' create periodic left window
  for (int i = nminwshalf; i < n; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  
  //' store maximum of temporal left neighbors
  res[0] = max_pair.second;
  
  // store maximum of temporal left neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i - 1;
    int out_range = i - windowsize_half - 1 + n;
    out_range %= n;
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first == out_range)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_right_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(windowsize_half, data[windowsize_half]);
  
  //' create periodic right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  
  //' store maximum of temporal right neighbors
  res[0] = max_pair.second;
  
  // store maximum of temporal right neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    max_range %= n;
    int out_range = i;
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first == out_range)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_both_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nminwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window_left;
  std::list< std::pair<int,double> > window_right;
  std::pair<int,double> max_pair_left = std::make_pair(nminwshalf, data[nminwshalf]);
  std::pair<int,double> max_pair_right = std::make_pair(1, data[1]);
  
  //' create periodic left window
  for (int i = nminwshalf; i < n; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window_left.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair_left.second)
    {
      max_pair_left.first = i;
      max_pair_left.second = data[i];
    }
  }
  //' create periodic right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window_right.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair_right.second)
    {
      max_pair_right.first = i;
      max_pair_right.second = data[i];
    }
  }
  
  //' store maximum of temporal neighbots
  res[0] = max_pair_left.second > max_pair_right.second ? max_pair_left.second : max_pair_right.second;
  
  // store maximum of temporal neighbors
  for (int i = 1; i < n; ++i)
  {
    //' left window
    int max_range_left = i - 1;
    int out_range_left = i - windowsize_half - 1 + n;
    out_range_left %= n;
    window_left.pop_front();
    std::pair<int,double> tmp_pair_left = std::make_pair(max_range_left,data[max_range_left]);
    window_left.push_back(tmp_pair_left);
    if (data[max_range_left] > max_pair_left.second)
    {
      max_pair_left.first = max_range_left;
      max_pair_left.second = data[max_range_left];
    }
    if (max_pair_left.first == out_range_left)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window_left);
      max_pair_left.first = tmp_max.first;
      max_pair_left.second = tmp_max.second;
    }
    //' right window
    int max_range_right = i + windowsize_half;
    max_range_right %= n;
    int out_range_right = i;
    window_right.pop_front();
    std::pair<int,double> tmp_pair_right = std::make_pair(max_range_right,data[max_range_right]);
    window_right.push_back(tmp_pair_right);
    if (data[max_range_right] > max_pair_right.second)
    {
      max_pair_right.first = max_range_right;
      max_pair_right.second = data[max_range_right];
    }
    if (max_pair_right.first == out_range_right)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window_right);
      max_pair_right.first = tmp_max.first;
      max_pair_right.second = tmp_max.second;
    }
    res[i] = max_pair_left.second >= max_pair_right.second ? max_pair_left.second : max_pair_right.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_all_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nminwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(0, data[0]);
  
  //' create periodic window
  for (int i = nminwshalf; i < n; ++i)
  {
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
        max_pair.first = i;
        max_pair.second = data[i];
    }
  }
  for (int i = 0; i <= windowsize_half; ++i)
  {
    if (i < n)
    {
      std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
      window.push_back(tmp_pair);
      //' check if data[i] is maximum in window
      if (data[i] > max_pair.second)
      {
        max_pair.first = i;
        max_pair.second = data[i];
      }
    }
  }
  
  //' store local maximum
  res[0] = max_pair.second;
  
  // store local maxima
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    max_range %= n;
    int out_range = i - windowsize_half - 1 + n;
    out_range %= n;
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first == out_range)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.second;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_left_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nminwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(0, data[0]);
  
  if (n >= windowsize)
  {
    //' create right window
    for (int i = 0; i < windowsize_half; ++i)
    {
      std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
      window.push_back(tmp_pair);
      //' check if data[i] is maximum in window
      if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  
  //' store maximum of temporal right neighbors
  res[windowsize_half] = max_pair.second;
  
  // store maximum of temporal right neighbors
  for (int i = windowsize_half + 1; i < nminwshalf; ++i)
  {
    int max_range = i - 1;
    int out_range = i - windowsize_half - 1;
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first == out_range)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.second;
  }
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_right_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nminwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL); //' data points in boundary are not available
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair;
  
  if (n >= windowsize)
  {
    //' create right window
    bool is_first = true;
    for (int i = windowsize_half + 1; i < windowsize; ++i)
    {
      std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
      window.push_back(tmp_pair);
      //' check if data[i] is maximum in window
      if (is_first)
      {
        max_pair.first = i;
        max_pair.second = data[i];
        is_first = false;
      } else 
      {
        if (data[i] > max_pair.second)
        {
          max_pair.first = i;
          max_pair.second = data[i];
        }
      }
    }
  
    //' store maximum of temporal right neighbors
    res[windowsize_half] = max_pair.second;
  
    // store maximum of temporal right neighbors
    for (int i = windowsize_half + 1; i < nminwshalf; ++i)
    {
      int max_range = i + windowsize_half;
      int out_range = i;
      window.pop_front();
      std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
      window.push_back(tmp_pair);
      if (data[max_range] > max_pair.second)
      {
        max_pair.first = max_range;
        max_pair.second = data[max_range];
      }
      if (max_pair.first == out_range)
      {
        std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
        max_pair.first = tmp_max.first;
        max_pair.second = tmp_max.second;
      }
      res[i] = max_pair.second;
    }
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_both_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nminwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL);
  std::list< std::pair<int,double> > window_left;
  std::list< std::pair<int,double> > window_right;
  std::pair<int,double> max_pair_left = std::make_pair(0, data[0]);
  std::pair<int,double> max_pair_right;
  
  if (n >= windowsize)
  {
    //' create left window
    for (int i = 0; i < windowsize_half; ++i)
    {
      std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
      window_left.push_back(tmp_pair);
      //' check if data[i] is maximum in window
      if (data[i] > max_pair_left.second)
      {
        max_pair_left.first = i;
        max_pair_left.second = data[i];
      }
    }
  
    //' create right window
    bool is_first = true;
    for (int i = windowsize_half + 1; i < windowsize; ++i)
    {
      std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
      window_right.push_back(tmp_pair);
      //' check if data[i] is maximum in window
      if (is_first)
      {
        max_pair_right.first = i;
        max_pair_right.second = data[i];
        is_first = false;
      } else 
      {
        if (data[i] > max_pair_right.second)
        {
          max_pair_right.first = i;
          max_pair_right.second = data[i];
        }
      }
    }

    //' store maximum of temporal neighbots
    res[windowsize_half] = max_pair_left.second > max_pair_right.second ? max_pair_left.second : max_pair_right.second;
    
    // store maximum of temporal neighbors
    for (int i = windowsize_half + 1; i < nminwshalf; ++i)
    {
      //' left window
      int max_range_left = i - 1;
      int out_range_left = i - windowsize_half - 1;
      window_left.pop_front();
      std::pair<int,double> tmp_pair_left = std::make_pair(max_range_left,data[max_range_left]);
      window_left.push_back(tmp_pair_left);
      if (data[max_range_left] > max_pair_left.second)
      {
        max_pair_left.first = max_range_left;
        max_pair_left.second = data[max_range_left];
      }
      if (max_pair_left.first == out_range_left)
      {
        std::pair<int,double> tmp_max =  max_temporal_neighbors(window_left);
        max_pair_left.first = tmp_max.first;
        max_pair_left.second = tmp_max.second;
      }
      //' right window
      int max_range_right = i + windowsize_half;
      int out_range_right = i;
      window_right.pop_front();
      std::pair<int,double> tmp_pair_right = std::make_pair(max_range_right,data[max_range_right]);
      window_right.push_back(tmp_pair_right);
      if (data[max_range_right] > max_pair_right.second)
      {
        max_pair_right.first = max_range_right;
        max_pair_right.second = data[max_range_right];
      }
      if (max_pair_right.first == out_range_right)
      {
        std::pair<int,double> tmp_max =  max_temporal_neighbors(window_right);
        max_pair_right.first = tmp_max.first;
        max_pair_right.second = tmp_max.second;
      }
      res[i] = max_pair_left.second >= max_pair_right.second ? max_pair_left.second : max_pair_right.second;
    }
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector max_neighbors_all_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nminwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL);
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(0, data[0]);
  
  if (n >= windowsize)
  {
    //' create discard window
    for (int i = 0; i < windowsize; ++i)
    {
      std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
      window.push_back(tmp_pair);
      //' check if data[i] is maximum in window
      if (data[i] > max_pair.second)
      {
        max_pair.first = i;
        max_pair.second = data[i];
      }
    }
  
    //' store local maximum
    res[windowsize_half] = max_pair.second;
  
    // check if elements are local maxima
    for (int i = windowsize_half + 1; i < nminwshalf; ++i)
    {
      int max_range = i + windowsize_half;
      int out_range = i - windowsize_half - 1;
      window.pop_front();
      std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
      window.push_back(tmp_pair);
      if (data[max_range] > max_pair.second)
      {
        max_pair.first = max_range;
        max_pair.second = data[max_range];
      }
      if (max_pair.first == out_range)
      {
        std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
        max_pair.first = tmp_max.first;
        max_pair.second = tmp_max.second;
      }
      res[i] = max_pair.second;
    }
  }
  return res;
}
