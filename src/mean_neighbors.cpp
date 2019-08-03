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
Rcpp::NumericVector mean_neighbors_left_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create reflecting left window
  for (int i = windowsize_half; i > 0; --i)
  {
    sum_window += data[i];
  }
  
  //' store mean of temporal left neighbors
  res[0] = sum_window / windowsize_half;
  
  // store meanimum of temporal left neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i - 1;
    int out_range = i - windowsize_half - 1;
    if (out_range < 0)
    {
      out_range = -out_range;
    }
    sum_window += data[max_range];
    sum_window -= data[out_range];
    res[i] = sum_window / windowsize_half;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_right_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create reflecting right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    sum_window += data[i];
  }
  
  //' store meanimum of temporal right neighbors
  res[0] = sum_window / windowsize_half;
  
  // store meanimum of temporal right neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    if (max_range >= n)
    {
      max_range = n - max_range + n - 2; // 2 * (n - 1) - max_range
    }
    int out_range = i;
    sum_window += data[max_range];
    sum_window -= data[out_range];
    res[i] = sum_window / windowsize_half;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_both_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int double_windowsize_half = windowsize_half + windowsize_half;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create reflecting left window
  for (int i = windowsize_half; i > 0; --i)
  {
    sum_window += data[i];
  }
  
  //' create reflecting right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    sum_window += data[i];
  }

  //' store meanimum of temporal neighbots
  res[0] = sum_window / double_windowsize_half;
  
  // store meanimum of temporal neighbors
  for (int i = 1; i < n; ++i)
  {
    //' left window
    int max_range_left = i - 1;
    int out_range_left = i - windowsize_half - 1;
    if (out_range_left < 0)
    {
      out_range_left = -out_range_left;
    }
    sum_window += data[max_range_left];
    sum_window -= data[out_range_left];
    //' right window
    int max_range_right = i + windowsize_half;
    int out_range_right = i;
    if (max_range_right >= n)
    {
      max_range_right = n - max_range_right + n - 2; // 2 * (n - 1) - max_range
    }
    sum_window += data[max_range_right];
    sum_window -= data[out_range_right];
    res[i] = sum_window / double_windowsize_half;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_all_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create reflecting window
  for (int i = windowsize_half; i > 0; --i)
  {
    sum_window += data[i];
  }
  for (int i = 0; i <= windowsize_half; ++i)
  {
    sum_window += data[i];
  }
  
  //' store local meanimum
  res[0] = sum_window / windowsize;
  
  //' store local meanimum
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
    sum_window += data[max_range];
    sum_window -= data[out_range];
    res[i] = sum_window / windowsize;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_left_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nmeanwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create periodic left window
  for (int i = nmeanwshalf; i < n; ++i)
  {
    sum_window += data[i];
  }
  
  //' store meanimum of temporal left neighbors
  res[0] = sum_window / windowsize_half;
  
  // store meanimum of temporal left neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i - 1;
    int out_range = i - windowsize_half - 1 + n;
    out_range %= n;
    sum_window += data[max_range];
    sum_window -= data[out_range];
    res[i] = sum_window / windowsize_half;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_right_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create periodic right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    sum_window += data[i];
  }
  
  //' store meanimum of temporal right neighbors
  res[0] = sum_window / windowsize_half;
  
  // store meanimum of temporal right neighbors
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    max_range %= n;
    int out_range = i;
    sum_window += data[max_range];
    sum_window -= data[out_range];
    res[i] = sum_window / windowsize_half;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_both_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int double_windowsize_half = windowsize_half + windowsize_half;
  const int n = data.size();
  const int nmeanwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create periodic left window
  for (int i = nmeanwshalf; i < n; ++i)
  {
    sum_window += data[i];
  }
  //' create periodic right window
  for (int i = 1; i <= windowsize_half; ++i)
  {
    sum_window += data[i];
  }
  
  //' store meanimum of temporal neighbots
  res[0] = sum_window / double_windowsize_half;
  
  // store meanimum of temporal neighbors
  for (int i = 1; i < n; ++i)
  {
    //' left window
    int max_range_left = i - 1;
    int out_range_left = i - windowsize_half - 1 + n;
    out_range_left %= n;
    sum_window += data[max_range_left];
    sum_window -= data[out_range_left];
    //' right window
    int max_range_right = i + windowsize_half;
    max_range_right %= n;
    int out_range_right = i;
    sum_window += data[max_range_right];
    sum_window -= data[out_range_right];
    res[i] = sum_window / double_windowsize_half;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_all_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nmeanwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n);
  double sum_window = 0.0;
  
  //' create periodic window
  for (int i = nmeanwshalf; i < n; ++i)
  {
    sum_window += data[i];
  }
  for (int i = 0; i <= windowsize_half; ++i)
  {
    sum_window += data[i];
  }
  
  //' store local meanimum
  res[0] = sum_window / windowsize;
  
  // store local meanima
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    max_range %= n;
    int out_range = i - windowsize_half - 1 + n;
    out_range %= n;
    sum_window += data[max_range];
    sum_window -= data[out_range];
    res[i] = sum_window / windowsize;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_left_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nmeanwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL);
  double sum_window = 0.0;
  
  if (n >= windowsize)
  {
    //' create left window
    for (int i = 0; i < windowsize_half; ++i)
    {
      sum_window += data[i];
    }
  
    //' store meanimum of temporal left neighbors
    res[windowsize_half] = sum_window / windowsize_half;
  
    // store meanimum of temporal left neighbors
    for (int i = windowsize_half + 1; i < nmeanwshalf; ++i)
    {
      int max_range = i - 1;
      int out_range = i - windowsize_half - 1;
      sum_window += data[max_range];
      sum_window -= data[out_range];
      res[i] = sum_window / windowsize_half;
    }
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_right_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nmeanwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL); //' data points in boundary are not available
  double sum_window = 0.0;
  
  if (n >= windowsize)
  {
    //' create right window
    for (int i = windowsize_half + 1; i < windowsize; ++i)
    {
      sum_window += data[i];
    }
  
    //' store meanimum of temporal right neighbors
    res[windowsize_half] = sum_window / windowsize_half;
  
    // store meanimum of temporal right neighbors
    for (int i = windowsize_half + 1; i < nmeanwshalf; ++i)
    {
      int max_range = i + windowsize_half;
      int out_range = i;
      sum_window += data[max_range];
      sum_window -= data[out_range];
      res[i] = sum_window / windowsize_half;
    }
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_both_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int double_windowsize_half = windowsize_half + windowsize_half;
  const int n = data.size();
  const int nmeanwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL);
  double sum_window = 0.0;
  
  if (n >= windowsize)
  {
    //' create left window
    for (int i = 0; i < windowsize_half; ++i)
    {
      sum_window += data[i];
    }
  
    //' create right window
    for (int i = windowsize_half + 1; i < windowsize; ++i)
    {
      sum_window += data[i];
    }

    //' store meanimum of temporal neighbots
    res[windowsize_half] = sum_window / double_windowsize_half;
    
    // store meanimum of temporal neighbors
    for (int i = windowsize_half + 1; i < nmeanwshalf; ++i)
    {
      //' left window
      int max_range_left = i - 1;
      int out_range_left = i - windowsize_half - 1;
      sum_window += data[max_range_left];
      sum_window -= data[out_range_left];
      //' right window
      int max_range_right = i + windowsize_half;
      int out_range_right = i;
      sum_window += data[max_range_right];
      sum_window -= data[out_range_right];
      res[i] = sum_window / double_windowsize_half;
    }
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::NumericVector mean_neighbors_all_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  const int nmeanwshalf = n - windowsize_half;
  Rcpp::NumericVector res(n, NA_REAL);
  double sum_window = 0.0;
  
  if (n >= windowsize)
  {
    //' create discard window
    for (int i = 0; i < windowsize; ++i)
    {
      sum_window += data[i];
    }
  
    //' store local meanimum
    res[windowsize_half] = sum_window / windowsize;
  
    // check if elements are local meanima
    for (int i = windowsize_half + 1; i < nmeanwshalf; ++i)
    {
      int max_range = i + windowsize_half;
      int out_range = i - windowsize_half - 1;
      sum_window += data[max_range];
      sum_window -= data[out_range];
      res[i] = sum_window / windowsize;
    }
  }
  return res;
}
