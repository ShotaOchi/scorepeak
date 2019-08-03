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
Rcpp::LogicalVector detect_localmaxima_reflecting_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::LogicalVector res(n);
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
  
  //' check if first element is a local maximum
  res[0] = max_pair.first == 0 ? true : false;
  
  // check if elements are local maxima
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    if (max_range >= n)
    {
      max_range = n + n - 2 - max_range; // n - 1 - (max_range - (n - 1))
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
    res[i] = max_pair.first == i ? true : false;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::LogicalVector detect_localmaxima_periodic_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::LogicalVector res(n);
  std::list< std::pair<int,double> > window;
  const int nminwshalf = n - windowsize_half;
  std::pair<int,double> max_pair = std::make_pair(nminwshalf, data[nminwshalf]);
  
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
    std::pair<int,double> tmp_pair = std::make_pair(i,data[i]);
    window.push_back(tmp_pair);
    //' check if data[i] is maximum in window
    if (data[i] > max_pair.second)
    {
      max_pair.first = i;
      max_pair.second = data[i];
    }
  }
  
  //' check if first element is a local maximum
  res[0] = max_pair.first == 0 ? true : false;
  
  // check if elements are local maxima
  for (int i = 1; i < n; ++i)
  {
    int max_range = i + windowsize_half;
    if (i >= nminwshalf)
    {
      max_range %= n;
    }
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
    res[i] = max_pair.first == i ? true : false;
  }
  return res;
}


// [[Rcpp::export]]
Rcpp::LogicalVector detect_localmaxima_discard_cpp(const Rcpp::NumericVector& data, const int& windowsize)
{
  const int windowsize_half = windowsize / 2;
  const int n = data.size();
  Rcpp::LogicalVector res(n);
  const int nminwshalf = n - windowsize_half;
  std::list< std::pair<int,double> > window;
  std::pair<int,double> max_pair = std::make_pair(0, data[0]);
  
  //' data points in boundary are not candidates of loal peak
  for (int i = 0; i < windowsize_half; ++i)
  {
    res[i] = false;
  }
  for (int i = nminwshalf; i < n; ++i)
  {
    res[i] = false;
  }
  
  //' create window
  for (int i = 0; i < windowsize; ++i)
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
  
  //' check if first element is a local maximum
  res[windowsize_half] = max_pair.first == windowsize_half ? true : false;
  
  // check if elements are local maxima
  for (int i = windowsize_half + 1; i < nminwshalf; ++i)
  {
    int max_range = i + windowsize_half;
    window.pop_front();
    std::pair<int,double> tmp_pair = std::make_pair(max_range,data[max_range]);
    window.push_back(tmp_pair);
    if (data[max_range] > max_pair.second)
    {
      max_pair.first = max_range;
      max_pair.second = data[max_range];
    }
    if (max_pair.first < i - windowsize_half)
    {
      std::pair<int,double> tmp_max =  max_temporal_neighbors(window);
      max_pair.first = tmp_max.first;
      max_pair.second = tmp_max.second;
    }
    res[i] = max_pair.first == i ? true : false;
  }
  return res;
}
