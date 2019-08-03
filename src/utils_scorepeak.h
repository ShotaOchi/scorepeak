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

#ifndef SCOREPEAK
#define SCOREPEAK

#include <Rcpp.h>

inline std::pair<int,double> max_temporal_neighbors(std::list< std::pair<int,double> >& window)
{
  std::pair<int,double> tmp_first = *(window.begin());
  std::pair<int,double> max_pair = std::make_pair(tmp_first.first, tmp_first.second);
  for (std::list< std::pair<int,double> >::iterator itr = window.begin(); itr != window.end(); ++itr)
  {
    std::pair<int,double> tmp_pair = *itr;
    if (tmp_pair.second > max_pair.second)
    {
      max_pair.first = tmp_pair.first;
      max_pair.second = tmp_pair.second;
    }
  }
  return max_pair;
}


inline std::pair<int,double> min_temporal_neighbors(std::list< std::pair<int,double> >& window)
{
  std::pair<int,double> tmp_first = *(window.begin());
  std::pair<int,double> min_pair = std::make_pair(tmp_first.first, tmp_first.second);
  for (std::list< std::pair<int,double> >::iterator itr = window.begin(); itr != window.end(); ++itr)
  {
    std::pair<int,double> tmp_pair = *itr;
    if (tmp_pair.second < min_pair.second)
    {
      min_pair.first = tmp_pair.first;
      min_pair.second = tmp_pair.second;
    }
  }
  return min_pair;
}


inline double sd_temporal_neighbors(double sum_window, double sum_squared_window, int size_window)
{
  return std::sqrt(sum_squared_window / size_window - sum_window * sum_window / ((double) size_window * size_window));
}

#endif // SCOREPEAK
