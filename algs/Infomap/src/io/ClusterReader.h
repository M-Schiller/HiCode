/**********************************************************************************

 Infomap software package for multi-level network clustering

 Copyright (c) 2013 Daniel Edler, Martin Rosvall

 For more information, see <http://www.mapequation.org>

 This file is part of Infomap software package.

 Infomap software package is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Infomap software package is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with Infomap software package.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/
#pragma once

#ifndef CLUSTERREADER_H_
#define CLUSTERREADER_H_

#include <vector>
#include <string>
using std::string;

class ClusterReader
{
public:
  ClusterReader(unsigned int numNodes);
  ~ClusterReader();

  void readData(std::string filename);

  const std::vector<unsigned int>& getClusterData() const
  {
    return m_clusterData;
  }

private:
  std::vector<unsigned int> m_clusterData;
};

#endif /* CLUSTERREADER_H_ */
