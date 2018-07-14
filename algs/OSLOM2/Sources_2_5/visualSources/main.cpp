/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *	This program is free software; you can redistribute it and/or modify         *
 *  it under the terms of the GNU General Public License as published by         *
 *  the Free Software Foundation; either version 2 of the License, or            *
 *  (at your option) any later version.                                          *
 *                                                                               *
 *  This program is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                *
 *  GNU General Public License for more details.                                 *
 *                                                                               *
 *  You should have received a copy of the GNU General Public License            *
 *  along with this program; if not, write to the Free Software                  *
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                               *
 *  Created by Andrea Lancichinetti on 5/5/09 (email: arg.lanci@gmail.com)       *
 *  Location: ISI foundation, Turin, Italy                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

 //./a.out -f network.dat

#include "visual_net.h"
#include "hier.h"
//#include "netvi.h"

int get_tp_files(std::string network_file)
{
  std::string s = network_file + "_oslo_files";
  std::string syst = "ls " + s + " > " + s + "_list";
  std::cout << syst << std::endl;

  char b[1000];
  cast_string_to_char(syst, b);
  systemCall(b);

  std::string f = s + "_list";
  cast_string_to_char(f, b);

  std::ifstream gin(b);
  int levels = 0;
  while (getline(gin, s))
  {
    if (!s.empty() && s[0] == 'n')
    {
      levels++;
    }
  }

  return levels;
}

int main(int argc, char * argv[])
{
  if (argc < 2)
  {
    std::cerr << argv[0] << " network_file" << std::endl;
    return -1;
  }

  srand_file();

  const int levels = get_tp_files(std::string(argv[1]));
  std::cout << "levels found: " << levels << std::endl;

  if (levels == 0)
  {
    std::cerr << "oslo directory not found..." << std::endl;
    return -1;
  }

  //return 0;
  all_levels(levels, std::string(argv[1]));

  //netgnu netsail;
  //netsail.set_networks(levels, string(argv[1]));
  //netsail.interface();

  return 0;
}
