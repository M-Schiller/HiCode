// File: main_hierarchy.cpp
// -- output community structure handling (number of levels, communities of one level)
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This program must not be distributed without agreement of the above mentionned authors.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

int display_level = -1;
char *filename = nullptr;

void print_usage(char *prog_name, const char *more)
{
  std::cerr << more;
  std::cerr << "print_usage: " << prog_name << " input_file [options]" << std::endl << std::endl;
  std::cerr << "input_file: read the community tree from this file." << std::endl;
  std::cerr << "-l xx\t display the community structure for the level xx." << std::endl;
  std::cerr << "\t outputs the community for each node." << std::endl;
  std::cerr << "\t xx must belong to [-1,N] if N is the number of levels." << std::endl;
  std::cerr << "-n\t displays the number of levels and the size of each level." << std::endl;
  std::cerr << "\t equivalent to -l -1." << std::endl;
  std::cerr << "-h\tshow this print_usage message." << std::endl;
}

void parse_args(int argc, char **argv)
{
  if (argc < 2)
  {
    print_usage(argv[0], "Bad arguments number\n");
  }

  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      switch (argv[i][1])
      {
      case 'l':
        display_level = std::atoi(argv[++i]);
        break;
      case 'n':
        display_level = -1;
        break;
      default:
        print_usage(argv[0], "Unknown option\n");
      }
    }
    else
    {
      if (filename == nullptr)
      {
        filename = argv[i];
      }
      else
      {
        print_usage(argv[0], "More than one filename\n");
      }
    }
  }
  if (filename == nullptr)
  {
    print_usage(argv[0], "No input file has been provided.\n");
  }
}
int main(int argc, char **argv) {
  parse_args(argc, argv);

  std::vector<std::vector<int>> levels;

  std::ifstream finput;
  finput.open(filename, std::fstream::in);

  int l = -1;
  while (!finput.eof())
  {
    int node, nodecomm;
    finput >> node >> nodecomm;

    if (finput)
    {
      if (node == 0)
      {
        l++;
        levels.resize(l + 1);
      }
      levels[l].push_back(nodecomm);
    }
  }

  if (display_level == -1)
  {
    std::cout << "Number of levels: " << levels.size() << std::endl;
    for (unsigned int i = 0; i < levels.size(); i++)
    {
      std::cout << "level " << i << ": " << levels[i].size() << " nodes" << std::endl;
    }
  }
  else if (display_level < 0 || (unsigned)display_level >= levels.size())
  {
    std::cerr << "Incorrect level\n";
  }
  else
  {
    std::vector<int> n2c(levels[0].size());

    for (unsigned int i = 0; i < levels[0].size(); ++i)
    {
      n2c[i] = i;
    }

    for (l = 0; l < display_level; l++)
    {
      for (unsigned int node = 0; node < levels[0].size(); node++)
      {
        n2c[node] = levels[l][n2c[node]];
      }
    }

    for (unsigned int node = 0; node < levels[0].size(); node++)
    {
      std::cout << node << " " << n2c[node] << std::endl;
    }
  }
}
