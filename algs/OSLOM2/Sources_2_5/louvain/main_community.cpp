// File: main_community.cpp
// -- community detection, sample main file
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
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <unistd.h>

#include "graph_binary.h"
#include "community.h"

char *filename = nullptr;
char *filename_w = nullptr;
char *filename_part = nullptr;
int type = UNWEIGHTED;
int nb_pass = 0;
double precision = 0.000001;
int display_level = -2;
int k1 = 16;

bool verbose = false;

void print_usage(char *prog_name, const char *more) {
  std::cerr << more;
  std::cerr << "print_usage: " << prog_name << " input_file [-w weight_file] [-p part_file] [-q epsilon] [-l display_level] [-v] [-h]" << std::endl << std::endl;
  std::cerr << "input_file: file containing the graph to decompose in communities." << std::endl;
  std::cerr << "-w file\tread the graph as a weighted one (weights are set to 1 otherwise)." << std::endl;
  std::cerr << "-p file\tstart the computation with a given partition instead of the trivial partition." << std::endl;
  std::cerr << "\tfile must contain lines \"node community\"." << std::endl;
  std::cerr << "-q eps\ta given pass stops when the modularity is increased by less than epsilon." << std::endl;
  std::cerr << "-l k\tdisplays the graph of level k rather than the hierachical structure." << std::endl;
  std::cerr << "\tif k=-1 then displays the hierarchical structure rather than the graph at a given level." << std::endl;
  std::cerr << "-v\tverbose mode: gives computation time, information about the hierarchy and modularity." << std::endl;
  std::cerr << "-h\tshow this print_usage message." << std::endl;
  std::exit(0);
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
      case 'w':
        type = WEIGHTED;
        filename_w = argv[++i];
        break;
      case 'p':
        filename_part = argv[++i];
        break;
      case 'q':
        precision = std::atof(argv[++i]);
        break;
      case 'l':
        display_level = std::atoi(argv[++i]);
        break;
      case 'k':
        k1 = std::atoi(argv[++i]);
        break;
      case 'v':
        verbose = true;
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
}

void
display_time(const char *str) {
  time_t rawtime;
  time(&rawtime);
  std::cerr << str << ": " << ctime(&rawtime);
}

int
main(int argc, char **argv) {
  srand(time(NULL) + getpid());

  parse_args(argc, argv);
  time_t time_begin, time_end;
  time(&time_begin);
  if (verbose)
    display_time("Begin");

  Community c(filename, filename_w, type, -1, precision);
  if (filename_part != nullptr)
  {
    c.init_partition(filename_part);
  }
  Graph g;
  bool improvement;
  double mod = c.modularity(), new_mod;
  int level = 0;

  do
  {
    if (verbose)
    {
      std::cerr << "level " << level << ":\n";
      display_time("  start computation");
      std::cerr << "  network size: "
        << c.m_graph.nb_nodes << " nodes, "
        << c.m_graph.nb_links << " links, "
        << c.m_graph.total_weight << " weight." << std::endl;
    }

    improvement = c.one_level();
    new_mod = c.modularity();
    if (++level == display_level)
    {
      g.display();
    }
    if (display_level == -1)
    {
      c.display_partition();
    }
    g = c.partition2graph_binary();
    c = Community(g, -1, precision);

    if (verbose)
    {
      std::cerr << "  modularity increased from " << mod << " to " << new_mod << std::endl;
    }

    mod = new_mod;
    if (verbose)
    {
      display_time("  end computation");
    }
    if (filename_part != nullptr
      && level == 1) // do at least one more computation if partition is provided
    {
      improvement = true;
    }
  } while (improvement);

  time(&time_end);
  if (verbose) {
    display_time("End");
    std::cerr << "Total duration: " << (time_end - time_begin) << " sec." << std::endl;
  }
  std::cerr << new_mod << std::endl;
}
