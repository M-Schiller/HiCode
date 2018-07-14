#pragma once

#ifndef NET_vi_included
#define NET_vi_included

#define UNFF -7412128

#include <map>

#include "static_network.h"
#include "hier.h"

class netvi : public static_network
{
public:

  netvi() :static_network() {};
  ~netvi() = default;

  int set_all(int level, const std::string& network_file);
  int print_positions(std::ostream & go, int color, double siz, const std::string& outfile, int level,
    const std::string& network_file, bool first);
  int print_internal(std::ostream & go, int color, const std::string& outfile, int level, const std::string
    & network_file, double edge_width);
  int print_external(std::ostream & go, int color, const std::string& outfile, int level, const std::string
    & network_file, double edge_width);
  int print_label(std::ostream & go, int color, const std::string& outfile, int level, const std::string
    & network_file, double siz);
  int print_module_id(int module);
  int print_module_names(int module, std::map<int, std::string> & names);
  int print_names_file(std::ostream & go, int color, const std::string& outfile, int level,
    const std::string& network_file, double siz, std::map<int, std::string> & names);
  int print_names_file(std::ostream & go, int color, const std::string& outfile, int level,
    const std::string& network_file, double siz, std::map<int, std::string> & names, std::deque<int> & group_ids);

  std::map<int, std::deque<int>> node_modules;		// this maps the id into its submodule

private:

  std::deque<double> label_x;		// label, not id!
  std::deque<double> label_y;

  std::deque<std::pair<double, double>> internal_links;
  std::deque<std::pair<double, double>> external_links;
};

int netvi::print_module_names(int module, std::map<int, std::string> & names)
{
  std::deque<int> & aa = node_modules[module];
  for (int i = 0; i < aa.size(); i++)
    std::cout << names[aa[i]] << " ";
  std::cout << std::endl;

  return 0;
}

inline int netvi::print_module_id(int module)
{
  /*cout<<"this network"<<endl;

  for(map<int, deque<int>> :: iterator itm_M = node_modules.begin(); itm_M!=node_modules.end(); itm_M++) {
    cout<<itm_M->first<<" ---------  ";
    prints(itm_M->second);
  }*/

  prints(node_modules[module]);

  return 0;
}

inline int netvi::print_positions(std::ostream & go, int color, double siz, const std::string&
  outfile, int level,
  const std::string& network_file, bool first)
{
  go << "set output \"" << outfile << ".eps\"" << std::endl;

  if (first)
    go << "p \"" << network_file << "_files/pos_" << level << "\" ps " << siz << " pt 4 lt " << color << std::endl;
  else
    go << "rep \"" << network_file << "_files/pos_" << level << "\" ps " << siz << " pt 4 lt " << color << std::endl;

  return 0;
}

inline int netvi::print_internal(std::ostream & go, int color, const std::string& outfile, int level,
  const std::string& network_file, double edge_width)
{
  //cout<<"internal ... "<<internal_links.size()<<endl;

  for (auto& internal_link : internal_links)
  {
    go << "set arrow from " << label_x[internal_link.first] << "," << label_y[internal_link.first] << " to " << label_x[
      internal_link.second] << "," << label_y[internal_link.second]
        << " nohead lw " << edge_width * log(1 + vertices[internal_link.first]->links->posweightof(
          internal_link.second).second) << std::endl;
  }

  return 0;
}

inline int netvi::print_external(std::ostream & go, int color, const std::string& outfile, int level,
  const std::string& network_file, double edge_width)
{
  for (int i = 0; i < external_links.size(); i++)
  {
    go << "set arrow from " << label_x[external_links[i].first] << "," << label_y[external_links[i].first] << " to " << label_x[external_links[i].second] << "," << label_y[external_links[i].second]
      << " nohead lw " << edge_width * log(vertices[external_links[i].first]->links->posweightof(external_links[i].second).second) << std::endl;
  }

  return 0;
}

inline int netvi::print_label(std::ostream & go, int color, const std::string& outfile, int level,
  const std::string& network_file, double siz)
{
  // color is not used here .........

  for (int i = 0; i < dim; i++)
    go << "set label \"" << id_of(i) << "\" at " << label_x[i] << "," << label_y[i] << " font \"Helvetica, " << siz << "\"" << std::endl;

  return 0;
}

inline int netvi::print_names_file(std::ostream & go, int color, const std::string& outfile, int level,
  const std::string& network_file, double siz, std::map<int, std::string> & names)
{
  // color is not used here .........

  for (int i = 0; i < dim; i++)
    go << "set label \"" << names[id_of(i)] << "\" at " << label_x[i] << "," << label_y[i] << " font \"Helvetica, " << siz << "\"" << std::endl;

  return 0;
}

inline int netvi::print_names_file(
  std::ostream & go,
  int color,
  const std::string& outfile,
  int level,
  const std::string& network_file,
  double siz,
  std::map<int, std::string> & names,
  std::deque<int> & group_ids)
{
  // color is not used here .........

  for (auto i = 0; i < group_ids.size(); i++)
  {
    go
      << "set label \"" << names[group_ids[i]] << "\" at " << label_x[i] << "," << label_y[i] << " font \"Helvetica, " << siz << "\"" << std::endl;
  }

  return 0;
}

inline int netvi::set_all(int level, const std::string& network_file)
{
  // I only need pos_%d from visual network

  char b1[1000];
  char b2[1000];

  //cast_string_to_char(network_file, b1);

  if (level > 0)
  {
    sprintf(b2, "%s_files/net%d", network_file.c_str(), level);
  }
  else
  {
    sprintf(b2, "%s", network_file.c_str());
  }

  std::cout << "setting graph... " << b2 << std::endl;

  set_graph(std::string(b2));

  std::map<int, int> A;
  get_id_label(A);

  if (level > 0)
  {
    std::map<int, int> sizes;

    if (level > 1)
    {
      sprintf(b2, "%s_files/tp%d", network_file.c_str(), level - 1);
    }
    else
    {
      sprintf(b2, "%s_files/tp", network_file.c_str());
    }

    //cout<<"sizes..."<<endl;

    //get_sizes("football.dat_files/tp1", sizes);
    get_sizes(std::string(b2), sizes);

    //prints(sizes);

    for (auto& size : sizes)
    {
      if (A.find(size.first) == A.end())
      {
        //std::cout<<"adding... "<<itm->first<<" "<<itm->second<<endl;

        add_isolated(size.first);
        //cherr();
      }
    }
  }

  get_id_label(A);

  sprintf(b2, "%s_files/pos_%d", network_file.c_str(), level);
  std::cout << "setting positions... " << b2 << std::endl;

  std::deque<double> xs;
  std::deque<double> ys;
  std::deque<int> labs;
  get_data_from_file(std::string(b2), xs, ys, 1, 2);
  get_data_from_file(std::string(b2), labs, 4);

  label_x.clear();
  label_y.clear();

  for (int i = 0; i < dim; i++)
  {
    label_x.push_back(0);
    label_y.push_back(0);
  }

  //prints(labs);

  for (int i = 0; i < labs.size(); i++)
  {
    label_x[A[labs[i]]] = xs[i];
    label_y[A[labs[i]]] = ys[i];
  }

  /*cout<<"A** "<<endl;
  prints(A);
  prints(label_x);*/

  if (level > 0)
  {
    sprintf(b2, "%s_files/short_tp%d", b1, level);
  }
  else
  {
    sprintf(b2, "%s_files/tp", b1);
  }

  std::cout << "setting partition... " << b2 << std::endl;
  std::deque<std::deque<int>> M;

  {
    std::ifstream gin(b2);
    if (gin.is_open()) {
      get_partition_from_file(std::string(b2), M);
      if (translate(M) == -1)
        return -1;
    }
    else {
      std::deque<int> aa;

      for (int i = 0; i < dim; i++) {
        aa.push_back(i);
      }

      M.push_back(aa);
    }
  }

  internal_links.clear();
  external_links.clear();

  for (int i = 0; i < M.size(); i++)
  {
    for (int ii = 0; ii < M[i].size(); ii++)
    {
      int & nodei = M[i][ii];

      for (int j = 0; j < vertices[nodei]->links->size(); j++) if (nodei < vertices[nodei]->links->l[j]) {
        //cout<<nodei<<" "<<vertices[nodei]->links->l[j]<<endl;

        if (binary_search(M[i].begin(), M[i].end(), vertices[nodei]->links->l[j]))
        {
          internal_links.emplace_back(nodei, vertices[nodei]->links->l[j]);		// there can be multiple links but
        }
        else
        {
          external_links.emplace_back(nodei, vertices[nodei]->links->l[j]);
        }
      }
    }
  }

  /*print_internal(cout, 1, "sss", level, "acc", 0.01);
  cherr();*/

  sort(internal_links.begin(), internal_links.end());
  sort(external_links.begin(), external_links.end());

  auto ituni = unique(internal_links.begin(), internal_links.end());
  internal_links.resize(ituni - internal_links.begin());
  ituni = unique(external_links.begin(), external_links.end());
  external_links.resize(ituni - external_links.begin());

  if (level > 1)
    sprintf(b2, "%s_files/tp%d", network_file.c_str(), level - 1);
  else if (level == 1)
    sprintf(b2, "%s_files/tp", network_file.c_str());
  else
    sprintf(b2, "%s_files/nodes", network_file.c_str());

  std::cout << "setting submodules... " << b2 << std::endl;

  if (level > 0)
  {
    //get_partition_from_file_tp_format(string(b2), node_modules);
    get_partition_from_file(std::string(b2), node_modules, 1);
    /*node_modules.clear();

    for(map<int, deque<int>> :: iterator itm_M = helpls.begin(); itm_M!=helpls.end(); itm_M++)
      node_modules[A[itm_M->first]]=itm_M->second;	*/
  }
  else
  {
    for (int i = 0; i < dim; i++)
    {
      std::deque<int> first;
      first.push_back(id_of(i));
      node_modules[id_of(i)] = first;
    }
  }

  return 0;
}

class netgnu
{
public:

  netgnu() = default;
  ~netgnu();

  int set_networks(int levels, const std::string& network_file);

  int set_position(int net, int color, double siz);
  int erase_position(int);

  int set_internal(int net, double edge_width);
  int erase_internal(int);

  int set_external(int net, double edge_width);
  int erase_external(int);

  int set_labs(int net, double siz);
  int erase_labs(int);

  int interface();
  int plot_commands();

  int read_names(const std::string& s);
  int print_id(int net, int module);
  int print_names(int net, int module);
  int print_names_on_graph(double siz);

private:

  std::deque<netvi *> networks;
  std::map<int, std::deque<double>> positions_commands;		// each row corresponds to a network - deque<double> are parameters
  std::map<int, std::deque<double>> internal_commands;		// each row corresponds to a network - deque<double> are parameters
  std::map<int, std::deque<double>> external_commands;		// each row corresponds to a network - deque<double> are parameters
  std::map<int, std::deque<double>> label_commands;		// each row corresponds to a network - deque<double> are parameters

  std::string network_file_saved;
  std::string outfile_saved;

  double minx;
  double maxx;

  double miny;
  double maxy;

  bool unset_key;
  double names_size;

  std::map<int, std::string> id_names;
};

netgnu::~netgnu()
{
  for (auto& network : networks)
  {
    delete network;
  }
}

inline int netgnu::set_networks(int levels, const std::string& network_file)
{
  for (int i = 0; i <= levels; i++)
  {
    auto* newnet = new netvi;
    newnet->set_all(i, network_file);
    networks.push_back(newnet);
  }

  network_file_saved = network_file;
  outfile_saved = "new_gnuplot_graph";

  minx = UNFF;
  miny = UNFF;

  unset_key = false;
  names_size = UNFF;

  return 0;
}

inline int netgnu::set_labs(int net, double siz)
{
  std::deque<double> h;
  h.push_back(siz);
  label_commands[net] = h;
  return 0;
}

inline int netgnu::erase_labs(int net)
{
  label_commands.erase(net);
  return 0;
}

inline int netgnu::set_position(int net, int color, double siz)
{
  std::deque<double> h;
  h.push_back(color);
  h.push_back(siz);
  positions_commands[net] = h;
  return 0;
}

inline int netgnu::erase_position(int net)
{
  positions_commands.erase(net);
  return 0;
}

inline int netgnu::set_internal(int net, double edge_width)
{
  std::deque<double> h;
  h.push_back(edge_width);
  internal_commands[net] = h;
  return 0;
}

inline int netgnu::erase_internal(int net)
{
  internal_commands.erase(net);
  return 0;
}

inline int netgnu::set_external(int net, double edge_width)
{
  std::deque<double> h;
  h.push_back(edge_width);
  external_commands[net] = h;
  return 0;
}

inline int netgnu::erase_external(int net)
{
  external_commands.erase(net);
  return 0;
}

inline int netgnu::print_id(int net, int module)
{
  networks[net]->print_module_id(module);
  return 0;
}

inline int netgnu::print_names(int net, int module)
{
  if (id_names.empty()) {
    std::cout << "names was not set" << std::endl;
    return -1;
  }

  networks[net]->print_module_names(module, id_names);
  return 0;
}

inline int netgnu::print_names_on_graph(double siz)
{
  if (id_names.empty())
  {
    std::cout << "names was not set" << std::endl;
    return -1;
  }

  names_size = siz;

  return 0;
}

inline int netgnu::interface()
{
  std::string s;
  std::cout << "press h for help" << std::endl;
  std::cout << "> ";

  while (std::cin >> s)
  {
    //prints(s);

    std::deque<std::string> sep;
    std::cout << "parsing ... ";
    separate_strings(s, sep);
    prints(sep);

    if (sep[0] == "quit"
      || sep[0] == "q")
    {
      break;
    }

    // ***************************************************************
    if (sep[0] == "pos")
    {
      double n;
      if (sep.size() < 2
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
        std::cout << "pos, <i>, (<color>,) (<size>,) : set positions of nodes for level i" << std::endl;
      }
      else if (sep.size() == 2)
      {
        set_position(cast_int(n), -1, 0.1);
      }
      else if (sep.size() == 3)
      {
        set_position(cast_int(n), cast_int(cast_string_to_double(sep[2])), 0.1);
      }
      else if (sep.size() == 4)
      {
        set_position(cast_int(n), cast_int(cast_string_to_double(sep[2])), cast_string_to_double(sep[3]));
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "print_names"
      || sep[0] == "g")
    {
      double n;
      if (sep.size() < 3
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
        std::cout << "print_names, <i>, <module_label> : print names which are in module <community_label> of label i" << std::endl;
      }
      else
      {
        print_names(cast_int(n), cast_int(cast_string_to_double(sep[2])));
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "print_id")
    {
      double n;
      if (sep.size() < 3
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
        std::cout << "print_id, <i>, <module_label> : print ids which are in module community_label of label i" << std::endl;
      }
      else
      {
        print_id(cast_int(n), cast_int(cast_string_to_double(sep[2])));
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "lab")
    {
      double n;
      if (sep.size() < 2 || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
        std::cout << "lab, <i>, (<size>,) : set labels of nodes for level i" << std::endl;
      }
      else if (sep.size() == 2)
      {
        set_labs(cast_int(n), 10);
      }
      else if (sep.size() == 3)
      {
        set_labs(cast_int(n), cast_int(cast_string_to_double(sep[2])));
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "nolab")
    {
      double n;
      if (sep.size() < 2 || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
        std::cout << "nolab,<i> : erase labels of nodes for level i" << std::endl;
      }
      else
        erase_labs(cast_int(n));
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "nop")
    {
      double n;
      if (sep.size() < 2 || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
        std::cout << "int,<i> : set internal links for level i" << std::endl;
      }
      else
      {
        erase_position(cast_int(n));
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "int")
    {
      double n;
      if (sep.size() < 2
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
      }
      else
      {
        set_internal(cast_int(n), 0.1);
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "noi")
    {
      double n;
      if (sep.size() < 2
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
      }
      else
      {
        erase_internal(cast_int(n));
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "ext")
    {
      double n;
      if (sep.size() < 2
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
      }
      else
      {
        set_external(cast_int(n), 0.1);
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "noe")
    {
      double n;
      if (sep.size() < 2
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
      }
      else
      {
        erase_external(cast_int(n));
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "plot"
      || sep[0] == "r"
      || sep[0] == "rep"
      || sep[0] == "p")
    {
      plot_commands();
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "unset_key")
    {
      unset_key = true;
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "zoom")
    {
      if (sep.size() < 4)
      {
        std::cout << "error" << std::endl;
        std::cout << "zoom,<x>,<y>,<delta> : zoom around (x, y) --  deltax=deltay=delta" << std::endl;
      }
      else
      {
        const double xx = cast_string_to_double(sep[1]);
        const double yy = cast_string_to_double(sep[2]);
        const double dd = cast_string_to_double(sep[3]);

        minx = xx - dd / 2;
        maxx = xx + dd / 2;
        miny = yy - dd / 2;
        maxy = yy + dd / 2;
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "names")
    {
      if (sep.size() < 2)
      {
        std::cout << "error" << std::endl;
        std::cout << "names,<file> : set names from file" << std::endl;
      }

      else {
        read_names(sep[1]);
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "lab_names")
    {
      double n;

      if (sep.size() < 2
        || !cast_string_to_double(sep[1], n))
      {
        std::cout << "error" << std::endl;
        std::cout << "lab_names, <size> : print names which are in module <community_label> of label i on the graph" << std::endl;
      }

      else
      {
        print_names_on_graph(n);
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "h"
      || sep[0] == "help")
    {
      std::cout << "lab_names, <i>, <module_label> : print names which are in module <community_label> of label i on the graph" << std::endl;
      std::cout << "print_names, <i>, <module_label> : print to screen names which are in module <community_label> of label i" << std::endl;
      std::cout << "print_id, <i>, <module_label> : print to screen ids which are in module community_label of label i" << std::endl;
      std::cout << "lab, <i>, (<size>,) : set labels of nodes for level i" << std::endl;
      std::cout << "pos, <i>, (<color>,) (<size>,) : set positions of nodes for level i" << std::endl;
      std::cout << "nolab,<i> : erase labels of nodes for level i" << std::endl;
      std::cout << "nop,<i> : erase positions of nodes for level i" << std::endl;
      std::cout << "int,<i> : set internal links for level i" << std::endl;
      std::cout << "noi,<i> : erase internal links for level i" << std::endl;
      std::cout << "ext,<i> : set external links for level i" << std::endl;
      std::cout << "noe,<i> : erase external links for level i" << std::endl;
      std::cout << "plot (r,p,rep) : plot" << std::endl;
      std::cout << "clear : clear" << std::endl;
      std::cout << "quit : quit the program" << std::endl;
      std::cout << "allint : draw all the internal links" << std::endl;
      std::cout << "open  : open " << outfile_saved << std::endl;
      std::cout << "unset_key  : unset key " << std::endl;
      std::cout << "zoom,<x>,<y>,<delta> : zoom around (x, y) --  deltax=deltay=delta" << std::endl;
      std::cout << "names,<file> : set names from file" << std::endl;

      std::cout << "there are " << networks.size() - 1 << " levels, so <i> can be between 0 and " << networks.size() - 1 << std::endl;
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "clear")
    {
      positions_commands.clear();
      internal_commands.clear();
      external_commands.clear();
      label_commands.clear();
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "allint")
    {
      for (auto i = 0; i < networks.size(); i++)
      {
        set_internal(i, 0.1);
      }
    }
    // ***************************************************************
    // ***************************************************************
    else if (sep[0] == "open")
    {
      systemCall("open " + outfile_saved);
    }
    // ***************************************************************

    else
      std::cout << "error" << std::endl;

    std::cout << "> ";
  }

  return 0;
}

inline int netgnu::plot_commands()
{
  if (positions_commands.empty())
  {
    std::cout << "no position set. there must be at least one" << std::endl;
    return -1;
  }

  std::ofstream go(outfile_saved);

  go << "set term post eps enh color \"Helvetica\" 20" << std::endl;
  go << "set output \"" << outfile_saved << ".eps\"" << std::endl;

  if (minx != UNFF) {
    go << "set xrange [" << minx << ":" << maxx << "]" << std::endl;
  }
  if (miny != UNFF) {
    go << "set yrange [" << miny << ":" << maxy << "]" << std::endl;
  }

  if (names_size != UNFF)
  {
    networks[0]->print_names_file(go, 0, outfile_saved, 0, network_file_saved, names_size, id_names);
  }

  if (unset_key)
    go << "unset key" << std::endl;

  for (auto& label_command : label_commands)
  {
    if (label_command.first >= 0
      && label_command.first < networks.size())
    {
      networks[label_command.first]->print_label(go, 0, outfile_saved, label_command.first, network_file_saved,
        label_command.second[0]);
    }
  }

  for (auto& internal_command : internal_commands)
  {
    if (internal_command.first >= 0
      && internal_command.first < networks.size())
    {
      networks[internal_command.first]->print_internal(go, 0, outfile_saved, internal_command.first, network_file_saved,
        internal_command.second[0]);
    }
  }
  for (auto& external_command : external_commands)
  {
    if (external_command.first >= 0
      && external_command.first < networks.size())
    {
      networks[external_command.first]->print_external(go, 0, outfile_saved, external_command.first, network_file_saved,
        external_command.second[0]);
    }
  }

  auto IT = positions_commands.end();
  bool first_to_print = true;
  while (true)
  {
    if (IT == positions_commands.begin())
      break;

    --IT;

    if (IT->first >= 0
      && IT->first < networks.size())
    {
      networks[IT->first]->print_positions(
        go, IT->second[0],
        IT->second[1],
        outfile_saved,
        IT->first,
        network_file_saved,
        first_to_print);
      first_to_print = false;
    }
  }

  /*
  char hb[1000];
  string hh= "gnuplot " + outfile_saved;
  //cout<<"plotting ... "<<hh<<endl;
  cast_string_to_char(hh, hb);
  //systemCall(hb);

  //hh= "open " + outfile_saved + ".eps";
  //cout<<hh<<endl;
  cast_string_to_char(hh, hb);
  systemCall(hb);*/

  return 0;
}

inline int netgnu::read_names(const std::string& s)
{
  std::string r;
  std::ifstream gin(s);
  while (gin >> r)
  {
    std::string w;
    gin >> w;
    id_names[cast_int(cast_string_to_double(w))] = r;
  }

  return 0;
}

#endif
