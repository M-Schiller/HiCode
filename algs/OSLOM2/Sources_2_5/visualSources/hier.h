#pragma once

#ifndef HIER_network_INCLUDED
#define HIER_network_INCLUDED

#include <utility>
#include "position.h"

#define LL 1.5

inline int check_position(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r)
{
  xypos Fplane;
  Fplane.set_all(lab_x, lab_y, lab_r, LL);
  Fplane.move_all();
  Fplane.get_new_pos(lab_x, lab_y);

  return 0;
}

inline int get_partition_from_file_tp_format(const std::string& S, std::map<int, std::deque<int>> & Mlab)
{
  Mlab.clear();

  std::ifstream inb(S);

  std::string st;
  while (getline(inb, st)) {
    std::deque<std::string>  vv;
    separate_strings(st, vv);
    if (!st.empty() && (vv[0] == "#module"
      || vv[0] == "#group"))
    {
      getline(inb, st);

      std::deque<int> v;
      cast_string_to_doubles(st, v);
      sort(v.begin(), v.end());
      if (!v.empty())
        Mlab[cast_int(cast_string_to_double(vv[1]))] = v;
    }
  }

  return 0;
}

inline int get_sizes(const std::string& file_for_sizes, std::map<int, int> & sizes)
{
  sizes.clear();
  std::ifstream gin(file_for_sizes);
  std::string st;
  while (getline(gin, st))
  {
    std::deque<std::string>  vv;
    separate_strings(st, vv);
    if (!st.empty()
      && (vv[0] == "#module" || vv[0] == "#group"))
    {
      sizes[cast_int(cast_string_to_double(vv[1]))] = cast_int(cast_string_to_double(vv[3]));
      getline(gin, st);
    }
  }
  return 0;
}

inline int insert_isolated(visual_network & luca, std::map<int, int> & sizes)
{
  std::map<int, int> A;
  luca.get_id_label(A);
  //prints(A);
  for (auto& size : sizes)
  {
    if (A.find(size.first) == A.end())
    {
      //cout<<"adding... "<<itm->first<<" "<<itm->second<<endl;

      luca.add_isolated(size.first);
      //cherr();
    }
  }

  return 0;
}

inline int o_level(
  double origin_x,
  double origin_y,
  double radius,
  visual_network & luca,
  std::map<int, int> & sizes,
  std::map<int, double> & lab_x,
  std::map<int, double> & lab_y,
  std::map<int, double> & lab_r)
{
  // this function takes a network and soem oher data and gives the position of the nodes and their radii

  lab_r.clear();
  lab_x.clear();
  lab_y.clear();

  /*visual_network luca;
  luca.set_graph("football.dat_files/net2");*/
  //cout<<"network:: "<<luca.size()<<" nodes"<<endl;

  std::deque<double> angles;
  std::deque<double> radii;
  const int seconds = cast_int(luca.size() * radius) + 1;
  luca.circles(seconds, sizes, angles, radii, radius);
  //luca.draw_pajek(angles, radii, radius, origin_x, origin_y);
  //cout<<"sizes.size() "<<sizes.size()<<endl;

  if (angles.size() == 1)
    radius = 0;

  for (int i = 0; i < angles.size(); i++) {
    lab_x[luca.id_of(i)] = radius * cos(angles[i]) + origin_x;
    lab_y[luca.id_of(i)] = radius * sin(angles[i]) + origin_y;
    lab_r[luca.id_of(i)] = radii[i];
  }
  return 0;
}

int another_level(
  std::string short_tpn,
  const std::string& netn,
  const std::string& tpn,
  std::map<int, double> & lab_x,
  std::map<int, double> & lab_y,
  std::map<int, double> & lab_r,
  std::map<int, double> & lab_x_next,
  std::map<int, double> & lab_y_next,
  std::map<int, double> & lab_r_next,
  int level,
  char * b1,
  std::deque<std::string> & edge_append)
{
  //cherr();

  lab_x_next.clear();
  lab_y_next.clear();
  lab_r_next.clear();

  std::map<int, std::deque<int>> Mlab;
  get_partition_from_file_tp_format(std::move(short_tpn), Mlab);

  //cout<<"Mlab... "<<Mlab.size()<<endl;

  visual_network luca;
  luca.set_graph(netn);
  std::cout << "network:: " << luca.size() << " nodes" << std::endl;

  std::map<int, int> sizes_hh;
  if (tpn == "_void_")
  {
    for (auto& itm : Mlab)
    {
      for (int i = 0; i < itm.second.size(); i++)
      {
        sizes_hh[itm.second[i]] = 1;
      }
    }
  }
  else
  {
    get_sizes(tpn, sizes_hh);
  }

  std::map<int, int> lucaAbef;
  luca.get_id_label(lucaAbef);

  //cout<<"luca: "<<luca.size()<<endl;
  //prints(lucaAbef);

  insert_isolated(luca, sizes_hh);

  std::map<int, int> id_occu;		// number of times a certain node appears

  std::map<int, int> lucaA;
  luca.get_id_label(lucaA);

  //cout<<"luca: "<<luca.size()<<endl;
  //prints(lucaA);

  std::map<int, int> sizes_h;
  for (auto itm = sizes_hh.begin(); itm != sizes_hh.end(); ++itm)
  {
    if (lucaA.find(itm->first) == lucaA.end())
      lucaA[itm->first] = -1;
    else
      sizes_h[lucaA[itm->first]] = itm->second;
  }

  /*cout<<"************  lucaA ************ "<<endl;
  prints(lucaA);

  cout<<"************ sizes hh ************ "<<endl;
  prints(sizes_hh);

  cout<<"************ sizes ************ "<<endl;
  prints(sizes_h);*/

  /*cout<<"************ sizes (lab luca -> sizes) ************ "<<endl;
  prints(sizes_h);*/

  std::multimap<int, int> internal_links;

  for (std::map<int, std::deque<int>> ::iterator itm_M = Mlab.begin(); itm_M != Mlab.end(); ++itm_M)
  {
    /*cout<<"submodule: "<<itm->first<<endl;		// this is the id of previous luca
    cout<<"x,y,r: "<<lab_x[itm->first]<<" "<<lab_y[itm->first]<<" "<<lab_r[itm->first]<<endl;*/

    std::deque<int> & group = itm_M->second;
    std::deque<std::deque<int>> link_per_node;
    std::deque<std::deque<double>> weights_per_node;

    /*cout<<"group before "<<endl;
    prints(group);*/

    /*for(int i=0; i<group.size(); i++) {
      if(lucaA.find(group[i])==lucaA.end())
        cout<<group[i]<<" not found"<<endl;
    }*/

    for (int& grp : group)
      grp = lucaA[grp];

    /*if(itm_M->first==6 && level==3)
      prints(group);*/

    luca.set_subgraph(group, link_per_node, weights_per_node);
    luca.internal_link_map(internal_links, group);

    //cout<<link_per_node.size()<<endl;

    visual_network giovanni;
    giovanni.set_graph(true, link_per_node, weights_per_node, group);

    //cout<<"network:: "<<giovanni.size()<<" nodes and "<<giovanni.edges()<<" edges;\t average degree = "<<2*giovanni.edges()/giovanni.size()<<endl;

    std::map<int, double> lab_x_h;
    std::map<int, double> lab_y_h;
    std::map<int, double> lab_r_h;

    o_level(lab_x[itm_M->first], lab_y[itm_M->first], lab_r[itm_M->first], giovanni, sizes_h, lab_x_h, lab_y_h, lab_r_h);

    // lab_x_h ha gli id di giovanni che sono gli i label di luca - non i suoi id !!!!!!!!!

    /*if(itm_M->first==6 && level==3) for(map<int, double> :: iterator itm= lab_x_h.begin(); itm!= lab_x_h.end(); itm++)
      cout<<"-------> "<<itm->second<<" "<<lab_y_h[itm->first]<<" "<<lab_r_h[itm->first]<<" "<<luca.id_of(itm->first)<<endl;*/

    for (auto itm = lab_x_h.begin(); itm != lab_x_h.end(); ++itm) {
      int_histogram(luca.id_of(itm->first), lab_x_next, itm->second);
      int_histogram(luca.id_of(itm->first), id_occu);
      /*if(level==3 && luca.id_of(itm->first)==28) {
        cout<<"ok "<<itm->first<<endl;
        cout<<"group"<<endl;
        prints(group);
      }*/
    }

    //cout<<"***<< >> "<<id_occu[28]<<endl;

    for (auto itm = lab_y_h.begin(); itm != lab_y_h.end(); ++itm)
      int_histogram(luca.id_of(itm->first), lab_y_next, itm->second);

    for (auto itm = lab_r_h.begin(); itm != lab_r_h.end(); ++itm)
      int_histogram(luca.id_of(itm->first), lab_r_next, itm->second);
  }

  /*
  prints(id_occu);
  cherr();
  //*/

  for (auto itm = lab_x_next.begin(); itm != lab_x_next.end(); ++itm)
    itm->second /= id_occu[itm->first];
  for (auto itm = lab_y_next.begin(); itm != lab_y_next.end(); ++itm)
    itm->second /= id_occu[itm->first];
  for (auto itm = lab_r_next.begin(); itm != lab_r_next.end(); ++itm)
    itm->second /= id_occu[itm->first];

  /*ofstream poso("posout");
  ofstream posop("posoutp");

  for(map<int, double> :: iterator itm= lab_x.begin(); itm!= lab_x.end(); itm++)
    posop<<itm->second<<" "<<lab_y[itm->first]<<" "<<lab_r[itm->first]<<" "<<itm->first<<endl;

  for(map<int, double> :: iterator itm= lab_x_next.begin(); itm!= lab_x_next.end(); itm++)
    poso<<itm->second<<" "<<lab_y_next[itm->first]<<" "<<lab_r_next[itm->first]<<" "<<itm->first<<endl;*/

    //  ******************** here you can check the overlaps between the nodes
  check_position(lab_x_next, lab_y_next, lab_r_next);
  //  ******************** here you can check the overlaps between the nodes

  // *********************** writing *************************************************

  char posf[1000];
  char graphf[1000];
  sprintf(posf, "%s_oslo_files/pos_%d", b1, level);
  sprintf(graphf, "%s_oslo_files/ggraph_%d", b1, level);

  std::cout << "writing..." << std::endl;
  luca.draw_gnuplot(lab_x_next, lab_y_next, lab_r_next, std::string(posf), std::string(graphf));
  luca.draw_gnuplot(lab_x_next, lab_y_next, lab_r_next, edge_append, internal_links);

  /*deque<string> justint;
  luca.draw_gnuplot(lab_x_next, lab_y_next, lab_r_next, justint, internal_links);
  sprintf(posf, "%s_oslo_files/justpos_%d", b1, level);
  sprintf(graphf, "%s_oslo_files/graph_just_%d", b1, level);
  luca.draw_gnuplot(lab_x_next, lab_y_next, lab_r_next, string(posf), string(graphf), justint);*/

  /*
  if(level==0) {
    sprintf(posf, "%s_oslo_files/specpos_%d", b1, level);
    sprintf(graphf, "%s_oslo_files/graph_spec_%d", b1, level);
    luca.draw_gnuplot(lab_x_next, lab_y_next, lab_r_next, string(posf), string(graphf), edge_append);
  }*/

  return 0;
}

int all_levels(int levels, std::string network_file) {
  std::deque<std::string>  edge_append;

  visual_network luca;
  char b1[1000];
  cast_string_to_char(std::move(network_file), b1);
  char b2[1000];
  sprintf(b2, "%s_oslo_files/net%d", b1, levels);

  luca.set_graph(std::string(b2));		// top hierarchy
  std::cout << "network:: " << luca.size() << " nodes" << std::endl;

  std::map<int, double> lab_x;
  std::map<int, double> lab_y;
  std::map<int, double> lab_r;

  // ************ initial center and radius
  double origin_x = 0.5;
  double origin_y = 0.5;
  double radius = 0.2;
  // ************ initial center and radius

  std::map<int, int> sizes;

  if (levels > 1)
    sprintf(b2, "%s_oslo_files/tp%d", b1, levels - 1);
  else
    sprintf(b2, "%s_oslo_files/tp", b1);

  //get_sizes("football.dat_files/tp1", sizes);
  get_sizes(std::string(b2), sizes);

  if (levels > 0)
    insert_isolated(luca, sizes);

  o_level(origin_x, origin_y, radius, luca, sizes, lab_x, lab_y, lab_r);
  check_position(lab_x, lab_y, lab_r);

  char posf[1000];
  char graphf[1000];
  sprintf(posf, "%s_oslo_files/pos_%d", b1, levels);
  sprintf(graphf, "%s_oslo_files/ggraph_%d", b1, levels);

  std::cout << "writing..." << std::endl;
  luca.draw_gnuplot(lab_x, lab_y, lab_r, std::string(posf), std::string(graphf));
  {
    std::multimap<int, int> internal_links;
    std::deque<int> group;
    for (int i = 0; i < luca.size(); i++)
      group.push_back(i);

    luca.internal_link_map(internal_links, group);
    //prints(internal_links);
    luca.draw_gnuplot(lab_x, lab_y, lab_r, edge_append, internal_links);
  }
  //*************************************************************************************************
  //*************************************************************************************************

  // now I need to get a partition. each node id of luca corresponds to a module here

  std::map<int, double> lab_x_next;		// id of luca (which is the right module number) -> position etc.
  std::map<int, double> lab_y_next;
  std::map<int, double> lab_r_next;

  // subgraphs					// network				// just for the sizes

  for (int i = levels - 1; i >= 0; i--)
  {
    if (i > 1)
      sprintf(b2, "%s_oslo_files/tp%d", b1, i - 1);
    else if (i == 1)
      sprintf(b2, "%s_oslo_files/tp", b1);
    else
      sprintf(b2, "_void_");

    char b3[1000];
    if (i > 0)
      sprintf(b3, "%s_oslo_files/net%d", b1, i);
    else
      sprintf(b3, "%s", b1);

    char b4[1000];
    if (i > 0)
      sprintf(b4, "%s_oslo_files/short_tp%d", b1, i);
    else
      sprintf(b4, "%s_oslo_files/tp", b1);

    std::cout << "reading ... " << std::string(b4) << " " << std::string(b3) << " " << std::string(b2) << std::endl;

    another_level(std::string(b4), std::string(b3), std::string(b2), lab_x, lab_y, lab_r, lab_x_next, lab_y_next, lab_r_next, i, b1, edge_append);

    lab_x = lab_x_next;
    lab_y = lab_y_next;
    lab_r = lab_r_next;
  }

  //another_level("football.dat_files/tp", "football.dat", "_void_", lab_x, lab_y, lab_r, lab_x_next, lab_y_next, lab_r_next, 2);

  // ************************ // ******************************************  //  **************************
}

#endif
