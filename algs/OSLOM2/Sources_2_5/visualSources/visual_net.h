#pragma once

#ifndef VISUAL_network_INCLUDED
#define VISUAL_network_INCLUDED

#include <utility>
#include "./standard_package/standard_include.cpp"
#include "static_network.h"
#include "wsarray.h"

class visual_network : public static_network
{
public:

  visual_network()
    : static_network()
    , label_size(0)
    , edge_width(0.1)
  {
  };
  ~visual_network() {};

  void internal_link_map(std::multimap<int, int> & int_edges, std::deque<int> & group);
  //double edge_list_map(map<int, int> & edge_list_app);
  double draw_gnuplot(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r,
    const std::string& position_file,
    const std::string& gfile, std::deque<std::string> & edge_append);
  double draw_gnuplot(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r, std::deque<std::string> & edge_append, std::multimap<int, int> & internal_links);
  double draw_gnuplot(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r,
    const std::string& position_file,
    const std::string& gfile);
  double draw_pajek(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r);
  double draw_pajek(std::deque<double> & angles, std::deque<double> & radii, double radius, double origin_x, double origin_y);
  double with_colors(std::deque<std::deque<int>>  & ten);
  int find_position_goniometric(int secs);
  int circles(int secs, std::map<int, int> & lab_sizes, std::deque<double> & angles, std::deque<double> & radii, double radius);

private:

  double label_size;
  double edge_width;

  double relabel(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r, std::map<int, double> & labx_h, std::map<int, double> & laby_h, std::map<int, double> & labr_h);

  //***************************** methods of find_position_goniometric **************************
  double distance_node(int node1);
  bool swap();
  std::deque<int> pos;
  std::deque<double> sin_tab;
  double theta;
  //***************************** methods of find_position_goniometric **************************
};

inline double visual_network::distance_node(int node1) {
  double distance = 0;
  for (int j = 0; j < vertices[node1]->links->size(); j++)
    distance += sin_tab[abs(pos[node1] - pos[vertices[node1]->links->l[j]])] * vertices[node1]->links->w[j];

  return distance;
}

inline bool visual_network::swap() {
  int node1 = irand(pos.size() - 1);
  int node2 = irand(pos.size() - 1);

  if (node1 == node2)
    return false;

  double op1 = distance_node(node1);
  double op2 = distance_node(node2);

  int op = pos[node1];
  pos[node1] = pos[node2];
  pos[node2] = op;

  double np1 = distance_node(node1);
  double np2 = distance_node(node2);

  if (np1 + np2 + 1e-7 < op1 + op2)
    return true;

  if (np1 + np2 > op1 + op2) {		// refuses
    op = pos[node1];
    pos[node1] = pos[node2];
    pos[node2] = op;

    return false;
  }
}

int visual_network::find_position_goniometric(int secs) {
  // this function is to position the nodes in a ring
  // secs is the maximum time you want to spend on it (seconds)

  sin_tab.clear();

  time_t rawtime = time(NULL);

  //cout<<ctime(& rawtime)<<endl;

  time_t rawtime2 = time(NULL);

  // here you need a table

  theta = 2 * M_PI / dim;
  double theta2 = M_PI / dim;
  //cout<<"theta: "<<theta<<endl;

  pos.clear();			// pos[i] is position of node i
  for (int i = 0; i < dim; i++) {
    pos.push_back(i);
    sin_tab.push_back(theta2 * i);
  }

  /*distance=0;

  for(int i=0; i<dim; i++)
    distance+=distance_node(i);*/
    //distance=distance/2;

  int stopper = 0;

  while (true) {
    if (swap() == false)
      stopper++;
    else
      stopper = 0;

    if (stopper == dim)
      break;

    if (ran4() < 1e-5) {
      //cout<<"stopper "<<stopper<<endl;

      rawtime2 = time(NULL);
      if (rawtime2 - rawtime > secs)		// secs maximum
        break;
    }
  }

  //prints(pos);

  return 0;
}

double visual_network::draw_pajek(std::deque<double> & angles, std::deque<double> & radii, double radius, double origin_x, double origin_y)
{
  std::ofstream subout("pajek_file");
  subout << "*Vertices " << dim << std::endl;
  for (int i = 0; i < dim; i++)
    subout << "   " << i + 1 << " \"" << id_of(i) << "\" " << radius * cos(angles[i]) + origin_x << " " << radius * sin(angles[i]) + origin_y << " " << radii[i] << " " << " box ic " << "red" << std::endl;

  subout << "*Edges" << std::endl;

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < vertices[i]->links->size(); j++)
      if (i < vertices[i]->links->l[j])
        subout << i + 1 << " " << vertices[i]->links->l[j] + 1 << " " << vertices[i]->links->w[j] << std::endl;
  }

  //systemCall("unixdos pajek_file pajek_file.net");
}

inline double visual_network::draw_pajek(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r) {
  std::map<int, int> A;
  get_id_label(A);

  std::map<int, double> labx_h;
  std::map<int, double> laby_h;
  std::map<int, double> labr_h;

  for (auto itm = lab_x.begin(); itm != lab_x.end(); ++itm)
    labx_h[A[itm->first]] = itm->second;

  for (auto itm = lab_y.begin(); itm != lab_y.end(); ++itm)
    laby_h[A[itm->first]] = itm->second;

  for (auto itm = lab_r.begin(); itm != lab_r.end(); ++itm)
    labr_h[A[itm->first]] = itm->second;

  std::ofstream subout("pajek_file");
  subout << "*Vertices " << dim << std::endl;
  for (int i = 0; i < dim; i++)
    subout << "   " << i + 1 << " \"" << id_of(i) << "\" " << labx_h[i] << " " << laby_h[i] << " " << labr_h[i] << " " << " box ic " << "red" << std::endl;

  subout << "*Edges" << std::endl;

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < vertices[i]->links->size(); j++)
      if (i < vertices[i]->links->l[j])
        subout << i + 1 << " " << vertices[i]->links->l[j] + 1 << " " << vertices[i]->links->w[j] << std::endl;
  }

  //systemCall("unixdos pajek_file pajek_file.net");

  return 0;
}

inline void visual_network::internal_link_map(std::multimap<int, int> & int_edges, std::deque<int> & group)
{
  std::sort(group.begin(), group.end());

  //prints(group);
  for (int ii = 0; ii < group.size(); ii++)
  {
    int nodei = group[ii];

    for (int j = 0; j < vertices[nodei]->links->size(); j++)
    {
      //cout<<nodei<<" "<<vertices[nodei]->links->l[j]<<endl;

      if (binary_search(group.begin(), group.end(), vertices[nodei]->links->l[j]) && nodei < vertices[nodei]->links->l[j])
        int_edges.emplace(nodei, vertices[nodei]->links->l[j]);
    }
  }

  //prints(int_edges);
}

/*
double visual_network::edge_list_map(map<int, int> & edge_list_app) {
  for(int j=0; j<vertices[i]->links->size(); j++)
    if(i < vertices[i]->links->l[j])
      edge_list_app[id_of(i)]=id_of(vertices[i]->links->l[j]);
}*/

double visual_network::draw_gnuplot(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r, std::deque<std::string> & edge_append, std::multimap<int, int> & internal_links) {
  std::map<int, double> labx_h;
  std::map<int, double> laby_h;
  std::map<int, double> labr_h;
  relabel(lab_x, lab_y, lab_r, labx_h, laby_h, labr_h);

  char b[1000];
  for (auto itm = internal_links.begin(); itm != internal_links.end(); ++itm)
  {
    //const int & i = itm->first;

    sprintf(b, "set arrow from %f,%f to %f,%f nohead lw %f",
      labx_h[itm->first],
      laby_h[itm->first],
      labx_h[itm->second],
      laby_h[itm->second],
      edge_width * log(1 + vertices[itm->first]->links->posweightof(itm->second).second));
    edge_append.push_back(std::string(b));
  }

  return 0;
}

inline double visual_network::draw_gnuplot(
  std::map<int, double> & lab_x,
  std::map<int, double> & lab_y,
  std::map<int, double> & lab_r,
  const std::string& position_file,
  const std::string& gfile)
{
  std::map<int, double> labx_h;
  std::map<int, double> laby_h;
  std::map<int, double> labr_h;

  relabel(lab_x, lab_y, lab_r, labx_h, laby_h, labr_h);

  std::ofstream suboutp(position_file);
  //for(int i=0; i<dim; i++)
    //suboutp<<labx_h[i]<<" "<<laby_h[i] <<" "<<labr_h[i]<<" "<<id_of(i)<<endl;

  for (int i = 0; i < dim; i++)
    suboutp << labx_h[i] << " " << laby_h[i] << " " << id_of(i) << std::endl;

  return 0;

  /*
  cast_string_to_char(gfile, b);
  ofstream suboutg(b);

  suboutg<<"set term post eps enh color \"Helvetica\" 20"<<endl;
  suboutg<<"set output \""<<gfile<<".eps\""<<endl;

  for(int i=0; i<dim; i++) if(label_size>0)
    suboutg<<"set label \""<<id_of(i)<<"\" at "<<labx_h[i]<<","<<laby_h[i]<<" font \"Helvetica, "<< min(10 *label_size*labr_h[i], label_size * 5. / dim)<<"\""<<endl;

  /*for(int i=0; i<dim; i++) if(labx_h[i]<0.05 && laby_h[i] <0.05)
    suboutg<<"set label \""<<id_of(i)<<"\" at "<<labx_h[i]<<","<<laby_h[i]<<" font \"Helvetica, "<<1<<"\""<<endl;* /

  for(int i=0; i<dim; i++) {
    for(int j=0; j<vertices[i]->links->size(); j++)
      if(i < vertices[i]->links->l[j])
        suboutg<<"set arrow from "<<labx_h[i]<<","<<laby_h[i]<<" to "<<labx_h[vertices[i]->links->l[j]]<<","<<laby_h[vertices[i]->links->l[j]]<<" nohead lw "<< edge_width * log(1+vertices[i]->links->w[j])<<endl;
  }

  suboutg<<"p \""<<position_file<<"\" ps 0.1 pt 4 lt -1"<<endl;

  char hb[1000];
  string hh= "gnuplot " + gfile;
  cout<<"plotting ... "<<hh<<endl;
  cast_string_to_char(hh, hb);
  systemCall(hb);

  /*
  hh="open " + gfile + ".eps";
  cast_string_to_char(hh, hb);
  systemCall(hb);
  */

  return 0;
}

inline double visual_network::draw_gnuplot(
  std::map<int, double> & lab_x,
  std::map<int, double> & lab_y,
  std::map<int, double> & lab_r,
  const std::string& position_file,
  const std::string& gfile,
  std::deque<std::string> & edge_append)
{
  std::map<int, double> labx_h;
  std::map<int, double> laby_h;
  std::map<int, double> labr_h;
  relabel(lab_x, lab_y, lab_r, labx_h, laby_h, labr_h);

  std::ofstream suboutp(position_file);
  //for(int i=0; i<dim; i++)
    //suboutp<<labx_h[i]<<" "<<laby_h[i] <<" "<<labr_h[i]<<" "<<id_of(i)<<endl;

  for (int i = 0; i < dim; i++)
    suboutp << labx_h[i] << " " << laby_h[i] << " " << id_of(i) << std::endl;

  return 0;

  /*

  ofstream suboutg(gfile);

  suboutg<<"set term post eps enh color \"Helvetica\" 20"<<endl;
  suboutg<<"set output \""<<gfile<<".eps\""<<endl;

  for(int i=0; i<dim; i++) if(label_size>0)
    suboutg<<"set label \""<<id_of(i)<<"\" at "<<labx_h[i]<<","<<laby_h[i]<<" font \"Helvetica, "<< min(10 *label_size*labr_h[i], label_size * 5. / dim)<<"\""<<endl;

  /*for(int i=0; i<dim; i++) if(labx_h[i]<0.05 && laby_h[i] <0.05)
    suboutg<<"set label \""<<id_of(i)<<"\" at "<<labx_h[i]<<","<<laby_h[i]<<" font \"Helvetica, "<<10<<"\""<<endl;* /

  for(int i=0; i<edge_append.size(); i++)
    suboutg<<edge_append[i]<<endl;

  suboutg<<"p \""<<position_file<<"\" ps 0.1 pt 4 lt -1"<<endl;

  char hb[1000];
  string hh= "gnuplot " + gfile;
  cout<<"plotting ... "<<hh<<endl;
  cast_string_to_char(hh, hb);
  systemCall(hb);

  /*
  hh="open " + gfile + ".eps";
  cast_string_to_char(hh, hb);
  systemCall(hb);
  return 0;
  */
}

inline double visual_network::relabel(
  std::map<int, double> & lab_x,
  std::map<int, double> & lab_y,
  std::map<int, double> & lab_r,
  std::map<int, double> & labx_h,
  std::map<int, double> & laby_h,
  std::map<int, double> & labr_h)
{
  labx_h.clear();
  laby_h.clear();
  labr_h.clear();

  std::map<int, int> A;
  get_id_label(A);

  /*for(map<int, double>:: iterator itm=lab_x.begin(); itm!=lab_x.end(); itm++) if(fabs(itm->second)<1e-5)
    cerr<<"node in x=0 "<<itm->first<<endl;

  for(map<int, double>:: iterator itm=lab_r.begin(); itm!=lab_r.end(); itm++) if(fabs(itm->second)<1e-5)
    cerr<<"node in r=0 "<<itm->first<<endl;

  for(map<int, double>:: iterator itm=lab_y.begin(); itm!=lab_y.end(); itm++) if(fabs(itm->second)<1e-5)
    cerr<<"node in y=0 "<<itm->first<<endl;*/

  for (auto itm = lab_x.begin(); itm != lab_x.end(); ++itm) if (A.find(itm->first) != A.end())
    labx_h[A[itm->first]] = itm->second;

  for (auto itm = lab_y.begin(); itm != lab_y.end(); ++itm) if (A.find(itm->first) != A.end())
    laby_h[A[itm->first]] = itm->second;

  for (auto itm = lab_r.begin(); itm != lab_r.end(); ++itm) if (A.find(itm->first) != A.end())
    labr_h[A[itm->first]] = itm->second;

  return 0;
}

double visual_network::with_colors(std::deque<std::deque<int>>  & ten) {
  std::deque<std::string> pajek_colors;
  pajek_colors.push_back("Red");
  pajek_colors.push_back("Green");
  pajek_colors.push_back("Black");
  pajek_colors.push_back("Yellow");
  pajek_colors.push_back("White");
  pajek_colors.push_back("Yellow");
  pajek_colors.push_back("Brown");
  pajek_colors.push_back("Mahogany");
  pajek_colors.push_back("Grey");
  pajek_colors.push_back("Blue");

  std::deque<int> mem(dim);
  for (int i = 0; i < ten.size(); i++)
    for (int j = 0; j < ten[i].size(); j++)
      mem[ten[i][j]] = i;

  std::ofstream subout("pajek_file");
  subout << "*Vertices " << dim << std::endl;
  for (int i = 0; i < dim; i++)
    subout << "   " << i + 1 << " \"" << id_of(i) << "\" " << (0.7*cos(theta * pos[i]) + 1) / 2.1 << " " << (0.7*sin(theta * pos[i]) + 1) / 2.1 << " 0.5" << " box ic " << pajek_colors[mem[i] % pajek_colors.size()] << std::endl;

  subout << "*Edges" << std::endl;

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < vertices[i]->links->size(); j++)
      if (i < vertices[i]->links->l[j])
        subout << i + 1 << " " << vertices[i]->links->l[j] + 1 << " " << vertices[i]->links->w[j] << std::endl;
  }

  //systemCall("unixdos pajek_file pajek_file.net");

  return 0;
}

inline int visual_network::circles(int secs, std::map<int, int> & lab_sizes, std::deque<double> & angles, std::deque<double> & radii, double radius) {
  // secs is the maximum time you want to spend on this
  // sizes is the number of nodes belonging to each module
  // radii are the radii of the circles

  angles.clear();
  radii.clear();

  if (dim == 0)
    return -1;

  if (dim == 1) {
    angles.push_back(0);
    radii.push_back(radius);
    return 0;
  }

  std::map<int, int> A;
  get_id_label(A);
  std::deque<int> sizes(dim);

  /*cout<<"***************** A giovanni *****"<<endl;
  prints(A);

  cout<<"***************** lucalab sizes giovanni *****"<<endl;
  prints(lab_sizes);*/

  for (auto itm = A.begin(); itm != A.end(); ++itm)
  {
    sizes[itm->second] = lab_sizes[itm->first];
    //cout<<"<?> "<<A[itm->first]<<" "<<itm->second<<endl;
  }

  /*cout<<"sizes... "<<endl;
  for(int i=0; i<dim; i++)
    cout<<id_of(i)<<" "<<sizes[i]<<endl;
  */

  find_position_goniometric(secs);

  int tots = accumulate(sizes.begin(), sizes.end(), 0);

  std::map<int, int> position_node;
  for (int i = 0; i < pos.size(); i++) {
    position_node[pos[i]] = i;
    angles.push_back(0);
    radii.push_back(0);
  }

  double current_angle = 0;

  for (auto itm = position_node.begin(); itm != position_node.end(); ++itm)
  {
    //cout<<id_of(itm->second)<<" "<<current_angle<<endl;

    current_angle += M_PI * double(sizes[itm->second]) / tots;
    angles[itm->second] = current_angle;
    current_angle += M_PI * double(sizes[itm->second]) / tots;
    radii[itm->second] = M_PI * double(sizes[itm->second]) / tots * radius / 2;
  }

  /*cout<<" "<<current_angle<<endl;
  cout<<"current angle"<<endl;*/

  return 0;
}

#endif
