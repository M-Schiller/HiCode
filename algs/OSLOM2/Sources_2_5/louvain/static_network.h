#pragma once

#ifndef STATIC_static_network_INCLUDED
#define STATIC_static_network_INCLUDED

#include <utility>
#include "./standard_package/standard_include.cpp"
#include "wsarray.h"

inline void statement()
{
  std::cout << "\nTo run the program type \n<program name> -f <filename> [-c(-l) <filename2>]\n" << std::endl;
}

inline int parse_command_line(bool & value, std::string & s, std::string & s2, int argc, char * argv[])
{
  value = false;

  int _arg_ = 1;
  if (argc <= 1)
  {
    statement();
    return -1;
  }

  s = argv[_arg_];

  if (s == "-f")
  {
    _arg_++;
    s = argv[_arg_];
  }

  else
  {
    std::cout << "ERROR" << std::endl << std::endl;
    statement();
    return -1;
  }

  _arg_++;

  if (argv[_arg_] != nullptr)
    s2 = argv[_arg_];
  else
    return 0;

  if (s2 == "-c"
    || s2 == "-l")
  {
    _arg_++;
    s2 = argv[_arg_];
    value = true;
  }
  else {
    std::cout << "ERROR" << std::endl << std::endl;
    statement();
    return -1;
  }
  return 0;
}

class static_network
{
public:
  static_network() = default;
  static_network(std::deque<std::deque<int>> &, std::deque<int> &);
  static_network(bool, std::deque<std::deque<int>> &, std::deque<std::deque<double>> &, std::deque<int> &);
  static_network(const std::string&);
  static_network(const std::string& file_name, bool & good_file);
  static_network(std::map<int, std::map<int, double>> & A);
  ~static_network();
  void add_isolated(int id_iso);

  int draw(std::string, bool);
  int draw(std::string);
  int draw_consecutive(const std::string&, const std::string&, bool);
  int draw_consecutive(std::string, std::string);

  void print_id(const std::deque<int> & a, std::ostream &);
  void print_id(const std::deque<std::deque<int>> &, std::ostream &);
  void print_id(const std::deque<std::set<int>> &, std::ostream &);
  void print_id(const std::set<int> &, std::ostream &);
  void deque_id(std::deque<int> &);

  int connected(const std::string&);
  void pajek_print_cluster(const std::string& file_name, std::deque<int> group, int number_of_shells);
  void print_subgraph(const std::string&, const std::deque<int> &);
  void print_component(const std::string& file_name, const std::deque<int> & group);
  void set_subgraph(std::deque<int> &, std::deque<std::deque<int>> &, std::deque<std::deque<double>> &);

  void print_random_version(const std::string&);
  void print_connected_components(std::ostream &);
  void same_component(int, std::set<int> &);
  void set_connected_components(std::deque<std::deque<int>> &);
  void set_connected_components(std::deque<std::set<int>> &);

  int translate(std::deque<std::deque<int>> &);
  void get_id_label(std::map <int, int> &);
  int get_degree(std::deque<int> &);

  int size() const { return dim; }
  double edges() const { return tstrength; }

  double kin(const std::deque<int> &);
  double kin(const std::set<int> &);
  double ktot(const std::deque<int> &);
  double ktot(const std::set<int> &);

  void erase_link(int, int);

  double newman_modularity(double kin_g, double ktot_g);
  double newman_modularity(std::set<int> & s);
  double newman_modularity(std::deque<std::set<int>> &);
  double newman_modularity(std::deque<int> & s);
  double newman_modularity(std::deque<std::deque<int>> &);

  int compute_betweeness_single_source(int, std::map<std::pair<int, int>, double> &);
  int component_betweeness(int, std::map<std::pair<int, int>, double> &, std::set<int> &);
  int all_betweeness(std::map<std::pair<int, int>, double> &);

  int id_of(int a) { return vertices[a]->id_num; };
  int set_mem_adj(std::deque<std::deque<int>> &);		// unweighted networks!

  int distances_from_i(int node, double &, std::deque<double> & ML, int);
  int propagate_distances(std::deque<int> & new_shell, std::set<int> & already_gone, std::deque<std::pair<int, int >> & distances_node, int shell, std::deque<double> & ML, int &, int);
  int diameter_and_asp(double & averageall, int & diameter, std::deque<double> & ML);

  int knn(std::map <int, double> &);
  int knn(std::map <int, std::deque<double >> & knn_hist);

  void clustering_coefficient(std::map <int, std::deque<double>> & c_of_k);
  double clustering_coefficient();

  void set_graph(std::map<int, std::map<int, double>> & A);
  void set_graph(std::multimap<int, int > & A);
  bool set_graph(const std::string& file_name);
  void set_graph(
    bool weighted_,
    std::deque<std::deque<int>> & link_per_node,
    std::deque<std::deque<double>> & weights_per_node,
    std::deque<int> & label_rows);
  void set_graph(std::deque<std::deque<int>> & link_per_node, std::deque<int> & label_rows);

  void clear();
  bool GN_bench(int nodes, int modules, double kout, double m_degree, std::deque<std::deque<int>> & ten);
  bool set_binary_random(const std::deque<int> & degrees);
  bool set_random_powlaw(int, double, double, int);
  bool set_random_powlaw_multiple(int num_nodes, double tau, double average_k, int max_degree, bool);
  bool set_multiple_random(const std::deque<int> & degrees, bool);

  void community_netknn(std::deque<std::deque<int>> &);
  void community_net(std::deque<std::deque<int>> & ten);
  double monte_carlo_asp();

protected:
  class  vertex
  {
  public:
    vertex(int, int, int);
    ~vertex();

    double kplus(const std::deque<int> &);
    double kplus(const std::set<int> &);

    int id_num;
    double strength;
    wsarray* links;
  };

  int dim;									// number of nodes
  double tstrength;					// number of links (weighted)
  bool weighted;
  std::deque<vertex*> vertices;

private:
  void propagate_bw(int*, int*, int, std::set<int> &, std::set<int> &, std::deque<int> &);
  void rewiring(std::deque<std::set<int>>& en);
  void next_shell(const std::set<int> & shell0, std::set<int> & shell1, std::set<int> & already);
};

inline static_network::vertex::vertex(int b, int c, int preall)
{
  id_num = b;
  links = new wsarray(preall);
}

inline static_network::vertex::~vertex()
{
  delete links;
  links = nullptr;
}

inline double static_network::vertex::kplus(const std::deque<int> &a)
{
  // computes the internal degree of the vertex respect with a
  double f = 0;
  for (int i : a)
  {
    f += links->posweightof(i).second;
  }
  return f;
}

inline double static_network::vertex::kplus(const std::set<int> &a)
{
  // computes the internal degree of the vertex respect with a (a is supposed to be sorted)
  double f = 0;
  for (int i = 0; i < links->size(); i++)
  {
    if (a.find(links->l[i]) != a.end())
    {
      f += links->w[i];
    }
  }
  return f;
}

inline static_network::static_network(const std::string& file_name, bool & good_file)
{
  good_file = set_graph(file_name);
}

inline static_network::static_network(const std::string& file_name)
{
  set_graph(file_name);
}

inline static_network::static_network(std::deque<std::deque<int>> & link_per_node, std::deque<int> & label_rows)
{
  set_graph(link_per_node, label_rows);
}

inline static_network::static_network(bool weighted_, std::deque<std::deque<int>> & link_per_node, std::deque<std::deque<double>> & weights_per_node, std::deque<int> & label_rows)
{
  set_graph(weighted_, link_per_node, weights_per_node, label_rows);
}

inline static_network::static_network(std::map<int, std::map<int, double>> & A)
{
  set_graph(A);
}

inline static_network::~static_network()
{
  clear();
}

inline void static_network::clear()
{
  for (auto& vertex : vertices)
  {
    delete vertex;
    vertex = nullptr;
  }

  vertices.clear();
  dim = 0;
  tstrength = 0;
}

inline void static_network::set_graph(std::multimap<int, int> & A)
{
  std::map<int, std::map<int, double>> B;

  for (auto& itm : A)
  {
    int m_name = itm.first;
    std::map<int, double> neigh_weight;
    B.emplace(m_name, neigh_weight);
  }

  for (auto &[i, j] : A)
  {
    //cout<<"link: "<<itm->first<<" "<<itm->second<< std::endl;

    B[i].emplace(j, 1.);
    B[j].emplace(i, 1.);
  }

  /*for(map<int, map<int, double>>:: iterator itm= B.begin(); itm!=B.end(); itm++)
    prints(itm->second);

  cout<<"ok"<< std::endl;*/
  set_graph(B);
}

inline void static_network::set_graph(std::map<int, std::map<int, double>> & A)
{
  // this maps the id into the usual stuff neighbors- weights
  std::deque<std::deque<int>> link_per_node;
  std::deque<std::deque<double>> weights_per_node;
  std::deque<int> label_rows;

  for (auto itm = A.begin(); itm != A.end(); ++itm)
  {
    label_rows.push_back(itm->first);
    std::deque<int> n;
    std::deque<double> w;

    for (auto itm2 = itm->second.begin(); itm2 != itm->second.end(); ++itm2)
    {
      if (itm2->second > 0)
      {
        n.push_back(itm2->first);
        w.push_back(itm2->second);
      }
    }

    link_per_node.push_back(n);
    weights_per_node.push_back(w);
  }

  /*
  prints(label_rows);
  printm(link_per_node);
  printm(weights_per_node);*/

  set_graph(true, link_per_node, weights_per_node, label_rows);
}

inline bool static_network::set_graph(const std::string& file_name)
{
  clear();

  weighted = false;

  std::map<int, int> newlabels;
  std::deque<int> link_i;

  bool good_file = true;

  {
    int label = 0;
    std::ifstream inb(file_name);
    std::string ins;

    while (getline(inb, ins))
    {
      if (!ins.empty()
        && ins[0] != '#')
      {
        std::deque<double> ds;
        cast_string_to_doubles(ins, ds);

        if (ds.size() < 2)
        {
          std::cerr << "From file " << file_name << ": string not readable " << ins << " " << std::endl;
          good_file = false;
          break;
        }

        int innum1 = cast_int(ds[0]);
        int innum2 = cast_int(ds[1]);

        auto itf = newlabels.find(innum1);
        if (itf == newlabels.end())
        {
          newlabels.emplace(innum1, label++);
          link_i.push_back(1);
        }
        else
        {
          link_i[itf->second]++;
        }

        itf = newlabels.find(innum2);
        if (itf == newlabels.end())
        {
          newlabels.emplace(innum2, label++);
          link_i.push_back(1);
        }
        else
          link_i[itf->second]++;
      }
    }
  }

  dim = newlabels.size();

  for (int i = 0; i < dim; i++)
    vertices.push_back(new vertex(0, 0, link_i[i]));

  for (auto itm = newlabels.begin(); itm != newlabels.end(); ++itm)
  {
    vertices[itm->second]->id_num = itm->first;
  }

  if (good_file)
  {
    std::ifstream inb(file_name);
    std::string ins;
    while (getline(inb, ins))
    {
      if (!ins.empty()
        && ins[0] != '#')
      {
        std::deque<double> ds;
        cast_string_to_doubles(ins, ds);

        int innum1 = cast_int(ds[0]);
        int innum2 = cast_int(ds[1]);

        double w = 1;
        if (ds.size() > 2)
        {
          if (ds[2] <= 0)
          {
            w = 1;
          }
          else
          {
            weighted = true;
            w = ds[2];
          }
        }
        int new1 = newlabels[innum1];
        int new2 = newlabels[innum2];

        if (new1 != new2)		// no self loops!
        {
          vertices[new1]->links->push_back(new2, w);
          vertices[new2]->links->push_back(new1, w);
        }
      }
    }

    tstrength = 0;

    for (int i = 0; i < dim; i++)
    {
      vertices[i]->links->freeze();
      double strength_i = 0;
      for (int j = 0; j < vertices[i]->links->size(); j++)
      {
        strength_i += vertices[i]->links->w[j];
      }
      vertices[i]->strength = strength_i;
      tstrength += strength_i;
    }
    tstrength = tstrength / 2.;
  }
  else
  {
    std::cerr << "File corrupted" << std::endl;
  }
  return good_file;
}

inline void static_network::set_graph(
  bool weighted_,
  std::deque<std::deque<int>>& link_per_node,
  std::deque<std::deque<double>>& weights_per_node,
  std::deque<int>& label_rows)
{
  clear();

  // there is no check between label_rows and link per node but they need to have the same labels
  // link_per_node and weights_per_node are the list of links and weights. label_rows[i] is the label corresponding to row i

  weighted = weighted_;

  std::map<int, int> newlabels;		// this maps the old labels with the new one
  for (int& label_row : label_rows)
  {
    newlabels.emplace(label_row, newlabels.size());
  }

  dim = newlabels.size();

  for (int i = 0; i < dim; i++)
  {
    vertices.push_back(new vertex(0, 0, link_per_node[i].size()));
  }

  for (auto &[new_label, old_label] : newlabels)
  {
    vertices[old_label]->id_num = new_label;
  }

  for (unsigned i = 0; i < link_per_node.size(); i++)
  {
    for (unsigned j = 0; j < link_per_node[i].size(); j++)
    {
      int new2 = newlabels[link_per_node[i][j]];
      vertices[i]->links->push_back(new2, weights_per_node[i][j]);
    }
  }

  tstrength = 0;

  for (int i = 0; i < dim; i++)
  {
    vertices[i]->links->freeze();

    double strength_i = 0;
    for (int j = 0; j < vertices[i]->links->size(); j++)
      strength_i += vertices[i]->links->w[j];

    vertices[i]->strength = strength_i;
    tstrength += strength_i;
  }
  tstrength = tstrength / 2.;
}

inline void static_network::set_graph(std::deque<std::deque<int>> & link_per_node, std::deque<int> & label_rows)
{
  // label rows contains the label of the subgraph. The very same labels must be in link_per_node

  // there is no check between label_rows and link per node but they need to have the same labels
  // link_per_node and weights_per_node are the list of links and weights. label_rows[i] is the label corresponding to row i
  clear();

  weighted = false;

  std::map<int, int> newlabels;		// this maps the old labels into the new ones
  for (int& label_row : label_rows)
  {
    newlabels.emplace(label_row, newlabels.size());
  }

  dim = newlabels.size();

  for (int i = 0; i < dim; i++)
  {
    vertices.push_back(new vertex(0, 0, link_per_node[i].size()));
  }

  for (auto &[old_label, new_label] : newlabels)
  {
    vertices[new_label]->id_num = old_label;
  }

  for (int i = 0; i < link_per_node.size(); i++)
  {
    for (int j = 0; j < link_per_node[i].size(); j++)
    {
      vertices[i]->links->push_back(newlabels[link_per_node[i][j]], 1.);
    }
  }

  tstrength = 0;
  for (int i = 0; i < dim; i++)
  {
    vertices[i]->links->freeze();

    double strength_i = 0;
    for (int j = 0; j < vertices[i]->links->size(); j++)
      strength_i += vertices[i]->links->w[j];

    vertices[i]->strength = strength_i;
    tstrength += strength_i;
  }
  tstrength = tstrength / 2.;
}

inline int static_network::get_degree(std::deque<int> &d)
{
  d.clear();
  for (int i = 0; i < dim; i++)
    d.push_back(vertices[i]->links->size());

  return 0;
}

inline double static_network::newman_modularity(const double kin_g, const double ktot_g)
{
  return ((kin_g) / (2.*tstrength) - pow((ktot_g) / (2.*tstrength), 2));
}

inline double static_network::newman_modularity(std::set<int> & s)
{
  return newman_modularity(kin(s), ktot(s));
}

inline double static_network::newman_modularity(std::deque<std::set<int>> & Comps)
{
  double mm = 0;
  for (auto& Comp : Comps)
  {
    mm += newman_modularity(Comp);
  }
  return mm;
}

inline double static_network::newman_modularity(std::deque<int> & s)
{
  return newman_modularity(kin(s), ktot(s));
}

inline double static_network::newman_modularity(std::deque<std::deque<int>> & Comps)
{
  double mm = 0;
  for (auto& Comp : Comps)
  {
    mm += newman_modularity(Comp);
  }
  return mm;
}

inline void static_network::erase_link(int a, int b)
{
  double sa = vertices[a]->links->posweightof(b).second;

  vertices[a]->links->erase(b);
  vertices[b]->links->erase(a);

  vertices[a]->strength -= sa;
  vertices[b]->strength -= sa;

  tstrength -= sa;
}

inline int static_network::connected(const std::string& str)
{
  int spannet = 0;
  std::deque <std::set<int>> partic;
  std::deque <int> present;
  present.assign(dim, 0);

  while (spannet != dim) {
    std::set <int> connected;
    std::set <int> newcon;

    for (int i = 0; i < dim; i++)
    {
      if (present[i] == 0) {
        connected.insert(i);
        newcon.insert(i);
        present[i] = 1;
        break;
      }
    }

    while (!newcon.empty())
    {
      std::set <int> nnewcon = newcon;
      newcon.clear();

      for (auto it = nnewcon.begin(); it != nnewcon.end(); ++it)
      {
        for (int near = 0; near != vertices[*it]->links->size(); ++near)
        {
          present[*it] = 1;
          if (connected.insert(vertices[*it]->links->l[near]).second)
          {
            newcon.insert(vertices[*it]->links->l[near]);
          }
        }
      }
    }

    partic.push_back(connected);
    spannet += connected.size();
  }

  std::ofstream con_out(str);

  //cout<<"number of connected components = "<<partic.size()<< std::endl;
  //cout<<"dimensions"<< std::endl;

  int max_pcon = 0;

  for (int i = 0; i < partic.size(); i++)
  {
    //cout<<partic[i].size()<<"   ";

    if (partic[i].size() >= partic[max_pcon].size())
      max_pcon = i;
  }

  //cout<< std::endl<< std::endl;

  for (auto it = partic[max_pcon].begin(); it != partic[max_pcon].end(); ++it)
  {
    for (int j = 0; j < vertices[*it]->links->size(); j++) if (vertices[*it]->id_num < vertices[vertices[*it]->links->l[j]]->id_num)
    {
      con_out << vertices[*it]->id_num << "\t" << vertices[vertices[*it]->links->l[j]]->id_num << "\t" << vertices[*it]->links->w[j] << std::endl;
    }
  }
  return (partic.size());
}

inline double static_network::kin(const std::deque<int> & seq)
{
  if (seq.size() > 2 * tstrength / dim)
  {
    std::set<int> H;
    deque_to_set(seq, H);
    return kin(H);
  }
  double k = 0;
  for (int i = 0; i < seq.size(); i++)
  {
    k += vertices[seq[i]]->kplus(seq);
  }

  return k;
}

inline double static_network::ktot(const std::deque<int> &seq)
{
  double k = 0;
  for (int i : seq)
    k += vertices[i]->strength;
  return k;
}

inline double static_network::ktot(const std::set <int> &s)
{
  double k = 0;
  for (auto it : s)
  {
    k += vertices[it]->strength;
  }
  return k;
}

inline double static_network::kin(const std::set <int> &s)
{
  double k = 0;
  for (auto it : s)
  {
    k += vertices[it]->kplus(s);
  }
  return k;
}

inline int static_network::draw(std::string file_name)
{
  return draw(std::move(file_name), weighted);
}

inline int static_network::draw_consecutive(std::string file_name1, std::string file_name2)
{
  return draw_consecutive(file_name1, file_name2, weighted);
}

inline int static_network::draw_consecutive(
  const std::string& file_name1,
  const std::string&
  file_name2,
  bool _weighted_)
{
  //cout<<"drawing in file "<<b<< std::endl;
  std::ofstream graph_out(file_name1);

  if (_weighted_)
  {
    for (int i = 0; i < vertices.size(); i++)
    {
      for (int j = 0; j < vertices[i]->links->size(); j++) if (vertices[i]->id_num <= vertices[vertices[i]->links->l[j]]->id_num)
      {
        graph_out << i << "\t" << vertices[i]->links->l[j] << "\t" << cast_int(vertices[i]->links->w[j]) << std::endl;
      }
    }
  }
  else
  {
    for (int i = 0; i < vertices.size(); i++)
    {
      for (int j = 0; j < vertices[i]->links->size(); j++) if (vertices[i]->id_num <= vertices[vertices[i]->links->l[j]]->id_num)
      {
        graph_out << i << "\t" << vertices[i]->links->l[j] << std::endl;
      }
    }
  }

  std::ofstream graph_out2(file_name2);
  for (int i = 0; i < vertices.size(); i++)
    graph_out2 << i << " " << vertices[i]->id_num << std::endl;

  return 0;
}

inline int static_network::draw(std::string file_name, bool _weighted_)
{
  int h = file_name.size();

  char b[h + 1];
  for (int i = 0; i < h; i++)
    b[i] = file_name[i];
  b[h] = '\0';

  std::ofstream graph_out(b);

  if (_weighted_)
  {
    for (auto& vertice : vertices)
    {
      for (int j = 0; j < vertice->links->size(); j++) /*if(vertices[i]->id_num <= vertices[vertices[i]->links->l[j]]->id_num)*/
      {
        graph_out << vertice->id_num << "\t" << vertices[vertice->links->l[j]]->id_num << "\t" <<
          vertice->links->w[j] << std::endl;
      }
    }
  }
  else
  {
    for (auto& vertice : vertices)
    {
      for (int j = 0; j < vertice->links->size(); j++) /*if(vertices[i]->id_num <= vertices[vertices[i]->links->l[j]]->id_num)*/
      {
        graph_out << vertice->id_num << "\t" << vertices[vertice->links->l[j]]->id_num << std::endl;
      }
    }
  }
  return 0;
}

inline void static_network::rewiring(std::deque<std::set<int>>& en)
{
  std::deque<std::pair<int, int>> Es;

  for (int i = 0; i < en.size(); i++) for (auto its = en[i].begin(); its != en[i].end(); ++its)
  {
    if (i < *its)
    {
      Es.emplace_back(i, *its);
    }
  }

  for (int kk = 0; kk < 10; kk++) for (int i = 0; i < Es.size(); i++)
  {
    int a = irand(Es.size() - 1);
    int b = irand(Es.size() - 1);

    while (true) {
      if (a == b)
        break;

      int a1 = Es[a].first;
      int a2 = Es[a].second;

      int b1, b2;

      if (ran4() < 0.5)
      {
        b1 = Es[b].first;
        b2 = Es[b].second;
      }
      else
      {
        b1 = Es[b].second;
        b2 = Es[b].first;
      }

      //cout<<a1<<" "<<a2<<"-- "<<b1<<" "<<b2<< std::endl;

      if (!(en[a1].find(b1) == en[a1].end() && en[a2].find(b2) == en[a2].end() && a1 != b1 && a2 != b2))
        break;

      en[a1].erase(a2);
      en[a2].erase(a1);
      en[b1].erase(b2);
      en[b2].erase(b1);

      en[a1].insert(b1);
      en[a2].insert(b2);
      en[b1].insert(a1);
      en[b2].insert(a2);

      Es[a].first = a1;
      Es[a].second = b1;
      Es[b].first = a2;
      Es[b].second = b2;

      break;
    }
  }
}

inline void static_network::print_random_version(const std::string& file_name)
{
  std::ofstream outt(file_name);

  std::deque<int> nodes;
  std::deque<int> degrees;

  for (int i = 0; i < dim; i++)
  {
    nodes.push_back(vertices[i]->id_num);
  }

  for (int i = 0; i < dim; i++)
  {
    degrees.push_back(cast_int(vertices[i]->strength));
  }
  /*
  std::cout<<"nodes"<< std::endl;
  prints(nodes);

  std::cout<<"degrees"<< std::endl;
  prints(degrees);
  //*/

  // this function is to build a network with the labels stored in nodes and the degree seq in degrees (correspondence is based on the vectorial index)

  // labels will be placed in the end
  std::deque<std::set<int>> en; // this is the E of the subgraph
  {
    std::set<int> first;
    for (int i = 0; i < nodes.size(); i++)
      en.push_back(first);
  }

  std::multimap <int, int> degree_node;

  for (int i = 0; i < degrees.size(); i++)
  {
    degree_node.insert(degree_node.end(), std::make_pair(degrees[i], i));
  }

  int var = 0;

  while (!degree_node.empty())
  {
    auto itlast = degree_node.end();
    --itlast;

    auto itit = itlast;
    std::deque <std::multimap<int, int>::iterator> erasenda;

    int inserted = 0;

    for (int i = 0; i < itlast->first; i++)
    {
      if (itit != degree_node.begin())
      {
        --itit;

        en[itlast->second].insert(itit->second);
        en[itit->second].insert(itlast->second);
        inserted++;

        erasenda.push_back(itit);
      }
      else
        break;
    }

    for (auto& i : erasenda)
    {
      if (i->first > 1)
        degree_node.emplace(i->first - 1, i->second);

      degree_node.erase(i);
    }
    var += itlast->first - inserted;
    degree_node.erase(itlast);
  }

  // this is to randomize the subgraph ------------------------------------------------------------------
  rewiring(en);

  for (int i = 0; i < en.size(); i++)
  {
    for (auto it = en[i].begin(); it != en[i].end(); ++it)
    {
      outt << vertices[i]->id_num << " " << vertices[*it]->id_num << std::endl;
    }
  }
}

inline void static_network::get_id_label(std::map <int, int> &a)
{
  for (int i = 0; i < dim; i++)
  {
    a.emplace(vertices[i]->id_num, i);
  }
}

inline void static_network::deque_id(std::deque<int> & a)
{
  for (int& i : a)
  {
    i = vertices[i]->id_num;
  }
}

inline void static_network::print_id(const std::deque<int> & a, std::ostream & pout)
{
  for (int i : a)
    pout << vertices[i]->id_num << "\t";
  pout << std::endl;
}

inline void static_network::print_id(const std::set<int> & a, std::ostream & pout)
{
  for (auto its = a.begin(); its != a.end(); ++its)
    pout << vertices[*its]->id_num << "\t";
  pout << std::endl;
}

void static_network::print_id(const std::deque<std::deque<int>> & a, std::ostream & pout)
{
  for (const auto& i : a)
    print_id(i, pout);
}

inline void static_network::print_id(const std::deque<std::set<int>> & a, std::ostream & pout)
{
  for (const auto& i : a)
    print_id(i, pout);
}

inline int static_network::translate(std::deque<std::deque<int>> & ten)
{
  std::map<int, int> A;
  get_id_label(A);

  for (auto& i : ten)
  {
    for (int j = 0; j < i.size(); j++)
    {
      const auto itf = A.find(i[j]);
      if (itf == A.end())
      {
        std::cerr << "ERROR: the nodes in the communities are different from those ones in the network!" << std::endl;
        return -1;
      }
      i[j] = itf->second;
    }
    std::sort(i.begin(), i.end());
  }
  return 0;
}

inline void static_network::print_component(const std::string& file_name, const std::deque<int> & group)
{
  std::ofstream subout(file_name);

  for (int nodei : group)
  {
    for (int j = 0; j < vertices[nodei]->links->size(); j++)
      subout << vertices[nodei]->id_num << " " << vertices[vertices[nodei]->links->l[j]]->id_num << " " << vertices[nodei]->links->w[j] << std::endl;
  }
}

inline void static_network::print_subgraph(const std::string& file_name, const std::deque<int> & group)
{
  std::ofstream subout(file_name);

  for (int i = 0; i < group.size(); i++)
  {
    int nodei = group[i];

    for (int j : group)
    {
      double wij = vertices[nodei]->links->posweightof(j).second;
      if (wij > 0 && vertices[nodei]->id_num < vertices[j]->id_num)
      {
        subout << vertices[nodei]->id_num << " " << vertices[j]->id_num << " " << wij << std::endl;
      }
    }
  }
}

inline void static_network::next_shell(const std::set<int> & shell0, std::set<int> & shell1, std::set<int> & already)
{
  shell1.clear();
  for (auto its = shell0.begin(); its != shell0.end(); ++its)
  {
    for (int j = 0; j < vertices[*its]->links->size(); j++) if (already.find(vertices[*its]->links->l[j]) == already.end())
    {
      shell1.insert(vertices[*its]->links->l[j]);
    }
  }
}

inline void static_network::pajek_print_cluster(const std::string& file_name, std::deque<int> group, int number_of_shells)
{
  std::ofstream subout(file_name);

  std::deque<std::string> pajek_colors;
  pajek_colors.emplace_back("Red");
  pajek_colors.emplace_back("Green");
  pajek_colors.emplace_back("Black");
  pajek_colors.emplace_back("Yellow");
  pajek_colors.emplace_back("White");
  pajek_colors.emplace_back("Yellow");
  pajek_colors.emplace_back("Brown");
  pajek_colors.emplace_back("Mahogany");
  pajek_colors.emplace_back("Grey");
  pajek_colors.emplace_back("Blue");

  std::map<int, int> shell;
  std::set<int> shell0;
  std::set<int> already;
  deque_to_set(group, shell0);

  int nsh = 0;
  for (auto its = shell0.begin(); its != shell0.end(); ++its)
  {
    shell[vertices[*its]->id_num] = nsh;
    already.insert(*its);
  }

  for (int i = 0; i < number_of_shells; i++)
  {
    std::set<int> shell1;
    next_shell(shell0, shell1, already);
    nsh++;

    shell0 = shell1;
    for (auto its = shell0.begin(); its != shell0.end(); ++its)
    {
      shell[vertices[*its]->id_num] = nsh;
      group.push_back(*its);
      already.insert(*its);
    }
  }

  std::map<int, int> node_plab;
  int pl = 0;
  for (int i : group)
  {
    if (node_plab.find(vertices[i]->id_num) == node_plab.end())
    {
      node_plab[vertices[i]->id_num] = ++pl;
    }
  }

  std::map<int, int> inverse;
  for (auto itm = node_plab.begin(); itm != node_plab.end(); ++itm)
    inverse[itm->second] = itm->first;

  subout << "*Vertices " << node_plab.size() << std::endl;
  for (auto itm = inverse.begin(); itm != inverse.end(); ++itm)
  {
    subout << "   " << itm->first << " \"" << itm->second << "\" box ic " << pajek_colors[shell[itm->second] % pajek_colors.size()] << std::endl;
  }

  subout << "*Edges" << std::endl;
  pl = group.size();

  for (int i : group)
  {
    for (int j = 0; j < vertices[i]->links->size(); j++)
    {
      if (vertices[i]->id_num < vertices[vertices[i]->links->l[j]]->id_num
        && node_plab.find(vertices[vertices[i]->links->l[j]]->id_num) != node_plab.end())
      {
        subout << node_plab[vertices[i]->id_num] << " "
          << node_plab[vertices[vertices[i]->links->l[j]]->id_num] << " "
          << vertices[i]->links->w[j] << std::endl;
      }
    }
  }

  char pbb[1000];
  sprintf(pbb, "unixdos %s %s.net", file_name.c_str(), file_name.c_str());
  systemCall(std::string(pbb));
}

inline void static_network::set_subgraph(std::deque<int> & group, std::deque<std::deque<int>> & link_per_node, std::deque<std::deque<double>> & weights_per_node)
{
  // in this function I'm not using id's... because I want to work with the same labels (don't want to translate)
  std::sort(group.begin(), group.end());

  weights_per_node.clear();
  link_per_node.clear();

  for (int node_i : group)
  {
    std::deque<int> link_i;
    std::deque<double> weight_i;

    for (int j = 0; j < vertices[node_i]->links->size(); j++) if (binary_search(group.begin(), group.end(), vertices[node_i]->links->l[j]))
    {
      link_i.push_back(vertices[node_i]->links->l[j]);
      weight_i.push_back(vertices[node_i]->links->w[j]);
    }

    link_per_node.push_back(link_i);
    weights_per_node.push_back(weight_i);
  }
}

inline void static_network::same_component(int source, std::set<int> & already_gone)
{
  already_gone.clear();

  int node = source;
  double  average;
  std::deque<double>  ML;
  int step;

  already_gone.insert(node);

  std::deque<int> new_shell;
  new_shell.push_back(node);

  std::deque<std::pair<int, int>> distances_node;
  int shell = 0;
  int reached = 1;

  if (shell >= ML.size())
  {
    ML.push_back(0);
  }

  ML[shell] += reached;

  propagate_distances(new_shell, already_gone, distances_node, shell, ML, reached, step);
}

inline void static_network::set_connected_components(std::deque<std::deque<int>> & comps)
{
  comps.clear();
  std::set<int> not_assigned;
  for (int i = 0; i < dim; i++)
    not_assigned.insert(i);

  while (!not_assigned.empty())
  {
    int source = *not_assigned.begin();

    std::set<int> mates;
    same_component(source, mates);

    std::deque<int> ccc;
    for (auto its = mates.begin(); its != mates.end(); ++its)
    {
      ccc.push_back(*its);
      not_assigned.erase(*its);
    }
    comps.push_back(ccc);
  }
}

inline void static_network::set_connected_components(std::deque<std::set<int>> & comps)
{
  comps.clear();
  std::set<int> not_assigned;
  for (int i = 0; i < dim; i++)
    not_assigned.insert(i);

  while (!not_assigned.empty())
  {
    int source = *not_assigned.begin();

    std::set<int> mates;
    same_component(source, mates);

    for (auto its = mates.begin(); its != mates.end(); ++its)
    {
      not_assigned.erase(*its);
    }

    comps.push_back(mates);
  }
}

inline void static_network::print_connected_components(std::ostream & outb)
{
  int spannet = 0;
  std::deque <std::set<int>> partic;
  std::deque <int> present;
  present.assign(dim, 0);

  while (spannet != dim)
  {
    std::set <int> connected;
    std::set <int> newcon;

    for (int i = 0; i < dim; i++)
    {
      if (present[i] == 0)
      {
        connected.insert(i);
        newcon.insert(i);
        present[i] = 1;
        break;
      }
    }

    while (!newcon.empty())
    {
      std::set <int> nnewcon = newcon;
      newcon.clear();

      for (auto it = nnewcon.begin(); it != nnewcon.end(); ++it)
      {
        for (int near = 0; near != vertices[*it]->links->size(); ++near)
        {
          present[*it] = 1;
          if (connected.insert(vertices[*it]->links->l[near]).second)
            newcon.insert(vertices[*it]->links->l[near]);
        }
      }
    }

    partic.push_back(connected);
    spannet += connected.size();
  }

  for (auto& i : partic)
  {
    for (auto its = i.begin(); its != i.end(); ++its)
    {
      outb << vertices[*its]->id_num << " ";
    }
    outb << std::endl;
  }
}

inline void static_network::propagate_bw(int* distances, int* weights, int source, std::set<int> & mates, std::set<int> & not_leaves, std::deque<int> & next_shell)
{
  mates.insert(source);

  for (int i = 0; i < vertices[source]->links->size(); i++)
  {
    int neigh = vertices[source]->links->l[i];

    if (mates.find(neigh) == mates.end())		// new vertex
    {
      //cout<<"new "<<vertices[neigh]->id_num<< std::endl;
      distances[neigh] = distances[source] + 1;
      weights[neigh] = weights[source];
      mates.insert(neigh);
      next_shell.push_back(neigh);
      not_leaves.insert(source);
    }
    else if (distances[neigh] == distances[source] + 1)
    {
      weights[neigh] += weights[source];
      not_leaves.insert(source);
    }
  }
}

inline int static_network::compute_betweeness_single_source(int source, std::map<std::pair<int, int>, double> & tot_edge_bw)
{
  // this function compute the edge betweenness of all the vertices in the component of source (only for source)
  std::set<int> mates;
  std::set<int> not_leaves;

  int distances[dim];
  int weights[dim];

  distances[source] = 0;
  weights[source] = 1;
  std::deque<int> present_shell;
  present_shell.push_back(source);

  while (true)
  {
    if (present_shell.empty())
      break;

    std::deque<int> next_shell;

    for (int i : present_shell)
      propagate_bw(distances, weights, i, mates, not_leaves, next_shell);

    //cout<<"------------ next "<< std::endl;
    //print_id(next_shell, std::cout);

    present_shell = next_shell;
  }

  /*
  prints(distances, dim);
  prints(weights, dim);
  */

  //print_id(mates, std::cout);

  std::deque<int> leaves;
  for (auto its = mates.begin(); its != mates.end(); ++its)
  {
    if (not_leaves.find(*its) == not_leaves.end())
    {
      leaves.push_back(*its);
    }
  }

  //cout<<"leaves"<< std::endl;
  //print_id(leaves, std::cout);

  std::map<std::pair<int, int>, double> edge_bw;	// map edge-betweenness
  for (auto mate : mates)
  {
    for (int i = 0; i < vertices[mate]->links->size(); i++)
    {
      if (mate < vertices[mate]->links->l[i])
      {
        edge_bw.emplace(std::make_pair(mate, vertices[mate]->links->l[i]), 0);
      }
    }
  }

  std::multimap<int, int> distance_not_leaves;
  for (auto not_leave : not_leaves)
  {
    distance_not_leaves.emplace(-distances[not_leave], not_leave);
  }

  for (auto its = leaves.begin(); its != leaves.end(); ++its)
  {
    for (int i = 0; i < vertices[*its]->links->size(); i++)
    {
      if (distances[*its] > distances[vertices[*its]->links->l[i]])
      {
        std::pair<int, int> ed;
        ed.first = std::min(*its, vertices[*its]->links->l[i]);
        ed.second = std::max(*its, vertices[*its]->links->l[i]);

        edge_bw[ed] = double(weights[vertices[*its]->links->l[i]]) / weights[*its];
        tot_edge_bw[ed] += double(weights[vertices[*its]->links->l[i]]) / weights[*its];
      }
    }
  }

  for (auto itm = distance_not_leaves.begin(); itm != distance_not_leaves.end(); ++itm)
  {
    //cout<<"node:::: "<<itm->second<< std::endl;

    double sum_of_weight = 0;

    for (int i = 0; i < vertices[itm->second]->links->size(); i++) {
      std::pair<int, int> ed;
      ed.first = std::min(itm->second, vertices[itm->second]->links->l[i]);
      ed.second = std::max(itm->second, vertices[itm->second]->links->l[i]);

      sum_of_weight += edge_bw[ed];
    }

    for (int i = 0; i < vertices[itm->second]->links->size(); i++)
    {
      if (distances[itm->second] > distances[vertices[itm->second]->links->l[i]])
      {
        std::pair<int, int> ed;
        ed.first = std::min(itm->second, vertices[itm->second]->links->l[i]);
        ed.second = std::max(itm->second, vertices[itm->second]->links->l[i]);

        edge_bw[ed] = double(weights[vertices[itm->second]->links->l[i]]) / weights[itm->second] * (1 + sum_of_weight);
        tot_edge_bw[ed] += double(weights[vertices[itm->second]->links->l[i]]) / weights[itm->second] * (1 + sum_of_weight);

        //cout<<"pred--> "<<vertices[itm->second]->links->l[i]<<" "<<double(weights[vertices[itm->second]->links->l[i]])/weights[itm->second]*(1 + sum_of_weight)<<" "<<sum_of_weight<< std::endl;
      }
    }
  }
  //cout<<"************************"<< std::endl;
  return 0;
}

inline int static_network::component_betweeness(int source, std::map<std::pair<int, int>, double> & tot_edge_bw, std::set<int> & mates)		// this compute the betweenness of the edges in the component of source
{
  mates.clear();

  same_component(source, mates);

  for (auto its = mates.begin(); its != mates.end(); ++its)
  {
    for (int j = 0; j < vertices[*its]->links->size(); j++) if (*its < vertices[*its]->links->l[j])
    {
      std::pair<int, int> ed;
      ed.first = *its;
      ed.second = vertices[*its]->links->l[j];

      tot_edge_bw[ed] = 0;
    }
  }

  for (auto its = mates.begin(); its != mates.end(); ++its)		// this must be made n times
  {
    compute_betweeness_single_source(*its, tot_edge_bw);						              // this requires m operations
  }

  //for(map<pair<int, int>, double>::iterator itm=tot_edge_bw.begin(); itm!=tot_edge_bw.end(); itm++)
    //cout<<vertices[itm->first.first]->id_num<<" "<<vertices[itm->first.second]->id_num<<" "<<itm->second<< std::endl;

  return 0;
}

inline int static_network::all_betweeness(std::map<std::pair<int, int>, double> & tot_edge_bw)
{
  // this compute the betweenness of all the edges
  std::set<int> not_assigned;
  for (int i = 0; i < dim; i++)
  {
    not_assigned.insert(i);
  }

  while (!not_assigned.empty())
  {
    int source = *not_assigned.begin();
    std::set<int> mates;

    component_betweeness(source, tot_edge_bw, mates);

    for (std::set<int>::iterator its = mates.begin(); its != mates.end(); ++its)
    {
      not_assigned.erase(*its);
    }
  }
  //for(map<pair<int, int>, double>::iterator itm=tot_edge_bw.begin(); itm!=tot_edge_bw.end(); itm++)
    //cout<<vertices[itm->first.first]->id_num<<" "<<vertices[itm->first.second]->id_num<<" "<<itm->second<< std::endl;
  return 0;
}

inline int static_network::set_mem_adj(std::deque<std::deque<int>> & mem_adj)		// unweighted networks!
{
  mem_adj.clear();

  for (int i = 0; i < dim; i++)
  {
    std::deque <int> first(dim);
    for (int j = 0; j < dim; j++)
    {
      first[j] = 0;
    }

    for (int j = 0; j < vertices[i]->links->size(); j++)
    {
      first[vertices[i]->links->l[j]] = 1;
    }

    mem_adj.push_back(first);
  }
  return 0;
}

inline int static_network::propagate_distances(
  std::deque<int> & new_shell,
  std::set<int> & already_gone,
  std::deque<std::pair<int, int>> & distances_node,
  int shell, std::deque<double> & ML,
  int & reached,
  int step)
{
  shell++;
  std::deque<int> next_shell;

  for (int i : new_shell)
    for (int j = 0; j < vertices[i]->links->size(); j++)
    {
      if (already_gone.insert(vertices[i]->links->l[j]).second)
      {
        distances_node.emplace_back(shell, vertices[i]->links->l[j]);
        next_shell.push_back(vertices[i]->links->l[j]);
      }
    }

  /*
  std::cout<<"new shell "<<shell<< std::endl;
  print_id(next_shell, std::cout);
  prints(ML);
  */

  if (!next_shell.empty())
  {
    if (shell >= ML.size())
      ML.push_back(dim*step);

    reached += next_shell.size();
    ML[shell] += reached;

    return propagate_distances(next_shell, already_gone, distances_node, shell, ML, reached, step);
  }
  return 0;
}

inline int static_network::distances_from_i(int node, double & average, std::deque<double> & ML, int step)
{
  std::set<int> already_gone;
  already_gone.insert(node);

  std::deque<int> new_shell;
  new_shell.push_back(node);

  std::deque<std::pair<int, int>> distances_node;
  int shell = 0;
  int reached = 1;

  if (shell >= ML.size())
    ML.push_back(0);

  ML[shell] += reached;

  propagate_distances(new_shell, already_gone, distances_node, shell, ML, reached, step);
  //sort(distances_node.begin(), distances_node.end());
  //cout<<"NODE "<<id_of(node)<< std::endl;

  average = 0;
  int maxd = 0;

  for (auto& i : distances_node)
  {
    average += i.first;

    if (i.first > maxd)
      maxd = i.first;

    //cout<<distances_node[i].first<<" "<<id_of(distances_node[i].second)<< std::endl;
  }

  for (int i = maxd + 1; i < ML.size(); i++)
  {
    ML[i] += dim;
  }

  //prints(ML);
  return maxd;
}

inline int static_network::diameter_and_asp(double & averageall, int & diameter, std::deque<double> & ML)
{
  ML.clear();
  averageall = 0;
  diameter = 0;

  for (int i = 0; i < dim; i++) {
    double average;

    //int previousMLsize=ML.size();

    int m = distances_from_i(i, average, ML, i);

    /*if(ML.size()>previousMLsize) {
      for(int j=previousMLsize; j<ML.size(); j++)
        ML[j]+=dim*(i);
    }*/

    if (m > diameter)
      diameter = m;

    averageall += average;
  }

  averageall /= (dim*(dim - 1));

  for (double& i : ML)
    i /= dim;
  //prints(ML);
  return 0;
}

inline int static_network::knn(std::map <int, std::deque<double>> & knn_hist)
{
  // this knn is not erased...
  std::deque<double> knn_nodo(dim);
  for (int i = 0; i < dim; i++)
  {
    knn_nodo[i] = 0;
  }

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < vertices[i]->links->size(); j++)
      knn_nodo[i] += vertices[vertices[i]->links->l[j]]->strength;
  }

  for (int i = 0; i < dim; i++) if (vertices[i]->strength > 0)
  {
    auto it = knn_hist.find(vertices[i]->links->size());
    if (it == knn_hist.end()) {
      std::deque<double> f;
      knn_hist[vertices[i]->links->size()] = f;
    }
    knn_hist[vertices[i]->links->size()].push_back(knn_nodo[i] / vertices[i]->strength);
  }
  return 0;
}

inline int static_network::knn(std::map <int, double> & knn_hist)
{
  double knn_nodo[dim];
  for (int i = 0; i < dim; i++)
  {
    knn_nodo[i] = 0;
  }

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < vertices[i]->links->size(); j++)
    {
      knn_nodo[i] += vertices[vertices[i]->links->l[j]]->strength;
    }
  }

  knn_hist.clear();
  std::map <int, int> degree_hist;

  for (int i = 0; i < dim; i++) if (vertices[i]->strength > 0)
  {
    std::map <int, double>::iterator itf = knn_hist.find(vertices[i]->links->size());
    std::map <int, int>::iterator itf2 = degree_hist.find(vertices[i]->links->size());
    if (itf != knn_hist.end())
    {
      itf->second += knn_nodo[i];
      itf2->second += vertices[i]->links->size();
    }
    else
    {
      knn_hist.emplace(vertices[i]->links->size(), knn_nodo[i]);
      degree_hist.emplace(vertices[i]->links->size(), vertices[i]->links->size());
    }
  }

  std::map <int, int>::iterator it2 = degree_hist.begin();

  for (auto& it1 : knn_hist)
  {
    it1.second = it1.second / (it2->second);
    ++it2;
  }
  /*
  double k2av=0;
  for(int i=0; i<dim; i++)
    k2av+=(vertices[i]->strength)*(vertices[i]->strength);

  std::cout<<"VALUE "<<k2av/(2*tstrength)<< std::endl;
  */
  return 0;
}

inline double static_network::clustering_coefficient()
{
  double number_of_triangles = 0;
  double possible = 0;

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < vertices[i]->links->size(); j++) for (int k = 0; k < vertices[vertices[i]->links->l[j]]->links->size(); k++)
    {
      number_of_triangles += std::min(vertices[i]->links->posweightof(vertices[vertices[i]->links->l[j]]->links->l[k]).first, 0) + 1;
    }
    possible = (vertices[i]->links->size())*(vertices[i]->links->size() - 1.) / 2.;
  }
  return number_of_triangles / possible;
}

inline void static_network::clustering_coefficient(std::map <int, std::deque<double>> & c_of_k)
{
  // c_of_k is a map kin -> clu. c. and it is not cleared

  for (int i = 0; i < dim; i++)
  {
    std::set<int> intneighs;		// neighbors of node c[i]
    for (int j = 0; j < vertices[i]->links->size(); j++)
      intneighs.insert(vertices[i]->links->l[j]);
    {
      const auto it = c_of_k.find(intneighs.size());
      if (it == c_of_k.end())
      {
        std::deque<double> f;
        c_of_k[intneighs.size()] = f;
      }
      if (intneighs.size() > 1)
        c_of_k[intneighs.size()].push_back(kin(intneighs) / (intneighs.size() * (intneighs.size() - 1)));
      else
        c_of_k[intneighs.size()].push_back(0);
    }
  }
}

inline bool static_network::GN_bench(int nodes, int modules, double kout, double m_degree, std::deque<std::deque<int>> & ten)
{
  ten.clear();

  int nm = nodes / modules;
  nodes = nm * modules;

  double pin = (m_degree - kout) / (nm - 1);
  double pout = kout / (nm*(modules - 1));

  std::cout << "pin " << pin << std::endl;
  if (modules > 1)
    std::cout << "pout " << pout << std::endl;

  int count = 1;
  std::deque<int> community;
  for (double i = 0; i < nodes; i++) {
    if (i >= nm * count) {
      ten.push_back(community);
      community.clear();
      count++;
    }

    community.push_back(cast_int(i));
  }

  ten.push_back(community);

  std::deque<std::deque<int>> e;
  std::deque<int> s1;

  for (double i = 0; i < nodes; i++)
    e.push_back(s1);

  count = 1;
  int edges = 0;

  for (int i = 0; i < nodes; i++)
  {
    if (i >= nm * count)
      count++;

    for (int j = i + 1; j < nodes; j++)
    {
      double p0;

      if (j < (nm*count))
        p0 = pin;
      else
        p0 = pout;

      if (ran4() <= p0)
      {
        e[i].push_back(j);
        e[j].push_back(i);
        edges++;
      }
    }
  }

  std::cout << "------------------------------------------------------------------" << std::endl;
  std::cout << "network di " << nodes << " nodi;\t" << edges << " lati;\t" << "grado medio: " << 2.*edges / nodes << std::endl;
  std::cout << "------------------------------------------------------------------" << std::endl;

  std::deque<int> labels;
  for (int i = 0; i < nodes; i++)
    labels.push_back(i);

  set_graph(e, labels);

  return true;
}

inline bool static_network::set_binary_random(const std::deque<int> & degrees)
{
  std::deque<std::set<int>> en;
  {
    std::set<int> first;
    for (int i = 0; i < degrees.size(); i++)
      en.push_back(first);
  }

  std::multimap <int, int> degree_node;

  for (int i = 0; i < degrees.size(); i++)
    degree_node.insert(degree_node.end(), std::make_pair(degrees[i], i));

  while (!degree_node.empty())
  {
    auto itlast = degree_node.end();
    --itlast;

    auto itit = itlast;
    std::deque <std::multimap<int, int>::iterator> erasenda;

    int inserted = 0;

    for (int i = 0; i < itlast->first; i++)
    {
      if (itit != degree_node.begin())
      {
        --itit;

        en[itlast->second].insert(itit->second);
        en[itit->second].insert(itlast->second);
        inserted++;

        erasenda.push_back(itit);
      }
      else
        break;
    }

    for (auto& i : erasenda)
    {
      if (i->first > 1)
        degree_node.emplace(i->first - 1, i->second);

      degree_node.erase(i);
    }
    degree_node.erase(itlast);
  }

  rewiring(en);

  std::deque<std::deque<int>> E;
  for (const auto& i : en)
  {
    std::deque<int> f;
    set_to_deque(i, f);
    E.push_back(f);
  }

  std::deque<int> labels;
  for (int i = 0; i < degrees.size(); i++)
    labels.push_back(i);

  set_graph(E, labels);

  return true;
}

// it computes the integral of a power law
inline double integral(double a, double b)
{
  if (fabs(a + 1.) > 1e-10)
    return (1. / (a + 1.)*pow(b, a + 1.));

  else
    return (log(b));
}

// it computes the correct (i.e. discrete) average of a power law
inline double integer_average(int n, int min, double tau)
{
  double a = 0;

  for (double h = min; h < n + 1; h++)
    a += pow((1. / h), tau);

  double pf = 0;
  for (double i = min; i < n + 1; i++)
    pf += 1 / a * pow((1. / (i)), tau)*i;

  return pf;
}

// it returns the average degree of a power law
double average_degree(const double &dmax, const double &dmin, const double &gamma)
{
  return (1. / (integral(gamma, dmax) - integral(gamma, dmin)))*(integral(gamma + 1, dmax) - integral(gamma + 1, dmin));
}

//bisection method to find the inferior limit, in order to have the expected average degree
double solve_dmin(const double& dmax, const double &dmed, const double &gamma)
{
  double dmin_l = 1;
  double dmin_r = dmax;
  double average_k1 = average_degree(dmin_r, dmin_l, gamma);
  double average_k2 = dmin_r;

  if ((average_k1 - dmed > 0)
    || (average_k2 - dmed < 0))
  {
    std::cerr << "\n***********************\nERROR: the average degree is out of range:";

    if (average_k1 - dmed > 0)
    {
      std::cerr << "\nyou should increase the average degree (bigger than " << average_k1 << ")" << std::endl;
      std::cerr << "(or decrease the maximum degree...)" << std::endl;
    }

    if (average_k2 - dmed < 0)
    {
      std::cerr << "\nyou should decrease the average degree (smaller than " << average_k2 << ")" << std::endl;
      std::cerr << "(or increase the maximum degree...)" << std::endl;
    }
    return -1;
  }

  while (fabs(average_k1 - dmed) > 1e-7)
  {
    double temp = average_degree(dmax, ((dmin_r + dmin_l) / 2.), gamma);
    if ((temp - dmed)*(average_k2 - dmed) > 0)
    {
      average_k2 = temp;
      dmin_r = ((dmin_r + dmin_l) / 2.);
    }
    else
    {
      average_k1 = temp;
      dmin_l = ((dmin_r + dmin_l) / 2.);
    }
  }
  return dmin_l;
}

inline bool static_network::set_random_powlaw(int num_nodes, double tau, double average_k, int max_degree)
{
  double dmin = solve_dmin(max_degree, average_k, -tau);
  if (dmin == -1)
    return false;

  int min_degree = int(dmin);

  double media1 = integer_average(max_degree, min_degree, tau);
  double media2 = integer_average(max_degree, min_degree + 1, tau);

  if (fabs(media1 - average_k) > fabs(media2 - average_k))
    min_degree++;

  std::deque <int> degree_seq;		//  degree sequence of the nodes
  std::deque <double> cumulative;
  powerlaw(max_degree, min_degree, tau, cumulative);

  for (int i = 0; i < num_nodes; i++)
  {
    int nn = lower_bound(cumulative.begin(), cumulative.end(), ran4()) - cumulative.begin() + min_degree;
    degree_seq.push_back(nn);
  }

  set_binary_random(degree_seq);
  return true;
}

inline bool static_network::set_multiple_random(const std::deque<int> & degrees, bool self_loops)
{
  // set self_loops = false to avoid self_loops

  int number_of_nodes = 0;
  int number_of_links = 0;
  for (int degree : degrees)
  {
    number_of_nodes++;
    number_of_links += degree;
  }

  std::deque<int> random_links(number_of_links);

  {
    std::deque<int> random_position(number_of_links);

    for (int i = 0; i < number_of_links; i++)
      random_position[i] = i;

    shuffle_s(random_position);

    int label = 0;
    int position = 0;
    for (int degree : degrees)
    {
      for (int k = 0; k < degree; k++)
        random_links[random_position[position++]] = label;
      label++;
    }
  }

  const int number_of_pairs = number_of_links / 2;

  if (number_of_links % 2 == 1)
    number_of_links--;

  //this is done to avoid self loops
  if (!self_loops)
  {
    for (int pos = 0; pos < number_of_links; pos += 2)
    {
      if (random_links[pos] == random_links[pos + 1])		// this is a self loop
      {
        bool flag = true;

        while (flag)
        {
          const int new_pair = irand(number_of_pairs - 1);

          if (random_links[2 * new_pair] != random_links[pos] && random_links[2 * new_pair + 1] != random_links[pos])
          {
            flag = false;
            random_links[pos + 1] = random_links[2 * new_pair + 1];
            random_links[2 * new_pair + 1] = random_links[pos];
          }
        }
      }
    }
  }

  std::map < std::pair<int, int>, int> links_and_weights;

  int position = 0;
  for (int i = 0; i < number_of_pairs; i++)
  {
    std::pair<int, int> p(random_links[position], random_links[position + 1]);

    if (random_links[position] > random_links[position + 1])
    {
      p.first = random_links[position + 1];
      p.second = random_links[position];
    }

    position += 2;

    auto itf = links_and_weights.find(p);

    if (itf == links_and_weights.end())
    {
      if (p.first == p.second)
        links_and_weights.emplace(p, 2);
      else
        links_and_weights.emplace(p, 1);
    }
    else
    {
      if (p.first == p.second)
      {
        itf->second += 2;
      }
      else
      {
        itf->second++;
      }
    }
  }

  std::deque<std::deque<int>> E;
  std::deque<std::deque<double>> Ew;
  std::deque<int> labels;

  for (unsigned i = 0; i < degrees.size(); i++)
  {
    E.emplace_back();
    Ew.emplace_back();
    labels.push_back(i);
  }

  for (auto &[link, weight] : links_and_weights)
  {
    E[link.first].push_back(link.second);
    E[link.second].push_back(link.first);

    Ew[link.first].push_back(weight);
    Ew[link.second].push_back(weight);
  }

  set_graph(true, E, Ew, labels);
  return true;
}

inline bool static_network::set_random_powlaw_multiple(int num_nodes, double tau, double average_k, int max_degree, bool self_loops)
{
  const double dmin = solve_dmin(max_degree, average_k, -tau);
  if (dmin == -1)
    return false;

  auto min_degree = int(dmin);

  const double media1 = integer_average(max_degree, min_degree, tau);
  const double media2 = integer_average(max_degree, min_degree + 1, tau);

  if (fabs(media1 - average_k) > fabs(media2 - average_k))
    min_degree++;

  std::deque <int> degree_seq;		//  degree sequence of the nodes
  std::deque <double> cumulative;
  powerlaw(max_degree, min_degree, tau, cumulative);

  for (int i = 0; i < num_nodes; i++)
  {
    int nn = lower_bound(cumulative.begin(), cumulative.end(), ran4()) - cumulative.begin() + min_degree;
    degree_seq.push_back(nn);
  }
  set_multiple_random(degree_seq, self_loops);

  return true;
}

inline void static_network::community_netknn(std::deque<std::deque<int>> & ten)
{
  // this function computes something...

  std::deque<int> members(dim);
  for (int i = 0; i < ten.size(); i++)
  {
    for (int j = 0; j < ten[i].size(); j++)
    {
      members[ten[i][j]] = i;
    }
  }

  prints(members);
  std::set<std::pair<int, int>> links;

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < vertices[i]->links->size(); j++)
    {
      const int ne = vertices[i]->links->l[j];
      if (members[i] != members[ne])
      {
        links.emplace(members[i], members[ne]);
      }
    }
  }

  std::map<int, std::deque<int>> size_size;
  for (const auto &[i, j] : links)
  {
    if (size_size.find(ten[i].size()) != size_size.end())
      size_size[ten[i].size()].push_back(ten[j].size());
    else
    {
      std::deque<int> f;
      size_size[ten[i].size()] = f;
      size_size[ten[i].size()].push_back(ten[j].size());
    }

    if (size_size.find(ten[i].size()) != size_size.end())
      size_size[ten[i].size()].push_back(ten[j].size());
    else
    {
      std::deque<int> f;
      size_size[ten[i].size()] = f;
      size_size[ten[i].size()].push_back(ten[j].size());
    }
  }

  std::ofstream sout("size_size");
  for (auto& it : size_size)
  {
    sout << it.first << " " << average_func(it.second) << std::endl;
  }
}

inline void static_network::community_net(std::deque<std::deque<int>> & ten)
{
  // this function computes something...
  std::deque<int> members(dim);
  for (unsigned i = 0; i < ten.size(); i++)
  {
    for (unsigned j = 0; j < ten[i].size(); j++)
    {
      members[ten[i][j]] = i;
    }
  }

  //prints(members);
  std::map<std::pair<int, int>, int> links;

  for (int i = 0; i < dim; i++)
  {
    for (int j = 0; j < vertices[i]->links->size(); j++)
    {
      const int ne = vertices[i]->links->l[j];
      if (members[i] != members[ne])
      {
        std::pair<int, int> P(members[i], members[ne]);

        if (links.find(P) == links.end())
          links[P] = 0;

        links[P]++;
      }
    }
  }

  std::ofstream sout("com_net");
  for (auto& link : links)
  {
    sout << link.first.first << " " << link.first.second << " " << link.second << std::endl;
  }
}

inline double static_network::monte_carlo_asp()
{
  if (dim < 1000)
  {
    std::deque<double> ML;
    double avt;
    int diam;
    diameter_and_asp(avt, diam, ML);

    return avt;
  }
  std::deque<double> asps;

  for (int i = 0; i < 1000; i++)
  {
    double average;

    std::deque<double> ML;
    int m = distances_from_i(irand(dim - 1), average, ML, irand(dim - 1));

    asps.push_back(average / dim);
  }
  return average_func(asps);
}

inline void static_network::add_isolated(int id_iso)
{
  vertices.push_back(new vertex(id_iso, 0, 0));
  dim++;
}

#endif
