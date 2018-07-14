#pragma once

#ifndef LOUVAIN_OSLOMNET_DIR_H
#define LOUVAIN_OSLOMNET_DIR_H

#include <deque>
#include <algorithm>
#include <map>

class oslom_module {
public:

  oslom_module(int a, int b) { nc = 1; kout_in = a; ktot_in = a; kout_out = b; ktot_out = b; };
  ~oslom_module() {};

  int nc;
  int kout_in;
  int kout_out;
  int ktot_in;
  int ktot_out;
};

class internal_links_weights
{
public:

  internal_links_weights(int a, double b, int c, double d)
    : k_in(a)
    , k_out(c)
    , win(b)
    , wout(d)
  {
  }
  ~internal_links_weights() = default;

  int k_in;
  int k_out;
  double win;
  double wout;
};

typedef std::map<int, internal_links_weights> mapip;
typedef std::map<int, oslom_module> map_int_om;

inline void int_histogram(const int & c, mapip  & hist, const int & w1, const double & w2, const int & w1_out, const double & w2_out)
{
  auto itf = hist.find(c);
  if (itf == hist.end())
  {
    internal_links_weights ILW(w1, w2, w1_out, w2_out);
    hist.emplace(c, ILW);
  }
  else
  {
    itf->second.k_in += w1;
    itf->second.win += w2;
    itf->second.k_out += w1_out;
    itf->second.wout += w2_out;
  }
}

inline void prints(map_int_om & M)
{
  for (auto& itm : M)
  {
    std::cout
      << "module: " << itm.first
      << "\t\t\t\tnc= " << itm.second.nc
      << "\t ktot_inout= " << itm.second.ktot_out + itm.second.ktot_in
      << "\t kout_in= " << itm.second.kout_in + itm.second.kout_out << std::endl;
  }
}

class oslomnet_louvain : public static_network
{
public:

  oslomnet_louvain() : static_network() {}
  ~oslomnet_louvain() = default;

  int collect_raw_groups_once(std::deque<std::deque<int>> &);

private:

  void weighted_favorite_of(const int & node, int & fi, int & kp_in, int & kop_in, int & kp_out, int & kop_out);
  void unweighted_favorite_of(const int & node, int & fi, int & kp_in, int & kop_in, int & kp_out, int & kop_out);

  void single_pass_unweighted();
  void single_pass_weighted();
  inline void update_modules(const int & i, const int & fi, const int & kp_in, const int & kop_in, const int & kp_out, const int & kop_out);

  void module_initializing();
  void set_partition_collected(std::deque<std::deque<int>> & ten2);

  //int check_all();

  std::map<int, oslom_module> label_module;
  std::deque<int> vertex_label;
  std::deque<int> vertex_order;
  std::deque<bool> vertex_to_check;
  std::deque<bool> vertex_to_check_next;
  int nodes_changed;
};

inline void oslomnet_louvain::module_initializing()
{
  for (int i = 0; i < dim; i++)
  {
    vertex_label.push_back(i);
    vertex_order.push_back(i);
    vertex_to_check.push_back(true);
    vertex_to_check_next.push_back(false);
    oslom_module newm(vertices[i]->instub_number, vertices[i]->outstub_number);
    label_module.emplace(i, newm);
  }
}

inline void oslomnet_louvain::unweighted_favorite_of(const int & node, int & fi, int & kp_in, int & kop_in, int & kp_out, int & kop_out)
{
  double min_fitness = 10;
  fi = vertex_label[node];
  kp_in = 0;
  kop_in = 0;
  kp_out = 0;
  kop_out = 0;

  std::map<int, std::pair<int, int>>  M;		// M is a map module_label -> kin_in, kin_out (internal stubs)

  for (int j = 0; j < vertices[node]->inlinks->size(); j++)
  {
    int_histogram(vertex_label[vertices[node]->inlinks->l[j]], M, vertices[node]->inlinks->w[j].first, 0);
  }

  for (int j = 0; j < vertices[node]->outlinks->size(); j++)
  {
    int_histogram(vertex_label[vertices[node]->outlinks->l[j]], M, 0, vertices[node]->outlinks->w[j].first);
  }

  for (auto &[label, internal_stubs] : M)
  {
    const auto itOM = label_module.find(label);

    double to_fit;

    if (label != vertex_label[node])
    {
      to_fit =
        compute_global_fitness_randomized(
          internal_stubs.first,
          itOM->second.kout_in,
          internal_stubs.second,
          itOM->second.kout_out,
          oneM - itOM->second.ktot_in,
          oneM - itOM->second.ktot_out,
          vertices[node]->instub_number,
          vertices[node]->outstub_number,
          0.,
          0.,
          1,
          1);
    }
    else
    {
      kop_in = internal_stubs.first;
      kop_out = internal_stubs.second;

      const int kout_prime_out = itOM->second.kout_out - vertices[node]->outstub_number + kop_in + kop_out;
      const int kout_prime_in = itOM->second.kout_in - vertices[node]->instub_number + kop_in + kop_out;

      to_fit = compute_global_fitness_randomized(
        internal_stubs.first,
        kout_prime_in,
        internal_stubs.second,
        kout_prime_out,
        oneM - itOM->second.ktot_in + vertices[node]->instub_number,
        oneM - itOM->second.ktot_out + vertices[node]->outstub_number,
        vertices[node]->instub_number, vertices[node]->outstub_number,
        0.,
        0.,
        1,
        1);

      to_fit *= 0.999;		// to break possible ties
    }

    if (to_fit < min_fitness)
    {
      kp_in = internal_stubs.first;
      kp_out = internal_stubs.second;
      min_fitness = to_fit;
      fi = label;
    }
  }
}

inline void oslomnet_louvain::weighted_favorite_of(const int & node, int & fi, int & kp_in, int & kop_in, int & kp_out, int & kop_out)
{
  double min_fitness = 10;
  fi = vertex_label[node];
  kp_in = 0;
  kop_in = 0;
  kp_out = 0;
  kop_out = 0;

  mapip M;		// M is a map module_label -> kin - win (internal stubs, internal weight)

  for (int j = 0; j < vertices[node]->inlinks->size(); j++)
  {
    int_histogram(vertex_label[vertices[node]->inlinks->l[j]], M, vertices[node]->inlinks->w[j].first, vertices[node]->inlinks->w[j].second, 0, 0.);
  }

  for (int j = 0; j < vertices[node]->outlinks->size(); j++)
  {
    int_histogram(vertex_label[vertices[node]->outlinks->l[j]], M, 0, 0., vertices[node]->outlinks->w[j].first, vertices[node]->outlinks->w[j].second);
  }

  for (auto &[label, internal_links] : M)
  {
    const auto itOM = label_module.find(label);

    double to_fit;

    if (label != vertex_label[node])
    {
      to_fit = compute_global_fitness_randomized(
        internal_links.k_in, itOM->second.kout_in,
        internal_links.k_out, itOM->second.kout_out,
        oneM - itOM->second.ktot_in, oneM - itOM->second.ktot_out, vertices[node]->instub_number, vertices[node]->outstub_number,
        internal_links.win,
        internal_links.wout,
        std::min(dim - itOM->second.nc, (itOM->second.kout_in + itOM->second.kout_out) / (internal_links.k_out + internal_links.k_in) + 1),
        dim - itOM->second.nc);
    }
    else
    {
      kop_in = internal_links.k_in;
      kop_out = internal_links.k_out;

      const int kout_prime_out = itOM->second.kout_out - vertices[node]->outstub_number + kop_in + kop_out;
      const int kout_prime_in = itOM->second.kout_in - vertices[node]->instub_number + kop_in + kop_out;

      to_fit = compute_global_fitness_randomized(
        internal_links.k_in, kout_prime_in,
        internal_links.k_out, kout_prime_out,
        oneM - itOM->second.ktot_in + vertices[node]->instub_number, oneM - itOM->second.ktot_out + vertices[node]->outstub_number,
        vertices[node]->instub_number, vertices[node]->outstub_number,
        internal_links.win,
        internal_links.wout,
        std::min(dim - itOM->second.nc, (kout_prime_out + kout_prime_in) / (kop_in + kop_out) + 1), dim - itOM->second.nc);

      to_fit *= 0.999;		// to break possible ties
    }

    if (to_fit < min_fitness)
    {
      kp_in = internal_links.k_in;
      kp_out = internal_links.k_out;
      min_fitness = to_fit;
      fi = label;
    }
  }
}

inline void oslomnet_louvain::update_modules(const int & i, const int & fi, const int & kp_in, const int & kop_in, const int & kp_out, const int & kop_out)
{
  if (fi != vertex_label[i])
  {
    nodes_changed++;

    for (int j = 0; j < vertices[i]->inlinks->size(); j++)
    {
      vertex_to_check_next[vertices[i]->inlinks->l[j]] = true;
    }
    for (int j = 0; j < vertices[i]->outlinks->size(); j++)
    {
      vertex_to_check_next[vertices[i]->outlinks->l[j]] = true;
    }

    auto itm = label_module.find(vertex_label[i]);
    --(itm->second.nc);
    if (itm->second.nc == 0)
      label_module.erase(itm);
    else {
      itm->second.kout_in -= vertices[i]->instub_number - kop_in - kop_out;
      itm->second.ktot_in -= vertices[i]->instub_number;
      itm->second.kout_out -= vertices[i]->outstub_number - kop_in - kop_out;
      itm->second.ktot_out -= vertices[i]->outstub_number;
    }

    itm = label_module.find(fi);
    ++(itm->second.nc);

    itm->second.kout_in += vertices[i]->instub_number - kp_in - kp_out;
    itm->second.ktot_in += vertices[i]->instub_number;
    itm->second.kout_out += vertices[i]->outstub_number - kp_in - kp_out;
    itm->second.ktot_out += vertices[i]->outstub_number;

    vertex_label[i] = fi;
  }

  //check_all();
}

inline void oslomnet_louvain::single_pass_unweighted()
{
  int fi, kp_in, kop_in, kp_out, kop_out;	// fi is the label node i likes most, kp is the number od internal stubs in module fi, kop is the same for vertex_label[i]

  for (auto itd = vertex_order.begin(); itd != vertex_order.end(); ++itd) {
    if (vertex_to_check[*itd]) {
      unweighted_favorite_of(*itd, fi, kp_in, kop_in, kp_out, kop_out);
      update_modules(*itd, fi, kp_in, kop_in, kp_out, kop_out);
    }
  }
}

inline void oslomnet_louvain::single_pass_weighted()
{
  int fi, kp_in, kop_in, kp_out, kop_out;	// fi is the label node i likes most, kp is the number od internal stubs in module fi, kop is the same for vertex_label[i]

  for (auto itd : vertex_order)
  {
    if (vertex_to_check[itd])
    {
      weighted_favorite_of(itd, fi, kp_in, kop_in, kp_out, kop_out);
      update_modules(itd, fi, kp_in, kop_in, kp_out, kop_out);
    }
  }
}

inline void oslomnet_louvain::set_partition_collected(std::deque<std::deque<int>> & ten2)
{
  ten2.clear();

  std::deque<std::deque<int>> M;

  // take partition from vertex_label  //******************************
  std::map<int, int> mems;
  for (int i = 0; i < dim; i++)
  {
    const std::pair<std::map<int, int>::iterator, bool>  itm_bool = mems.emplace(vertex_label[i], mems.size());

    if (itm_bool.second)
    {
      M.push_back(std::deque<int>());
    }

    M[itm_bool.first->second].push_back(i);
  }

  // check if subgraphs are connected  //******************************
  for (auto& i : M)
  {
    std::deque<std::deque<int>> link_per_node;
    std::deque<std::deque<std::pair<int, double>> > weights_per_node;
    set_subgraph(i, link_per_node, weights_per_node);
    static_network giovanni;
    giovanni.set_graph(link_per_node, weights_per_node, i);
    std::deque<std::deque<int>> gM;
    giovanni.set_connected_components(gM);

    if (gM.size() == 1)
    {
      ten2.push_back(i);
    }
    else
    {
      for (auto& j : gM)
      {
        giovanni.deque_id(j);
        ten2.push_back(j);
      }
    }
  }
}

inline int oslomnet_louvain::collect_raw_groups_once(std::deque<std::deque<int>> & P)
{
  module_initializing();

  int stopper = 0;
  int previous_nodes_changed = dim;
  int iteration = 0;

  while (true)
  {
    nodes_changed = 0;

    for (int i = 0; i<int(vertex_to_check.size()); i++)
      vertex_to_check_next[i] = false;
    shuffle_s(vertex_order);

    if (paras.weighted)
    {
      single_pass_weighted();
    }
    else
    {
      single_pass_unweighted();
    }

    if (paras.print_flag_subgraph && iteration % 20 == 0)
    {
      std::cout << "iteration: " << iteration << " number of modules: " << label_module.size() << std::endl;
    }

    ++iteration;

    /* this conditions means that at least max_iteration_convergence iterations are done
     * and, after that, the number of nodes changed has to decrease (up to a facto 1e-3) */
    if (double(nodes_changed - previous_nodes_changed) > 1e-3 * previous_nodes_changed
      && iteration > paras.max_iteration_convergence)
    {
      stopper++;
    }

    if (stopper == paras.max_iteration_convergence
      || nodes_changed == 0
      || iteration == dim)
    {
      break;
    }

    vertex_to_check = vertex_to_check_next;
    previous_nodes_changed = nodes_changed;
  }

  set_partition_collected(P);

  if (paras.print_flag_subgraph)
  {
    std::cout << "collection done " << std::endl << std::endl << std::endl;
  }

  label_module.clear();
  vertex_label.clear();
  vertex_order.clear();
  vertex_to_check.clear();
  vertex_to_check_next.clear();
  nodes_changed = 0;

  return 0;
}

/*
int oslomnet_louvain::check_all() {
  bool print_stuff_oslom_local= false;

  std::deque<std::deque<int>> M;
  set_partition_collected(M);

  if(print_stuff_oslom_local)
    prints(label_module);

  double one=0;
  double two=0;

  for(map_int_om :: iterator itm = label_module.begin(); itm!=label_module.end(); itm++)
    one+=sin(itm->second.nc) + log(itm->second.ktot_in + itm->second.ktot_out ) + cos(itm->second.kout_in + itm->second.kout_out);

  for(int i=0; i<int(M.size()); i++) {
    int _ktot= ktot_m(M[i]).first + ktot_m(M[i]).second ;
    int _kin= kin_m(M[i]);
    if(print_stuff_oslom_local)
      cout<<"nc: "<<M[i].size()<<" kout_inout: "<<_ktot - 2* _kin<<" ktot_inout: "<<_ktot<< std::endl;
    two+=sin(M[i].size()) + log( _ktot) + cos(_ktot - 2*_kin);
  }

  cherr(one-two, 1e-6);
  if(print_stuff_oslom_local)
    cout<<"one, two: "<<one<<" "<<two<< std::endl;

  cout<<"check passed"<< std::endl;
  return 0;
}
//*/
#endif // LOUVAIN_OSLOMNET_DIR_H
