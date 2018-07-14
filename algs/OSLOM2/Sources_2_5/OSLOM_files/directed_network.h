#pragma once

#ifndef STATIC_static_network_INCLUDED
#define STATIC_static_network_INCLUDED

#include "./standard_package/standard_include.cpp"
#include "wsarray.h"

class static_network
{
public:
  static_network() = default;
  ~static_network();

  //int draw(std::string, bool);
  int draw(const std::string&);
  int draw_consecutive(const std::string& file_name1, const std::string& file_name2);
  int draw_with_weight_probability(const std::string& file_name);

  void print_id(const std::deque<int> & a, std::ostream &);
  void print_id(const std::deque<std::deque<int>> &, std::ostream &);
  void print_id(const std::deque<std::set<int>> &, std::ostream &);
  void print_id(const std::set<int> &, std::ostream &);
  void deque_id(std::deque<int> &);

  void set_subgraph(std::deque<int> & group, std::deque<std::deque<int>> & link_per_node, std::deque<std::deque<std::pair<int, double>> > & weights_per_node);

  int translate(std::deque<int> & ten);
  int translate(std::deque<std::deque<int>> &);
  int translate_anyway(std::deque<std::deque<int>> & ten);

  void get_id_label(std::map <int, int> &);
  int id_of(int a) { return vertices[a]->id_num; };

  int size() const
  {
    return dim;
  };

  double stubs() const
  {
    return oneM;
  }

  int kin_m(const std::deque<int> &);
  int kin_m(const std::set<int> &);
  std::pair<int, int> ktot_m(const std::deque<int> &);
  std::pair<int, int> ktot_m(const std::set<int> &);

  void set_graph(std::map<int, std::map<int, std::pair<int, double>> > & A);
  bool set_graph(const std::string& file_name);
  void set_graph(std::deque<std::deque<int>> & link_per_node, std::deque<std::deque<std::pair<int, double>> > & weights_per_node, std::deque<int> & label_rows);
  void clear();

  void set_proper_weights();

  void set_connected_components(std::deque<std::deque<int>> &);
  int propagate_distances(std::deque<int> & new_shell, std::set<int> & already_gone, std::deque<std::pair<int, int >> & distances_node, int shell, std::deque<double> & ML, int &, int);
  void same_component(int, std::set<int> &);

  int set_upper_network(std::map<int, std::map<int, std::pair<int, double>> > & neigh_weight_f, module_collection & module_coll);

  void print_degree_of_homeless(std::deque<int> & homel, std::ostream & outt);

protected:
  class vertex
  {
  public:
    vertex(int b, int c, int preall_i, int preall_o);
    ~vertex();

    std::pair<int, int> kplus_m(const std::deque<int> &);
    std::pair<double, double> kplus_w(const std::deque<int> &) const;

    std::pair<int, int> kplus_m(const std::set<int> &) const;

    int id_num;							// id

    double instrength;					// sum of the weights
    int instub_number;					// number of stubs
    double outstrength;					// sum of the weights
    int outstub_number;					// number of stubs

    wsarray* inlinks;					// array with label of neighbor, multiple links, sum of the weights towards it
    wsarray* outlinks;					// array with label of neighbor, multiple links, sum of the weights towards it
    //deque<double> in_original_weights;
    std::deque<double> out_original_weights;
  };

  int dim;									// number of nodes
  int oneM;									// number of in(out)-stubs

  std::deque <vertex*> vertices;

  void set_oneM_etc();
};

inline static_network::vertex::vertex(int b, int c, int preall_o, int preall_i)
  : id_num(b)
{
  outlinks = new wsarray(preall_o);
  inlinks = new wsarray(preall_i);
}

inline static_network::vertex::~vertex()
{
  delete inlinks;
  inlinks = nullptr;

  delete outlinks;
  outlinks = nullptr;
}

inline std::pair<int, int> static_network::vertex::kplus_m(const std::deque<int> &a)
{
  // computes the internal degree of the vertex respect with a
  int ins = 0;
  //double f=0;
  for (int i : a)
  {
    std::pair<int, double> A = inlinks->posweightof(i).second;
    ins += A.first;
  }

  int outs = 0;
  //double f=0;
  for (int i : a)
  {
    std::pair<int, double> A = outlinks->posweightof(i).second;
    outs += A.first;
  }

  return std::make_pair(ins, outs);
}

inline std::pair<double, double> static_network::vertex::kplus_w(const std::deque<int> &a) const
{
  // computes the internal degree of the vertex respect with a

  double ins = 0;
  //double f=0;
  for (int i : a)
  {
    std::pair<int, double> A = inlinks->posweightof(i).second;
    ins += A.second;
  }

  double outs = 0;
  //double f=0;
  for (int i : a)
  {
    std::pair<int, double> A = outlinks->posweightof(i).second;
    //cout<<a[i]<<" -*-*-* "<<A.first<<" "<<A.second<< std::endl;

    outs += A.second;
  }
  return std::make_pair(ins, outs);
}

inline std::pair<int, int> static_network::vertex::kplus_m(const std::set<int> &a) const
{
  // computes the internal degree of the vertex respect with a (a is supposed to be sorted)

  int ins = 0;
  //double f=0;

  for (int i = 0; i < inlinks->size(); i++)
  {
    if (a.find(inlinks->l[i]) != a.end())
    {
      ins += inlinks->w[i].first;
    }
  }

  int outs = 0;
  for (int i = 0; i < outlinks->size(); i++)
  {
    if (a.find(outlinks->l[i]) != a.end())
    {
      outs += outlinks->w[i].first;
    }
  }
  return std::make_pair(ins, outs);
}

inline static_network::~static_network()
{
  clear();
}

inline void static_network::clear()
{
  for (auto& vertice : vertices)
  {
    delete vertice;
    vertice = nullptr;
  }

  vertices.clear();
  dim = 0;
  oneM = 0;
}

inline void static_network::set_graph(std::map<int, std::map<int, std::pair<int, double>>> & A)
{
  // this maps the id into the usual stuff neighbors-weights
  // A gives the out-links				 ****IMPORTANT****

  std::deque<std::deque<int>> link_per_node;
  std::deque<std::deque<std::pair<int, double>>> weights_per_node;
  std::deque<int> label_rows;

  for (auto& itm : A)
  {
    label_rows.push_back(itm.first);
    std::deque<int> n;
    std::deque<std::pair<int, double>> w;

    for (auto& itm2 : itm.second)
    {
      if (itm2.second.first > 0)
      {
        n.push_back(itm2.first);
        w.push_back(itm2.second);
      }
    }
    link_per_node.push_back(n);
    weights_per_node.push_back(w);
  }
  /*
  prints(label_rows);
  printm(link_per_node);
  printm(weights_per_node);//*/
  set_graph(link_per_node, weights_per_node, label_rows);
}

//   all this stuff here should be improved
//   all this stuff here should be improved
//   all this stuff here should be improved

inline void static_network::set_oneM_etc()
{
  oneM = 0;

  for (int i = 0; i < dim; i++)
  {
    vertices[i]->inlinks->freeze();
    vertices[i]->outlinks->freeze();

    int stub_number_i = 0;
    double strength_i = 0;

    for (int j = 0; j < vertices[i]->inlinks->size(); j++)
    {
      stub_number_i += vertices[i]->inlinks->w[j].first;
      strength_i += vertices[i]->inlinks->w[j].second;
      //cout<<"-> "<<vertices[i]->links->w[j].second<< std::endl;
    }

    vertices[i]->instub_number = stub_number_i;
    vertices[i]->instrength = strength_i;

    stub_number_i = 0;
    strength_i = 0;

    for (int j = 0; j < vertices[i]->outlinks->size(); j++)
    {
      stub_number_i += vertices[i]->outlinks->w[j].first;
      strength_i += vertices[i]->outlinks->w[j].second;
    }

    vertices[i]->outstub_number = stub_number_i;
    vertices[i]->outstrength = strength_i;

    oneM += stub_number_i;
  }
}

inline bool static_network::set_graph(const std::string& file_name)
{
  clear();

  std::map<int, int> newlabels;
  std::deque<std::deque<int>>  link_per_node;
  std::deque<std::deque<std::pair<int, double>> >  weights_per_node;
  std::deque<int> label_rows;

  int label = 0;

  std::ifstream inb(file_name);
  std::string ins;

  while (getline(inb, ins)) if (!ins.empty() && ins[0] != '#') {
    std::deque<double> ds;
    cast_string_to_doubles(ins, ds);

    if (ds.size() < 2)
    {
      std::cerr << "From file " << file_name << ": string not readable " << ins << " " << std::endl;
      return false;
    }

    //prints(ds);
    //cout<<"-------------------------------------"<< std::endl;

    int innum1 = cast_int(ds[0]);
    int innum2 = cast_int(ds[1]);

    if (innum1 != innum2)
    {
      if (newlabels.find(innum1) == newlabels.end())
      {
        newlabels.emplace(innum1, label++);
        label_rows.push_back(innum1);
        std::deque<int> first;
        link_per_node.push_back(first);
        std::deque<std::pair<int, double>> sec;
        weights_per_node.push_back(sec);
      }

      if (newlabels.find(innum2) == newlabels.end())
      {
        newlabels.emplace(innum2, label++);
        label_rows.push_back(innum2);

        std::deque<int> first;
        link_per_node.push_back(first);
        std::deque<std::pair<int, double>> sec;
        weights_per_node.push_back(sec);
      }

      int pos = newlabels[innum1];
      link_per_node[pos].push_back(innum2);

      double w = 1;
      int multiple_l = 1;

      //--------------------------------------------------------------
      //--------------------------------------------------------------
      if (ds.size() >= 4)
      {
        if (paras.weighted)
        {
          //---------------------------------------------
          if (ds[2] > 0)
          {
            w = ds[2];
          }
          else
          {
            std::cerr << "error: not positive weights" << std::endl;
            return false;
          }
          //---------------------------------------------

          if (ds[3] > 0.99) {
            //cout<<ds[2]<<"<- "<< std::endl;
            multiple_l = cast_int(ds[3]);
          }
          //---------------------------------------------
        }
        else
        {
          if (ds[2] > 0.99) {
            //cout<<ds[2]<<"<- "<< std::endl;
            multiple_l = cast_int(ds[2]);
          }
        }
      }

      if (ds.size() == 3)
      {
        if (paras.weighted)
        {
          if (ds[2] > 0)
            w = ds[2];
          else {
            std::cerr << "error: not positive weights" << std::endl;
            return false;
          }
        }
        else {
          if (ds[2] > 0.99)
            multiple_l = cast_int(ds[2]);
        }
      }
      //--------------------------------------------------------------
      //--------------------------------------------------------------

      weights_per_node[pos].emplace_back(multiple_l, w);
    }
  }

  /*
  prints(label_rows);
  printm(link_per_node);
  //printm(weights_per_node);
  //*/

  set_graph(link_per_node, weights_per_node, label_rows);
  return true;
}

inline void static_network::set_graph(
  std::deque<std::deque<int>> & link_per_node,
  std::deque<std::deque<std::pair<int, double>>> & weights_per_node,
  std::deque<int> & label_rows)
{
  clear();

  // there is no check between label_rows and link per node but they need to have the same labels
  // link_per_node and weights_per_node are the list of out-links and weights. label_rows[i] is the label corresponding to row i

  std::map<int, int> newlabels;		// this maps the old labels with the new one
  for (int& label_row : label_rows)
  {
    newlabels.emplace(label_row, newlabels.size());
  }

  dim = newlabels.size();

  std::deque<int> inlinks_pernode;		// this counts how many in-links a given nodes receives
  for (int i = 0; i < dim; i++)
    inlinks_pernode.push_back(0);

  for (auto& i : link_per_node)
  {
    for (unsigned j = 0; j < i.size(); j++)
    {
      inlinks_pernode[newlabels[i[j]]]++;
    }
  }

  for (int i = 0; i < dim; i++)
    vertices.push_back(new vertex(0, 0, link_per_node[i].size(), inlinks_pernode[i]));

  for (auto& newlabel : newlabels)
  {
    vertices[newlabel.second]->id_num = newlabel.first;
  }

  for (unsigned i = 0; i < link_per_node.size(); i++)
  {
    for (unsigned j = 0; j < link_per_node[i].size(); j++)
    {
      int new2 = newlabels[link_per_node[i][j]];
      vertices[i]->outlinks->push_back(new2, weights_per_node[i][j].first, weights_per_node[i][j].second);
      vertices[new2]->inlinks->push_back(i, weights_per_node[i][j].first, weights_per_node[i][j].second);
    }
  }

  set_oneM_etc();

  //draw_consecutive();
  //draw("GIIO");
  if (paras.weighted)
    set_proper_weights();
}

inline int static_network::kin_m(const std::deque<int> & seq)
{
  if (seq.size() > double(oneM) / dim)
  {
    std::set<int> H;
    deque_to_set(seq, H);
    return kin_m(H);
  }

  int ki = 0;
  //int ko=0;
  for (unsigned i = 0; i < seq.size(); i++)
  {
    ki += vertices[seq[i]]->kplus_m(seq).first;
    //ko+=vertices[seq[i]]->kplus_m(seq).second;
  }

  return ki;
}

inline std::pair<int, int> static_network::ktot_m(const std::deque<int> &seq)
{
  int ki = 0;
  int ko = 0;

  for (int i : seq)
  {
    ki += vertices[i]->instub_number;
    ko += vertices[i]->outstub_number;
  }
  return std::make_pair(ki, ko);
}

inline std::pair<int, int> static_network::ktot_m(const std::set<int> &s)
{
  int ki = 0;
  int ko = 0;

  for (auto it : s)
  {
    ki += vertices[it]->instub_number;
    ko += vertices[it]->outstub_number;
  }

  return std::make_pair(ki, ko);
}

inline int static_network::kin_m(const std::set <int> &s)
{
  int ki = 0;

  for (auto it : s)
  {
    ki += vertices[it]->kplus_m(s).first;
  }

  return ki;
}

inline int static_network::draw_consecutive(const std::string& file_name1, const std::string&
  file_name2)
{
  //cout<<"drawing in file "<<b<< std::endl;
  std::ofstream graph_out(file_name1);

  if (paras.weighted)
  {
    for (unsigned i = 0; i < vertices.size(); i++)
    {
      for (int j = 0; j < vertices[i]->outlinks->size(); j++)
      {
        graph_out << i
          << "\t" << vertices[i]->outlinks->l[j]
          << "\t" << cast_int(vertices[i]->out_original_weights[j]) << std::endl;
      }
    }
  }
  else
  {
    for (unsigned i = 0; i < vertices.size(); i++)
    {
      for (int j = 0; j < vertices[i]->outlinks->size(); j++)
      {
        graph_out << i
          << "\t" << vertices[i]->outlinks->l[j]
          << "\t" << vertices[i]->outlinks->w[j].first << std::endl;
      }
    }
  }

  std::ofstream graph_out2(file_name2);
  for (unsigned i = 0; i < vertices.size(); i++)
    graph_out2 << i << " " << vertices[i]->id_num << std::endl;

  return 0;
}

inline int static_network::draw(const std::string& file_name)
{
  std::ofstream graph_out(file_name);

  if (paras.weighted)
  {
    for (auto& vertex : vertices)
    {
      for (int j = 0; j < vertex->outlinks->size(); j++)
      {
        graph_out << vertex->id_num
          << "\t" << vertices[vertex->outlinks->l[j]]->id_num
          << "\t" << vertex->out_original_weights[j]
          << "\t" << vertex->outlinks->w[j].first << std::endl;
      }
    }
  }
  else
  {
    for (auto& vertex : vertices)
    {
      for (int j = 0; j < vertex->outlinks->size(); j++)
      {
        graph_out << vertex->id_num
          << "\t" << vertices[vertex->outlinks->l[j]]->id_num
          << "\t" << vertex->outlinks->w[j].first << std::endl;
      }
    }
    //"\t"<<vertices[i]->out_original_weights[j]<<"\t"<<vertices[i]->outlinks->w[j].second<< std::endl;
  }

  /*
  for (int i=0; i<vertices.size(); i++)
    for (int j=0; j<vertices[i]->inlinks->size(); j++)
      cout<<vertices[i]->id_num<<"\t"<<vertices[vertices[i]->inlinks->l[j]]->id_num<<"\t"<<vertices[i]->inlinks->w[j].first<<
      "\t"<<77<<"\t"<<vertices[i]->inlinks->w[j].second<< std::endl;

  */

  return 0;
}

inline int static_network::draw_with_weight_probability(const std::string& file_name)
{
  std::ofstream graph_out(file_name);

  if (paras.weighted)
  {
    for (auto& vertex : vertices)
    {
      for (int j = 0; j < vertex->outlinks->size(); j++)
      {
        graph_out << vertex->id_num
          << "\t" << vertices[vertex->outlinks->l[j]]->id_num
          << "\t" << vertex->out_original_weights[j]
          << "\t" << vertex->outlinks->w[j].first
          << " " << vertex->outlinks->w[j].second << std::endl;
      }
    }
  }
  else
  {
    for (auto& vertex : vertices)
    {
      for (int j = 0; j < vertex->outlinks->size(); j++)
      {
        graph_out << vertex->id_num
          << "\t" << vertices[vertex->outlinks->l[j]]->id_num
          << "\t" << vertex->outlinks->w[j].first << std::endl;
      }
    }
    //"\t"<<vertices[i]->out_original_weights[j]<<"\t"<<vertices[i]->outlinks->w[j].second<< std::endl;
  }

  /*
  for (int i=0; i<vertices.size(); i++)
    for (int j=0; j<vertices[i]->inlinks->size(); j++)
      cout<<vertices[i]->id_num<<"\t"<<vertices[vertices[i]->inlinks->l[j]]->id_num<<"\t"<<vertices[i]->inlinks->w[j].first<<
      "\t"<<77<<"\t"<<vertices[i]->inlinks->w[j].second<< std::endl;

  */

  return 0;
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
    i = vertices[i]->id_num;
}

inline void static_network::print_id(const std::deque<int> & a, std::ostream & pout)
{
  for (int i : a)
    pout << vertices[i]->id_num << "\t";
  pout << std::endl;
}

inline void static_network::print_id(const std::set<int> & a, std::ostream & pout)
{
  for (auto its : a)
  {
    pout << vertices[its]->id_num << "\t";
  }
  pout << std::endl;
}

inline void static_network::print_id(const std::deque<std::deque<int>> & a, std::ostream & pout)
{
  for (const auto& i : a)
    print_id(i, pout);
}

inline void static_network::print_id(const std::deque<std::set<int>> & a, std::ostream & pout)
{
  for (const auto& i : a)
    print_id(i, pout);
}

inline int static_network::translate(std::deque<std::deque<int>> & ten) {
  std::map<int, int> A;
  get_id_label(A);

  std::deque<std::deque<int>> ten2;

  for (unsigned i = 0; i < ten.size(); i++)
  {
    std::deque<int> ff;

    for (unsigned j = 0; j < ten[i].size(); j++)
    {
      auto itf = A.find(ten[i][j]);
      if (itf == A.end())
      {
        std::cerr << "warning: the nodes in the communities are different from those ones in the network!" << std::endl;
        //return -1;
      }
      else
        ff.push_back(itf->second);
    }

    if (!ff.empty())
      ten2.push_back(ff);
  }

  ten = ten2;

  return 0;
}

inline int static_network::translate(std::deque<int> & ten) {
  std::map<int, int> A;
  get_id_label(A);

  for (unsigned i = 0; i < ten.size(); i++)
  {
    auto itf = A.find(ten[i]);
    if (itf == A.end())
    {
      std::cerr << "warning: the nodes in the communities are different from those ones in the network!" << std::endl;
      //return -1;
    }
    else
      ten[i] = itf->second;
  }

  return 0;
}

inline void static_network::set_subgraph(
  std::deque<int> & group,
  std::deque<std::deque<int>> & link_per_node,
  std::deque<std::deque<std::pair<int, double>>> & weights_per_node)
{
  // in this function I'm not using id's... because I want to work with the same labels (don't want to translate)

  std::sort(group.begin(), group.end());

  weights_per_node.clear();
  link_per_node.clear();

  for (unsigned i = 0; i < group.size(); i++)
  {
    int nodei = group[i];

    std::deque<int> link_i;
    std::deque<std::pair<int, double>> weight_i;

    for (int j = 0; j < vertices[nodei]->outlinks->size(); j++)
    {
      if (binary_search(group.begin(), group.end(), vertices[nodei]->outlinks->l[j]))
      {
        link_i.push_back(vertices[nodei]->outlinks->l[j]);
        if (paras.weighted)
          weight_i.emplace_back(vertices[nodei]->outlinks->w[j].first, vertices[nodei]->out_original_weights[j]);
        else
          weight_i.emplace_back(vertices[nodei]->outlinks->w[j].first, 1);

        //weight_i.push_back(std::make_pair(vertices[nodei]->links->w[j].first, vertices[nodei]->links->w[j].second));
      }
    }

    link_per_node.push_back(link_i);
    weights_per_node.push_back(weight_i);
  }
}

inline void static_network::set_proper_weights() {
  // this function id to normalize the weights in order to have the -log(prob(Weiight>w)) which is simply w[i].second / <w_ij>
  // <w_ij> is s_i * s_j / k_i / k_j / eta
  // eta is <s_i/k_i>

  if (dim == 0)
  {
    //cout<<"network empty"<< std::endl;
    //cherr();
  }
  else
  {
    for (int i = 0; i < dim; i++) {
      vertices[i]->out_original_weights.clear();

      for (int j = 0; j < vertices[i]->outlinks->size(); j++)
      {
        vertices[i]->out_original_weights.push_back(vertices[i]->outlinks->w[j].second);
      }
    }

    for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < vertices[i]->outlinks->size(); j++)
      {
        int & neigh = vertices[i]->outlinks->l[j];

        double w1 = (vertices[i]->outstrength / vertices[i]->outstub_number) * vertices[i]->outlinks->w[j].first;
        double w2 = (vertices[neigh]->instrength / vertices[neigh]->instub_number) * vertices[i]->outlinks->w[j].first;

        vertices[i]->outlinks->w[j].second /= 2. / (1. / w1 + 1. / w2);
        int posneigh = vertices[neigh]->inlinks->find(i);
        vertices[neigh]->inlinks->w[posneigh].second /= 2. / (1. / w1 + 1. / w2);
      }
    }
  }
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

inline void static_network::same_component(int source, std::set<int> & already_gone)
{
  already_gone.clear();
  already_gone.insert(source);

  std::deque<int> this_shell;
  this_shell.push_back(source);

  while (!this_shell.empty()) {
    std::deque<int> next_shell;

    for (int i = 0; i<int(this_shell.size()); i++)
    {
      for (int j = 0; j < vertices[this_shell[i]]->inlinks->size(); j++)
      {
        if (already_gone.insert(vertices[this_shell[i]]->inlinks->l[j]).second)
          next_shell.push_back(vertices[this_shell[i]]->inlinks->l[j]);
      }

      for (int j = 0; j < vertices[this_shell[i]]->outlinks->size(); j++)
      {
        if (already_gone.insert(vertices[this_shell[i]]->outlinks->l[j]).second)
          next_shell.push_back(vertices[this_shell[i]]->outlinks->l[j]);
      }
    }

    this_shell = next_shell;
  }
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

  for (unsigned i = 0; i < new_shell.size(); i++)
  {
    for (int j = 0; j < vertices[new_shell[i]]->outlinks->size(); j++)
    {
      if (already_gone.insert(vertices[new_shell[i]]->outlinks->l[j]).second)
      {
        distances_node.emplace_back(shell, vertices[new_shell[i]]->outlinks->l[j]);
        next_shell.push_back(vertices[new_shell[i]]->outlinks->l[j]);
      }
    }

    for (int j = 0; j < vertices[new_shell[i]]->inlinks->size(); j++)
    {
      if (already_gone.insert(vertices[new_shell[i]]->inlinks->l[j]).second)
      {
        distances_node.emplace_back(shell, vertices[new_shell[i]]->inlinks->l[j]);
        next_shell.push_back(vertices[new_shell[i]]->inlinks->l[j]);
      }
    }
  }

  /*
  cout<<"new shell "<<shell<< std::endl;
  print_id(next_shell, std::cout);
  prints(ML);
  */

  if (!next_shell.empty())
  {
    if (shell >= int(ML.size()))
      ML.push_back(dim*step);

    reached += next_shell.size();
    ML[shell] += reached;

    return propagate_distances(next_shell, already_gone, distances_node, shell, ML, reached, step);
  }

  return 0;
}

inline int static_network::translate_anyway(std::deque<std::deque<int>> & ten)
{
  std::map<int, int> A;
  get_id_label(A);

  std::deque<std::deque<int>> ten2;

  for (auto& i : ten)
  {
    std::deque<int> ff;

    for (unsigned j = 0; j < i.size(); j++)
    {
      auto itf = A.find(i[j]);
      if (itf != A.end())
      {
        ff.push_back(itf->second);
      }
    }

    if (!ff.empty())
      ten2.push_back(ff);
  }

  ten = ten2;

  return 0;
}

inline int static_network::set_upper_network(std::map<int, std::map<int, std::pair<int, double>>> & neigh_weight_f, module_collection & module_coll)
{
  // loop on all the edges of the network...

  neigh_weight_f.clear();

  if (module_coll.empty())
    return 0;

  std::map<int, std::map<int, std::pair<double, double>> >neigh_weight_s;

  for (auto& module_b : module_coll.module_bs)
  {
    //cout<<"NAMES:: "<<m_name<< std::endl;
    std::map<int, std::pair<double, double>> neigh_weight;
    std::map<int, std::pair<int, double>> ooo;
    neigh_weight_s.emplace(module_b.first, neigh_weight);
    neigh_weight_f.emplace(module_b.first, ooo);
  }

  for (int i = 0; i < dim; i++)
  {
    std::set<int> & mem1 = module_coll.memberships[i];

    for (int j = 0; j < vertices[i]->outlinks->size(); j++) {
      int & neigh = vertices[i]->outlinks->l[j];
      std::set<int> & mem2 = module_coll.memberships[neigh];

      double denominator = mem1.size() * mem2.size();
      // I add a link between all different modules

      //cout<<"denomi "<<denominator<< std::endl;

      //**************************************************************************************************
      if (paras.weighted)
      {
        for (auto itk = mem1.begin(); itk != mem1.end(); ++itk)
        {
          for (auto itkk = mem2.begin(); itkk != mem2.end(); ++itkk)
          {
            if (*itk != *itkk)
            {
              int_histogram(
                *itkk,
                neigh_weight_s[*itk],
                double(vertices[i]->outlinks->w[j].first) / denominator,
                vertices[i]->out_original_weights[j] / denominator);
            }
          }
        }
      }
      else
      {
        for (auto itk = mem1.begin(); itk != mem1.end(); ++itk)
        {
          for (auto itkk = mem2.begin(); itkk != mem2.end(); ++itkk)
          {
            if (*itk != *itkk)
            {
              int_histogram(*itkk, neigh_weight_s[*itk], double(vertices[i]->outlinks->w[j].first) / denominator, double(vertices[i]->outlinks->w[j].first) / denominator);
            }
          }
        }
      }
      //**************************************************************************************************
    }
  }

  for (auto &[neighbor, weights] : neigh_weight_s)
  {
    for (auto &[w1, w2] : weights)
    {
      int intml = cast_int(w2.first);
      if (intml > 0)
      {
        neigh_weight_f[neighbor].emplace(w1, std::make_pair(intml, w2.second));
      }
    }
  }

  return 0;
}

inline void static_network::print_degree_of_homeless(std::deque<int> & homel, std::ostream & outt)
{
  std::deque<int> degree_of_homeless;
  for (int i : homel)
  {
    degree_of_homeless.push_back(vertices[i]->instub_number);
  }
  outt << "average in-degree of homeless nodes: " << average_func(degree_of_homeless)
    << " dev: " << sqrt(variance_func(degree_of_homeless)) << std::endl;

  degree_of_homeless.clear();
  for (int i : homel)
  {
    degree_of_homeless.push_back(vertices[i]->outstub_number);
  }
  outt << "average out-degree of homeless nodes: " << average_func(degree_of_homeless)
    << " dev: " << sqrt(variance_func(degree_of_homeless)) << std::endl;

  degree_of_homeless.clear();
  for (int i : homel)
  {
    degree_of_homeless.push_back(vertices[i]->instub_number + vertices[i]->outstub_number);
  }
  outt << "average in+out-degree of homeless nodes: " << average_func(degree_of_homeless)
    << " dev: " << sqrt(variance_func(degree_of_homeless)) << std::endl;
}

#endif
