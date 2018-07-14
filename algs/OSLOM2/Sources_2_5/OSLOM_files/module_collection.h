#pragma once

#ifndef MODULE_COLLECTION_H
#define MODULE_COLLECTION_H

#include "./standard_package/standard_include.cpp"

class module_collection
{
  /* all the labels refer to the index in int_matrix modules */
public:
  module_collection(int dim);
  ~module_collection() = default;

  int size() const { return module_bs.size(); }

  bool empty() const { return  module_bs.empty(); };

  bool insert(std::deque<int> & c, double bs, int & new_name);
  bool insert(std::deque<int> & c, double bs);
  bool erase(int);
  void print(std::ostream & outt, std::deque<int> & netlabels, bool);

  void fill_gaps();
  void put_gaps();
  void homeless(std::deque<int> & h);
  int coverage();
  int effective_groups();

  void set_partition(std::deque<std::deque<int>> & A);
  void set_partition(std::deque<std::deque<int>> & A, std::deque<double> & b);

  void compute_inclusions();
  void erase_included();
  bool almost_equal(int module_id, std::deque<int> & smaller);

  void compact();
  void sort_modules(std::deque<int> &);
  void merge(std::deque<int> & c);

  /*************************** DATA ***************************/

  std::deque<std::set<int>> memberships;
  int_matrix modules;
  std::map<int, double> module_bs;						/* it maps the module id into the b-score */

  /***********************************************************/

private:
  void _set_(int dim);
  bool check_already(const std::deque<int> & c);
  bool erase_first_shell(std::map<int, std::deque<int>> & erase_net);
  bool egomodules_to_merge(std::deque<int> & egom, std::deque<int> & smaller);
};

inline module_collection::module_collection(int dim)
{
  _set_(dim);
}

inline void module_collection::_set_(int dim)
{
  const std::set<int> first;
  for (int i = 0; i < dim; i++)
  {
    memberships.push_back(first);
  }
}

inline bool module_collection::insert(std::deque<int> & c, double bs)
{
  int new_name;
  return insert(c, bs, new_name);
}

inline bool module_collection::insert(std::deque<int> & c, double bs, int & new_name)
{
  if (bs == 0)
  {
    bs = ran4() * 1e-100;
  }

  std::sort(c.begin(), c.end());
  new_name = -1;

  if (check_already(c)) {
    new_name = modules.size();
    for (int i : c)
    {
      memberships[i].insert(new_name);
    }
    modules.push_back(c);
    module_bs[new_name] = bs;
    return true;
  }
  return false;
}

inline bool module_collection::erase(int a)
{
  // it erases module a
  if (module_bs.find(a) == module_bs.end())		// it only erases not empty modules
  {
    return false;
  }
  std::deque<int> & nodes_a = modules[a];
  for (int i : nodes_a)
  {
    memberships[i].erase(a);
  }

  modules[a].clear();
  module_bs.erase(a);
  return true;
}

inline void module_collection::print(std::ostream & outt, std::deque<int> & netlabels, bool not_homeless)
{
  int nmod = 0;
  for (auto& module_b : module_bs)
  {
    if (!not_homeless
      || modules[module_b.first].size() > 1)
    {
      nmod++;
      std::deque<int> & module_nodes = modules[module_b.first];
      outt
        << "#module " << module_b.first
        << " size: " << modules[module_b.first].size()
        << " bs: " << module_bs[module_b.first] << std::endl;

      std::deque<int> labseq;
      for (int module_node : module_nodes)
      {
        labseq.push_back(netlabels[module_node]);
      }

      std::sort(labseq.begin(), labseq.end());

      for (int i : labseq)
      {
        outt << i << " ";
      }
      outt << std::endl;
    }
  }
}

inline void module_collection::fill_gaps()
{
  for (int i = 0; i<int(memberships.size()); i++)
  {
    if (memberships[i].empty())
    {
      std::deque<int> new_d;
      new_d.push_back(i);
      insert(new_d, 1.);
    }
  }
}

inline void module_collection::put_gaps()
{
  std::deque<int> to_erase;
  for (int i = 0; i<int(modules.size()); i++)
  {
    if (modules[i].size() == 1)
    {
      to_erase.push_back(i);
    }
  }

  for (int i : to_erase)
  {
    erase(i);
  }
}

//*/
inline void module_collection::homeless(std::deque<int> & h)
{
  h.clear();

  for (int i = 0; i < int(memberships.size()); i++)
  {
    if (memberships[i].empty())
    {
      h.push_back(i);
    }
  }

  for (auto& module : modules)
  {
    if (module.size() == 1)
    {
      h.push_back(module[0]);
    }
  }
  std::sort(h.begin(), h.end());
}

inline int module_collection::coverage()
{
  // this function returns the number of nodes which are covered by at least one module
  int cov = 0;
  for (auto& membership : memberships)
  {
    if (!membership.empty())
    {
      cov++;
    }
  }
  return cov;
}

inline int module_collection::effective_groups()
{
  int nmod = 0;
  for (auto& module_b : module_bs)
  {
    if (modules[module_b.first].size() > 1)
    {
      nmod++;
    }
  }
  return nmod;
}

inline void module_collection::set_partition(std::deque<std::deque<int>> & A)
{
  A.clear();
  for (auto& module_b : module_bs)
  {
    if (modules[module_b.first].size() > 1)
    {
      A.push_back(modules[module_b.first]);
    }
  }
}

inline void module_collection::set_partition(std::deque<std::deque<int>> & A, std::deque<double> & b)
{
  A.clear();
  b.clear();

  for (auto& module_b : module_bs)
  {
    if (modules[module_b.first].size() > 1)
    {
      A.push_back(modules[module_b.first]);
      b.push_back(module_bs[module_b.first]);
    }
  }
}

inline bool module_collection::check_already(const std::deque<int> & c)
{
  // returns false if the module is already present
  std::map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)

  for (int i : c)
  {
    for (auto membership : memberships[i])
    {
      int_histogram(membership, com_ol);
    }
  }

  for (auto& item : com_ol)
  {
    if (item.second == int(c.size())
      && item.second == int(modules[item.first].size()))
    {
      return false;
    }
  }
  return true;
}

inline void module_collection::compute_inclusions()
{
  put_gaps();
  erase_included();
  compact();
}

inline void module_collection::erase_included()
{
  std::map<int, std::deque<int>> erase_net;
  for (auto& module_b : module_bs)
  {
    std::deque<int> smaller;
    almost_equal(module_b.first, smaller);
    erase_net[module_b.first] = smaller;
  }

  while (true)
  {
    if (!erase_first_shell(erase_net))
      break;
  }
}

inline bool module_collection::erase_first_shell(std::map<int, std::deque<int>> & erase_net)
{
  bool again = false;
  std::set<int> roots;

  for (auto& module_b : module_bs)
  {
    roots.insert(module_b.first);
  }

  for (auto& itm : erase_net)
  {
    std::deque<int> & smaller = itm.second;
    for (int i : smaller)
    {
      roots.erase(i);
    }
  }

  //cout<<"roots:"<< std::endl;
  //prints(roots);

  for (auto root : roots)
  {
    std::deque<int> & smaller = erase_net[root];
    for (int i : smaller)
    {
      if (module_bs.find(i) != module_bs.end())
      {
        erase(i);
        erase_net.erase(i);
        again = true;
      }
    }
  }
  return again;
}

inline bool module_collection::almost_equal(int module_id, std::deque<int> & smaller)
{
  // c is the module you want to know about
  // smaller is set to contain the module ids contained by module_id
  smaller.clear();
  std::deque<int> & c = modules[module_id];
  std::map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)
  for (int i : c)
  {
    for (auto itj : memberships[i])
    {
      int_histogram(itj, com_ol);
    }
  }

  for (auto& itm : com_ol)
  {
    if (itm.first != module_id
      && modules[itm.first].size() <= c.size())
    {
      const unsigned& other_size = modules[itm.first].size();

      if (double(itm.second) / other_size >= paras.coverage_inclusion_module_collection)
      {
        if (c.size() > other_size)
        {
          smaller.push_back(itm.first);
        }
        else if (c.size() == other_size && module_bs[module_id] < module_bs[itm.first])
        {
          smaller.push_back(itm.first);
        }
      }
    }
  }
  return true;
}

inline void module_collection::compact()
{
  /* this function is used to have continuos ids */
  put_gaps();

  std::map<int, int> from_old_index_to_new;
  {
    std::deque<std::deque<int>> modules2;
    std::map<int, double> module_bs2;

    for (auto& module_b : module_bs)
    {
      from_old_index_to_new.emplace(module_b.first, from_old_index_to_new.size());
      modules2.push_back(modules[module_b.first]);
      module_bs2[from_old_index_to_new.size() - 1] = module_b.second;
    }
    modules = modules2;
    module_bs = module_bs2;
  }

  for (auto& membership : memberships)
  {
    std::set<int> first;
    for (auto its : membership)
    {
      first.insert(from_old_index_to_new[its]);
    }
    membership = first;
  }
}

inline void module_collection::sort_modules(std::deque<int> & module_order)
{
  module_order.clear();
  std::multimap<double, int> rank_id;		/* modules are sorted from the biggest to the smallest. if they have equal size, we look at the score */

  for (auto& module_b : module_bs)
  {
    //std::cout<<modules[itm->first].size()<<" ... "<< std::endl;
    rank_id.insert(
      std::make_pair(-double(modules[module_b.first].size()) + 1e-2 * module_b.second, module_b.first));
  }
  //std::cout<<"rank_id"<< std::endl;
  //prints(rank_id);
  for (auto& itm : rank_id)
  {
    module_order.push_back(itm.second);
  }
}

inline bool module_collection::egomodules_to_merge(std::deque<int> & egom, std::deque<int> & smaller)
{
  // egom is the module you want to know about
  // smaller is set to contain the module ids to merge with egom
  smaller.clear();
  std::map<int, int> com_ol;		// it maps the index of the modules into the overlap (overlap=numeber of overlapping nodes)

  for (int i : egom)
  {
    for (auto itj : memberships[i])
    {
      int_histogram(itj, com_ol);
    }
  }

  //cout<<"egomodules_to_merge"<< std::endl;
  //prints(egom);

  for (auto& itm : com_ol)
  {
    //cout<<" other group "<<itm->second<< std::endl;
    //prints(modules[itm->first]);

    const unsigned& other_size = std::min(modules[itm.first].size(), egom.size());
    if (double(itm.second) / other_size >= paras.coverage_inclusion_module_collection)
    {
      smaller.push_back(itm.first);
    }
  }
  return true;
}

inline void module_collection::merge(std::deque<int>& c)
{
  std::deque<int> to_merge;
  egomodules_to_merge(c, to_merge);

  //cout<<"module c: "<< std::endl;
  //prints(c);

  if (to_merge.empty())
  {
    insert(c, 1.);
  }
  else
  {
    for (int i : to_merge)
    {
      //cout<<"to_merge"<< std::endl;
      //prints(modules[to_merge[i]]);

      std::set<int> si;
      deque_to_set_app(modules[to_merge[i]], si);
      deque_to_set_app(c, si);
      erase(i);

      //prints(si);

      std::deque<int> to_insert;
      set_to_deque(si, to_insert);
      insert(to_insert, 1.);
    }
  }
  erase_included();
}
#endif // MODULE_COLLECTION_H
