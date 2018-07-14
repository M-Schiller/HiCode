#pragma once

#ifndef DIRECTED_OSLOMNET_EVALUATE_H
#define DIRECTED_OSLOMNET_EVALUATE_H

#include <deque>
#include <string>

inline double log_zero(double a)
{
  if (a <= 0)
    return -1e20;
  return log(a);
}

class oslomnet_evaluate : public oslomnet_louvain
{
public:
  oslomnet_evaluate(std::deque<std::deque<int>> & b, std::deque<std::deque<std::pair<int, double>>> & c, std::deque<int> & d) : oslomnet_louvain()
  {
    set_graph(b, c, d);
    set_maxbord();
    set_changendi_cum();
  }

  oslomnet_evaluate(std::string a) : oslomnet_louvain()
  {
    set_graph(a);
    set_maxbord();
    set_changendi_cum();
  }

  oslomnet_evaluate(std::map<int, std::map<int, std::pair<int, double>> > & A) : oslomnet_louvain()
  {
    set_graph(A);
    set_maxbord();
    set_changendi_cum();
  }

  ~oslomnet_evaluate() {};

  double CUP_both(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int);
  double CUP_check(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int);
  double group_inflation(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int);
  int try_to_assign_homeless_help(module_collection & module_coll, std::map<int, std::deque<int>> & to_check);

private:
  double CUP_iterative(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int);
  double CUP_search(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int);
  void erase_cgroup(int wnode);
  void insert_cgroup(int wnode);
  bool erase_the_worst(int & wnode);

  int set_maxbord();
  void set_cgroup_and_neighs(const std::deque<int> & G);
  double all_external_test(int kout_g_in, int tmin, int kout_g_out, int tmout,
    int Nstar, int nneighs, const double & max_r_one, const double & maxr_two, std::deque<int> & gr_cleaned, bool only_c, weighted_tabdeg & previous_tab_c);
  double cup_on_list(cup_data_struct & a, std::deque<int> & gr_cleaned);
  void get_external_scores(weighted_tabdeg & neighs, cup_data_struct & fitness_label_to_sort, int kout_g_in, int tmin, int kout_g_out, int tmout,
    int Nstar, int nneighs, const double & max_r, bool only_c, weighted_tabdeg & previous_tab_c);

  double CUP_runs(weighted_tabdeg & previous_tab_c, weighted_tabdeg & previous_tab_n, int kin_cgroup_prev, int ktot_cgroup_prev_in, int ktot_cgroup_prev_out, std::deque<int> & gr_cleaned, bool only_c, int number_of_runs);
  void initialize_for_evaluation(const std::deque<int> & _c_, weighted_tabdeg & previous_tab_c, weighted_tabdeg & previous_tab_n, int & kin_cgroup_prev, int & ktot_cgroup_prev_in, int & ktot_cgroup_prev_out);
  void initialize_for_evaluation(weighted_tabdeg & previous_tab_c, weighted_tabdeg & previous_tab_n, int & kin_cgroup_prev, int & ktot_cgroup_prev_in, int & ktot_cgroup_prev_out);
  double partial_CUP(weighted_tabdeg & previous_tab_c, weighted_tabdeg & previous_tab_n, int kin_cgroup_prev, int ktot_cgroup_prev_in, int ktot_cgroup_prev_out, std::deque<int> & border_group, bool only_c);
  void set_changendi_cum();

  void insertion(int changendi);
  bool insert_the_best();

  /* DATA ***************************************************/

  double  max_r_bord;								// this is the maximum r allowed for the external nodes (we don't want to look at all the graph, it would take too long)
  int maxb_nodes;									// this is the maximum number of nodes allowed in the border (similar as above)
  std::deque<double> changendi_cum;					// this is the cumulative distribution of the number of nodes to add to the cluster in the group_inflation function

  // ************* things to update *************************
  weighted_tabdeg cgroup;									//*
  weighted_tabdeg neighs;									//*
                              //*
  int kin_cgroup;											//*
  int ktot_cgroup_in;										//*
  int ktot_cgroup_out;									//*
  /*********************************************************/
};

inline void oslomnet_evaluate::set_changendi_cum()
{
  if (dim != 0
    && oneM != 0)
  {
    int flat_until = cast_int(oneM / dim * 3);
    flat_until = std::min(dim / 2, flat_until);

    int max_p = std::max(paras.CUT_Off, flat_until);		// this is something which might be optimized
    max_p = std::min(dim / 2, max_p);

    powerlaw(max_p, flat_until + 1, 3, changendi_cum);
    std::deque<double> distr;
    distribution_from_cumulative(changendi_cum, distr);
    double ac = 1;

    if (!distr.empty())
      ac = distr[0];

    for (int i = 0; i < flat_until; i++)
      distr.push_front(ac);

    normalize_one(distr);
    cumulative_from_distribution(changendi_cum, distr);
  }
}

inline int oslomnet_evaluate::set_maxbord()
{
  max_r_bord = paras.maxbg_ordinary;
  maxb_nodes = paras.maxborder_nodes;
  return 0;
}

inline void oslomnet_evaluate::erase_cgroup(int wnode)
{
  std::map<int, facts>::iterator itm = cgroup.lab_facts.find(wnode);

  if (itm != cgroup.lab_facts.end()) {
    int kpin = itm->second.internal_indegree;
    int ktin = itm->second.indegree;
    int kpout = itm->second.internal_outdegree;
    int ktout = itm->second.outdegree;
    double mtlwin = itm->second.minus_log_total_wrin;
    double mtlwout = itm->second.minus_log_total_wrout;

    kin_cgroup -= kpin + kpout;
    ktot_cgroup_out -= ktout;
    ktot_cgroup_in -= ktin;

    int kout_g_in = ktot_cgroup_in - kin_cgroup;
    int kout_g_out = ktot_cgroup_out - kin_cgroup;

    int tmin = oneM - ktot_cgroup_in;
    int tmout = oneM - ktot_cgroup_out;

    double fi = compute_global_fitness_ofive(kpin, kout_g_in, kpout, kout_g_out, tmin, tmout, ktin, ktout, mtlwin, mtlwout, neighs.size() + 1, dim - cgroup.size() + 1);
    neighs.edinsert(wnode, kpin, kpout, ktin, ktout, mtlwin, mtlwout, fi);

    //cout<<"node erased: "<<id_of(wnode)<< std::endl;
    cgroup.erase(wnode);

    std::deque<int> tobe;
    std::pair <int, std::pair<int, double>> OPA;

    for (int i = 0; i < vertices[wnode]->inlinks->size(); i++)
    {
      OPA = vertices[wnode]->outlinks->posweightof(vertices[wnode]->inlinks->l[i]);

      if (!cgroup.update_group(
        vertices[wnode]->inlinks->l[i],
        -vertices[wnode]->inlinks->w[i].first,
        -OPA.second.first,
        -vertices[wnode]->inlinks->w[i].second,
        -OPA.second.second,
        dim - cgroup.size(),
        neighs.size(),
        kout_g_in,
        kout_g_out,
        tmin,
        tmout,
        vertices[vertices[wnode]->inlinks->l[i]]->instub_number,
        vertices[vertices[wnode]->inlinks->l[i]]->outstub_number,
        tobe))
        neighs.update_neighs(vertices[wnode]->inlinks->l[i], -vertices[wnode]->inlinks->w[i].first, -OPA.second.first, -vertices[wnode]->inlinks->w[i].second, -OPA.second.second,
          dim - cgroup.size(), kout_g_in, kout_g_out, tmin, tmout, vertices[vertices[wnode]->inlinks->l[i]]->instub_number, vertices[vertices[wnode]->inlinks->l[i]]->outstub_number);
    }

    for (int i = 0; i < vertices[wnode]->outlinks->size(); i++)
    {
      OPA = vertices[wnode]->inlinks->posweightof(vertices[wnode]->outlinks->l[i]);

      if (OPA.first == -1)
      {
        if (!cgroup.update_group(
          vertices[wnode]->outlinks->l[i],
          0,
          -vertices[wnode]->outlinks->w[i].first,
          0,
          -vertices[wnode]->outlinks->w[i].second,
          dim - cgroup.size(),
          neighs.size(),
          kout_g_in, kout_g_out,
          tmin,
          tmout,
          vertices[vertices[wnode]->outlinks->l[i]]->instub_number,
          vertices[vertices[wnode]->outlinks->l[i]]->outstub_number,
          tobe))
        {
          neighs.update_neighs(
            vertices[wnode]->outlinks->l[i],
            0,
            -vertices[wnode]->outlinks->w[i].first,
            0,
            -vertices[wnode]->outlinks->w[i].second,
            dim - cgroup.size(),
            kout_g_in,
            kout_g_out,
            tmin,
            tmout,
            vertices[vertices[wnode]->outlinks->l[i]]->instub_number,
            vertices[vertices[wnode]->outlinks->l[i]]->outstub_number);
        }
      }
    }
    for (int i : tobe)
      erase_cgroup(i);
  }
}

inline bool oslomnet_evaluate::erase_the_worst(int & wnode)
{
  // this function is to look for the worst node in cgroup and to erase it
  int Nstar = dim - cgroup.size();
  int nn = neighs.size();
  int kout_g_in = ktot_cgroup_in - kin_cgroup;
  int kout_g_out = ktot_cgroup_out - kin_cgroup;
  int tmin = oneM - ktot_cgroup_in;
  int tmout = oneM - ktot_cgroup_out;

  double wf;

  cgroup.worst_node(wnode, wf, kout_g_in, kout_g_out, Nstar, nn, tmin, tmout);

  if (cgroup.empty())
  {
    return false;
  }
  erase_cgroup(wnode);

  return true;
}

inline void oslomnet_evaluate::insert_cgroup(int wnode)
{
  /* this function is to insert benode into cgroup  updating all the system, neighs - kin_cgroup - ktot_cgroup */

  int kpin, ktin, kpout, ktout;
  double mtlwin, mtlwout;
  {
    const auto itm = neighs.lab_facts.find(wnode);
    if (itm != neighs.lab_facts.end())
    {
      kpin = itm->second.internal_indegree;
      ktin = itm->second.indegree;

      kpout = itm->second.internal_outdegree;
      ktout = itm->second.outdegree;

      mtlwin = itm->second.minus_log_total_wrin;
      mtlwout = itm->second.minus_log_total_wrout;
    }
    else
    {
      kpin = 0;
      kpout = 0;
      ktin = vertices[wnode]->instub_number;
      ktout = vertices[wnode]->outstub_number;
      mtlwin = 0;
      mtlwout = 0;
    }
  }

  int kout_g_in = ktot_cgroup_in - kin_cgroup;
  int kout_g_out = ktot_cgroup_out - kin_cgroup;
  int tmin = oneM - ktot_cgroup_in;
  int tmout = oneM - ktot_cgroup_out;

  double fi = compute_global_fitness_ofive(kpin, kout_g_in, kpout, kout_g_out, tmin, tmout, ktin, ktout, mtlwin, mtlwout, neighs.size(), dim - cgroup.size());

  kin_cgroup += kpin + kpout;
  ktot_cgroup_in += ktin;
  ktot_cgroup_out += ktout;
  kout_g_in = ktot_cgroup_in - kin_cgroup;
  kout_g_out = ktot_cgroup_out - kin_cgroup;
  tmin = oneM - ktot_cgroup_in;
  tmout = oneM - ktot_cgroup_out;

  cgroup.edinsert(wnode, kpin, kpout, ktin, ktout, mtlwin, mtlwout, fi);
  neighs.erase(wnode);

  std::deque<int> tobe;
  std::pair <int, std::pair<int, double>> OPA;

  for (int i = 0; i < vertices[wnode]->inlinks->size(); i++)
  {
    OPA = vertices[wnode]->outlinks->posweightof(vertices[wnode]->inlinks->l[i]);

    //cout<<"------> "<<vertices[wnode]->inlinks->w[i].second<<" "<<OPA.second.second<< std::endl;

    if (!cgroup.update_group(
      vertices[wnode]->inlinks->l[i],
      vertices[wnode]->inlinks->w[i].first,
      OPA.second.first,
      vertices[wnode]->inlinks->w[i].second,
      OPA.second.second,
      dim - cgroup.size(),
      neighs.size(),
      kout_g_in,
      kout_g_out,
      tmin,
      tmout,
      vertices[vertices[wnode]->inlinks->l[i]]->instub_number,
      vertices[vertices[wnode]->inlinks->l[i]]->outstub_number,
      tobe))
    {
      neighs.update_neighs(
        vertices[wnode]->inlinks->l[i],
        vertices[wnode]->inlinks->w[i].first,
        OPA.second.first,
        vertices[wnode]->inlinks->w[i].second,
        OPA.second.second,
        dim - cgroup.size(),
        kout_g_in,
        kout_g_out,
        tmin,
        tmout,
        vertices[vertices[wnode]->inlinks->l[i]]->instub_number,
        vertices[vertices[wnode]->inlinks->l[i]]->outstub_number);
    }
  }

  for (auto i = 0; i < vertices[wnode]->outlinks->size(); i++)
  {
    OPA = vertices[wnode]->inlinks->posweightof(vertices[wnode]->outlinks->l[i]);

    if (OPA.first == -1)
    {
      if (!cgroup.update_group(
        vertices[wnode]->outlinks->l[i],
        0,
        vertices[wnode]->outlinks->w[i].first,
        0,
        vertices[wnode]->outlinks->w[i].second,
        dim - cgroup.size(),
        neighs.size(),
        kout_g_in,
        kout_g_out,
        tmin,
        tmout,
        vertices[vertices[wnode]->outlinks->l[i]]->instub_number,
        vertices[vertices[wnode]->outlinks->l[i]]->outstub_number,
        tobe))
      {
        neighs.update_neighs(
          vertices[wnode]->outlinks->l[i],
          0,
          vertices[wnode]->outlinks->w[i].first,
          0, vertices[wnode]->outlinks->w[i].second,
          dim - cgroup.size(),
          kout_g_in,
          kout_g_out,
          tmin,
          tmout,
          vertices[vertices[wnode]->outlinks->l[i]]->instub_number,
          vertices[vertices[wnode]->outlinks->l[i]]->outstub_number);
      }
    }
  }
}

inline void oslomnet_evaluate::set_cgroup_and_neighs(const std::deque<int> & G)
{
  kin_cgroup = 0;
  ktot_cgroup_in = 0;
  ktot_cgroup_out = 0;
  cgroup.clear();
  neighs.clear();
  for (int group : G)
  {
    insert_cgroup(group);
  }
}

inline double oslomnet_evaluate::cup_on_list(cup_data_struct & a, std::deque<int> & gr_cleaned)
{
  int Nstar;
  if (paras.weighted)
    Nstar = neighs.size();
  else
    Nstar = dim - cgroup.size();

  const double critical_xi = -log(1 - paras.threshold) / fitted_exponent(Nstar);
  int pos = Nstar;

  int until = -1;								// until tells how many nodes should be included into the cluster - actually the number of good nodes are (until +1)
  double probability_a, probability_b;		// these are the two extremes of a possible good node I could have found
  double c_min = 1;								// this is the score we give to the border we are evaluating here
  //cout<<"critical_xi: "<<critical_xi<<" --------------------------------------- "<<neighs.size()<<" cgroup "<<cgroup.size()<< std::endl<< std::endl<< std::endl;

  cup_data_struct::iterator itl;

  for (itl = a.begin(); itl != a.end(); ++itl)
  {
    double c_pos = order_statistics_left_cumulative(Nstar, pos, itl->first);

    //cout<<"position .... "<<pos<<" "<<order_statistics_left_cumulative(Nstar, pos, itl->first + itl->second.second)<<" >>AAA<< "<<Nstar<<" +++ "<<itl->first + itl->second.second<<" "<<itl->first - itl->second.second<< std::endl;
    c_min = std::min(c_pos, c_min);

    if (c_pos < critical_xi)
    {
      /*
      this is the basic condition of the order statistics test
      it's saying: look, this guy (itl->second.first) has an average fitness (itl->first)
      whose order_statistics_left_cumulative is below the threshold
      */
      if (until == -1)			// this node is the first node to be below the threshold
      {
        until = Nstar - pos;
        c_min = c_pos;
        probability_a = itl->first - itl->second.second;
        probability_b = itl->first + itl->second.second;
      }
      else
      {
        /*
        the previous node was already below the threshold.
        In this case I need to know if I should stop now or go on.
        The condition is related to the probability_to_overtake the previous guy
        */

        const double probability_to_overtake = compare_r_variables(probability_a, probability_b, itl->first - itl->second.second, itl->first + itl->second.second);

        //cout<<"probability_to_overtake: "<<probability_to_overtake<< std::endl;

        if (probability_to_overtake > 0.4999)/*preliminary check: this node is basically equivalent to the previous guy, I consider it good*/
        {
          until = Nstar - pos;
          c_min = c_pos;
          probability_a = itl->first - itl->second.second;
          probability_b = itl->first + itl->second.second;
        }
        else
        {
          /*now I need to compute the bootstrap probability that the previous guy would have stopped the process*/
          if ((probability_to_overtake == 0)
            || ((1. - probability_to_overtake) * compute_probability_to_stop(probability_a, probability_b, critical_xi, Nstar, pos + 1) > 0.5001))
          {
            if (equivalent_check_gather(a, until, probability_a, probability_b, Nstar, critical_xi))
            {
              break;
            }
          }
          until = Nstar - pos;
          c_min = c_pos;
          probability_a = itl->first - itl->second.second;
          probability_b = itl->first + itl->second.second;
        }
      }
    }
    else		/* this node is not below the threshold */
    {
      if (until != -1)		/* this means that this node is not good and the previous one was good. So, I stop here */
      {
        if (equivalent_check_gather(a, until, probability_a, probability_b, Nstar, critical_xi))
          break;
      }
    }
    --pos;
  }

  // equalizer check
  // this check is important to see if the procedure stopped just because there were a lot of equivalents nodes
  if (until != -1
    && itl == a.end())
  {
    equivalent_check_gather(a, until, probability_a, probability_b, Nstar, critical_xi);
  }
  // inserting nodes in gr_cleaned
  int nodes_added = -1;
  itl = a.begin();

  while (itl != a.end())
  {
    if (nodes_added == until)
    {
      break;
    }
    gr_cleaned.push_back(itl->second.first);
    ++itl;
    ++nodes_added;
  }
  return pron_min_exp(Nstar, c_min);
}

inline double oslomnet_evaluate::all_external_test(
  int kout_g_in,
  int tmin,
  int kout_g_out,
  int tmout,
  int Nstar,
  int nneighs,
  const double & max_r_one,
  const double & maxr_two,
  std::deque<int> & gr_cleaned,
  bool only_c,
  weighted_tabdeg & previous_tab_c)
{
  const double max_r = std::min(max_r_one, maxr_two);

  cup_data_struct fitness_label_to_sort;

  get_external_scores(neighs, fitness_label_to_sort, kout_g_in, tmin, kout_g_out, tmout, Nstar, nneighs, max_r, only_c, previous_tab_c);

  return cup_on_list(fitness_label_to_sort, gr_cleaned);
}

inline void oslomnet_evaluate::initialize_for_evaluation(
  weighted_tabdeg & previous_tab_c,
  weighted_tabdeg & previous_tab_n,
  int & kin_cgroup_prev,
  int & ktot_cgroup_prev_in,
  int & ktot_cgroup_prev_out)
{
  int Nstar = dim - cgroup.size();
  int nn = neighs.size();
  int kout_g_in = ktot_cgroup_in - kin_cgroup;
  int kout_g_out = ktot_cgroup_out - kin_cgroup;

  int tmin = oneM - ktot_cgroup_in;
  int tmout = oneM - ktot_cgroup_out;

  previous_tab_c.set_and_update_group(Nstar, nn, kout_g_out, tmout, kout_g_in, tmin, cgroup);
  previous_tab_n.set_and_update_neighs(Nstar, nn, kout_g_out, tmout, kout_g_in, tmin, neighs);

  kin_cgroup_prev = kin_cgroup;
  ktot_cgroup_prev_in = ktot_cgroup_in;
  ktot_cgroup_prev_out = ktot_cgroup_out;
}

inline void oslomnet_evaluate::initialize_for_evaluation(
  const std::deque<int> & _c_,
  weighted_tabdeg & previous_tab_c,
  weighted_tabdeg & previous_tab_n,
  int & kin_cgroup_prev, int & ktot_cgroup_prev_in,
  int & ktot_cgroup_prev_out)
{
  set_cgroup_and_neighs(_c_);

  int Nstar = dim - cgroup.size();
  int nn = neighs.size();
  int kout_g_in = ktot_cgroup_in - kin_cgroup;
  int kout_g_out = ktot_cgroup_out - kin_cgroup;

  int tmin = oneM - ktot_cgroup_in;
  int tmout = oneM - ktot_cgroup_out;

  previous_tab_c.set_and_update_group(Nstar, nn, kout_g_out, tmout, kout_g_in, tmin, cgroup);
  previous_tab_n.set_and_update_neighs(Nstar, nn, kout_g_out, tmout, kout_g_in, tmin, neighs);

  kin_cgroup_prev = kin_cgroup;
  ktot_cgroup_prev_in = ktot_cgroup_in;
  ktot_cgroup_prev_out = ktot_cgroup_out;
}

inline double oslomnet_evaluate::partial_CUP(
  weighted_tabdeg & previous_tab_c,
  weighted_tabdeg & previous_tab_n,
  int kin_cgroup_prev,
  int ktot_cgroup_prev_in,
  int ktot_cgroup_prev_out,
  std::deque<int> & border_group,
  bool only_c) {
  //draw_with_weight_probability("graph____");

  // still there is some stochasticity due to possible ties

  /*	previous_stuff is the module-stuff before the CUP (Clean Up Procedure)
    cgroup + border_group is the module cleaned					*/

  border_group.clear();

  //previous_tab_c.print_nodes(cout);

  cgroup._set_(previous_tab_c);
  neighs._set_(previous_tab_n);
  kin_cgroup = kin_cgroup_prev;
  ktot_cgroup_in = ktot_cgroup_prev_in;
  ktot_cgroup_out = ktot_cgroup_prev_out;

  if (cgroup.size() == dim)
  {
    return 1;
  }

  double bscore = 1;
  while (true)
  {
    bscore = all_external_test(ktot_cgroup_in - kin_cgroup, oneM - ktot_cgroup_in, ktot_cgroup_out - kin_cgroup, oneM - ktot_cgroup_out,
      dim - cgroup.size(), neighs.size(), maxb_nodes / double(dim - cgroup.size()), max_r_bord, border_group, only_c, previous_tab_c);

    //cout<<cgroup.size()<< std::endl;

    if (!border_group.empty())
      break;

    if (cgroup.empty())
      break;

    int wnode;
    erase_the_worst(wnode);
  }
  return bscore;
}

inline double oslomnet_evaluate::CUP_runs(
  weighted_tabdeg & previous_tab_c,
  weighted_tabdeg & previous_tab_n,
  int kin_cgroup_prev,
  int ktot_cgroup_prev_in,
  int ktot_cgroup_prev_out,
  std::deque<int> & gr_cleaned,
  bool only_c,
  int number_of_runs)
{
  /* this if statemets are here to speed up the program if there are big clusters */
  if (previous_tab_c.size() > 100000)
    number_of_runs = 3;
  else if (previous_tab_c.size() > 10000)
    number_of_runs = 5;
  else if (previous_tab_c.size() > 1000)
    number_of_runs = 10;

  gr_cleaned.clear();

  int max_gr_size = 0;
  double bscore = 1;

  int good_runs = 0;

  if (previous_tab_c.empty())
    return 1;

  for (int i = 0; i < number_of_runs; i++)
  {
    std::deque<int> gr_run_i;
    const double score_i = partial_CUP(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out, gr_run_i, only_c);

    if (cgroup.size() + int(gr_run_i.size()) > max_gr_size)
    {
      bscore = score_i;
      cgroup.set_deque(gr_cleaned);
      for (int j : gr_run_i)
      {
        gr_cleaned.push_back(j);
      }

      max_gr_size = gr_cleaned.size();
      sort(gr_cleaned.begin(), gr_cleaned.end());
    }

    if (!gr_run_i.empty())
    {
      ++good_runs;
      if (good_runs >= 0.55 * number_of_runs)
        return bscore;
    }
  }

  if (good_runs < 0.55 * number_of_runs)
  {
    gr_cleaned.clear();
    bscore += paras.threshold;
    bscore = std::min(1., bscore);
  }
  return bscore;
}

inline double oslomnet_evaluate::CUP_check(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int number_of_runs = paras.clean_up_runs)
{
  /*_c_ is the module to clean up and gr_cleaned is the result */

  weighted_tabdeg previous_tab_c;
  weighted_tabdeg previous_tab_n;
  int kin_cgroup_prev;
  int ktot_cgroup_prev_in, ktot_cgroup_prev_out;

  initialize_for_evaluation(_c_, previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out);

  return CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out, gr_cleaned, true, number_of_runs);;
}

inline double oslomnet_evaluate::CUP_search(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int number_of_runs = paras.clean_up_runs)
{
  /*_c_ is the module to clean up and gr_cleaned is the result */

  weighted_tabdeg previous_tab_c;
  weighted_tabdeg previous_tab_n;
  int kin_cgroup_prev;
  int ktot_cgroup_prev_in, ktot_cgroup_prev_out;

  initialize_for_evaluation(_c_, previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out);

  return CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out, gr_cleaned, false, number_of_runs);
}

inline double oslomnet_evaluate::CUP_both(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int number_of_runs = paras.clean_up_runs)
{
  /*_c_ is the module to clean up and gr_cleaned is the result */

  weighted_tabdeg previous_tab_c;
  weighted_tabdeg previous_tab_n;
  int kin_cgroup_prev;
  int ktot_cgroup_prev_in, ktot_cgroup_prev_out;

  initialize_for_evaluation(_c_, previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out);

  CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out, gr_cleaned, false, number_of_runs);

  initialize_for_evaluation(_c_, previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out);

  /* this "truely" means I can only look at nodes in previous_tab_c */
  return CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out, gr_cleaned, true, number_of_runs);
}

inline double oslomnet_evaluate::CUP_iterative(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int number_of_runs = paras.clean_up_runs)
{
  double bs = CUP_both(_c_, gr_cleaned, number_of_runs);
  int stopp = 0;

  do
  {
    std::deque<int> _c_temp = gr_cleaned;
    bs = CUP_search(_c_temp, gr_cleaned, number_of_runs);
    ++stopp;

    if (stopp == paras.iterative_stopper)
      break;
  } while (gr_cleaned.size() > _c_.size());
  return bs;
}

inline bool oslomnet_evaluate::insert_the_best()
{
  int Nstar = dim - cgroup.size();
  int nn = neighs.size();

  int kout_g_in = ktot_cgroup_in - kin_cgroup;
  int kout_g_out = ktot_cgroup_out - kin_cgroup;
  int tmin = oneM - ktot_cgroup_in;
  int tmout = oneM - ktot_cgroup_out;

  double lowest_r;
  int benode;
  neighs.best_node(benode, lowest_r, kout_g_in, kout_g_out, Nstar, nn, tmin, tmout);

  if (benode == -1)
    return false;
  insert_cgroup(benode);
  return true;
}

inline void oslomnet_evaluate::insertion(int changendi)
{
  for (int i = 0; i < changendi; i++)
  {
    insert_the_best();
  }
}

inline double oslomnet_evaluate::group_inflation(const std::deque<int> & _c_, std::deque<int> & gr_cleaned, int number_of_runs = paras.inflate_runs)
{
  /* preliminary check */
  //cout<<"evaluating group of "<<_c_.size()<<" nodes"<< std::endl;
  double bscore = CUP_iterative(_c_, gr_cleaned);

  if (!gr_cleaned.empty()) {
    //cout<<"test passed, bscore: "<<bscore<<"  size: "<<gr_cleaned.size()<< std::endl;
    return bscore;
  }
  /* preliminary check */

  weighted_tabdeg _c_tab_c;
  weighted_tabdeg _c_tab_n;
  int kin_cgroup_c;
  int ktot_cgroup_c_in, ktot_cgroup_c_out;

  //print_id(_c_, std::cout);
  initialize_for_evaluation(_c_, _c_tab_c, _c_tab_n, kin_cgroup_c, ktot_cgroup_c_in, ktot_cgroup_c_out);

  weighted_tabdeg previous_tab_c;
  weighted_tabdeg previous_tab_n;
  int kin_cgroup_prev;
  int ktot_cgroup_prev_in, ktot_cgroup_prev_out;

  int stopper = 0;

  while (true)
  {
    cgroup._set_(_c_tab_c);
    neighs._set_(_c_tab_n);
    kin_cgroup = kin_cgroup_c;
    ktot_cgroup_in = ktot_cgroup_c_in;
    ktot_cgroup_out = ktot_cgroup_c_out;

    /*
    deque<int> hhh;
    cgroup.set_deque(hhh);
    print_id(hhh, std::cout);
    //*/

    int changendi = lower_bound(changendi_cum.begin(), changendi_cum.end(), ran4()) - changendi_cum.begin() + 1;
    changendi = std::min(changendi, neighs.size());
    insertion(changendi);

    /*
    std::cout<<"... after "<<cgroup.size()<< std::endl;
    cgroup.set_deque(hhh);
    print_id(hhh, std::cout);//*/

    if (cgroup.size() == dim)
      return 1;

    /*here it make a CUP_search using c_group with the nodes added*/
    initialize_for_evaluation(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out);

    CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out, gr_cleaned, false, number_of_runs);

    if (!gr_cleaned.empty())
    {
      /*the first clean up passed. now it makes the CUP_check*/

      initialize_for_evaluation(gr_cleaned, previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out);

      bscore = CUP_runs(previous_tab_c, previous_tab_n, kin_cgroup_prev, ktot_cgroup_prev_in, ktot_cgroup_prev_out, gr_cleaned, true, paras.clean_up_runs);

      if (!gr_cleaned.empty())
      {
        //cout<<"test passed, bscore: "<<bscore<<"  size: "<<gr_cleaned.size()<< std::endl;
        return bscore;
      }
    }
    ++stopper;
    if (stopper == paras.inflate_stopper)
      break;
  }
  return 1;
}

inline void oslomnet_evaluate::get_external_scores(
  weighted_tabdeg & neighs,
  cup_data_struct & fitness_label_to_sort,
  int kout_g_in,
  int tmin,
  int kout_g_out,
  int tmout,
  int Nstar,
  int nneighs,
  const double & max_r,
  bool only_c,
  weighted_tabdeg & previous_tab_c)
{
  int counter = 0;
  for (auto bit = neighs.fitness_lab.begin(); bit != neighs.fitness_lab.end(); ++bit)
  {
    const auto itm = neighs.lab_facts.find(bit->second);
    double interval;

    /***************************************************************************************
    std::cout<<"boring check but somebody ...."<< std::endl;
    std::cout<<"node: "<<vertices[itm->first]->id_num<< std::endl;
    deque<int> hhh;
    cgroup.set_deque(hhh);
    std::cout<<"group"<< std::endl;
    print_id(hhh, std::cout);

    std::cout<<"ktot_cgroup_in: "<<ktot_cgroup_in<<" ktot_cgroup_out: "<<ktot_cgroup_out<<" kin_cgroup: "<<kin_cgroup<< std::endl;
    std::cout<<"internal in out  "<<itm->second.internal_indegree<<" "<<itm->second.internal_outdegree<< std::endl;

    *-**************************************************************************************/

    double F = compute_global_fitness(itm->second.internal_indegree, kout_g_in, itm->second.internal_outdegree, kout_g_out, tmin, tmout,
      itm->second.indegree, itm->second.outdegree, itm->second.minus_log_total_wrin, itm->second.minus_log_total_wrout, nneighs, Nstar, interval);

    //cout<<"tmin: "<<tmin<<" tmout "<<tmout<<" "<<oneM<<" ... "<<oneM - ktot_cgroup_out<< std::endl;
    //cout<<"------------------------------------------------------------------------------------"<< std::endl;
    if (F > max_r) {
      counter++;
      if (counter > num_up_to)
        break;
    }
    else
    {
      if (!only_c
        || previous_tab_c.is_internal(itm->first))
        fitness_label_to_sort.emplace(F, std::make_pair(itm->first, interval));
    }
  }
}

inline int oslomnet_evaluate::try_to_assign_homeless_help(module_collection & module_coll, std::map<int, std::deque<int>> & to_check)
{
  to_check.clear();
  module_coll.put_gaps();

  //if(paras.print_cbs)
    //cout<<"checking homeless nodes "<< std::endl;

  std::deque<int> homel;
  module_coll.homeless(homel);

  const int before_procedure = homel.size();

  if (homel.empty())
    return before_procedure;

  /*cout<<"homel"<< std::endl;
  print_id(homel, std::cout);*/

  std::set<int> called;						// modules connected to homeless nodes
  std::map<int, std::set<int>> homel_module;		// maps the homeless node with the modules it's connected to

  for (unsigned i = 0; i < homel.size(); i++)
  {
    std::set<int> thish;
    for (int j = 0; j < vertices[homel[i]]->inlinks->size(); j++)
    {
      int & neigh = vertices[homel[i]]->inlinks->l[j];

      for (auto neighbor : module_coll.memberships[neigh])
      {
        called.insert(neighbor);
        thish.insert(neighbor);
      }
    }

    for (int j = 0; j < vertices[homel[i]]->outlinks->size(); j++)
    {
      int & neigh = vertices[homel[i]]->outlinks->l[j];

      for (auto neighbor : module_coll.memberships[neigh])
      {
        called.insert(neighbor);
        thish.insert(neighbor);
      }
    }

    if (!thish.empty())
      homel_module[homel[i]] = thish;
  }

  std::map<int, int> module_kin;
  std::map<int, int> module_ktotin;
  std::map<int, int> module_ktotout;

  for (auto its : called)
  {
    module_kin[its] = cast_int(kin_m(module_coll.modules[its]));
    module_ktotin[its] = cast_int(ktot_m(module_coll.modules[its]).first);
    module_ktotout[its] = cast_int(ktot_m(module_coll.modules[its]).second);
  }

  for (auto& itm : homel_module)
  {
    double cmin = 1.1;
    int belongs_to = -1;

    //cout<<"homeless node: "<<id_of(itm->first)<< std::endl;

    for (auto its = itm.second.begin(); its != itm.second.end(); ++its)
    {
      int kin_node_in = cast_int(vertices[itm.first]->kplus_m(module_coll.modules[*its]).first);
      int kin_node_out = cast_int(vertices[itm.first]->kplus_m(module_coll.modules[*its]).second);

      int kout_g_in = module_ktotin[*its] - module_kin[*its];
      int tmin = cast_int(oneM - module_ktotin[*its]);
      double kinw_in = vertices[itm.first]->kplus_w(module_coll.modules[*its]).first;

      int kout_g_out = module_ktotout[*its] - module_kin[*its];
      int tmout = cast_int(oneM - module_ktotout[*its]);
      double kinw_out = vertices[itm.first]->kplus_w(module_coll.modules[*its]).second;

      double rh1 = compute_global_fitness_randomized_short(
        kin_node_in,
        kout_g_out,
        tmin,
        vertices[itm.first]->instub_number,
        kinw_in);

      double rh2 = compute_global_fitness_randomized_short(
        kin_node_out,
        kout_g_in,
        tmout, vertices[itm.first]->outstub_number,
        kinw_out);

      double rh = std::min(rh1, rh2);

      if (rh < cmin)
      {
        cmin = rh;
        belongs_to = *its;
      }
    }

    if (belongs_to != -1)
    {
      if (to_check.find(belongs_to) == to_check.end())
      {
        std::deque<int> void_d;
        to_check[belongs_to] = void_d;
      }
      to_check[belongs_to].push_back(itm.first);
    }
    //if(paras.print_cbs)
      //cout<<"homeless node: "<<id_of(itm->first)<<" belongs_to "<<belongs_to<<" cmin... "<<cmin<< std::endl;
      //cherr();
  }
  //if(paras.print_cbs)
    //cout<<"homeless node: "<<homel.size()<<" try_to_assign: "<<homel_module.size()<<" modules to check: "<<to_check.size()<< std::endl;
  return before_procedure;
}

#endif // DIRECTED_OSLOMNET_EVALUATE_H
