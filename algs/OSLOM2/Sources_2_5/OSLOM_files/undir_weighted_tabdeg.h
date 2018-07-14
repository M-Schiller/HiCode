#pragma once

#ifndef UNDIR_WEIGHTED_TABDEG_H
#define UNDIR_WEIGHTED_TABDEG_H

#include "./standard_package/standard_include.cpp"

# define sqrt_two 1.41421356237
# define num_up_to 5
# define bisection_precision 1e-2

typedef std::multimap<double, std::pair<int, double>> cup_data_struct;

inline double compare_r_variables(double a, double b, double c, double d) {
  // r1 \in (a,b)
  // r2 \in (c,d)
  // compute the probability p(r2<r1)

  //cout<<a<<" "<<b<<" "<<c<<" "<<d<< std::endl;

  if (c < a)
    return (1 - compare_r_variables(c, d, a, b));

  if (c > b)
    return 0;

  if (d > b)
    return 0.5*(b - c)*(b - c) / ((d - c)*(b - a));
  return (b - 0.5*(d + c)) / (b - a);
}

inline double right_error_function(double x)
{
  return 0.5*erfc(x / sqrt_two);
}

inline double log_together(double minus_log_total, int number) {
  if (number < 11) {
    double fa = 1;
    double zum = 1;

    for (int i = 1; i < number; i++) {
      fa *= minus_log_total / i;
      zum += fa;
    }

    return std::max(zum*exp(-minus_log_total), 1e-100);
  }

  double mu = number;
  return std::max(right_error_function((minus_log_total - mu) / sqrt(mu)), 1e-100);
}

inline double fitted_exponent(int N) {
  double l = log(double(N));

  if (N > 100)
    return 4.2*l - 8.5;

  if (N > 30)
    return 3.5*l - 5.5;

  if (N > 7)
    return 2.5*l - 2;

  if (N > 1)
    return 1.3*l + 0.1;

  return 1;
}

inline double order_statistics_left_cumulative(int N, int pos, double x) {
  // the routine computes the probality c_pos=  p(X_pos <= x)
  // N is the total number of variables, pos is from 1 to N. N is the smallest.

  return LOG_TABLE.cum_binomial_right(N - pos + 1, N, x);
}

inline double inverse_order_statistics(int sample_dim, int pos, const double  & zerof, double lo, double hi)
{
  // anyway it's possible that I will tabulated this stuff for N small like sample_dim < 10000 which can already help a lot.
  // moreover, i should not use bisection but something better
  // finally there is the gaussian approx which could be studied more...
  //file: /Users/admin/Desktop/ambiente/right_binomial_cum/approx.cpp

  //double zerof= - log(1-threshold)/fitted_exponent(sample_dim);

  double mid = (hi + lo) / 2;

  while ((mid != lo) && (mid != hi)) {
    double fmid = order_statistics_left_cumulative(sample_dim, pos, mid);

    if (fabs(fmid - zerof) < bisection_precision*zerof)
      break;

    if ((fmid - zerof) <= 0)
      lo = mid;
    else
      hi = mid;

    mid = (hi + lo) / 2;
  }

  return mid;
}

inline double pron_min_exp(int N, double xi)
{
  // this should return the probability that the minimum of the quantiles p(c_min<=xi)
  //cout<<"-> "<<fitted_exponent(N)<< std::endl;

  return 1 - exp(-fitted_exponent(N)*xi);
}

inline double compute_probability_to_stop(const double & a, const double & b, const double & critical_xi, int Nstar, int pos)
{
  /*	this routine is to compute the bootstrap probaility that the node with extremes a and b will be below threshold
    when I already know that the average is below critical_xi */

  if (order_statistics_left_cumulative(Nstar, pos, b) <= critical_xi)
    return 1;

  return (inverse_order_statistics(Nstar, pos, critical_xi, (a + b)*0.5, b) - a) / (b - a);
}

inline bool equivalent_check(
  int pos_first,
  int pos_last,
  double & A_average,
  double & B_average,
  int equivalents,
  int Nstar,
  const double & critical_xi)
{
  // returns true is the test was passed

  /*small_simulation(pos_first, pos_last, A_average, B_average, equivalents, Nstar, critical_xi);
  cout<<pos_first<<" "<<pos_last<<" A, B "<<A_average<<" "<<B_average<<" "<<equivalents<<" "<<Nstar<<" "<<critical_xi<< std::endl;*/

  int pos = pos_first;
  double cr_previous = A_average;

  for (int i = equivalents; i >= 1; i--) {
    // loop which starts from the best node

    //cout<<"i... "<<i<< std::endl;

    if (order_statistics_left_cumulative(Nstar, pos, cr_previous) <= critical_xi) {
      if (order_statistics_left_cumulative(Nstar, pos, B_average) <= critical_xi)
        return true;

      double cr = inverse_order_statistics(Nstar, pos, critical_xi, cr_previous, B_average);
      //cout<<i<<" cr: "<<cr<<" "<<order_statistics_left_cumulative(equivalents, i, (cr-A_average)/(B_average-A_average))<< std::endl;
      //cout<<" this,  "<<order_statistics_left_cumulative(Nstar, pos, cr)<<" "<<Nstar<<" "<<pos<< std::endl;

      if (order_statistics_left_cumulative(equivalents, i, (cr - A_average) / (B_average - A_average)) > 0.5)
        return true;

      cr_previous = cr;
    }

    --pos;
  }

  return false;
}

inline bool equivalent_check_gather(
  cup_data_struct & a,
  int & until,
  const double & probability_a,
  const double & probability_b,
  int Nstar,
  const double & critical_xi)
{
  int nodes_added = -1;

  double A_average = 0;
  double B_average = 0;
  int pos_first = -1;
  int pos_last = -1;

  for (auto& itl : a)
  {
    if (nodes_added == until)
      break;

    ++nodes_added;

    if (compare_r_variables(
      probability_a,
      probability_b,
      itl.first - itl.second.second,
      itl.first + itl.second.second) > paras.equivalence_parameter)
    {
      if (pos_first == -1)
        pos_first = Nstar - nodes_added;

      A_average += itl.first - itl.second.second;
      B_average += itl.first + itl.second.second;
      pos_last = Nstar - nodes_added;
      //cout<<"pos_first: "<<pos_first<<" ... "<<pos_last<< std::endl;
    }
  }

  const int equivalents = pos_first - pos_last + 1;
  A_average /= equivalents;
  B_average /= equivalents;

  if (equivalents == 1)
    return true;

  if (equivalent_check(pos_first, pos_last, A_average, B_average, equivalents, Nstar, critical_xi))
  {
    //cout<<"check passed"<< std::endl;
    return true;
  }
  until = -1;
  //cout<<"check not passed"<< std::endl;
  return false;
}

inline double hyper_table(int kin_node, int kout_g, int tm, int degree_node)
{
  return LOG_TABLE.hyper(kin_node, kout_g, tm, degree_node);
}

inline double topological_05(int kin_node, int kout_g, int tm, int degree_node)
{
  return  LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node + 1) + ran4() * (hyper_table(kin_node, kout_g, tm, degree_node));
}

inline double compute_global_fitness(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total, int number_of_neighs, int Nstar, double & boot_interval) {
  /*  kin_node is referred to the node and not to the module */
  /* boot_interval is half the domain of the fitness. We assume that also in the weighted case the fitness can be linearized */
  /* which is true if the boot_interval is small enough	*/

  boot_interval = (0.5 + 1e-6*(ran4() - 0.5)) * hyper_table(kin_node, kout_g, tm, degree_node);
  double topologic = LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node + 1) + boot_interval;

  //cout<<"----------> "<<LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node+1)<<"  boot_interval "<<boot_interval<<" kin_node: "<<kin_node<<" / "<<degree_node<< std::endl;
  //cout<<"<> "<<LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node)<<"  boot_interval "<<boot_interval<<" kin_node: "<<kin_node<<" / "<<degree_node<< std::endl;

  if (!paras.weighted)
  {
    if (topologic > 1)
    {
      topologic = 1;
    }

    if (1. - topologic < boot_interval)
    {
      boot_interval = 1. - topologic;
    }
    if (topologic < boot_interval)
    {
      boot_interval = topologic;
    }

    return std::max(topologic, 1e-100);
  }

  topologic *= (Nstar + 1.) / (number_of_neighs + 1.);
  boot_interval *= (Nstar + 1.) / (number_of_neighs + 1.);

  if (topologic > 1) {
    topologic = 1;
    boot_interval = 1e-100;
  }

  //cout<<"to: "<<topologic<< std::endl;
  const double weight_part = log_together(minus_log_total, kin_node);
  //cout<<"we: "<<weight_part<<" - "<<minus_log_total<< std::endl;

  if (topologic <= 1e-100 || weight_part <= 1e-100) {
    boot_interval = 1e-100;
    return 1e-100;
  }

  const double _log_sum_ = -log(topologic*weight_part);
  boot_interval *= weight_part * _log_sum_;

  const double global_score = log_together(_log_sum_, 2);

  if (1. - global_score < boot_interval)
    boot_interval = 1. - global_score;
  if (global_score < boot_interval)
    boot_interval = global_score;

  return global_score;
}

inline double compute_global_fitness_step(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total, int number_of_neighs, int Nstar, double _step_)
{
  double topologic = LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node + 1) + _step_ * (hyper_table(kin_node, kout_g, tm, degree_node));

  if (!paras.weighted)
  {
    return std::max(topologic, 1e-100);
  }

  topologic *= (Nstar + 1.) / (number_of_neighs + 1.);

  if (topologic > 1)
    topologic = 1;

  const double weight_part = log_together(minus_log_total, kin_node);

  if (topologic <= 1e-100
    || weight_part <= 1e-100)
  {
    return 1e-100;
  }

  return log_together(-log(topologic * weight_part), 2);
}

inline double compute_global_fitness_ofive(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total, int number_of_neighs, int Nstar)
{
  return compute_global_fitness_step(kin_node, kout_g, tm, degree_node, minus_log_total, number_of_neighs, Nstar, 0.5);
}

inline double compute_global_fitness_randomized(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total, int number_of_neighs, int Nstar)
{
  return compute_global_fitness_step(kin_node, kout_g, tm, degree_node, minus_log_total, number_of_neighs, Nstar, ran4());
}

inline double compute_global_fitness_randomized_short(int kin_node, int kout_g, int tm, int degree_node, double minus_log_total)
{
  // this function is used in try_to_assign_homeless.
  // the usual problem is that we don't know the number of neighbors of the module.
  // this could be taken in account with some more thinking...

  double b2 = LOG_TABLE.right_cumulative_function(degree_node, kout_g, tm, kin_node + 1);

  const double topologic = b2 + 0.5 * (hyper_table(kin_node, kout_g, tm, degree_node));

  if (!paras.weighted)
    return std::max(topologic, 1e-100);

  const double weight_part = log_together(minus_log_total, kin_node);

  if (topologic <= 1e-100
    || weight_part <= 1e-100)
  {
    return 1e-100;
  }

  return std::min(topologic, weight_part);
}

class facts
{
public:

  facts(int a, double b, std::multimap<double, int>::iterator c, int d) { internal_degree = a; minus_log_total_wr = b; fitness_iterator = c; degree = d; };
  ~facts() = default;

  int degree;
  int internal_degree;
  double minus_log_total_wr;								// wr is the right part of the exponential for the weights, this is the sum over the internal stubs of that
  std::multimap<double, int>::iterator fitness_iterator;
};

class weighted_tabdeg
{
public:

  weighted_tabdeg() = default;
  ~weighted_tabdeg() = default;

  void _set_(weighted_tabdeg &);

  void clear();
  void edinsert(int a, int kp, int kt, double mtlw, double fit);
  bool erase(int a);
  void set_deque(std::deque<int> &);

  int size() const
  {
    return lab_facts.size();
  }

  bool empty() const
  {
    return lab_facts.empty();
  }

  void print_nodes(std::ostream &, std::deque<int> &);

  void set_and_update_group(int nstar, int nn, int kout_g, int tm, weighted_tabdeg & one);
  void set_and_update_neighs(int nstar, int nn, int kout_g, int tm, weighted_tabdeg & one);
  bool update_group(int a, int delta_degree, double delta_mtlw, int nstar, int nn, int kout_g, int tm, int kt, std::deque<int> & to_be_erased);
  bool update_neighs(int a, int delta_degree, double delta_mtlw, int nstar, int kout_g, int tm, int kt);

  int best_node(int & lab, double & best_fitness, int kout_g, int Nstar, int nneighs, int tm);
  int worst_node(int & lab, double & worst_fitness, int kout_g, int Nstar, int nneighs, int tm);
  bool is_internal(int a);

  std::map<int, facts> lab_facts;					// maps the label into the facts
  std::multimap<double, int> fitness_lab;			// maps the fitness into the label  (this can be optimized)
};

inline void weighted_tabdeg::clear()
{
  lab_facts.clear();
  fitness_lab.clear();
}

inline void weighted_tabdeg::edinsert(int a, int kp, int kt, double mtlw, double fit)		// this function inserts element a (or edit it if it was already inserted)
{
  erase(a);

  const auto fiit = fitness_lab.emplace(fit, a);
  facts F(kp, mtlw, fiit, kt);

  lab_facts.emplace(a, F);
}

inline bool weighted_tabdeg::erase(int a)		// this function erases element a if exists (and returns true)
{
  const auto itm = lab_facts.find(a);
  if (itm != lab_facts.end())
  {
    fitness_lab.erase(itm->second.fitness_iterator);
    lab_facts.erase(itm);

    return true;
  }

  return false;
}

inline bool weighted_tabdeg::is_internal(int a)
{
  return lab_facts.find(a) != lab_facts.end();
}

inline void weighted_tabdeg::set_deque(std::deque<int> & vv)
{
  vv.clear();

  for (auto& lab_fact : lab_facts)
    vv.push_back(lab_fact.first);
}

inline void weighted_tabdeg::print_nodes(std::ostream & outb, std::deque<int> & lab_id)
{
  std::cout << "printing nodes:.. (lab intk mtlw fitness degree) " << size() << std::endl;

  for (auto &[label, fact] : lab_facts)
    std::cout << lab_id[label] << " "
    << fact.internal_degree << " "
    << fact.minus_log_total_wr << " "
    << (fact.fitness_iterator)->first << " "
    << fact.degree << std::endl;
}

inline int weighted_tabdeg::worst_node(int & lab, double & worst_fitness, int kout_g, int Nstar, int nneighs, int tm)
{
  //cout<<"worst_node fitness - lab - (cgroup)"<< std::endl;
  //prints(fitness_lab);

  lab = -1;
  worst_fitness = -1;

  if (fitness_lab.begin() == fitness_lab.end())
    return -1;

  int stopper = 0;

  for (auto bit = fitness_lab.rbegin(); bit != fitness_lab.rend(); ++bit)
  {
    const std::map<int, facts> ::iterator itm = lab_facts.find(bit->second);

    const double F = compute_global_fitness_randomized(
      itm->second.internal_degree,
      kout_g + 2 * itm->second.internal_degree - itm->second.degree,
      tm + itm->second.degree,
      itm->second.degree,
      itm->second.minus_log_total_wr,
      nneighs + 1,
      Nstar + 1);

    if (F > worst_fitness)
    {
      worst_fitness = F;
      lab = itm->first;
    }

    stopper++;
    if (stopper == num_up_to)
      break;
  }

  return 0;
}

inline int weighted_tabdeg::best_node(int & lab, double & best_fitness, int kout_g, int Nstar, int nneighs, int tm)
{
  // I can try to compute the fitness here

  /*cout<<"NE BEST NODE "<< std::endl;

  cout<<"fitness_lab  "<< std::endl;
  prints(fitness_lab);*/

  lab = -1;
  best_fitness = 1;

  if (fitness_lab.begin() == fitness_lab.end())
  {
    return -1;
  }

  int stopper = 0;
  for (auto &[fitness, label] : fitness_lab)
  {
    const std::map<int, facts> ::iterator itm = lab_facts.find(label);

    const double F = compute_global_fitness_randomized(itm->second.internal_degree, kout_g, tm, itm->second.degree, itm->second.minus_log_total_wr, nneighs, Nstar);

    //cout<<itm->first<<" "<<F<<" ... node-fit"<< std::endl;

    if (F < best_fitness)
    {
      best_fitness = F;
      lab = itm->first;
    }

    stopper++;
    if (stopper == num_up_to)
      break;
  }

  return 0;
}

inline void weighted_tabdeg::_set_(weighted_tabdeg & one)
{
  clear();
  for (auto &[label, fact] : one.lab_facts)
  {
    edinsert(label, fact.internal_degree, fact.degree, fact.minus_log_total_wr, (fact.fitness_iterator)->first);
  }
}

inline bool weighted_tabdeg::update_group(int a, int delta_degree, double delta_mtlw, int nstar, int nn, int kout_g, int tm, int kt, std::deque<int> & to_be_erased)
{
  // this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)

  auto itm = lab_facts.find(a);
  if (itm == lab_facts.end())
    return false;

  itm->second.minus_log_total_wr += delta_mtlw;
  itm->second.internal_degree += delta_degree;

  if (itm->second.internal_degree == 0 && size() > 1)
  {
    to_be_erased.push_back(a);
    return true;
  }

  //cout<<"UPdating... group "<<a<< std::endl;

  double fit = compute_global_fitness_ofive(
    itm->second.internal_degree,
    kout_g + 2 * itm->second.internal_degree - itm->second.degree,
    tm + itm->second.degree,
    itm->second.degree,
    itm->second.minus_log_total_wr,
    nn + 1,
    nstar + 1);

  fitness_lab.erase(itm->second.fitness_iterator);
  const auto fiit = fitness_lab.emplace(fit, a);
  itm->second.fitness_iterator = fiit;

  return true;
}

inline bool weighted_tabdeg::update_neighs(int a, int delta_degree, double delta_mtlw, int nstar, int kout_g, int tm, int kt)
{
  // this function is to change the internal degree and mtlw of a certain node (to insert it or erase if necessary)

  //cout<<"UPdating... neighs "<<a<<" "<<kt<< std::endl;

  auto itm = lab_facts.find(a);
  if (itm == lab_facts.end())
  {
    edinsert(a, 0, kt, 0, 1);
    itm = lab_facts.find(a);
  }

  itm->second.internal_degree += delta_degree;
  if (itm->second.internal_degree == 0) {
    //cout<<"erased from neigh update "<<a<< std::endl;
    erase(a);

    return true;
  }

  itm->second.minus_log_total_wr += delta_mtlw;

  double fit = compute_global_fitness_ofive(itm->second.internal_degree, kout_g, tm, itm->second.degree, itm->second.minus_log_total_wr, size(), nstar);

  fitness_lab.erase(itm->second.fitness_iterator);
  itm->second.fitness_iterator = fitness_lab.emplace(fit, a);

  return true;
}

inline void weighted_tabdeg::set_and_update_group(int nstar, int nn, int kout_g, int tm, weighted_tabdeg & one)
{
  /*this function is to set and update the fitnesses of all the nodes in cgroup*/

  clear();
  for (auto &[label, fact] : one.lab_facts)
  {
    double fit =
      compute_global_fitness_ofive(
        fact.internal_degree,
        kout_g + 2 * fact.internal_degree - fact.degree,
        tm + fact.degree,
        fact.degree,
        fact.minus_log_total_wr,
        nn + 1,
        nstar + 1);

    edinsert(label, fact.internal_degree, fact.degree, fact.minus_log_total_wr, fit);
  }
}

inline void weighted_tabdeg::set_and_update_neighs(int nstar, int nn, int kout_g, int tm, weighted_tabdeg & one)
{
  /*this function is to set and update the fitnesses of all the nodes in neighs*/

  clear();
  for (auto &[label, fact] : one.lab_facts)
  {
    const double fit = compute_global_fitness_ofive(fact.internal_degree, kout_g, tm, fact.degree, fact.minus_log_total_wr, nn, nstar);
    edinsert(label, fact.internal_degree, fact.degree, fact.minus_log_total_wr, fit);
  }
}
#endif // UNDIR_WEIGHTED_TABDEG_H
