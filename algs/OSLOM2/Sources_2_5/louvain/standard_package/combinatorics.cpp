#pragma once

#ifndef COMBINATORICS_INCLUDED
#define COMBINATORICS_INCLUDED

#include <deque>
#include <map>
#include <set>
#include <iostream>
#include <algorithm>

template <typename Seq>
double average_func(Seq &sq)
{
  if (sq.empty())
    return 0;

  double av = 0;
  auto it = sq.begin();
  while (it != sq.end())
    av += *(it++);

  av = av / sq.size();

  return av;
}

template <typename Seq>
double variance_func(Seq &sq)
{
  if (sq.empty())
    return 0;

  double av = 0;
  double var = 0;

  typename Seq::iterator it = sq.begin();
  while (it != sq.end()) {
    av += *(it);
    var += (*(it))*(*(it));
    ++it;
  }

  av = av / sq.size();
  var = var / sq.size();
  var -= av * av;

  if (var < 1e-7)
    return 0;

  return var;
}

// this returns the average of the discrete probability function stored in Seq
template <typename Seq>
double average_pf(Seq &sq)
{
  double av = 0;
  int h = 0;

  typename Seq::iterator it = sq.begin();
  while (it != sq.end()) {
    av += *(it)*h;
    ++it;
    h++;
  }
  return av;
}

template <typename Seq>
double variance_pf(Seq &sq)
{
  double av = 0;
  double var = 0;
  int h = 0;

  typename Seq::iterator it = sq.begin();
  while (it != sq.end())
  {
    av += *(it)* h;
    var += (*(it)) * h * h;
    ++it;
    h++;
  }

  var -= av * av;

  if (var < 1e-7)
    return 0;

  return var;
}

double log_factorial(int num)
{
  double log_result = 0;
  for (int i = 1; i <= num; i++)
    log_result += log(i);

  return (log_result);
}

double log_combination(int n, int k)
{
  if (k == 0)
    return 0;

  if (n < k)
    return 0;

  if (n - k < k)
    k = n - k;

  double log_c = 0;
  for (int i = n - k + 1; i <= n; i++)
    log_c += log(i);

  for (int i = 1; i <= k; i++)
    log_c -= log(i);

  return log_c;
}

double binomial(int n, int x, double p)		//	returns the binomial distribution, n trials, x successes, p probability{
{
  if (p == 0)
  {
    if (x == 0)
    {
      return 1;
    }
    return 0;
  }

  if (p >= 1)
  {
    if (x == n)
    {
      return 1;
    }
    return 0;
  }

  double log_b = 0;
  log_b += log_combination(n, x) + x * log(p) + (n - x)*log(1 - p);
  return (exp(log_b));
}

//to draw a number:

/*

deque <double> cumulative;
binomial_cumulative(10, 0.5, cumulative);
int nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin();

*/

int binomial_cumulative(int n, double p, std::deque<double> &cum)
{
  cum.clear();

  double c = 0;
  for (int i = 0; i <= n; i++)
  {
    c += binomial(n, i, p);
    cum.push_back(c);
  }
  return 0;
}

// this function sets "cumulative" as the cumulative function of (1/x)^tau, with range= [min, n]
//to draw a number:
//int nn=lower_bound(cumulative.begin(), cumulative.end(), ran4())-cumulative.begin()+min_degree;

int powerlaw(int n, int min, double tau, std::deque<double> &cumulative)
{
  cumulative.clear();
  double a = 0;

  for (double h = min; h < n + 1; h++)
    a += pow((1. / h), tau);

  double pf = 0;
  for (double i = min; i < n + 1; i++)
  {
    pf += 1 / a * pow((1. / (i)), tau);
    cumulative.push_back(pf);
  }
  return 0;
}

int distribution_from_cumulative(const std::deque<double> &cum, std::deque<double> &distr)		// cum is the cumulative, distr is set equal to the distribution
{
  distr.clear();
  double previous = 0;
  for (double i : cum)
  {
    distr.push_back(i - previous);
    previous = i;
  }
  return 0;
}

int cumulative_from_distribution(std::deque<double> &cum, const std::deque<double> &distr)		// cum is set equal to the cumulative, distr is the distribution
{
  cum.clear();
  double sum = 0;
  for (int i = 0; i < distr.size(); i++) {
    sum += distr[i];
    cum.push_back(sum);
  }

  return 0;
}

double poisson(int x, double mu)
{
  return (exp(-mu + x * log(mu) - log_factorial(x)));
}

int shuffle_and_set(int *due, int dim)		// it sets due as a random sequence of integers from 0 to dim-1
{
  std::multimap <double, int> uno;
  for (int i = 0; i < dim; i++)
  {
    uno.emplace(ran4(), i);
  }

  int h = 0;
  for (auto& it : uno)
  {
    due[h++] = it.second;
  }

  return 0;
}

template<typename T>
int shuffle_s(std::deque<T> &sq)
{
  int siz = sq.size();
  if (siz == 0)
    return -1;

  for (int i = 0; i < sq.size(); i++)
  {
    int random_pos = irand(siz - 1);

    T random_card_ = sq[random_pos];

    sq[random_pos] = sq[siz - 1];
    sq[siz - 1] = random_card_;
    siz--;
  }
  return 0;
}

template <typename T>
int shuffle_s(T *a, int b)
{
  int siz = b;
  if (siz == 0)
    return -1;

  for (int i = 0; i < b; i++)
  {
    int random_pos = irand(siz - 1);

    T random_card_ = a[random_pos];

    a[random_pos] = a[siz - 1];
    a[siz - 1] = random_card_;
    siz--;
  }
  return 0;
}

double compute_r(int x, int k, int kout, int m)
{
  double r = 0;

  for (int i = x; i <= k; i++)
    r += binomial(k, i, double(kout) / double(m));

  return r;
}

int add_factors(std::deque<double> & num, std::deque<double> &den, int  n, int k)
{
  if (n < k)
    return -1;

  if (n - k < k)
    k = n - k;

  if (k == 0)
    return 0;

  for (int i = n - k + 1; i <= n; i++)
    num.push_back(double(i));

  for (int i = 1; i <= k; i++)
    den.push_back(double(i));

  return 0;
}

double compute_hypergeometric(int i, int k, int kout, int m) {
  if (i > k || i > kout || k > m || kout > m)
    return 0;

  double prod = 1;
  std::deque <double> num;
  std::deque <double> den;

  if (add_factors(num, den, kout, i) == -1)
    return 0;

  if (add_factors(num, den, m - kout, k - i) == -1)
    return 0;

  if (add_factors(den, num, m, k) == -1)
    return 0;

  std::sort(num.begin(), num.end());
  std::sort(den.begin(), den.end());

  //prints(den);

  for (double h : den)
  {
    if (h <= 0)
    {
      std::cerr << "denominator has zero or less (in the hypergeometric)" << std::endl;
      return 0;
    }
  }

  for (double h : num)
  {
    if (num[h] <= 0)
    {
      std::cerr << "numerator has zero or less (in the hypergeometric)" << std::endl;
      return 0;
    }
  }

  //cout<<"sizes: "<<num.size()<<" "<<den.size()<< std::endl;
  for (int i = 0; i < num.size(); i++)
  {
    prod = prod * num[i] / den[i];
  }
  return prod;
}

/*
double compute_self_links(int k, int n, int x) {
  if (2*x > k)
    return 0;

  double prod= log_combination(n/2, k-x) + log_combination(k-x, x) + (k-2*x) * log(2) - log_combination(n, k);

  return exp(prod);
}

//*/

int random_from_set(std::set<int> & s)
{
  auto it1 = s.begin();
  std::advance(it1, irand(s.size() - 1));
  return *it1;
}

#endif
