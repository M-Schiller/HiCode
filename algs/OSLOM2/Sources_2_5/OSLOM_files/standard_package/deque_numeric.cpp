#pragma once

#ifndef DEQUE_NUMERIC
#define DEQUE_NUMERIC

#include <deque>
#include <set>
#include <numeric>

bool compare(std::deque<double> & one, std::deque<double> & two)
{
  if (one.size() != two.size())
  {
    return false;
  }

  for (int i = 0; i<int(one.size()); i++)
  {
    if (abs(one[i] - two[i]) > 1e-7)
    {
      return false;
    }
  }
  return true;
}

double Euclidean_norm(const std::deque<double> & a)
{
  return sqrt(std::inner_product(a.begin(), a.end(), a.begin(), 0.));
}

int Euclidean_normalize(std::deque<double> & a)
{
  double norm = Euclidean_norm(a);
  for (double& i : a)
  {
    i /= norm;
  }
  return 0;
}

double scalar_product(std::deque<double> & a, std::deque<double> & b)
{
  return std::inner_product(a.begin(), a.end(), b.begin(), 0.);
}

int orthogonalize(std::deque<double> & a, std::deque<std::deque<double>> & M)
{
  Euclidean_normalize(a);

  for (auto& i : M)
  {
    double w = scalar_product(a, i);

    for (int j = 0; j<int(a.size()); j++)
    {
      a[j] -= w * i[j];
    }
  }
  Euclidean_normalize(a);
  return 0;
}

int matrix_time_vector(std::deque<std::deque<double>> & Q, std::deque<double> & v, std::deque<double> & new_s)
{
  new_s.clear();

  for (auto& i : Q)
  {
    double n = 0;
    for (int j = 0; j<int(i.size()); j++)
    {
      n += i[j] * v[j];
    }

    new_s.push_back(n);
  }
  return 0;
}

template<typename T>
void set_to_deque(const std::set<T> & s, std::deque<T> & a)
{
  a.clear();
  a.insert(a.begin(), s.begin(), s.end());
}

template<typename T>
void deque_to_set(const std::deque<T> & a, std::set<T> & s)
{
  s.clear();
  s.insert(a.begin(), a.end());
}

template<typename T>
void deque_to_set_app(const std::deque<T> & a, std::set<T> & s)
{
  s.insert(a.begin(), a.end());
}

template<typename T>
T norm_one(const std::deque<T> & a)
{
  return std::accumulate(a.begin(), a.end(), T());
}

int normalize_one(std::deque<double> & a)
{
  double norm = norm_one(a);
  for (double& i : a)
  {
    i /= norm;
  }

  return 0;
}

double jaccard(std::set<int> & a1, std::set<int> & a2)
{
  std::deque<int> group_intsec;
  std::deque<int> group_union;
  set_intersection(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group_intsec));
  set_union(a1.begin(), a1.end(), a2.begin(), a2.end(), back_inserter(group_union));

  return double(group_intsec.size()) / group_union.size();
}

#endif
