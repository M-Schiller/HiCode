#pragma once

#ifndef PRINT_INCLUDED
#define PRINT_INCLUDED

#include <iostream>
#include <cmath>

void cherr()
{
  std::cerr << "the check failed" << std::endl;
  int e;
  std::cin >> e;
}

void cherr(double a)
{
  std::cerr << "the check failed because of " << a << std::endl;
  int e;
  std::cin >> e;
}

void cherr(double a, double ee)
{
  if (abs(a) > ee)
  {
    std::cerr << "the check failed because of " << a << std::endl;
    int e;
    std::cin >> e;
  }
}

template <typename T, typename U>
void prints(std::pair <T, U> &pair, std::ostream &out)
{
  out << pair.first << "\t" << pair.second << std::endl;
}

template <typename T, typename U>
void prints(std::pair <T, U> &sq)
{
  std::cout << sq.first << "\t" << sq.second << std::endl;
}

template <typename T, typename U>
void prints(std::map <T, U> &sq, std::ostream &out)
{
  for (auto it = sq.begin(); it != sq.end(); ++it)
  {
    out << it->first << "\t" << it->second << std::endl;
  }
  out << std::endl;
}

template <typename T, typename U>
void prints(std::multimap <T, U> &sq, std::ostream &out)
{
  for (typename std::multimap<T, U>::iterator it = sq.begin(); it != sq.end(); ++it)
  {
    out << it->first << "\t" << it->second << std::endl;
  }
  out << std::endl;
}

template <typename Seq>
void prints(Seq &sq, std::ostream &out)
{
  for (auto it = sq.begin(); it != sq.end(); ++it)
  {
    out << *it << "\t";
  }
  out << std::endl;
}

template <typename T>
void prints(T *a, int b)
{
  for (int i = 0; i < b; i++)
  {
    std::cout << a[i] << " ";
  }
  std::cout << std::endl;
}

template<typename T, template<typename> class C>
void printm(C<T>& c, std::ostream &out)
{
  for (typename C<T>::iterator it = c.begin(); it != c.end(); ++it)
  {
    prints(*it, out);
  }
  out << std::endl;
}

template <typename T, typename U>
void prints(std::map <T, U> &sq)
{
  for (typename std::map <T, U>::iterator it = sq.begin(); it != sq.end(); ++it)
  {
    std::cout << it->first << "\t" << it->second << std::endl;
  }
  std::cout << std::endl;
}

template <typename T, typename U>
void prints(std::multimap <T, U> &sq)
{
  for (typename std::map <T, U>::iterator it = sq.begin(); it != sq.end(); ++it)
  {
    std::cout << it->first << "\t" << it->second << std::endl;
  }
  std::cout << std::endl;
}

template <typename T, typename U>
void prints(std::deque<T> & a, std::deque<U> &b)
{
  for (int i = 0; i < a.size(); i++)
  {
    std::cout << a[i] << " " << b[i] << std::endl;
  }
}

template <typename T, typename U>
void prints(std::deque<T> & a, std::deque<U> &b, std::ostream &out)
{
  for (int i = 0; i < a.size(); i++)
  {
    out << a[i] << " " << b[i] << std::endl;
  }
}

template <typename Seq>
void prints(Seq &sq) {
  for (typename Seq::iterator it = sq.begin(); it != sq.end(); ++it)
  {
    std::cout << *it << "\t";
  }
  std::cout << std::endl;
}

template <typename T>
void prints(const std::deque<T> & sq)
{
  for (auto& elem : sq)
  {
    std::cout << elem << "\t";
  }
  std::cout << std::endl;
}

template <typename T>
void prints(const std::vector<T> & sq)
{
  for (auto& elem : sq)
  {
    std::cout << elem << "\t";
  }
  std::cout << std::endl;
}

template <typename T>
void printm(std::deque<T> & M)
{
  for (auto& elem : M)
  {
    prints(elem);
  }
  std::cout << std::endl;
}

template <typename T>
void printm(std::vector<T> & M)
{
  for (auto& elem : M)
  {
    prints(elem);
  }
  std::cout << std::endl;
}

template <typename T>
void get_data_from_file(std::string s, std::deque<T> & a1, int col)
{
  // default will be col=1

  std::ifstream lin(s);

  a1.clear();
  col--;

  std::string sas;
  while (getline(lin, sas))
  {
    std::deque<double> s;
    cast_string_to_doubles(sas, s);
    if (s.size() > col) {
      a1.push_back(s[col]);
    }
  }
}

template <typename T>
void get_data_from_file(std::string s, std::deque<T> & a1)
{
  get_data_from_file(s, a1, 1);
}

template <typename T>
void get_data_from_file(std::string s, std::deque<T> & a1, std::deque<T> & a2, int col1, int col2)
{
  // default will be col1=1, col2=2

  std::ifstream lin(s);

  a1.clear();
  a2.clear();
  col1--;
  col2--;

  std::string sas;
  while (getline(lin, sas))
  {
    std::deque<double> s;
    cast_string_to_doubles(sas, s);
    if (s.size() > col2)
    {
      a1.push_back(s[col1]);
      a2.push_back(s[col2]);
    }
  }
}

template <typename T>
void get_data_from_file(std::string s, std::deque<T> & a1, std::deque<T> & a2)
{
  get_data_from_file(s, a1, a2, 1, 2);
}

#endif
