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
  if (std::abs(a) > ee)
  {
    std::cerr << "the check failed because of " << a << std::endl;
    int e;
    std::cin >> e;
  }
}

template <typename uno, typename due>
void prints(std::pair <uno, due> &sq, std::ostream &out)
{
  out << sq.first << "\t" << sq.second << std::endl;
}

template <typename uno, typename due>
void prints(std::pair <uno, due> &sq)
{
  std::cout << sq.first << "\t" << sq.second << std::endl;
}

template <typename uno, typename due>
void prints(std::map <uno, due> &sq, std::ostream &out)
{
  typename std::map <uno, due>::iterator it = sq.begin();
  while (it != sq.end())
  {
    out << it->first << "\t" << it->second << std::endl;
    ++it;
  }

  out << std::endl;
}

template <typename uno, typename due>
void prints(std::multimap <uno, due> &sq, std::ostream &out)
{
  typename std::map <uno, due>::iterator it = sq.begin();
  while (it != sq.end())
  {
    out << it->first << "\t" << it->second << std::endl;
    ++it;
  }

  out << std::endl;
}

template <typename Seq>
void prints(Seq &sq, std::ostream &out)
{
  typename Seq::iterator it = sq.begin();
  while (it != sq.end())
  {
    out << *(it++) << "\t";
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

template <typename uno, typename due>
void prints(std::map <uno, due> &sq)
{
  for (typename std::map <uno, due>::iterator it = sq.begin(); it != sq.end(); ++it)
  {
    std::cout << it->first << "\t" << it->second << std::endl;
  }
  std::cout << std::endl;
}

template <typename uno, typename due>
void prints(std::multimap <uno, due> &sq)
{
  for (typename std::map<uno, due>::iterator it = sq.begin(); it != sq.end(); ++it)
  {
    std::cout << it->first << "\t" << it->second << std::endl;
  }
  std::cout << std::endl;
}

template <typename uno, typename due>
void prints(std::deque<uno> & a, std::deque<due> &b, std::ostream &out = std::cout) {
  for (int i = 0; i < a.size(); i++)
  {
    out << a[i] << " " << b[i] << std::endl;
  }
}

template <typename Seq>
void prints(Seq &sq)
{
  for (typename Seq::iterator it = sq.begin(); it != sq.end(); ++it)
  {
    std::cout << *it << "\t";
  }
  std::cout << std::endl;
}

template <typename T>
void prints(const std::deque<T> & sq) {
  for (int i = 0; i < sq.size(); i++)
    std::cout << sq[i] << "\t";
  std::cout << std::endl;
}

template <typename T>
void prints(const std::vector<T> & sq) {
  for (int i = 0; i < sq.size(); ++i)
  {
    std::cout << sq[i] << "\t";
  }
  std::cout << std::endl;
}

template <typename T>
void printm(std::deque<T> & M)
{
  for (int i = 0; i < M.size(); ++i)
    prints(M[i]);
  std::cout << std::endl;
}

template <typename T>
void printm(std::vector<T> & M)
{
  for (int i = 0; i < M.size(); ++i)
  {
    prints(M[i]);
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
  while (getline(lin, sas)) {
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
  while (getline(lin, sas)) {
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

void get_data_from_file_string(std::string s, std::deque<std::string> & a1, int col)
{
  // default will be col=1

  std::ifstream lin(s);

  a1.clear();
  col--;

  std::string sas;
  while (getline(lin, sas))
  {
    std::deque<std::string> v;
    separate_strings(sas, v);

    //prints(v);
    if (int(s.size()) > col)
    {
      a1.push_back(v[col]);
    }
  }
}

#endif
