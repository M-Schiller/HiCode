#pragma once

#ifndef CAST_INCLUDED
#define CAST_INCLUDED

#include <string>
#include <deque>
#include <cmath>

bool cast_string_to_double(std::string &str, double &h)
{
  // set h= the number written in b[];
  // return false if there is an error

  h = 0;
  int num_exp_markers = 0;

  for (char c : str)
  {
    if (c == 'e' || c == 'E')
    {
      num_exp_markers++;
    }
  }

  if (num_exp_markers > 1)
  {
    return false;
  }

  if (num_exp_markers == 1)
  {
    std::string sbefe;
    std::string safte;

    num_exp_markers = 0;

    for (char c : str)
    {
      if (num_exp_markers == 0
        && (c != 'e' && c != 'E'))
      {
        sbefe.push_back(c);
      }

      if (num_exp_markers == 1 && (c != 'e' && c != 'E'))
      {
        safte.push_back(c);
      }

      if (c == 'e' || c == 'E')
      {
        num_exp_markers++;
      }
    }

    double number;
    bool bone = cast_string_to_double(sbefe, number);

    double exp;
    bool btwo = cast_string_to_double(safte, exp);

    if (bone && btwo) {
      h = number * pow(10, exp);
      return true;
    }
    return false;
  }

  if (str.empty())
    return false;

  int sign = 1;

  if (str[0] == '-')
  {
    str[0] = '0';
    sign = -1;
  }

  if (str[0] == '+')
    str[0] = '0';

  int digits_before = 0;
  for (int i = 0; i < str.size(); i++)
    if (str[i] != '.')
      digits_before++;
    else
      break;

  int j = 0;

  while (j != digits_before) {
    int number = (int(str[j]) - 48);
    h += number * pow(10, digits_before - j - 1);

    if (number < 0 || number>9) {
      //cout<<"err "<<number<< std::endl;
      return false;
    }

    j++;
  }

  j = digits_before + 1;

  while (j < str.size()) {
    int number = (int(str[j]) - 48);
    h += number * pow(10, digits_before - j);

    if (number < 0 || number>9)
      return false;

    j++;
  }
  h = sign * h;
  return true;
}

double cast_string_to_double(std::string &b)
{
  return std::stod(b);
}

inline int cast_int(double u)
{
  return lround(u);
}

int cast_string_to_char(std::string file_name, char *b)
{
  for (int i = 0; i < file_name.size(); i++)
    b[i] = file_name[i];
  b[file_name.size()] = '\0';

  return 0;
}

bool cast_string_to_doubles(std::string &b, std::deque<double> & v)
{
  v.clear();
  std::string s1;

  for (unsigned int i = 0; i < b.size(); ++i)
  {
    if (!b.empty()
      && b[0] != '#')
    {
      if (b[i] != '1'
        && b[i] != '2'
        && b[i] != '3'
        && b[i] != '4'
        && b[i] != '5'
        && b[i] != '6'
        && b[i] != '7'
        && b[i] != '8'
        && b[i] != '9'
        && b[i] != '0'
        && b[i] != '.'
        && b[i] != 'e'
        && b[i] != '+'
        && b[i] != '-'
        && b[i] != 'E')
      {
        double num;
        if (cast_string_to_double(s1, num))
          v.push_back(num);

        s1.clear();
      }
      else
        s1.push_back(b[i]);

      if (i == b.size() - 1)
      {
        double num;
        if (cast_string_to_double(s1, num))
          v.push_back(num);

        s1.clear();
      }
    }
    return true;
  }
  return false;
}

bool cast_string_to_doubles(std::string &b, std::deque<int> & v)
{
  v.clear();
  std::deque<double> d;
  cast_string_to_doubles(b, d);
  for (double i : d)
  {
    v.push_back(cast_int(i));
  }
  return true;
}

bool separate_strings(std::string &b, std::deque<std::string> & v)
{
  v.clear();
  std::string s1;

  for (unsigned i = 0; i < b.size(); i++)
  {
    if (b[i] == ' '
      || b[i] == '\t'
      || b[i] == '\n'
      || b[i] == ',')
    {
      if (!s1.empty())
        v.push_back(s1);

      s1.clear();
    }
    else
      s1.push_back(b[i]);

    if (i == b.size() - 1)
    {
      if (!s1.empty())
        v.push_back(s1);
      s1.clear();
    }
  }
  return true;
}

double approx(double a, int digits)
{
  digits--;
  bool neg = false;
  if (a < 0)
  {
    neg = true;
    a = -a;
  }

  //cout<<a<< std::endl;

  int tpow = 0;
  while (a < pow(10, digits))
  {
    tpow++;
    a *= 10;
  }
  while (a > pow(10, digits + 1))
  {
    tpow--;
    a /= 10;
  }

  if (!neg)
    return cast_int(a) / pow(10, tpow);
  return -cast_int(a) / pow(10, tpow);
}

#endif
