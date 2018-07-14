#pragma once

#ifndef MUTUAL_INCLUDED
#define MUTUAL_INCLUDED

int overlap_grouping(std::deque<std::deque<int>> ten, int unique)		//hrepiguhpueh
{
  std::set <int> conta;
  int all = 0;

  for (auto& i : ten)
  {
    for (int j : i)
    {
      conta.insert(j);
    }

    all += i.size();
  }

  unique = conta.size();

  const int overlap = all - unique;

  return overlap;
}

double mutual(std::deque<std::deque<int>> en, std::deque<std::deque<int>> ten)
{
  // en e ten are two partitions of integer numbers
  int dim;

  {
    std::set <int> conta;
    //set <int> ten_;
    //set <int> en_;

    for (auto& i : ten)
    {
      std::sort(i.begin(), i.end());
      for (int j : i)
      {
        conta.insert(j);
        //ten_.insert(ten[i][j]);
      }
    }

    for (auto& i : en)
    {
      std::sort(i.begin(), i.end());
      for (int j : i)
      {
        conta.insert(j);
        //en_.insert(en[i][j]);
      }
    }

    dim = conta.size();

    /*
    for (set<int>::iterator its=conta.begin(); its!=conta.end(); its++) {
      if(ten_.find(*its)==ten_.end()) {
        deque <int> first;
        first.push_back(*its);
        ten.push_back(first);
      }

      if(en_.find(*its)==en_.end()) {
        deque <int> first;
        first.push_back(*its);
        en.push_back(first);
      }
    }
    */
  }

  //cout<<"dim:\t"<<dim<<endl;

  std::deque<std::deque<double>> N;
  std::deque<double> first(en.size());

  for (int i = 0; i < ten.size(); i++)
  {
    N.emplace_back(first);
  }

  std::deque <int> s(dim);
  for (int i = 0; i < ten.size(); i++)
  {
    for (int j = 0; j < en.size(); j++)
    {
      N[i][j] = set_intersection(ten[i].begin(), ten[i].end(), en[j].begin(), en[j].end(), s.begin()) - s.begin();
    }
  }
  //printv(N);

  /*
  cout<<"one:"<<endl;
  printm(ten);

  cout<<"two:"<<endl;
  printm(en);

  cout<<"confusion matrix"<<endl;
  printm(N, cout);
  */

  //cout<<"confusion matrix"<<endl;
  //printm(N, cout);

  std::deque <double> NR(ten.size());
  std::deque <double> NC(en.size());
  double NTOT = dim;

  for (int i = 0; i < ten.size(); i++)
  {
    for (int j = 0; j < en.size(); j++)
    {
      NR[i] += N[i][j];
      NC[j] += N[i][j];
    }
  }

  double IN = 0;
  double ID1 = 0;
  double ID2 = 0;

  for (int i = 0; i < ten.size(); i++)
  {
    for (int j = 0; j < en.size(); j++)
    {
      if (N[i][j] != 0)
      {
        IN += N[i][j] * log(N[i][j] * NTOT / (NR[i] * NC[j]));
      }
    }
  }

  IN = -2.*IN;

  for (int i = 0; i < ten.size(); i++)
  {
    if (NR[i] != 0)
    {
      ID1 += NR[i] * log(NR[i] / (NTOT));
    }
  }

  for (unsigned j = 0; j < en.size(); j++)
  {
    if (NC[j] != 0)
    {
      ID2 += NC[j] * log(NC[j] / NTOT);
    }
  }

  double I = IN / (ID1 + ID2);

  if ((ID1 + ID2) == 0)
  {
    I = -2;
  }

  return I;
}

double H(double a)
{
  if (a <= 0)
    return 0;
  return (-a * log(a));
}

double H(std::deque <double> &p)
{
  double h = 0;
  for (auto it = p.begin(); it != p.end(); ++it)
  {
    if (*it != 0)
    {
      h += (*it)*log(*it);
    }
  }

  return (-h);
}

double H_x_given_y(std::deque<std::deque<int>> &en, std::deque<std::deque<int>> &ten, int dim)
{
  // you know y and you want to find x according to a certain index labelling.
  // so, for each x you look for the best y.

  double H_x_y = 0;
  double H2 = 0;

  for (auto& j : en)
  {
    std::deque <double> p;
    auto I2 = double(j.size());
    double O2 = (dim - I2);
    p.push_back(I2 / dim);
    p.push_back(O2 / dim);
    double H2_ = H(p);
    p.clear();

    H2 += H2_;

    double diff = H2_;

    for (int i = 0; i < ten.size(); i++)
    {
      auto I1 = double(ten[i].size());
      double O1 = (dim - I1);

      //cout<<"I1 "<<I1<<" O1 "<<O1<<endl;

      p.push_back(I1 / dim);
      p.push_back(O1 / dim);
      double H1_ = H(p);
      p.clear();

      //prints(ten[i]);
      //cout<<"H1_: "<<H1_<<"\t";

      std::deque <int> s(dim);
      double I1_I2 = set_intersection(ten[i].begin(), ten[i].end(), j.begin(), j.end(), s.begin()) - s.begin();	// common
      double I1_02 = set_difference(ten[i].begin(), ten[i].end(), j.begin(), j.end(), s.begin()) - s.begin();
      double O1_I2 = set_difference(j.begin(), j.end(), ten[i].begin(), ten[i].end(), s.begin()) - s.begin();
      double O1_02 = dim - I1_I2 - I1_02 - O1_I2;

      p.push_back(I1_I2 / dim);
      p.push_back(O1_02 / dim);

      double H12_positive = H(p);

      p.clear();
      p.push_back(I1_02 / dim);
      p.push_back(O1_I2 / dim);

      double H12_negative = H(p);

      double H12_ = H12_negative + H12_positive;

      p.clear();

      if (H12_negative > H12_positive)
      {
        H12_ = H1_ + H2_;
        //cout<<"the negative part is bigger"<<endl;
        //prints(en[j]);
        //prints(ten[i]);
      }

      /*

      cout<<"worst case "<<H1_+H2_<<"\ttotal "<<H12_negative+H12_positive<<"\tnegative part "<<H12_negative<<"\tpositive part "<<H12_positive<<endl;
      prints(en[j]);
      prints(ten[i]);

      */

      //cout<<"entropies "<<H1_<<" "<<H2_<<" "<<H12_<<endl;

      if ((H12_ - H1_) < diff)
      {
        diff = (H12_ - H1_);
      }
    }

    //H_x_y+=diff;
    if (H2_ == 0)
    {
      H_x_y += 1;
    }
    else
    {
      H_x_y += (diff / H2_);
    }
  }

  //if (H2==0)
  //	return 1;

  return (H_x_y / (en.size()));
}

double mutual2(std::deque<std::deque<int>> en, std::deque<std::deque<int>> ten)
{
  if (en.empty()
    || ten.empty())
  {
    return 0;
  }

  // en e ten are two partitions of integer numbers
  int dim;

  //printm(ten);
  //printm(en);

  {
    std::map <int, int> all;		// node, index
    //set <int> ten_;
    //set <int> en_;

    for (auto& i : ten)
    {
      for (int& j : i)
      {
        all.emplace(j, all.size());
        //ten_.insert(ten[i][j]);
      }
    }

    for (auto& i : en)
    {
      for (int& j : i)
      {
        all.emplace(j, all.size());
        //en_.insert(en[i][j]);
      }
    }

    dim = all.size();

    /*
    for (map<int, int>::iterator its=all.begin(); its!=all.end(); its++) {
      if(ten_.find(its->first)==ten_.end()) {
        deque <int> first;
        first.push_back(its->first);
        ten.push_back(first);
      }

      if(en_.find(its->first)==en_.end()) {
        deque <int> first;
        first.push_back(its->first);
        en.push_back(first);
      }
    }
    */

    for (auto& i : ten)
    {
      for (int& j : i)
      {
        j = all[j];
      }

      sort(i.begin(), i.end());
    }

    for (auto& i : en)
    {
      for (int& j : i)
      {
        j = all[j];
      }

      std::sort(i.begin(), i.end());
    }
  }

  return (0.5*(2. - H_x_given_y(ten, en, dim) - H_x_given_y(en, ten, dim)));
}

double H_x_given_y3(std::deque<std::deque<int>> &en, std::deque<std::deque<int>> &ten, int dim)
{
  // you know y and you want to find x according to a certain index labelling.
  // so, for each x you look for the best y.

  std::deque<std::deque<int>> mems;

  for (int i = 0; i < dim; i++)
  {
    mems.emplace_back(std::deque<int>());
  }

  for (int ii = 0; ii < ten.size(); ii++)
  {
    for (int i = 0; i < ten[ii].size(); i++)
    {
      mems[ten[ii][i]].push_back(ii);
    }
  }

  double H_x_y = 0;
  double H2 = 0;

  for (auto& c : en)
  {
    std::deque <double> p;
    auto I2 = double(c.size());
    double O2 = (dim - I2);
    p.push_back(I2 / dim);
    p.push_back(O2 / dim);
    double H2_ = H(p);
    p.clear();

    double diff = H2_;

    // I need to know all the group wuth share nodes with en[k]

    std::map<int, int> com_ol;		// it maps the index of the ten into the overlap with en[k]

    for (int i : c)
    {
      for (int j = 0; j < mems[i].size(); j++)
      {
        int_histogram(mems[c[i]][j], com_ol);
      }
    }

    for (auto& itm : com_ol)
    {
      auto I1 = double(ten[itm.first].size());
      double O1 = (dim - I1);

      p.push_back(I1 / dim);
      p.push_back(O1 / dim);
      double H1_ = H(p);
      p.clear();

      double I1_I2 = itm.second;
      double I1_02 = ten[itm.first].size() - I1_I2;
      double O1_I2 = c.size() - I1_I2;
      double O1_02 = dim - I1_I2 - I1_02 - O1_I2;

      p.push_back(I1_I2 / dim);
      p.push_back(O1_02 / dim);

      double H12_positive = H(p);

      p.clear();
      p.push_back(I1_02 / dim);
      p.push_back(O1_I2 / dim);

      double H12_negative = H(p);

      double H12_ = H12_negative + H12_positive;

      p.clear();

      if (H12_negative > H12_positive)
      {
        H12_ = H1_ + H2_;
      }

      if ((H12_ - H1_) < diff)
      {
        diff = (H12_ - H1_);
      }
    }

    if (H2_ == 0)
    {
      H_x_y += 1;
    }
    else
    {
      H_x_y += (diff / H2_);
    }
  }

  return (H_x_y / (en.size()));
}

double mutual3(std::deque<std::deque<int>> en, std::deque<std::deque<int>> ten)
{
  if (en.empty()
    || ten.empty())
  {
    return 0;
  }

  // en e ten are two partitions of integer numbers
  int dim;

  {
    std::map <int, int> all;		// node, index
    //set <int> ten_;
    //set <int> en_;

    for (auto& i : ten)
    {
      for (int& j : i)
      {
        all.emplace(j, all.size());
        //ten_.insert(ten[i][j]);
      }
    }

    for (auto& i : en)
    {
      for (int& j : i)
      {
        all.emplace(j, all.size());
      }
    }

    dim = all.size();

    for (auto& i : ten)
    {
      for (int& j : i)
      {
        j = all[j];
      }

      std::sort(i.begin(), i.end());
    }

    for (auto& i : en)
    {
      for (int& j : i)
      {
        j = all[j];
      }

      std::sort(i.begin(), i.end());
    }
  }

  return (0.5 * (2. - H_x_given_y3(ten, en, dim) - H_x_given_y3(en, ten, dim)));
}

#endif
