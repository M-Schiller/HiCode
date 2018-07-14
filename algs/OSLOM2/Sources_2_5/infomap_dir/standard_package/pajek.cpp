#include <deque>
#include <map>
#include <iostream>
#include <fstream>

int number_together(std::deque<int> & a)
{
  std::string s;

  //cout<<"a"<< std::endl;
  //prints(a);

  char b[100];

  for (int i = 0; i < a.size(); i++)
  {
    sprintf(b, "%d", a[i]);
    s += b;
  }

  //cout<<s<< std::endl;
  int l = cast_int(cast_string_to_double(s));
  //cout<<"l::::  "<<l<< std::endl;
  return l;
}

int get_partition_from_file_list_pajek(const std::string& s, std::deque<std::deque<int>> & ten, std::deque<int> & oldlabels)
{
  ten.clear();

  std::ifstream inb(s);
  std::string sst;
  getline(inb, sst);

  int a = 1;
  int bb;

  std::map<int, int> id_value;
  std::map <int, int> value_com;
  while (inb >> bb)
  {
    value_com.insert(std::make_pair(bb, value_com.size()));
    id_value.insert(std::make_pair(a, bb));
    a++;
  }

  int com = 0;
  for (auto it = value_com.begin(); it != value_com.end(); ++it)
  {
    it->second = com++;
  }

  std::deque <int> first;
  for (int i = 0; i < value_com.size(); i++)
  {
    ten.push_back(first);
  }

  for (auto itm = id_value.begin(); itm != id_value.end(); ++itm)
  {
    ten[value_com[itm->second]].push_back(oldlabels[itm->first - 1]);
  }

  return 0;
}

int set_partition_from_list(std::deque<int> & mems, std::deque<std::deque<int>> & ten)
{
  ten.clear();

  //cout<<"mems"<< std::endl;
  //prints(mems);

  int a = 1;

  std::map<int, int> id_value;
  std::map <int, int> value_com;
  for (int i = 0; i < mems.size(); i++)
  {
    value_com.insert(std::make_pair(mems[i], value_com.size()));
    id_value.insert(std::make_pair(a, mems[i]));
    ++a;
  }

  int com = 0;
  for (auto it = value_com.begin(); it != value_com.end(); ++it)
  {
    it->second = com++;
  }

  std::deque <int> first;
  for (int i = 0; i < value_com.size(); i++)
  {
    ten.push_back(first);
  }

  for (auto itm = id_value.begin(); itm != id_value.end(); ++itm)
  {
    ten[value_com[itm->second]].push_back(itm->first);
  }
}

int get_partition_from_file_list_pajek_tree(const std::string& s, std::deque<std::deque<std::deque<int>>> & Ten)
{
  Ten.clear();

  std::ifstream inb(s);
  std::string sst;
  std::getline(inb, sst);

  std::deque<std::deque<int>> TREE;
  while (getline(inb, sst))
  {
    std::deque<int> sss;
    cast_string_to_doubles(sst, sss);
    TREE.push_back(sss);
  }

  if (TREE.empty())
    return 0;

  std::deque<int> mems(TREE.size());

  //int tree_length=TREE[0].size();

  for (int ih = 0; ih < 2; ih++)		// I want the first two partitions
  {
    for (int i = 0; i < TREE.size(); i++)
    {
      int finest_level = TREE[i].size() - 3;
      int coarse_level = finest_level - 1;

      std::deque<int> history;
      int this_l = 0;
      if (ih == 0)
        this_l = coarse_level;
      else
        this_l = finest_level;

      for (int k = 0; k < this_l; k++)
        history.push_back(TREE[i][k]);

      //prints(history);

      //cout<<y<<"-<<<"<< std::endl;
      //cout<<(TREE[i][tree_length-1]) -1<<" "<<tree_length<<" "<<TREE[i][TREE[i].size()-1]<<" "<<TREE[i].size()<< std::endl;
      mems[(TREE[i][TREE[i].size() - 1]) - 1] = number_together(history);
      //cout<<"oks"<< std::endl;
    }
    std::deque<std::deque<int>> ten;
    set_partition_from_list(mems, ten);
    Ten.push_back(ten);
  }

  std::cout << "number of levels: " << Ten.size() << std::endl;

  /*for(int i=0; i<Ten.size(); i++) {
    std::cout<<"-----------------------"<< std::endl;
    printm(Ten[i]);
  }
  */
  return 0;
}

int pajek_format(const std::string& filename, bool directed)
{
  // filename is the edge list, directed=true if directed

  std::map<int, int> newlabels;
  std::deque<int> old_labels;

  {
    int label = 0;
    std::ifstream inb(filename);
    std::string ins;

    while (std::getline(inb, ins))
    {
      if (!ins.empty() && ins[0] != '#')
      {
        std::deque<double> ds;
        cast_string_to_doubles(ins, ds);

        if (ds.size() < 2)
        {
          std::cout << "From file " << filename << ": string not readable:" << ins << std::endl;
        }
        else
        {
          int innum1 = cast_int(ds[0]);
          int innum2 = cast_int(ds[1]);

          auto itf = newlabels.find(innum1);
          if (itf == newlabels.end())
          {
            newlabels.insert(std::make_pair(innum1, label++));
            old_labels.push_back(innum1);
          }

          itf = newlabels.find(innum2);
          if (itf == newlabels.end())
          {
            newlabels.insert(std::make_pair(innum2, label++));
            old_labels.push_back(innum2);
          }
        }
      }
    }
  }

  int dim = newlabels.size();
  std::ofstream out("net.paj");
  out << "*Vertices " << dim << std::endl;

  for (int i = 0; i < dim; i++) {
    out << i + 1 << " \"" << old_labels[i] << "\"" << std::endl;
  }

  if (directed == false)
    out << "*Edges " << std::endl;
  else
    out << "*Arcs " << std::endl;

  std::set<std::pair<int, int>> _pairs_;
  {
    std::ifstream inb(filename);
    std::string ins;
    while (std::getline(inb, ins))
    {
      if (!ins.empty()
        && ins[0] != '#')
      {
        std::deque<std::string> ds;
        separate_strings(ins, ds);

        if (ds.size() > 1)
        {
          //prints(ds);

          int innum1 = cast_int(cast_string_to_double(ds[0]));
          int innum2 = cast_int(cast_string_to_double(ds[1]));

          double w = 1;
          if (ds.size() > 2)
          {
            if (!cast_string_to_double(ds[2], w))
            {
              w = 1;
            }
          }

          int new1 = newlabels[innum1];
          int new2 = newlabels[innum2];

          //cout<<"***********************************"<< std::endl;
          //cout<<innum1<<" "<<innum2<<" "<<new1<<" "<<new2<< std::endl;

          if (!directed)
          {
            int an1 = std::min(new1, new2);
            int an2 = std::max(new1, new2);
            new1 = an1;
            new2 = an2;
          }

          std::pair<int, int> pp(new1, new2);
          if (_pairs_.find(pp) == _pairs_.end() && new1 != new2)
          {
            out << new1 + 1 << " " << new2 + 1 << " " << w << std::endl;
            _pairs_.insert(pp);
          }
        }
      }
    }
  }

  std::ofstream lout("labels.dat");
  prints(old_labels, lout);
  return 0;
}
