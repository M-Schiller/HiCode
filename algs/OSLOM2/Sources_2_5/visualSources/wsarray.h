#pragma once

#ifndef WSARRAY_INCLUDED
#define WSARRAY_INCLUDED

// this class is a container of int and double;
// it needs a preallocation when constructing the class;
// then you need to use the function push_back() to insert the numbers; this is the only way
// once you did, you can freeze the numbers; multiple entries are neglected; the preallocation can be reduced; you cannot use push_back anymore
// find returns the position of the integer you are looking for;
// posweightof returns a pair position-weight
// to access the integers use l[] and w[]
// you can erase elements (after you freezed them) keeping the wsarray sorted and right located

class wsarray {
public:

  wsarray(int a);
  ~wsarray();

  int find(int);
  std::pair <int, double> posweightof(int x);
  int size() const
  {
    return position;
  };
  void freeze();
  void push_back(int, double);
  bool erase(int);		// erase this int and its double

  double* w;
  int* l;

private:

  int _size_;
  int position;
};

inline wsarray::wsarray(int a)
{
  position = 0;
  _size_ = a;

  l = new int[_size_];
  w = new double[_size_];
}

inline wsarray::~wsarray()
{
  delete[] l;
  l = nullptr;

  delete[] w;
  w = nullptr;
}

inline std::pair <int, double> wsarray::posweightof(int x)
{
  int i = find(x);
  if (i == -1)
    return (std::make_pair(-1, 0));

  return (std::make_pair(i, w[i]));
};

int wsarray::find(int a) {
  int one = 0;
  int two = position - 1;

  if (position == 0)
    return -1;

  if (a<l[one] || a>l[two])
    return -1;

  if (a == l[one])
    return one;

  if (a == l[two])
    return two;

  while (two - one > 1) {
    int middle = (two - one) / 2 + one;

    if (a == l[middle])
      return middle;

    if (a > l[middle])
      one = middle;
    else
      two = middle;
  }

  return -1;
}

inline void wsarray::push_back(int a, double b)
{
  l[position] = a;
  w[position++] = b;
}

inline bool wsarray::erase(int y)
{
  int pos_of_y = find(y);
  if (pos_of_y == -1)
    return false;

  //cout<<"you want to erase "<<y<<endl;

  _size_--;
  position = _size_;

  //cout<<"_size_ "<<_size_<<endl;

  int ll[_size_];
  double ww[_size_];

  int poi = 0;
  for (int i = 0; i < _size_ + 1; i++) if (l[i] != y) {
    ll[poi] = l[i];
    ww[poi] = w[i];
    poi++;
  }

  delete[] l;
  l = nullptr;

  delete[] w;
  w = nullptr;

  l = new int[_size_];
  w = new double[_size_];

  for (int i = 0; i < _size_; i++)
  {
    l[i] = ll[i];
    w[i] = ww[i];
  }

  return true;
}

inline void wsarray::freeze()
{
  std::map<int, double> M;
  for (int i = 0; i < position; i++)
  {
    /*	//this is to sum up multiple entries
    map<int, double>::iterator itf=M.find(l[i]);
    if (itf==M.end())
      M.insert(make_pair(l[i], w[i]));
    else
      itf->second+=w[i];
    //*/

    M.emplace(l[i], w[i]);
  }

  if (_size_ != M.size())
  {
    delete[] l;
    l = nullptr;

    delete[] w;
    w = nullptr;

    _size_ = M.size();
    position = M.size();

    l = new int[_size_];
    w = new double[_size_];
  }

  int poi = 0;
  for (auto itm = M.begin(); itm != M.end(); ++itm) {
    l[poi] = itm->first;
    w[poi] = itm->second;
    poi++;
  }
}

inline void prints(wsarray &a) {
  for (int i = 0; i < a.size(); i++)
    std::cout << a.l[i] << "\t" << a.w[i] << std::endl;
  std::cout << std::endl;
}

inline void prints(wsarray &a, std::ostream &out) {
  for (int i = 0; i < a.size(); i++)
    std::cout << a.l[i] << "\t" << a.w[i] << std::endl;
  out << std::endl;
}

inline void prints(wsarray *a, std::ostream & out) {
  for (int i = 0; i < a->size(); i++)
    std::cout << a->l[i] << "\t" << a->w[i] << std::endl;
  out << std::endl;
}

void prints(wsarray *a) {
  prints(a, std::cout);
}

#endif
