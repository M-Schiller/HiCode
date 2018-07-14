#pragma once

#ifndef POSITION_H
#define POSITION_H

#include "./standard_package/standard_include.cpp"
#define bseparator 1.3

typedef std::multimap<double, int>::iterator map_idi;
typedef std::pair<map_idi, map_idi> pmap_idi;

class xypos
{
public:

  xypos() = default;
  ~xypos() = default;

  bool erase(int);
  void edinsert(int, double, double);
  int size() const
  {
    return xlabel.size();
  };
  void print_points(std::ostream &);
  void clear();
  int neighbors(int a, std::deque<int> & group_intsec);
  void set_all(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r, double);
  int move(int a);
  int move_all();
  void get_new_pos(std::map<int, double> & lab_x, std::map<int, double> & lab_y);

private:

  int xy_close(std::multimap<double, int> & m, double x, double dx, std::set<int> & close);

  std::map<int, double> label_radius;
  std::multimap <double, int> xlabel;
  std::multimap <double, int> ylabel;
  std::map<int, pmap_idi >  label_points;
  double center_x;
  double center_y;
};

inline void xypos::set_all(std::map<int, double> & lab_x, std::map<int, double> & lab_y, std::map<int, double> & lab_r, double L)
{
  // L is how far apart you want the module (1.2 should be ok)

  for (auto& itm : lab_x)
  {
    edinsert(itm.first, itm.second, lab_y[itm.first]);
  }

  for (auto& itm : lab_r)
  {
    label_radius[itm.first] = L * itm.second;
  }

  center_x = 0;
  center_y = 0;

  for (auto &[label, point] : label_points)
  {
    center_x += point.first->first;
    center_y += point.second->first;
  }

  center_x /= label_points.size();
  center_y /= label_points.size();
}

inline void xypos::clear()
{
  xlabel.clear();
  ylabel.clear();
  label_points.clear();
}

inline bool xypos::erase(int a)		// this function erases element a if exists (and returns true)
{
  auto itm = label_points.find(a);
  if (itm != label_points.end())
  {
    xlabel.erase(itm->second.first);
    ylabel.erase(itm->second.second);
    label_points.erase(itm);
    return true;
  }

  return false;
}

inline void xypos::edinsert(int a, double x, double y)		// this function inserts element a (or edit it if it was already inserted)
{
  erase(a);

  auto itfx = xlabel.emplace(x, a);
  auto itfy = ylabel.emplace(y, a);

  label_points.emplace(a, make_pair(itfx, itfy));
}

inline void xypos::print_points(std::ostream & outb) {
  outb << "#x, y, label" << std::endl;
  for (auto &[label, point] : label_points)
  {
    outb << point.first->first << "\t" << point.second->first << "\t" << label << std::endl;
  }
}

inline void xypos::get_new_pos(std::map<int, double> & lab_x, std::map<int, double> & lab_y)
{
  lab_x.clear();
  lab_y.clear();

  for (auto &[label, point] : label_points)
  {
    lab_x[label] = point.first->first;
    lab_y[label] = point.second->first;
  }
}

inline int xypos::neighbors(int a, std::deque<int> & group_intsec)
{
  auto itm = label_points.find(a);
  double x = (itm->second.first)->first;
  double y = (itm->second.second)->first;
  double delta = 2 * label_radius[a];

  // returns all the guys which are delta close to (x,y)

  //cout<<"neighs of "<<a<<endl;

  std::set<int> xclose;
  std::set<int> yclose;

  xy_close(xlabel, x, delta, xclose);
  //prints(xclose);
  //cout<<">>>>>>>>>>><<<<<<<<<<<<<"<<endl;

  xy_close(ylabel, y, delta, yclose);
  //prints(yclose);
  //cout<<">>>>>>>>>>><<<<<<<<<<<<<"<<endl;

  xclose.erase(a);
  yclose.erase(a);

  group_intsec.clear();
  set_intersection(xclose.begin(), xclose.end(), yclose.begin(), yclose.end(), back_inserter(group_intsec));

  return 0;
}

inline int xypos::xy_close(std::multimap<double, int> & m, double x, double dx, std::set<int> & close)
{
  //prints(m);

  //cout<<"looking for "<<x - 1e-2<<endl;
  auto xit = m.lower_bound(x - 1e-10);
  auto xit_m = xit;

  while (true)
  {
    if (xit == m.end())
      break;

    //cout<<"xit: "<<xit->first<<" "<<xit->second<<endl;

    if (xit->first - x <= dx)
    {
      close.insert(xit->second);
      //cout<<"-> "<<xit->first - x<<endl;
    }
    else
      break;

    ++xit;
  }

  while (true)
  {
    if (xit_m == m.begin())
      break;

    --xit_m;

    //cout<<"xitm: "<<xit_m->first<<" "<<xit_m->second<<endl;

    if (x - xit_m->first <= dx) {
      close.insert(xit_m->second);
      //cout<<"-> "<<x - xit_m->first<<endl;
    }
    else
      break;
  }

  /*cout<<"close:"<<endl;
  prints(close);*/
}

inline int xypos::move(int a)
{
  auto itm = label_points.find(a);
  double x = (itm->second.first)->first;
  double y = (itm->second.second)->first;

  /*cout<<"moving::: "<<a<<endl;
  cout<<"x... "<<x<<" "<<y<<endl;*/

  std::deque<int> neighs;
  neighbors(a, neighs);

  std::deque<int> eff;

  //cout<<"neighbors ************** "<<endl;
  for (int neigh : neighs)
  {
    auto itm2 = label_points.find(neigh);
    double x2 = (itm2->second.first)->first;
    double y2 = (itm2->second.second)->first;

    //cout<<neighs[i]<<" "<<x2<<" "<<y2<<endl;

    double l = (x - x2)*(x - x2) + (y - y2)*(y - y2);

    //cout<<"l: "<<sqrt(l)<<" "<<label_radius[a] + label_radius[neighs[i]]<<endl;
    if (sqrt(l) < label_radius[a] + label_radius[neigh])
    {
      eff.push_back(neigh);
    }
  }

  //cout<<"*****************"<<endl;

  if (eff.empty())
    return -1;

  for (int rn : eff)
  {
    const std::map<int, pmap_idi>::iterator itm2 = label_points.find(rn);
    double x2 = (itm2->second.first)->first;
    double y2 = (itm2->second.second)->first;

    double di = sqrt((x - x2)*(x - x2) + (y - y2)*(y - y2));
    double bi = bseparator * (label_radius[a] + label_radius[rn] - di);

    double o1 = (x - center_x)*(x - center_x) + (y - center_y)*(y - center_y);
    double o2 = (x2 - center_x)*(x2 - center_x) + (y2 - center_y)*(y2 - center_y);

    //cout<<di<<"<----"<<endl;

    if (o1 > o2)
    {
      double cos_, sin_;
      if (di > 0)
      {
        cos_ = (x - x2) / di;
        sin_ = (y - y2) / di;
      }
      else
      {
        cos_ = 1.;
        sin_ = 0;
      }

      const double nx = x + bi * cos_;
      const double ny = y + bi * sin_;
      edinsert(a, nx, ny);
    }

    else
    {
      double cos_, sin_;
      if (di > 0) {
        cos_ = (x2 - x) / di;
        sin_ = (y2 - y) / di;
      }
      else
      {
        cos_ = 1.;
        sin_ = 0;
      }

      const double nx = x2 + bi * cos_;
      const double ny = y2 + bi * sin_;
      edinsert(rn, nx, ny);
    }
  }
  //cout<<center_x<<" "<<center_y<<" "<<o1<<" "<<rn<<endl;

  return 0;
}

inline int xypos::move_all()
{
  std::deque<int> points;
  for (auto& label_point : label_points)
  {
    points.push_back(label_point.first);
  }

  int stopper = 0;

  while (stopper < 100)
  {
    stopper++;
    std::cout << "checking overlaps... " << stopper << std::endl;

    bool again = false;
    for (int point : points)
    {
      int h = move(point);
      if (h == 0)
        again = true;
    }

    if (!again)
    {
      break;
    }
  }

  return 0;
}
#endif // POSITION_H
