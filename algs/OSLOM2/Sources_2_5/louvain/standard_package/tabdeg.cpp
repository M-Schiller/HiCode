#include <cmath>
#include <iostream>
#include <deque>
#include <set>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <ctime>
#include <iterator>
#include <algorithm>

typedef std::multiset<std::pair<double, int>>  muspi;

class tabdeg
{
public:
  tabdeg() = default;
  ~tabdeg() = default;

  bool is_internal(int);
  void edinsert(int, double);
  bool erase(int);
  double indegof(int);
  int size() const { return nodes_indeg.size(); };
  void print_nodes(std::ostream &);
  int best_node();
  int worst_node();
private:
  std::map <int, muspi::iterator> nodes_indeg;				// for each node belonging to the cluster, here you have where to find the internal degree
  muspi number_label;
};

bool tabdeg::is_internal(int a)
{
  return nodes_indeg.find(a) != nodes_indeg.end();
}

void tabdeg::edinsert(int a, double kp)		// this function inserts element a (or edit it if it was already inserted)
{
  erase(a);
  auto itms = number_label.insert(std::make_pair(kp, a));
  nodes_indeg.insert(std::make_pair(a, itms));
}

bool tabdeg::erase(int a)		// this function erases element a if exists (and returns true)
{
  const auto itm = nodes_indeg.find(a);
  if (itm != nodes_indeg.end())
  {
    number_label.erase(itm->second);
    nodes_indeg.erase(itm);
    return true;
  }
  return false;
}

double tabdeg::indegof(int a)		// return the internal degree of a, 0 if it's not internal
{
  const auto itm = nodes_indeg.find(a);
  if (itm != nodes_indeg.end())
  {
    return itm->second->first;
  }
  return 0;
}

void tabdeg::print_nodes(std::ostream & outb)
{
  for (auto& itm : nodes_indeg)
  {
    outb << itm.first << "\t" << itm.second->first << std::endl;
  }
}

int tabdeg::best_node()
{
  auto itm = number_label.end();
  if (number_label.empty())
    return -1;
  --itm;
  return itm->second;
}

int tabdeg::worst_node()
{
  const auto itm = number_label.begin();
  if (number_label.empty())
    return -1;
  return itm->second;
}

int main()
{
  tabdeg C;
  C.edinsert(2, 14);
  C.edinsert(7, 24);
  C.edinsert(10, -1.2);
  C.edinsert(11, 12.8);

  std::cout << "--------------------" << std::endl;
  C.print_nodes(std::cout);

  C.erase(11);

  std::cout << "--------------------" << std::endl;
  C.print_nodes(std::cout);

  std::cout << "indegof " << C.indegof(10) << std::endl;

  const int best_node = C.worst_node();
  std::cout << "best node " << best_node << std::endl;
  C.erase(best_node);
  std::cout << "--------------------" << std::endl;
  C.print_nodes(std::cout);
  return 0;
}
