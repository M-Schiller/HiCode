#include "static_network.h"

int main(int argc, char * argv[])
{
  //cout<<"the result is going to be written in blondel.part"<< std::endl;

  bool value;

  std::string s;
  std::string s2;

  if (parse_command_line(value, s, s2, argc, argv) == -1)
    return -1;

  static_network luca(s);
  //cout<<"network:: "<<luca.size()<<" nodes and "<<luca.edges()<<" edges;\t average degree = "<<2*luca.edges()/luca.size()<< std::endl;

  luca.draw_consecutive("netconsec.dat", "labelconsec.dat", true);

  systemCall("./convert -i netconsec.dat -o network.bin -w");
  systemCall("./community network.bin -l -1 -w > network.tree");
  systemCall("./hierarchy network.tree -l 1 > leve1.dat");

  std::deque<std::deque<int>> one_consec;
  get_partition_from_file_list("leve1.dat", one_consec);

  std::ofstream outt("louvain.part");
  luca.print_id(one_consec, outt);

  return 0;
}
