#include "standard_package/standard_include.cpp"

int main(int argc, char * argv[])
{
  if (argc < 4) {
    std::cerr << argv[0] << " " << "[filename] [seed] [number_of_runs]" << " -- output is infomap_net.net infomap.part" << std::endl;
    return 0;
  }

  pajek_format(std::string(argv[1]), false);
  systemCall("mv net.paj infomap_net.net");

  char bbs[200];
  sprintf(bbs, "./infomap_undir %s infomap_net.net %s", argv[2], argv[3]);

  std::cout << "running: " << bbs << std::endl;
  systemCall(bbs);

  std::deque<int> oldlabs;
  std::ifstream lin("labels.dat");
  int a;
  while (lin >> a)
    oldlabs.push_back(a);

  std::deque<std::deque<int>> one;
  get_partition_from_file_list_pajek("infomap_net.clu", one, oldlabs);
  std::ofstream oneout("infomap.part");
  for (auto& i : one)
  {
    prints(i, oneout);
  }

  return 0;
}
