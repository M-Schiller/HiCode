#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <utility>   // for pairs
#include <algorithm> // for swap
#include<unistd.h>
#include <sys/time.h>

#include "clusterJaccsFile.h"

int main(int argc, char const *argv[]) {
  //************* make sure args are present:
  if (argc != 5) {
    std::cout << "ERROR: something wrong with the inputs" << std::endl;
    std::cout << "usage:\n    " << argv[0] << " network.pairs network.jaccs network.clusters network.cluster_stats threshold" << std::endl;
    exit(1);
  }

  //for saving the time
  struct timeval start, end;
  gettimeofday(&start, NULL);

  double maxDm = -1000000.0;
  double Dmns = 0.0;
  double counti = 0.1;
  for (int i = 1; i <= 9; ++i)
  {
    double tmpd = clusterJaccsFile(argc, argv, maxDm, counti);
    //cout<<"the tmpd:"<<tmpd<< std::endl;
    if (tmpd > maxDm)
      maxDm = tmpd;
    else
      break;
    counti += 0.1;
  }

  gettimeofday(&end, NULL);
  std::ifstream inTime;  inTime.open("./LCtime.txt");
  double all;
  inTime >> all;
  inTime.close();

  double span = end.tv_sec - start.tv_sec + (end.tv_usec - start.tv_usec) / 1000000.0;

  std::ofstream otime;
  otime.open("./LCtime.txt", std::ios::out);
  otime << all + span << std::endl;
  otime.close();

  return 0;
}
