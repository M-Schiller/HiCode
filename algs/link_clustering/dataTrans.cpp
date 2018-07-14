#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <algorithm> // for swap

int main(int argc, char const *argv[])
{
  //printf("%f", pow(2, -3));

  std::ifstream inFile;
  inFile.open(argv[1]);
  if (!inFile)
  {
    std::cout << "ERROR: unable to open input file" << std::endl;
    exit(1); // terminate with error
  }
  std::ofstream outfile;
  outfile.open(argv[2], std::ios::out);

  int ni, nj, wij, index = 0;
  int tempV; //inFile >> tempV;
  while (inFile >> ni >> nj)
  {
    ni--;
    nj--;
    outfile << ni << "\t" << nj << std::endl;
  }
  inFile.close(); inFile.clear();
  outfile.close();
  return 0;
}
