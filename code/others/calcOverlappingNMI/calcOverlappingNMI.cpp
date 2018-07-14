#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <sstream>
#include <cctype>
#include <fstream>

std::map<std::string, std::string> config;

#include "basic.h"
#include "Graph.h"
#include "Community.h"
#include "Metrics.h"

int main(int argc, char** argv)
{
  std::cout << calcNMI(argv[1], argv[2]) << std::endl;
}
