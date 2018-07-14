#include "static_network.h"

int main(int argc, char * argv[])
{
  systemCall("make clean");
  systemCall("make");
  systemCall("ls > list_out");

  std::string s;
  std::ifstream gin;
  gin.open("list_out", std::ios::in);
  bool compiled = false;

  while (gin >> s)
  {
    if (s == "convert")
      compiled = true;
  }
  gin.close();
  if (compiled)
    return 0;

  std::cout << "error in compiling" << std::endl;
  return -1;
}
