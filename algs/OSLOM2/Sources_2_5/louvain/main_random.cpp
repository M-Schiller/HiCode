#include <cstddef>
#include <iostream>
#include <cstdlib>

char *outfile = nullptr;

int main(int argc, char **argv)
{
  srand(time(NULL) + getpid());

  int n = std::atoi(argv[1]);
  int degree = std::atoi(argv[2]);

  for (unsigned int i = 0; i < n; i++)
  {
    for (unsigned int j = 0; j < n; j++)
    {
      int r = rand() % n;
      if (r < degree)
      {
        std::cout << i << " " << j << std::endl;
      }
    }
  }
}
