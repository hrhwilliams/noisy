#include <iostream>

#include "encoder.hpp"
#include "noise.hpp"

int main()
{
  auto noise = Noise::SimplexNoise(5);
  std::cout << noise.generator(0.1, 0.1) << ", "
    << noise.generator(0.101, 0.101) << std::endl;
  return 0;
}
