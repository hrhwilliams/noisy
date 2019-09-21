#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "encoder.hpp"
#include "noise.hpp"

std::vector<uint8_t> convert_vec(const std::vector<float>& data)
{
  std::vector<uint8_t> bytes;

  for (auto &&x: data) {
    bytes.push_back(std::round(x * 255));
  }

  return bytes;
}

int main(int argc, char *argv[])
{
  std::random_device rd;
  std::stringstream seed_str;
  uint64_t seed;

  if (argc < 2) {
    // Generate a 64-bit seed
    seed = rd();
    seed <<= 32;
    seed += rd();
  } else {
    seed = std::atol(argv[1]);
  }
  seed_str << seed << ".png";

  auto noise = Noise::SimplexNoise(seed);
  auto encoder = PNG::PngEncoder();

  int x_dim, y_dim, octaves;
  float scale, persistence, lacunarity;

  std::cout << "Image width: ";
  std::cin >> x_dim;
  std::cout << "Image height: ";
  std::cin >> y_dim;
  std::cout << "Scale: ";
  std::cin >> scale;
  std::cout << "Octaves: ";
  std::cin >> octaves;
  std::cout << "Persistence: ";
  std::cin >> persistence;
  std::cout << "Lacunarity: ";
  std::cin >> lacunarity;

  auto image = noise.generate(x_dim, y_dim, scale, octaves, persistence,
    lacunarity);

  encoder.open(seed_str.str());
  encoder.write_ihdr(x_dim, y_dim, 8, PNG::ColorType::Grayscale);
  encoder.write_idat(convert_vec(image));
  encoder.write_iend();
  encoder.close();

  std::cout << "Wrote out " << seed_str.str() << std::endl;

  return 0;
}
