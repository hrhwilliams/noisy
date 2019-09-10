#include <iostream>
#include <cmath>

#include "encoder.hpp"
#include "noise.hpp"

std::vector<uint8_t> convert_vec(const std::vector<double>& data)
{
  std::vector<uint8_t> bytes;

  for (auto &&x: data) {
    bytes.push_back(std::round(x * 255));
  }

  return bytes;
}

int main()
{
  auto noise = Noise::SimplexNoise(5);
  auto encoder = PNG::PngEncoder();
  auto image = noise.generate(256, 256, 100, 4, 0.333, 2.5);

  encoder.open("noise.png");
  encoder.write_ihdr(256, 256, 8, PNG::ColorType::Grayscale);
  encoder.write_idat(convert_vec(image));
  encoder.write_iend();
  encoder.close();
  return 0;
}
