#include <iostream>
#include <algorithm>
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

std::vector<double> normalize(const std::vector<double>& data, double a,
  double b)
{
  std::vector<double> normalized;
  double min = *std::min_element(data.begin(), data.end());
  double max = *std::max_element(data.begin(), data.end());

  // copy[y][x] = a + ((copy[y][x] - data_min) * (b - a)) / (data_max - data_min)
  for (auto&& x : data) {
    normalized.push_back(a + ((x - min) * (b - a)) / (max - min));
  }

  return normalized;
}

int main()
{
  auto noise = Noise::SimplexNoise(5);
  auto encoder = PNG::PngEncoder();
  auto image = noise.generate(256, 256, 100, 4, 0.333, 2.5);

  encoder.open("noise.png");
  encoder.write_ihdr(256, 256, 8, PNG::ColorType::Grayscale);
  encoder.write_idat(convert_vec(normalize(image, 0, 1)));
  encoder.write_iend();
  encoder.close();
  return 0;
}
