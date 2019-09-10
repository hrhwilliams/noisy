#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "noise.hpp"

using namespace Noise;

// std::random_device dev;
// auto seed = dev();

// def generate(self, x_dim, y_dim, scale, octaves, persistence, lacunarity):
//     data = [[0 for _ in range(x_dim)] for _ in range(y_dim)]
//
//     for y in range(y_dim):
//         for x in range(x_dim):
//             amplitude = 1
//             frequency = 1
//             noise_height = 0
//
//             for _ in range(octaves):
//                 sample_x = (x / scale) * frequency
//                 sample_y = (y / scale) * frequency
//                 noise = self._generator(sample_x, sample_y)
//                 noise_height += noise * amplitude
//
//                 amplitude *= persistence
//                 frequency *= lacunarity
//
//             data[y][x] = noise_height
//     return normalize(data, 0, 1)

NoiseGenerator::NoiseGenerator(uint64_t seed)
{
	mt(seed);
}

std::vector<double> NoiseGenerator::generate(int x_dim, int y_dim, float scale,
  int octaves, float persistence, float lacunarity);
{
  std::vector<double> data;
  double amplitude = 1.0;
  double frequency = 1.0;
  double noise_height = 0.0;

  for (int y = 0; y < y_dim; y++) {
    for (int x = 0; x < x_dim; x++) {
      amplitude = 1.0;
      frequency = 1.0;
      noise_height = 0.0;

      for (int i = 0; i < octaves; i++) {
        sample_x = (x / scale) * frequency;
        sample_y = (y / scale) * frequency;
        noise = this->generator(sample_x, sample_y);
        noise_height += noise * amplitude

        amplitude *= persistence;
        frequency *= lacunarity;
      }
      data.push_back(noise_height);
    }
  }

  return data;
}

SimplexNoise::SimplexNoise(uint64_t seed) : NoiseGenerator(seed)
{
	for (int i = 0; i < perm_table.size(); i++) {
		perm_table[i] = rand_byte(mt);
		perm_mod12[i] = perm_table[i] % 12;
	}
}

double SimplexNoise::generator(double x, double y)
{
  double n0, n1, n2;
  double s = (x + y) * skewing_factor;
  int i = (int) (x + s);
  int j = (int) (y + s);
  double t = (i + j) * unskewing_factor;
  double x_origin = i - t;
  double y_origin = j - t;
  double x0 = x - x_origin;
  double y0 = y - y_origin;

  int i1, j1;
  if(x0 > y0) {
    i1 = 1;
    j1 = 0;
  } else {
    i1 = 0;
    j1 = 1;
  }

  double x1 = x0 - i1 + unskewing_factor;
  double y1 = y0 - j1 + unskewing_factor;
  double x2 = x0 - 1.0 + 2.0 * unskewing_factor;
  double y2 = y0 - 1.0 + 2.0 * unskewing_factor;

  int ii = i & 255;
  int jj = j & 255;
  int gi0 = perm_mod12[ii + perm_table[jj]];
  int gi1 = perm_mod12[ii + i1 + perm_table[jj + j1]];
  int gi2 = perm_mod12[ii + 1 + perm_table[jj + 1]];

  double t0 = 0.5 - x0 * x0 - y0 * y0;
  if(t0 < 0) n0 = 0.0;
  else {
    t0 *= t0;
    n0 = t0 * t0 * dot(gradients[gi0], Grad(x0, y0, 0));
  }
  double t1 = 0.5 - x1 * x1 - y1 * y1;
  if(t1 < 0) n1 = 0.0;
  else {
    t1 *= t1;
    n1 = t1 * t1 * dot(gradients[gi1], Grad(x1, y1, 0));
  }
  double t2 = 0.5 - x2*x2-y2*y2;
  if(t2 < 0) n2 = 0.0;
  else {
    t2 *= t2;
    n2 = t2 * t2 * dot(gradients[gi2], Grad(x2, y2, 0));
  }
  // Add contributions from each corner to get the final noise value.
  // The result is scaled to return values in the interval [-1,1].
  return 70.0 * (n0 + n1 + n2);
}
