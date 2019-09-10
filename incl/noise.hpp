#ifndef __NOISE_NOISE_HPP_
#define __NOISE_NOISE_HPP_

#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

namespace Noise {
  class NoiseGenerator {
  protected:
    uint64_t seed;
    std::mt19937_64 mt;
    std::uniform_int_distribution<uint8_t> rand_byte{0,255};
  public:
    NoiseGenerator(uint64_t s);
    virtual double generator(double x, double y) = 0;
    std::vector<double> generate(int x_dim, int y_dim, double scale, int octaves,
      double persistence, double lacunarity);
  };

  // class PerlinNoise : public NoiseGenerator {
  //
  // };

  class SimplexNoise : public NoiseGenerator {
  private:
    struct Grad {
      double x, y, z;
      Grad(double x, double y, double z) : x(x), y(y), z(z) {}
    };
    double dot(Grad g1, Grad g2) { return g1.x * g2.x + g1.y * g2.y + g1.z * g2.z; }
    const std::array<Grad, 12> gradients =
      { Grad(1,1,0), Grad(-1,1,0), Grad(1,-1,0), Grad(-1,-1,0),
        Grad(1,0,1), Grad(-1,0,1), Grad(1,0,-1), Grad(-1,0,-1),
        Grad(0,1,1), Grad(0,-1,1), Grad(0,1,-1), Grad(0,-1,-1) };
    std::array<uint8_t, 512> perm_table;
    std::array<uint8_t, 512> perm_mod12;
    double skewing_factor = 0.5 * (std::sqrt(3.0) - 1.0);
    double unskewing_factor = (3.0 - std::sqrt(3.0)) / 6.0;
  public:
    SimplexNoise(uint64_t s);
    double generator(double x, double y);
  };

  // class OpenSimplexNoise : public NoiseGenerator {
  //
  // };
};

#endif // __NOISE_NOISE_HPP_
