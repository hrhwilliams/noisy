#ifndef __NOISE_NOISE_HPP_
#define __NOISE_NOISE_HPP_

#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

namespace Noise {
  class NoiseGenerator {
  private:
    uint64_t seed;
    std::mt19937 mt19937(seed);
    std::uniform_int_distribution<uint8_t> rand_byte(0,255);
  public:
    NoiseGenerator(uint32_t seed) : seed(seed) {};
    virtual float generator(float x, float y) const = 0;
    std::vector<double> generate(int x_dim, int y_dim, float scale, int octaves,
      float persistence, float lacunarity);
  };

  // class PerlinNoise : public NoiseGenerator {
  //
  // };

  class SimplexNoise : NoiseGenerator {
  private:
    uint64_t seed;
    std::mt19937 mt19937();
    std::uniform_int_distribution<uint8_t> rand_byte(0,255);
    struct Grad {
      float x, y, z;
      Grad(float x, float y, float z) : x(x), y(y), z(z) {}
    };
    float dot(Grad g1, Grad g2) { return g1.x * g2.x + g1.y * g2.y + g1.z * g2.z; }
    const std::array<Grad, 12> gradients =
      { Grad(1,1,0), Grad(-1,1,0), Grad(1,-1,0), Grad(-1,-1,0),
        Grad(1,0,1), Grad(-1,0,1), Grad(1,0,-1), Grad(-1,0,-1),
        Grad(0,1,1), Grad(0,-1,1), Grad(0,1,-1), Grad(0,-1,-1) };
    std::array<uint8_t, 512> perm_table;
    std::array<uint8_t, 512> perm_mod12;
    static const double skewing_factor = 0.5 * (sqrt(3.0) - 1.0);
    static const double unskewing_factor = (3.0 - sqrt(3.0)) / 6.0;
  public:
    SimplexNoise(uint32_t seed);
    double generator(double x, double y);
  };

  // class OpenSimplexNoise : public NoiseGenerator {
  //
  // };
};

#endif // __NOISE_NOISE_HPP_
