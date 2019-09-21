#ifndef __NOISE_NOISE_HPP_
#define __NOISE_NOISE_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

// Noise Generator Classes
//
// NoiseGenerator is the base class for each noise generation algorithm.
// generate returns a 1D vector which can be interpreted as a 2D array of
// size x_dim * y_dim. Each derived class implements generator which is that
// class's noise generation algorithm.
namespace Noise {
    struct Grad {
        float x, y, z;
        Grad() : x(0), y(0), z(0) {}
        Grad(float x, float y, float z) : x(x), y(y), z(z) {}
        float &operator[](std::size_t i)
        {
            // if (!(i > 0 && i < 3)) return nullptr;
            return (i == 0) ? x : (i == 1) ? y : z;
        }
    };

    class NoiseGenerator {
    protected:
        uint64_t seed;
        std::mt19937_64 mt;
        std::uniform_int_distribution<uint8_t> rand_byte{0,255};
        std::vector<double> normalize(const std::vector<double>& data, double a,
            double b);
    public:
        NoiseGenerator(uint64_t s);
        virtual double generator(double x, double y) = 0;
        std::vector<double> generate(int x_dim, int y_dim, double scale, int octaves,
            double persistence, double lacunarity);
    };

    class PerlinNoise : public NoiseGenerator {
    private:
        Grad normalize2(const Grad &g);
        float s_curve(const float& t);
        std::uniform_real_distribution<float> rand_float{-1.0,1.0};
        std::uniform_int_distribution<unsigned int> rand_uint{0,511};
        std::array<uint8_t, 514> p;
        std::array<float, 514> g1;
        std::array<Grad, 514> g2;
        std::array<Grad, 32> points;
    public:
        PerlinNoise(uint64_t s);
        double generator(double x, double y);
    };

    class SimplexNoise : public NoiseGenerator {
    private:
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
