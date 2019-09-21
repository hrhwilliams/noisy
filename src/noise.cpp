#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

#include "noise.hpp"

using namespace Noise;

float lerp(float t, float x, float y)
{
    return x + t * (y - x);
}

NoiseGenerator::NoiseGenerator(uint64_t s) : seed(s)
{
    mt = std::mt19937_64(seed);
}

std::vector<float> NoiseGenerator::normalize(const std::vector<float>& data,
  float a, float b)
{
    std::vector<float> normalized;
    float min = *std::min_element(data.begin(), data.end());
    float max = *std::max_element(data.begin(), data.end());

    for (auto&& x : data) {
        normalized.push_back(a + ((x - min) * (b - a)) / (max - min));
    }

    return normalized;
}

std::vector<float> NoiseGenerator::generate(int x_dim, int y_dim, float scale,
  int octaves, float persistence, float lacunarity)
{
    std::vector<float> data;
    float amplitude = 1.0;
    float frequency = 1.0;
    float noise_height = 0.0;

    for (int y = 0; y < y_dim; y++) {
        for (int x = 0; x < x_dim; x++) {
            float amplitude = 1.0;
            float frequency = 1.0;
            float noise_height = 0.0;

            for (int i = 0; i < octaves; i++) {
                float sample_x = (x / scale) * frequency;
                float sample_y = (y / scale) * frequency;
                float noise = this->generator(sample_x, sample_y);
                noise_height += noise * amplitude;

                amplitude *= persistence;
                frequency *= lacunarity;
            }
        data.push_back(noise_height);
        }
    }

    return normalize(data, 0.0, 1.0);
}

// std::vector<std::vector<float>> NoiseGenerator::generate_2d_array()

// ----- //

Grad PerlinNoise::normalize2(const Grad &g)
{
    float norm = std::sqrt(g.x * g.x + g.y * g.y + g.z * g.z);
    return Grad(g.x / norm, g.y / norm, g.z / norm);
}

float PerlinNoise::s_curve(const float& t)
{
    return t * t * (3 - t - t);
}

PerlinNoise::PerlinNoise(uint64_t s) : NoiseGenerator(s)
{
    int i, j, k;
    float u, v, w, U, V, W, Hi, Lo;
    for (i = 0; i < 256; i++) {
        p[i] = i;
        g1[i] = rand_float(mt);

        do {
            u = rand_float(mt);
            v = rand_float(mt);
        } while (u * u + v * v > 1 || std::abs(u) > 2.5 * std::abs(v)
            || std::abs(v) > 2.5 * std::abs(u)
            || std::abs(std::abs(u) - std::abs(v)) < .4);
        g2[i] = normalize2(Grad(u, v, 0));

        do {
            u = rand_float(mt);
            v = rand_float(mt);
            w = rand_float(mt);
            U = std::abs(u);
            V = std::abs(v);
            W = std::abs(w);
            Lo = std::min(U, std::min(V, W));
            Hi = std::max(U, std::max(V, W));
        } while (u * u + v * v + w * w > 1 || Hi > 4 * Lo
            || std::min(std::abs(U - V), std::min(std::abs(U - W), std::abs(V - W))) < .2);
    }

    while (--i > 0) {
        k = p[i];
        j = (int) (rand_byte(mt) & 255);
        p[i] = p[j];
        p[j] = k;
    }
    for (i = 0; i < 256 + 2; i++) {
        p[256 + i] = p[i];
        g1[256 + i] = g1[i];
        for (j = 0; j < 2; j++) {
            g2[256 + i] = g2[i];
        }
    }

    points[3] = Grad(std::sqrt(1.0f / 3.0f), std::sqrt(1.0f / 3.0f), std::sqrt(1.0f / 3.0f));
    float r2 = std::sqrt(1.0f / 2.0f);
    float sqr = std::sqrt(2.0f + r2 + r2);

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            points[i][j] = (i == j ? 1 + r2 + r2 : r2) / sqr;
        }
    }
    for (i = 0; i <= 1; i++) {
        for (j = 0; j <= 1; j++) {
            for (k = 0; k <= 1; k++) {
                int n = i + j * 2 + k * 4;
                if (n > 0) {
                    for (int m = 0; m < 4; m++) {
                        points[4 * n + m] = Grad((i == 0 ? 1 : -1) * points[m].x,
                            (j == 0 ? 1 : -1) * points[m].y, (k == 0 ? 1 : -1) * points[m].z);
                    }
                }
            }
        }
    }
}

float PerlinNoise::generator(float x, float y)
{
    int bx0, bx1, by0, by1, b00, b10, b01, b11;
    float rx0, rx1, ry0, ry1, sx, sy, a, b, t, u, v;
    Grad q;
    int i, j;

    t = x + 256;
    bx0 = ((int) t) & 255;
    bx1 = (bx0 + 1) & 255;
    rx0 = t - (int) t;
    rx1 = rx0 - 1;

    t = y + 256;
    by0 = ((int) t) & 255;
    by1 = (by0 + 1) & 255;
    ry0 = t - (int) t;
    ry1 = ry0 - 1;

    i = p[bx0];
    j = p[bx1];

    b00 = p[i + by0];
    b10 = p[j + by0];
    b01 = p[i + by1];
    b11 = p[j + by1];

    sx = s_curve(rx0);
    sy = s_curve(ry0);

    q = g2[b00];
    u = rx0 * q.x + ry0 * q.y;
    q = g2[b10];
    v = rx1 * q.x + ry0 * q.y;
    a = lerp(sx, u, v);

    q = g2[b01];
    u = rx0 * q.x + ry1 * q.y;
    q = g2[b11];
    v = rx1 * q.x + ry1 * q.y;
    b = lerp(sx, u, v);

    return lerp(sy, a, b);
}

// ----- //

SimplexNoise::SimplexNoise(uint64_t s) : NoiseGenerator(s)
{
    for (int i = 0; i < perm_table.size(); i++) {
        perm_table[i] = rand_byte(mt);
        perm_mod12[i] = perm_table[i] % 12;
    }
}

// 2D Simplex Noise algorithm adapted from lines 81-132 of:
// https://github.com/SRombauts/SimplexNoise/blob/master/references/SimplexNoise.java
// Based on example code by Stefan Gustavson (stegu@itn.liu.se) and Peter Eastman
// (peastman@drizzle.stanford.edu).
float SimplexNoise::generator(float x, float y)
{
    float n0, n1, n2;
    float s = (x + y) * skewing_factor;
    int i = (int) (x + s);
    int j = (int) (y + s);
    float t = (i + j) * unskewing_factor;
    float x_origin = i - t;
    float y_origin = j - t;
    float x0 = x - x_origin;
    float y0 = y - y_origin;

    int i1, j1;
    if (x0 > y0) {
        i1 = 1;
        j1 = 0;
    } else {
        i1 = 0;
        j1 = 1;
    }

    float x1 = x0 - i1 + unskewing_factor;
    float y1 = y0 - j1 + unskewing_factor;
    float x2 = x0 - 1.0 + 2.0 * unskewing_factor;
    float y2 = y0 - 1.0 + 2.0 * unskewing_factor;

    int ii = i & 255;
    int jj = j & 255;
    int gi0 = perm_mod12[ii + perm_table[jj]];
    int gi1 = perm_mod12[ii + i1 + perm_table[jj + j1]];
    int gi2 = perm_mod12[ii + 1 + perm_table[jj + 1]];

    float t0 = 0.5 - x0 * x0 - y0 * y0;
    if (t0 < 0)
        n0 = 0.0;
    else {
        t0 *= t0;
        n0 = t0 * t0 * dot(gradients[gi0], Grad(x0, y0, 0));
    }
    float t1 = 0.5 - x1 * x1 - y1 * y1;
    if (t1 < 0)
        n1 = 0.0;
    else {
        t1 *= t1;
        n1 = t1 * t1 * dot(gradients[gi1], Grad(x1, y1, 0));
    }
    float t2 = 0.5 - x2*x2-y2*y2;
    if (t2 < 0)
        n2 = 0.0;
    else {
        t2 *= t2;
        n2 = t2 * t2 * dot(gradients[gi2], Grad(x2, y2, 0));
    }

    return (n0 + n1 + n2);
}
