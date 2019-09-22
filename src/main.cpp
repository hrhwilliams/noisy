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
    uint64_t seed;
    std::string generator_name;
    bool generator_input_flag = false;
    bool need_generator = true;
    bool seed_input_flag = false;
    bool need_seed = true;

    for (int i = 1; i < argc; i++) {
        std::string str = std::string(argv[i]);
        if (str == "--generator" || str == "-g") {
            generator_input_flag = true;
            continue;
        }

        if (str == "--seed" || str == "-s") {
            seed_input_flag = true;
            continue;
        }

        if (generator_input_flag) {
            generator_name = str;
            need_generator = false;
            generator_input_flag = false;
        }

        if (seed_input_flag) {
            seed = std::atol(argv[i]);
            need_seed = false;
            seed_input_flag = false;
        }
    }

    if (need_generator) {
        std::cerr << argv[0] << ": need generator." << std::endl;
        return -1;
    }

    std::stringstream seed_str;

    if (need_seed) {
        // Generate a 64-bit seed
        seed = rd();
        seed <<= 32;
        seed += rd();
    } else {
        seed = std::atol(argv[1]);
    }
    seed_str << seed << ".png";

    Noise::NoiseGenerator *noise;
    auto *perlin = new Noise::PerlinNoise(seed);
    auto *simplex = new Noise::SimplexNoise(seed);

    if (generator_name == "perlin") {
        noise = perlin;
    } else if (generator_name == "simplex") {
        noise = simplex;
    } else {
        std::cerr << argv[0] << ": no match." << std::endl;
        return -1;
    }


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

    auto image = noise->generate(x_dim, y_dim, scale, octaves, persistence,
        lacunarity);

    encoder.open(seed_str.str());
    encoder.write_ihdr(x_dim, y_dim, 8, PNG::ColorType::Grayscale);
    encoder.write_idat(convert_vec(image));
    encoder.write_iend();
    encoder.close();

    std::cout << "Wrote out " << seed_str.str() << std::endl;

    return 0;
}
