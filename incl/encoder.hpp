#ifndef _PNG_ENCODER_HPP_
#define _PNG_ENCODER_HPP_

#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>

#define DEFLATE_BUF_SZ 65536

// PNG Encoder
//
// Encoder that writes PNG files. Reference:
// https://www.w3.org/TR/2003/REC-PNG-20031110/
namespace PNG {
  enum class ColorType {
    Grayscale = 0,
    RGBTriple = 2,
    PLTEIndex = 3,
    GrayscaleAlpha = 4,
    RGBTripleAlpha = 6
  };

  struct PngChunk {
    uint32_t length;
    uint32_t crc;
    std::string type_code;
    std::vector<uint8_t> data;

    PngChunk(std::string tc, std::vector<uint8_t> data);
  };

  class PngEncoder {
  private:
    const std::array<uint8_t, 8> png_sig = {137, 80, 78, 71, 13, 10, 26, 10};
    std::fstream fs;
    bool fsig;
    bool ihdr;
    bool iend;
    uint32_t chunks;
    uint32_t img_width;
    uint32_t img_height;
    uint8_t bit_depth;
  public:
    PngEncoder();
    ~PngEncoder();
    void open(std::string fp);
    void close();

    void write(PngChunk& chunk);
    void operator<<(PngChunk& chunk);

    void write_ihdr(uint32_t width, uint32_t height, uint8_t bit_depth,
        ColorType color_type);
    void write_idat(const std::vector<uint8_t>& data);
    void write_iend();
  };
};

template <typename T>
T swap_endianness(const T& x)
{
    uint8_t *byte = (uint8_t*)&x + sizeof(T);
    uint8_t byte_array[sizeof(T)];
    for (std::size_t i = 0; i < sizeof(T); i++)
        byte_array[i] = *(--byte);
    return *reinterpret_cast<T*>(byte_array);
}

#endif // _PNG_ENCODER_HPP_
