#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "encoder.hpp"

// PNG encoder.

struct crc_table {
private:
  std::array<uint32_t, 256> arr;
public:
  constexpr crc_table() : arr()
  {
    for (int n = 0; n < 256; n++)
    {
      uint32_t c = (uint32_t) n;
      for (int k = 0; k < 8; k++)
      {
        if (c & 1)
          c = 0xEDB88320L ^ (c >> 1);
        else
          c >>= 1;
      }
      arr[n] = c;
    }
  }
  const uint32_t &operator[](std::size_t i) const { return arr[i]; }
};

constexpr auto crctable = crc_table();

uint32_t update_crc(uint32_t crc, const std::vector<uint8_t> &buf)
{
  uint32_t c = crc;

  for (int n = 0; n < buf.size(); n++) {
    c = crctable[(c ^ buf[n]) & 0xff] ^ (c >> 8);
  }

  return c;
}

uint32_t create_crc(const std::vector<uint8_t> &buf)
{
  return update_crc(0xFFFFFFFFL, buf) ^ 0xFFFFFFFFL;
}

uint32_t lsbtomsb(const uint32_t& val)
{
  return ((val >> 24) & 0xFF) | ((val >> 8) & 0xFF00) | ((val << 8) & 0xFF0000)
    | ((val << 24) & 0xFF000000);
}

// https://zlib.net/zlib_how.html
std::vector<uint8_t> compress_vector(const std::vector<uint8_t>& vec,
  int compressionlevel = Z_BEST_COMPRESSION)
{
    z_stream zs;
    memset(&zs, 0, sizeof(zs));

    if (deflateInit(&zs, compressionlevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    int ret;
    std::array<uint8_t, DEFLATE_BUF_SZ> buf;
    std::vector<uint8_t> out;

    zs.next_in = (Bytef*)vec.data();
    zs.avail_in = vec.size();
    // zs.next_out = (Bytef*)buf.data();
    // zs.avail_out = buf.size();

    do {
        zs.next_out = (Bytef*)(buf.data());
        zs.avail_out = buf.size();

        ret = deflate(&zs, Z_FINISH);

        std::copy(buf.begin(), buf.begin() + zs.total_out, std::back_inserter(out));
    } while (ret == Z_OK);
    deflateEnd(&zs);
    return out;
}

// ------------------------------ //
using namespace PNG;

PngChunk::PngChunk(std::string tc, std::vector<uint8_t> vec) : type_code(tc)
{
  length = vec.size();
  uint32_t msb_len = lsbtomsb(length);
  std::copy((uint8_t*)&msb_len, (uint8_t*)&msb_len + sizeof(msb_len), std::back_inserter(data));
  std::copy(tc.begin(), tc.end(), std::back_inserter(data));
  std::copy(vec.begin(), vec.end(), std::back_inserter(data));

  std::vector<uint8_t> crc_data;
  std::copy(data.begin() + 4, data.end(), std::back_inserter(crc_data));

  crc = create_crc(crc_data);
  uint32_t msb_crc = lsbtomsb(crc);

  std::copy((uint8_t*)&msb_crc, (uint8_t*)&msb_crc + sizeof(msb_crc), std::back_inserter(data));
}

PngEncoder::PngEncoder()
{

}

PngEncoder::~PngEncoder()
{
  this->close();
}

void PngEncoder::open(std::string fp)
{
  fs.open(fp, std::ios::out | std::ios::trunc | std::ios::binary);
  fsig = false;
  ihdr = false;
  iend = false;
  chunks = 0;
  img_width = 0;
  img_height = 0;
  bit_depth = 0;
  fs.write((char*)(png_sig.begin()), 8);
  fsig = true;

  fs.flush();
}

void PngEncoder::close()
{
  if (fs.is_open()) {
    if (!iend) this->write_iend();
    fs.close();
  }
}

void PngEncoder::write(PngChunk& chunk)
{
  if (!fs.is_open()) return;

  fs.write((char*)&chunk.data[0], chunk.data.size());

  if (chunk.type_code == "IHDR")
    ihdr = true;
  if (chunk.type_code == "IEND")
    iend = true;
}

void PngEncoder::operator<<(PngChunk& chunk)
{
  this->write(chunk);
}

void PngEncoder::write_ihdr(uint32_t width, uint32_t height, uint8_t bitd,
  ColorType color_type)
{
  std::vector<uint8_t> data;
  img_width = width;
  img_height = height;
  bit_depth = bitd;
  uint32_t msb_width = lsbtomsb(width);
  uint32_t msb_height = lsbtomsb(height);
  std::copy((uint8_t*)&msb_width, (uint8_t*)&msb_width + sizeof(msb_width),
    std::back_inserter(data));
  std::copy((uint8_t*)&msb_height, (uint8_t*)&msb_height + sizeof(msb_height),
    std::back_inserter(data));
  data.push_back(bit_depth);
  data.push_back((uint8_t)color_type);

  for (int i = 0; i < 3; i++)
    data.push_back(0);

  auto chunk = PngChunk("IHDR", data);
  this->write(chunk);
}

void PngEncoder::write_idat(const std::vector<uint8_t>& data)
{
  std::vector<uint8_t> filtered_data;
  for (int i = 0; i < data.size(); i++)
  {
    if (i % (img_width * (bit_depth / 8)) == 0)
      filtered_data.push_back(0);
    filtered_data.push_back(data[i]);
  }
  auto compressed_data = compress_vector(filtered_data);
  auto chunk = PngChunk("IDAT", compressed_data);
  this->write(chunk);
}

void PngEncoder::write_iend()
{
  auto chunk = PngChunk("IEND", {});
  this->write(chunk);
}
