#include <array>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <zlib.h>

#include "encoder.hpp"

// --------------------------- Utility Functions --------------------------- //
// crc_table stores the CRC code for every value of an 8-bit number into a table at compile time
// using constexpr. The table is used when calculating the CRC of a PNG chunk as to avoid
// recalculating the CRC code for each byte thousands of times.
struct crc_table {
private:
    uint32_t polynomial = 0xEDB88320L; // magic number which yields the results of
    uint32_t tableLength = 256;        // polynomial division when XOR'd
    std::array<uint32_t, 256> tableArray;
public:
    constexpr crc_table() : tableArray() {
        for (uint32_t i = 0; i < tableLength; i++) {
            uint32_t coefficient = i;
            for (int j = 0; j < 8; j++) {
                coefficient = (coefficient & 1) ? (polynomial ^ (coefficient >> 1)) : (coefficient >> 1);
            }
            tableArray[i] = coefficient;
        }
    }
    const uint32_t &operator[](std::size_t ind) const { return tableArray[ind]; }
};

constexpr auto crcTable = crc_table();

uint32_t update_crc(uint32_t crcCode, const std::vector<uint8_t> &vectorBytes)
{
    for (int i = 0; i < buffer.size(); i++) {
        crcCode = crcTable[(crcCode ^ vectorBytes[i]) & 0xFF] ^ (crcCode >> 8);
    }

    return crcCode;
}

uint32_t create_crc(const std::vector<uint8_t> &vectorBytes)
{
    uint32_t startingValue = 0xFFFFFFFF;
    return update_crc(startingValue, vectorBytes) ^ startingValue;
}

// Reference: https://zlib.net/zlib_how.html
std::vector<uint8_t> compress_vector(const std::vector<uint8_t>& vectorBytes,
    int compressionLevel = Z_BEST_COMPRESSION)
{
    using zLibBytePtr = Bytef*;
    z_stream zLibStream;
    std::memset(&zLibStream, 0, sizeof(zLibStream));

    if (deflateInit(&zLibStream, compressionLevel) != Z_OK)
        throw(std::runtime_error("deflateInit failed while compressing."));

    int retValue;
    uint32_t prevOutputLength = 0;
    std::array<uint8_t, DEFLATE_BUF_SZ> buffer;
    std::vector<uint8_t> deflatedBuffer;

    zLibStream.next_in = (zLibBytePtr)(vectorBytes.data());
    zLibStream.avail_in = vectorBytes.size();

    do {
        zLibStream.next_out = (zLibBytePtr)(buffer.data());
        zLibStream.avail_out = buffer.size();

        retValue = deflate(&zLibStream, Z_FINISH);
        prevOutputLength = zLibStream.total_out - prevOutputLength;

        std::copy(buffer.begin(), buffer.begin() + prevOutputLength,
            std::back_inserter(deflatedBuffer));
        prevOutputLength = zLibStream.total_out;
    } while (retValue == Z_OK);

    deflateEnd(&zLibStream);
    return deflatedBuffer;
}

// --------------------------- Utility Functions --------------------------- //

// --------------------------- Class Definitions --------------------------- //
using namespace PNG;

PngChunk::PngChunk(std::string tc, std::vector<uint8_t> vec) : type_code(tc)
{
    length = vec.size();
    uint32_t msb_len = swap_endianness(length);
    std::copy((uint8_t*)&msb_len, (uint8_t*)&msb_len + sizeof(msb_len), std::back_inserter(data));
    std::copy(tc.begin(), tc.end(), std::back_inserter(data));
    std::copy(vec.begin(), vec.end(), std::back_inserter(data));

    std::vector<uint8_t> crc_data;
    std::copy(data.begin() + 4, data.end(), std::back_inserter(crc_data));

    crc = create_crc(crc_data);
    uint32_t msb_crc = swap_endianness(crc);

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
    uint32_t msb_width = swap_endianness(width);
    uint32_t msb_height = swap_endianness(height);
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
    for (int i = 0; i < data.size(); i++) {
        if (i % (img_width * (bit_depth / 8)) == 0) {
            filtered_data.push_back(0);
        }
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
