#ifndef GELEX_DATA_DECODE_LUT_H
#define GELEX_DATA_DECODE_LUT_H

#include <array>
#include <cstdint>

namespace gelex
{

constexpr std::array<double, 4> generate_lut_entry(uint8_t byte, bool reverse)
{
    std::array<double, 4> entry{};

    // 00 -> Homozygote 1
    // 01 -> Missing
    // 10 -> Heterozygote
    // 11 -> Homozygote 2

    // IsReverse=false:
    // 00(0)->2.0, 01(1)->0.0, 10(2)->1.0, 11(3)->0.0 we set the missing to
    // major allele for simplicity following PLINK convention(major allele is
    // A2).
    constexpr std::array<double, 4> std_map = {2.0, 0.0, 1.0, 0.0};

    // IsReverse=true, swap A1 A2:
    // 00->0.0, 01->NaN, 10->1.0, 11->2.0
    constexpr std::array<double, 4> rev_map = {0.0, 0.0, 1.0, 2.0};

    const std::array<double, 4> map = reverse ? rev_map : std_map;

    for (int i = 0; i < 4; ++i)
    {
        entry[i] = map[(byte >> (2 * i)) & 0x03];
    }
    return entry;
}

template <bool Reverse>
consteval std::array<std::array<double, 4>, 256> generate_full_lut()
{
    std::array<std::array<double, 4>, 256> table{};
    for (std::size_t i = 0; i < 256; ++i)
    {
        table[i] = generate_lut_entry(static_cast<uint8_t>(i), Reverse);
    }
    return table;
}

using LutEntry = std::array<double, 4>;
static constexpr std::size_t kLutSize = 256;

// Forward lookup table (standard mapping)
inline const std::array<LutEntry, kLutSize> kDecodeLut
    = generate_full_lut<false>();

}  // namespace gelex

#endif  // GELEX_DATA_DECODE_LUT_H
