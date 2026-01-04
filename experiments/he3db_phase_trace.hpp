#pragma once

// HE3DB phase 计算与三段格式打印工具（仅供 experiments 使用）
// 约束：不修改 HE3DB/TFHEpp 原有源码；只在 experiments/ 新文件中实现。

#include <tfhe++.hpp>

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>

namespace HE3DBTrace {

template <class T>
constexpr uint32_t torus_bits_v =
    static_cast<uint32_t>(std::numeric_limits<std::make_unsigned_t<T>>::digits);

template <class T>
using WideT =
    std::conditional_t<(torus_bits_v<T> <= 32), uint64_t, unsigned __int128>;

// =========================
// A1) 固定宽度二进制字符串
// =========================
inline std::string bits_fixed_u64(uint64_t x, uint32_t width)
{
    std::string s;
    s.reserve(width);
    for (uint32_t i = 0; i < width; ++i) {
        const uint32_t bitpos = width - 1 - i;
        s.push_back(((x >> bitpos) & 1ULL) ? '1' : '0');
    }
    return s;
}

inline uint64_t mask_u64(uint32_t width)
{
    if (width == 0) return 0;
    if (width >= 64) return ~uint64_t{0};
    return (uint64_t{1} << width) - 1;
}

// =========================
// A2) 三段格式：m区|Δ间隔区|误差区
// =========================
inline std::string format_phase_triplet(uint64_t phase_u, uint32_t W,
                                        uint32_t plain_bits,
                                        uint32_t err_bits)
{
    const uint32_t gap_bits = W - plain_bits - err_bits;
    const uint64_t m_part = phase_u >> (gap_bits + err_bits);
    const uint64_t gap_part = (phase_u >> err_bits) & mask_u64(gap_bits);
    const uint64_t err_part = phase_u & mask_u64(err_bits);
    return bits_fixed_u64(m_part, plain_bits) + " | " +
           bits_fixed_u64(gap_part, gap_bits) + " | " +
           bits_fixed_u64(err_part, err_bits);
}

template <class T>
std::string format_phase_triplet(T phase_u, uint32_t plain_bits,
                                 uint32_t err_bits)
{
    using U = std::make_unsigned_t<T>;
    constexpr uint32_t W = torus_bits_v<U>;
    return format_phase_triplet(static_cast<uint64_t>(static_cast<U>(phase_u)), W,
                                plain_bits, err_bits);
}

// =========================
// A3) 计算 TLWE 的 phase：b - Σ(a_i*s_i) (mod 2^W)
// =========================
template <class P>
typename P::T tlwe_phase(const TFHEpp::TLWE<P> &ct, const TFHEpp::Key<P> &sk)
{
    using T = typename P::T;
    using Wide = WideT<T>;
    constexpr uint32_t W = torus_bits_v<T>;
    const Wide mask = (Wide{1} << W) - 1;

    Wide acc = static_cast<Wide>(ct[P::k * P::n]) & mask;
    for (size_t i = 0; i < P::k * P::n; ++i) {
        acc = (acc - (static_cast<Wide>(ct[i]) * static_cast<Wide>(sk[i]))) &
              mask;
    }
    return static_cast<T>(acc);
}

template <class T>
uint64_t ideal_u64(uint64_t m, uint32_t plain_bits)
{
    constexpr uint32_t W = torus_bits_v<T>;
    const uint32_t shift = W - plain_bits;
    if constexpr (W == 64) {
        return (m << shift);
    }
    else {
        return (m << shift) & mask_u64(W);
    }
}

// =========================
// A4) 专用打印：三段布局 + ideal + error
// =========================
template <class T>
void print_phase_layout(std::ostream &os, std::string_view tag, T phase_u,
                        uint64_t m, uint32_t plain_bits, uint32_t err_bits,
                        uint32_t W)
{
    using U = std::make_unsigned_t<T>;
    const uint32_t gap_bits = W - plain_bits - err_bits;
    const uint64_t m_mask = mask_u64(plain_bits);
    const uint64_t m_bin_u = m & m_mask;

    const uint64_t phase64 = static_cast<uint64_t>(static_cast<U>(phase_u));
    const uint64_t ideal64 = ideal_u64<U>(m_bin_u, plain_bits);
    const uint64_t error_mod =
        (W == 64) ? (phase64 - ideal64) : ((phase64 - ideal64) & mask_u64(W));
    const uint64_t error_low = error_mod & mask_u64(err_bits);

    const uint64_t m_top = phase64 >> (W - plain_bits);
    const uint64_t m_top_msb = (plain_bits == 0) ? 0 : ((m_top >> (plain_bits - 1)) & 1ULL);
    const uint64_t m_in_msb = (plain_bits == 0) ? 0 : ((m_bin_u >> (plain_bits - 1)) & 1ULL);

    os << "\n【" << tag << "】\n";
    os << "W=" << W << "，plain_bits=" << plain_bits << "，err_bits=" << err_bits
       << "，gap_bits=" << gap_bits << "\n";
    os << "m(十进制)=" << m << "，m(二进制t位)=" << bits_fixed_u64(m_bin_u, plain_bits)
       << "\n";
    os << "Δ间隔位数=gap_bits=" << gap_bits << "，Δ间隔=2^" << gap_bits << "\n";

    os << "ideal(十进制)=" << ideal64
       << "，ideal(三段)=" << format_phase_triplet(ideal64, W, plain_bits, err_bits)
       << "\n";
    os << "phase(十进制)=" << phase64
       << "，phase(三段)=" << format_phase_triplet(phase64, W, plain_bits, err_bits)
       << "\n";

    const unsigned __int128 mod128 = (static_cast<unsigned __int128>(1) << W);
    const unsigned __int128 half128 =
        (static_cast<unsigned __int128>(1) << (W - 1));
    const __int128 centered128 =
        (static_cast<unsigned __int128>(error_mod) >= half128)
            ? (static_cast<__int128>(static_cast<unsigned __int128>(error_mod)) -
               static_cast<__int128>(mod128))
            : static_cast<__int128>(static_cast<unsigned __int128>(error_mod));

    os << "error_mod=(phase-ideal) mod 2^W = " << error_mod
       << "，二进制(W位)=" << bits_fixed_u64(error_mod, W) << "\n";
    os << "噪声e_centered(带符号)=" << static_cast<int64_t>(centered128)
       << "（即 phase = ideal + e_centered mod 2^W）\n";
    os << "error_mod(三段)="
       << format_phase_triplet(error_mod, W, plain_bits, err_bits) << "\n";
    os << "error_low(低err_bits)=" << bits_fixed_u64(error_low, err_bits)
       << "（十进制=" << error_low << "）\n";

    os << "从phase解码 m_top=phase>>(W-plain_bits) = " << m_top
       << "，二进制=" << bits_fixed_u64(m_top, plain_bits) << "\n";
    os << "m区MSB：输入m=" << m_in_msb << "，解码m_top=" << m_top_msb << "\n";
    if (m_top_msb != m_in_msb) {
        os << "提示：m区最高位发生反转（可能与连续进位/回绕有关）\n";
    }
}

}  // namespace HE3DBTrace
