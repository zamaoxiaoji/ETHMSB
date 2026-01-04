#pragma once

// Experiment-only MSB6 implementation.
//
// Hard constraints:
// - Do NOT modify any existing HE3DB/TFHEpp sources.
// - Implement by copying/renaming logic in experiments/ only.
//
// Idea:
// - HE3DB 的 ExtractMSB5/MSBGateBootstrapping 设计上对应 plain_bits=5 的整数编码：
//     Δ = 2^(W-plain_bits)，offset=Δ/2
// - 当输入实际是 plain_bits=6 时，原实现仍用固定 offset=2^(W-6)，这会变成 Δ 而非 Δ/2，
//   在 worst-case（m=2^t-1）时容易把相位推到回绕临界点，从而出现大量错误。
// - 本文件提供一个 MSB6 版本：offset = 2^(W-plain_bits-1) = Δ/2（plain_bits=6）。

#include "HEDB/comparison/tfhepp_utils.h"
#include "HEDB/utils/types.h"

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace HE3DBExperiment {

// Copied from: src/HEDB/comparison/tfhepp_utils.cpp (MSBGateBootstrapping)
// Change: offset depends on plain_bits: offset = 2^(W - plain_bits - 1) = Δ/2
inline void MSBGateBootstrappingDeltaHalf(HEDB::TLWELvl1 &res,
                                          const HEDB::TLWELvl1 &tlwe,
                                          const HEDB::TFHEEvalKey &ek,
                                          bool result_type,
                                          uint32_t plain_bits)
{
    using P = TFHEpp::lvl1param;
    using T = typename P::T;
    constexpr uint32_t W = std::numeric_limits<T>::digits;  // e.g. 32
    if (plain_bits == 0 || plain_bits + 1 >= W) {
        throw std::runtime_error(
            "plain_bits 非法：需要满足 1 <= plain_bits <= W-2");
    }

    uint32_t μ = 1U << 29;
    if (IS_ARITHMETIC(result_type)) μ = μ << 1;

    const uint64_t offset = 1ULL << (W - plain_bits - 1);  // Δ/2

    TFHEpp::TLWE<P> tlweoffset = tlwe;
    tlweoffset[P::k * P::n] += static_cast<T>(offset);

    TFHEpp::TLWE<TFHEpp::lvl0param> tlwelvl0;
    TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(tlwelvl0, tlweoffset,
                                                  *ek.iksklvl10);
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<TFHEpp::lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<P>(μ));

    if (IS_ARITHMETIC(result_type)) res[P::k * P::n] += static_cast<T>(μ);
}

// MSB6: for a ciphertext encoding a 6-bit integer (plain_bits=6).
inline void ExtractMSB(HEDB::TLWELvl1 &res, const HEDB::TLWELvl1 &tlwe,
                        const HEDB::TFHEEvalKey &ek, bool result_type,)
{
    MSBGateBootstrappingDeltaHalf(res, tlwe, ek, result_type, /*plain_bits=*/);
}

}  // namespace HE3DBExperiment

