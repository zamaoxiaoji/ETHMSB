#pragma once

// experiments/he3db_msb6.hpp
//
// 实验版“通用 MSB 门”（只改 offset=Δ/2），并提供可控 trace。
//
// 背景：
// - HE3DB 的 MSBGateBootstrapping/ExtractMSB5 内部固定加 offset=2^(W-6)=1/64；
//   这等价于假设输入明文位宽 plain_bits=5 时的 Δ/2。
// - 当输入位宽 plain_bits>5（例如 6 位、8 位…），固定 offset 会“侵入 message 区”、
//   在 worst-case m=2^plain_bits-1 时把相位推到 2^W 回绕边界附近甚至越界，导致 MSB 失败。
// - 自然修正：对任意 plain_bits=t，使用 offset=Δ/2=2^(W-t-1)；LUT 仍复用 TFHEpp 的
//   μ_polygen（即同一个 MSB 门 LUT）。
//
// 约束（硬性）：
// - 不修改任何 HE3DB/TFHEpp 既有源码（src/thirdparty/include/...）。
// - 仅在 experiments/ 下“复制一份实现到新命名空间”并可插桩打印。
// - 不在大循环里无条件打印；trace 由调用者在失败样本复跑时开启。

#include "experiments/he3db_phase_trace.hpp"
#include "experiments/tfhepp_instrumented_pbs.hpp"

#include "HEDB/comparison/tfhepp_utils.h"
#include "HEDB/utils/types.h"

#include <tfhe++.hpp>

#include <cstdint>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <stdexcept>

namespace HE3DBExperiment {

template <class TorusT>
constexpr uint32_t torus_bits_v =
    HE3DBTrace::torus_bits_v<std::make_unsigned_t<TorusT>>;

template <class TorusT>
uint64_t delta_half_offset_u64(uint32_t plain_bits)
{
    constexpr uint32_t W = torus_bits_v<TorusT>;
    if (plain_bits == 0 || plain_bits >= W) {
        throw std::runtime_error("plain_bits 非法：需要满足 1 <= plain_bits <= W-1");
    }
    const uint32_t shift = W - plain_bits - 1;
    return (uint64_t{1} << shift);
}

template <class PhaseT>
void print_phase_layout_safe(std::ostream &os, std::string_view tag, PhaseT phase_u,
                             uint64_t m, uint32_t plain_bits, uint32_t err_bits)
{
    using U = std::make_unsigned_t<PhaseT>;
    constexpr uint32_t W = torus_bits_v<U>;
    if (plain_bits + err_bits > W) {
        const uint64_t phase64 = static_cast<uint64_t>(static_cast<U>(phase_u));
        os << "\n【" << tag << "】\n";
        os << "W=" << W << "，plain_bits=" << plain_bits << "，err_bits=" << err_bits
           << "（plain_bits+err_bits>W，无法按三段切分，改为 W 位二进制输出）\n";
        os << "phase(十进制)=" << phase64
           << "，phase(W位)=" << HE3DBTrace::bits_fixed_u64(phase64, W) << "\n";
        return;
    }
    HE3DBTrace::print_phase_layout(os, tag, phase_u, m, plain_bits, err_bits, W);
}

// ============================================================
// 复制版 MSB 门（只改 offset=Δ/2），并在 trace_enable 时打印：
// - 输入 phase
// - 加 offset 后 phase
// - IKS 后 lvl0 phase（这是进入 BlindRotate 之前的 TLWE）
// - BlindRotate 量化细节（qb/bbar/q_phase_direct/q_phase_parts）
// - SampleExtract 后输出 bit TLWE 的 phase(raw) + 解密 bit
// ============================================================

// Copied from: src/HEDB/comparison/tfhepp_utils.cpp (MSBGateBootstrapping)
// Change:
//   - offset: 固定 2^(W-6) -> 2^(W-plain_bits-1)（Δ/2）
//   - bootstrapping: 使用 experiments/tfhepp_instrumented_pbs.hpp 的复制版接口
//     以便在失败样本复跑时打印 BlindRotate 的量化细节
inline void MSBGateBootstrappingDeltaHalfTrace(
    HEDB::TLWELvl1 &res, const HEDB::TLWELvl1 &tlwe, const HEDB::TFHEEvalKey &ek,
    bool result_type, uint32_t plain_bits, uint32_t err_bits, uint64_t m,
    uint64_t seed, const HEDB::TFHESecretKey *sk_dbg, bool trace_enable,
    uint32_t max_nonzero_a_print)
{
    using InP = TFHEpp::lvl1param;
    using OutP = TFHEpp::lvl1param;
    using InT = typename InP::T;
    using OutT = typename OutP::T;
    constexpr uint32_t W = torus_bits_v<InT>;

    const uint32_t mu_base = 1U << 29;
    uint32_t mu = mu_base;
    if (IS_ARITHMETIC(result_type)) mu = mu << 1;

    const uint64_t offset_u = delta_half_offset_u64<InT>(plain_bits);

    if (trace_enable && sk_dbg != nullptr) {
        std::cout << "\n（通用MSB门 Δ/2）参数：W=" << W << "，plain_bits=" << plain_bits
                  << "，Δ=2^" << (W - plain_bits) << "，offset=Δ/2=2^"
                  << (W - plain_bits - 1) << "，seed=" << seed << "，m=" << m
                  << "\n";
        const auto phase_in =
            HE3DBTrace::tlwe_phase<InP>(tlwe, sk_dbg->key.lvl1);
        print_phase_layout_safe(std::cout, "通用MSB门：输入phase", phase_in, m,
                                plain_bits, err_bits);
    }

    TFHEpp::TLWE<InP> tlweoffset = tlwe;
    tlweoffset[InP::k * InP::n] += static_cast<InT>(offset_u);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_off =
            HE3DBTrace::tlwe_phase<InP>(tlweoffset, sk_dbg->key.lvl1);
        print_phase_layout_safe(std::cout, "通用MSB门：加offset后phase", phase_off,
                                m, plain_bits, err_bits);
    }

    TFHEpp::TLWE<TFHEpp::lvl0param> tlwelvl0;
    TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(tlwelvl0, tlweoffset,
                                                  *ek.iksklvl10);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_lvl0 =
            HE3DBTrace::tlwe_phase<TFHEpp::lvl0param>(tlwelvl0, sk_dbg->key.lvl0);
        print_phase_layout_safe(std::cout, "通用MSB门：IKS后(lvl0) phase（进入BlindRotate前）",
                                phase_lvl0, m, plain_bits, err_bits);
    }

    TFHEppInstrumented::TraceConfig cfg;
    TFHEppInstrumented::TraceRecord rec;
    TFHEppInstrumented::TraceConfig *cfg_p = nullptr;
    TFHEppInstrumented::TraceRecord *rec_p = nullptr;
    const TFHEpp::Key<TFHEpp::lvl0param> *sk0_p = nullptr;
    if (trace_enable && sk_dbg != nullptr) {
        cfg.enable = true;
        cfg.max_nonzero_a_print = max_nonzero_a_print;
        cfg.os = &std::cout;
        rec.t = plain_bits;
        rec.m = m;
        rec.seed = seed;
        cfg_p = &cfg;
        rec_p = &rec;
        sk0_p = &sk_dbg->key.lvl0;
    }

    TFHEppInstrumented::GateBootstrappingTLWE2TLWEFFTTrace<TFHEpp::lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<OutP>(mu), sk0_p,
        cfg_p, rec_p);

    if (trace_enable && sk_dbg != nullptr) {
        // res 是 bit 密文（语义已不是“整数 m”），这里只打印 raw phase + 解密结果。
        const auto phase_out =
            HE3DBTrace::tlwe_phase<OutP>(res, sk_dbg->key.lvl1);
        const uint64_t phase64 = static_cast<uint64_t>(static_cast<OutT>(phase_out));
        std::cout << "\n【通用MSB门：SampleExtract后输出bit phase(raw)】\n";
        std::cout << "phase(十进制)=" << phase64 << "，phase(W位)="
                  << HE3DBTrace::bits_fixed_u64(phase64, torus_bits_v<OutT>)
                  << "\n";
        const bool got = TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(res, sk_dbg->key.lvl1);
        std::cout << "解密得到 bit=" << (got ? 1 : 0)
                  << "（result_type=" << (IS_ARITHMETIC(result_type) ? "ARITHMETIC" : "LOGIC") << "）\n";
    }

    if (IS_ARITHMETIC(result_type)) {
        res[OutP::k * OutP::n] += static_cast<OutT>(mu);
    }
}

inline void MSBGateBootstrappingDeltaHalfTrace(
    HEDB::TLWELvl1 &res, const HEDB::TLWELvl2 &tlwe, const HEDB::TFHEEvalKey &ek,
    bool result_type, uint32_t plain_bits, uint32_t err_bits, uint64_t m,
    uint64_t seed, const HEDB::TFHESecretKey *sk_dbg, bool trace_enable,
    uint32_t max_nonzero_a_print)
{
    using InP = TFHEpp::lvl2param;
    using OutP = TFHEpp::lvl1param;
    using InT = typename InP::T;
    using OutT = typename OutP::T;
    constexpr uint32_t W = torus_bits_v<InT>;

    const uint32_t mu_base = 1U << 29;
    uint32_t mu = mu_base;
    if (IS_ARITHMETIC(result_type)) mu = mu << 1;

    const uint64_t offset_u = delta_half_offset_u64<InT>(plain_bits);

    if (trace_enable && sk_dbg != nullptr) {
        std::cout << "\n（通用MSB门 Δ/2，输入lvl2->输出lvl1）参数：W=" << W
                  << "，plain_bits=" << plain_bits << "，offset=2^"
                  << (W - plain_bits - 1) << "，seed=" << seed << "，m=" << m
                  << "\n";
        const auto phase_in =
            HE3DBTrace::tlwe_phase<InP>(tlwe, sk_dbg->key.lvl2);
        print_phase_layout_safe(std::cout, "通用MSB门：输入phase(lvl2)", phase_in,
                                m, plain_bits, err_bits);
    }

    TFHEpp::TLWE<InP> tlweoffset = tlwe;
    tlweoffset[InP::k * InP::n] += static_cast<InT>(offset_u);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_off =
            HE3DBTrace::tlwe_phase<InP>(tlweoffset, sk_dbg->key.lvl2);
        print_phase_layout_safe(std::cout, "通用MSB门：加offset后phase(lvl2)",
                                phase_off, m, plain_bits, err_bits);
    }

    TFHEpp::TLWE<TFHEpp::lvl0param> tlwelvl0;
    TFHEpp::IdentityKeySwitch<TFHEpp::lvl20param>(tlwelvl0, tlweoffset,
                                                  *ek.iksklvl20);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_lvl0 =
            HE3DBTrace::tlwe_phase<TFHEpp::lvl0param>(tlwelvl0, sk_dbg->key.lvl0);
        print_phase_layout_safe(std::cout, "通用MSB门：IKS后(lvl0) phase（进入BlindRotate前）",
                                phase_lvl0, m, plain_bits, err_bits);
    }

    TFHEppInstrumented::TraceConfig cfg;
    TFHEppInstrumented::TraceRecord rec;
    TFHEppInstrumented::TraceConfig *cfg_p = nullptr;
    TFHEppInstrumented::TraceRecord *rec_p = nullptr;
    const TFHEpp::Key<TFHEpp::lvl0param> *sk0_p = nullptr;
    if (trace_enable && sk_dbg != nullptr) {
        cfg.enable = true;
        cfg.max_nonzero_a_print = max_nonzero_a_print;
        cfg.os = &std::cout;
        rec.t = plain_bits;
        rec.m = m;
        rec.seed = seed;
        cfg_p = &cfg;
        rec_p = &rec;
        sk0_p = &sk_dbg->key.lvl0;
    }

    TFHEppInstrumented::GateBootstrappingTLWE2TLWEFFTTrace<TFHEpp::lvl01param>(
        res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<OutP>(mu), sk0_p,
        cfg_p, rec_p);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_out =
            HE3DBTrace::tlwe_phase<OutP>(res, sk_dbg->key.lvl1);
        const uint64_t phase64 = static_cast<uint64_t>(static_cast<OutT>(phase_out));
        std::cout << "\n【通用MSB门：SampleExtract后输出bit phase(raw)】\n";
        std::cout << "phase(十进制)=" << phase64 << "，phase(W位)="
                  << HE3DBTrace::bits_fixed_u64(phase64, torus_bits_v<OutT>)
                  << "\n";
        const bool got = TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(res, sk_dbg->key.lvl1);
        std::cout << "解密得到 bit=" << (got ? 1 : 0) << "\n";
    }

    if (IS_ARITHMETIC(result_type)) {
        res[OutP::k * OutP::n] += static_cast<OutT>(mu);
    }
}

inline void MSBGateBootstrappingDeltaHalfTrace(
    HEDB::TLWELvl2 &res, const HEDB::TLWELvl2 &tlwe, const HEDB::TFHEEvalKey &ek,
    bool result_type, uint32_t plain_bits, uint32_t err_bits, uint64_t m,
    uint64_t seed, const HEDB::TFHESecretKey *sk_dbg, bool trace_enable,
    uint32_t max_nonzero_a_print)
{
    using InP = TFHEpp::lvl2param;
    using OutP = TFHEpp::lvl2param;
    using InT = typename InP::T;
    using OutT = typename OutP::T;
    constexpr uint32_t W = torus_bits_v<InT>;

    uint64_t mu = 1ULL << 61;
    if (IS_ARITHMETIC(result_type)) mu = mu << 1;

    const uint64_t offset_u = delta_half_offset_u64<InT>(plain_bits);

    if (trace_enable && sk_dbg != nullptr) {
        std::cout << "\n（通用MSB门 Δ/2，输入lvl2->输出lvl2）参数：W=" << W
                  << "，plain_bits=" << plain_bits << "，offset=2^"
                  << (W - plain_bits - 1) << "，seed=" << seed << "，m=" << m
                  << "\n";
        const auto phase_in =
            HE3DBTrace::tlwe_phase<InP>(tlwe, sk_dbg->key.lvl2);
        print_phase_layout_safe(std::cout, "通用MSB门：输入phase(lvl2)", phase_in,
                                m, plain_bits, err_bits);
    }

    TFHEpp::TLWE<InP> tlweoffset = tlwe;
    tlweoffset[InP::k * InP::n] += static_cast<InT>(offset_u);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_off =
            HE3DBTrace::tlwe_phase<InP>(tlweoffset, sk_dbg->key.lvl2);
        print_phase_layout_safe(std::cout, "通用MSB门：加offset后phase(lvl2)",
                                phase_off, m, plain_bits, err_bits);
    }

    TFHEpp::TLWE<TFHEpp::lvl0param> tlwelvl0;
    TFHEpp::IdentityKeySwitch<TFHEpp::lvl20param>(tlwelvl0, tlweoffset,
                                                  *ek.iksklvl20);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_lvl0 =
            HE3DBTrace::tlwe_phase<TFHEpp::lvl0param>(tlwelvl0, sk_dbg->key.lvl0);
        print_phase_layout_safe(std::cout, "通用MSB门：IKS后(lvl0) phase（进入BlindRotate前）",
                                phase_lvl0, m, plain_bits, err_bits);
    }

    TFHEppInstrumented::TraceConfig cfg;
    TFHEppInstrumented::TraceRecord rec;
    TFHEppInstrumented::TraceConfig *cfg_p = nullptr;
    TFHEppInstrumented::TraceRecord *rec_p = nullptr;
    const TFHEpp::Key<TFHEpp::lvl0param> *sk0_p = nullptr;
    if (trace_enable && sk_dbg != nullptr) {
        cfg.enable = true;
        cfg.max_nonzero_a_print = max_nonzero_a_print;
        cfg.os = &std::cout;
        rec.t = plain_bits;
        rec.m = m;
        rec.seed = seed;
        cfg_p = &cfg;
        rec_p = &rec;
        sk0_p = &sk_dbg->key.lvl0;
    }

    TFHEppInstrumented::GateBootstrappingTLWE2TLWEFFTTrace<TFHEpp::lvl02param>(
        res, tlwelvl0, *ek.bkfftlvl02, TFHEpp::μ_polygen<OutP>(mu), sk0_p,
        cfg_p, rec_p);

    if (trace_enable && sk_dbg != nullptr) {
        const auto phase_out =
            HE3DBTrace::tlwe_phase<OutP>(res, sk_dbg->key.lvl2);
        const uint64_t phase64 = static_cast<uint64_t>(static_cast<OutT>(phase_out));
        std::cout << "\n【通用MSB门：SampleExtract后输出bit phase(raw)】\n";
        std::cout << "phase(十进制)=" << phase64 << "，phase(W位)="
                  << HE3DBTrace::bits_fixed_u64(phase64, torus_bits_v<OutT>)
                  << "\n";
        const bool got = TFHEpp::tlweSymDecrypt<HEDB::Lvl2>(res, sk_dbg->key.lvl2);
        std::cout << "解密得到 bit=" << (got ? 1 : 0) << "\n";
    }

    if (IS_ARITHMETIC(result_type)) {
        res[OutP::k * OutP::n] += static_cast<OutT>(mu);
    }
}

// 便捷封装：仅做 MSB 门（bit 输出）
inline void ExtractMSBDeltaHalf(HEDB::TLWELvl1 &res, const HEDB::TLWELvl1 &tlwe,
                                const HEDB::TFHEEvalKey &ek, bool result_type,
                                uint32_t plain_bits)
{
    MSBGateBootstrappingDeltaHalfTrace(res, tlwe, ek, result_type, plain_bits,
                                       /*err_bits=*/0, /*m=*/0, /*seed=*/0,
                                       /*sk_dbg=*/nullptr,
                                       /*trace_enable=*/false,
                                       /*max_nonzero_a_print=*/0);
}

// MSB6: 兼容旧实验程序接口（plain_bits 固定为 6）。
inline void ExtractMSB6(HEDB::TLWELvl1 &res, const HEDB::TLWELvl1 &tlwe,
                        const HEDB::TFHEEvalKey &ek, bool result_type)
{
    ExtractMSBDeltaHalf(res, tlwe, ek, result_type, /*plain_bits=*/6);
}

}  // namespace HE3DBExperiment
