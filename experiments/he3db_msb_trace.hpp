#pragma once

// HE3DB MSB(5) trace helpers.
// Hard constraints:
// - Do NOT modify any existing HE3DB/TFHEpp sources.
// - Instrumentation is done by copying the original implementations into this
//   header under a new namespace, and adding optional trace prints only.

#include "experiments/he3db_phase_trace.hpp"

#include "HEDB/comparison/tfhepp_utils.h"
#include "HEDB/utils/types.h"

#include <tfhe++.hpp>

#include <cstdint>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>

namespace HE3DBTrace {

// ============================
// Copied HE3DB implementations
// ============================

// Copied from: src/HEDB/comparison/tfhepp_utils.cpp
// Original symbol: TFHEpp::MSBGateBootstrapping(TLWE<lvl1param>&, ...)
// 行为保持一致；仅在 trace_enable=true 时打印关键 phase（不改计算顺序）。
inline void
MSBGateBootstrappingTrace(HEDB::TLWELvl1 &res, const HEDB::TLWELvl1 &tlwe,
                          const HEDB::TFHEEvalKey &ek, bool result_type,
                          uint32_t plain_bits, uint32_t err_bits, uint64_t m,
                          const HEDB::TFHESecretKey *sk_dbg, bool trace_enable,
                          uint32_t /*max_print*/) {
  constexpr uint32_t W = torus_bits_v<TFHEpp::lvl1param::T>;

  // 这里打印的是“MSB门内部中间态”的 phase：
  // - 输入 tlwe 的 phase 仍然承载整数编码语义
  // - 加 offset 后的 tlweoffset phase 仍然承载整数编码语义
  // 最终输出 res 是 bit 密文，不再适用三段格式，所以不打印 res 的三段。
  if (trace_enable && sk_dbg != nullptr) {
    TFHEpp::TLWE<TFHEpp::lvl1param> tlwe_in = tlwe; // 保守：确保类型一致
    const auto phase_in =
        tlwe_phase<TFHEpp::lvl1param>(tlwe_in, sk_dbg->key.lvl1);
    print_phase_layout(std::cout, "MSB门：输入phase", phase_in, m, plain_bits,
                       err_bits, W);
  }

  uint32_t μ = 1U << 29;
  if (IS_ARITHMETIC(result_type))
    μ = μ << 1;
  constexpr uint64_t offset =
      1ULL << (std::numeric_limits<TFHEpp::lvl1param::T>::digits - 6);
  TFHEpp::TLWE<TFHEpp::lvl1param> tlweoffset = tlwe;
  tlweoffset[TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] += offset;

  if (trace_enable && sk_dbg != nullptr) {
    const auto phase_offset =
        tlwe_phase<TFHEpp::lvl1param>(tlweoffset, sk_dbg->key.lvl1);
    print_phase_layout(std::cout, "MSB门：加offset后phase", phase_offset, m,
                       plain_bits, err_bits, W);
    std::cout << "（注）MSB门最终输出为bit密文(res)"
                 "，不再承载整数编码语义，因此不对res做三段解析。\n";
  }

  TFHEpp::TLWE<TFHEpp::lvl0param> tlwelvl0;
  TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(tlwelvl0, tlweoffset,
                                                *ek.iksklvl10);
  if (trace_enable && sk_dbg != nullptr) {
    constexpr uint32_t W0 = torus_bits_v<TFHEpp::lvl0param::T>;
    const auto phase_lvl0 =
        tlwe_phase<TFHEpp::lvl0param>(tlwelvl0, sk_dbg->key.lvl0);
    print_phase_layout(std::cout, "MSB门：IKS后lvl0 phase", phase_lvl0, m,
                       plain_bits, err_bits, W0);
  }
  TFHEpp::GateBootstrappingTLWE2TLWEFFT<TFHEpp::lvl01param>(
      res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<TFHEpp::lvl1param>(μ));
  if (IS_ARITHMETIC(result_type))
    res[TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] += μ;
}

// Copied from: src/HEDB/comparison/extract_msb.cpp
// Original symbol: HEDB::ExtractMSB5(...)
// 行为保持一致；trace 仅用于输出 got_msb5（不打印中间 phase）。
inline void ExtractMSB5Trace(HEDB::TLWELvl1 &res, const HEDB::TLWELvl1 &tlwe,
                             const HEDB::TFHEEvalKey &ek, bool result_type,
                             uint32_t plain_bits, uint32_t err_bits, uint64_t m,
                             const HEDB::TFHESecretKey *sk_dbg,
                             bool trace_enable, uint32_t max_print) {
  if (trace_enable && sk_dbg != nullptr) {
    std::cout << "\n（失败链路）ExtractMSB5Trace：内部 MSB 门相位\n";
  }
  MSBGateBootstrappingTrace(res, tlwe, ek, result_type, plain_bits, err_bits, m,
                            sk_dbg, trace_enable, max_print);

  if (trace_enable && sk_dbg != nullptr) {
    const bool got = TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(res, sk_dbg->key.lvl1);
    std::cout << "复跑得到 got_msb5=" << (got ? 1 : 0) << "（仅供确认）\n";
  }
}

// Copied from: src/HEDB/comparison/tfhepp_utils.cpp
// Original symbol: TFHEpp::IdeGateBootstrapping(TLWE<lvl1param>&, ...)
// 说明：原实现里 plain_bits 是常量；这里为了实验“刷新但不改
// m”，我们将其参数化。
inline void
IdeGateBootstrappingTrace(HEDB::TLWELvl1 &res, const HEDB::TLWELvl1 &tlwe,
                          uint32_t scale_bits, const HEDB::TFHEEvalKey &ek,
                          uint32_t plain_bits, uint32_t err_bits, uint64_t m,
                          const HEDB::TFHESecretKey *sk_dbg, bool trace_enable,
                          uint32_t /*max_print*/) {
  constexpr uint32_t W = torus_bits_v<HEDB::Lvl1::T>;
  if (trace_enable && sk_dbg != nullptr) {
    const auto phase_in = tlwe_phase<HEDB::Lvl1>(tlwe, sk_dbg->key.lvl1);
    print_phase_layout(std::cout, "加密后phase", phase_in, m, plain_bits,
                       err_bits, W);
  }

  constexpr uint64_t offset =
      1ULL << (std::numeric_limits<TFHEpp::lvl1param::T>::digits - 6);
  TFHEpp::TLWE<TFHEpp::lvl1param> tlweoffset = tlwe;
  tlweoffset[TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] += offset;

  // 额外打印：刷新流程里“加offset后”的相位，帮助区分：
  // - 是否 offset 本身就已经造成回绕/跨区
  // - 还是 PBS/量化误差导致跨区
  if (trace_enable && sk_dbg != nullptr) {
    const auto phase_off =
        tlwe_phase<TFHEpp::lvl1param>(tlweoffset, sk_dbg->key.lvl1);
    print_phase_layout(std::cout, "刷新流程：加offset后phase", phase_off, m,
                       plain_bits, err_bits, W);
  }

  TFHEpp::TLWE<TFHEpp::lvl0param> tlwelvl0;
  TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(tlwelvl0, tlweoffset,
                                                *ek.iksklvl10);

  TFHEpp::GateBootstrappingTLWE2TLWEFFT<TFHEpp::lvl01param>(
      res, tlwelvl0, *ek.bkfftlvl01,
      TFHEpp::gpolygen<TFHEpp::lvl1param>(plain_bits, scale_bits));

  if (trace_enable && sk_dbg != nullptr) {
    const auto phase_out = tlwe_phase<HEDB::Lvl1>(res, sk_dbg->key.lvl1);
    print_phase_layout(std::cout, "刷新后phase", phase_out, m, plain_bits,
                       err_bits, W);
  }
}

} // namespace HE3DBTrace
