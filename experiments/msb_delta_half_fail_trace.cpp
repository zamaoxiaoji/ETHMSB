// Experiment: "universal MSB gate" with offset=Δ/2 for arbitrary plain_bits.
//
// Goal:
// - 验证：当 plain_bits 增大（Δ 变小）且输入选择 worst-case m=2^plain_bits-1 时，
//   即便 offset 按 Δ/2 自适应，PBS/量化误差也可能侵入 message 区并触发回绕，导致 MSB 失败。
// - 在失败样本上用同 seed 复跑并开启 trace，打印：
//   - 输入 phase
//   - 加 offset 后 phase
//   - IKS 后 lvl0 phase（进入 BlindRotate 前）
//   - BlindRotateTrace 内部的 q_phase_direct/q_phase_parts 对照
//   - SampleExtract 后输出 bit 的 raw phase + 解密结果
//
// Hard constraints:
// - Do NOT modify any existing HE3DB/TFHEpp sources.
// - Only add/modify files under experiments/.
// - No unconditional printing inside trial loops: only on failure replay.

#include "experiments/he3db_msb6.hpp"

#include "HEDB/utils/types.h"

#include <tfhe++.hpp>

#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

#ifndef USE_RANDEN
#error "This experiment requires -DUSE_RANDEN so TFHEpp::generator is seedable."
#endif

struct Args {
    uint32_t level = 1;  // 1: 输入 TLWELvl1；2: 输入 TLWELvl2（仍输出 lvl1 bit）
    uint32_t trials = 1000;
    uint32_t plain_bits_min = 6;
    uint32_t plain_bits_max = 24;
    uint32_t err_bits = 8;
    uint32_t max_fail_traces_per_bits = 5;
    uint32_t max_nonzero_a_print = 10;  // BlindRotateTrace 内最多打印多少个非零 ai_q
    uint64_t base_seed = 0xC0FFEEULL;
    bool run_control = true;
    bool stop_on_first_bits_with_failure = false;
    bool collect_qphase = false;  // 仅收集 q_phase_* 统计（不运行 PBS 解密比对）
    bool csv = false;
    enum class CaseMode { msb0, msb1 } case_mode = CaseMode::msb1;
};

std::optional<uint64_t> parse_u64(std::string_view s)
{
    if (s.size() >= 2 && s[0] == '0' && (s[1] == 'x' || s[1] == 'X')) {
        uint64_t v = 0;
        for (size_t i = 2; i < s.size(); ++i) {
            const char c = s[i];
            uint8_t digit = 0;
            if (c >= '0' && c <= '9')
                digit = static_cast<uint8_t>(c - '0');
            else if (c >= 'a' && c <= 'f')
                digit = static_cast<uint8_t>(10 + (c - 'a'));
            else if (c >= 'A' && c <= 'F')
                digit = static_cast<uint8_t>(10 + (c - 'A'));
            else
                return std::nullopt;
            if (v > (std::numeric_limits<uint64_t>::max() >> 4)) return std::nullopt;
            v = (v << 4) | digit;
        }
        return v;
    }

    uint64_t v = 0;
    for (char c : s) {
        if (c < '0' || c > '9') return std::nullopt;
        const uint64_t digit = static_cast<uint64_t>(c - '0');
        if (v > (std::numeric_limits<uint64_t>::max() - digit) / 10)
            return std::nullopt;
        v = v * 10 + digit;
    }
    return v;
}

Args parse_args(int argc, char **argv)
{
    Args args;
    auto require = [&](int &i) -> std::string_view {
        if (i + 1 >= argc) throw std::runtime_error("missing value");
        return argv[++i];
    };

    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "--level") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --level");
            args.level = static_cast<uint32_t>(*v);
        }
        else if (a == "--trials") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --trials");
            args.trials = static_cast<uint32_t>(*v);
        }
        else if (a == "--min") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --min");
            args.plain_bits_min = static_cast<uint32_t>(*v);
        }
        else if (a == "--max") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --max");
            args.plain_bits_max = static_cast<uint32_t>(*v);
        }
        else if (a == "--err-bits") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --err-bits");
            args.err_bits = static_cast<uint32_t>(*v);
        }
        else if (a == "--K") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --K");
            args.max_fail_traces_per_bits = static_cast<uint32_t>(*v);
        }
        else if (a == "--max-nonzero-a") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --max-nonzero-a");
            args.max_nonzero_a_print = static_cast<uint32_t>(*v);
        }
        else if (a == "--case") {
            const std::string_view v = require(i);
            if (v == "msb0")
                args.case_mode = Args::CaseMode::msb0;
            else if (v == "msb1")
                args.case_mode = Args::CaseMode::msb1;
            else
                throw std::runtime_error("invalid --case (use msb0/msb1)");
        }
        else if (a == "--base-seed") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --base-seed");
            args.base_seed = *v;
        }
        else if (a == "--no-control") {
            args.run_control = false;
        }
        else if (a == "--stop") {
            args.stop_on_first_bits_with_failure = true;
        }
        else if (a == "--collect-qphase") {
            args.collect_qphase = true;
        }
        else if (a == "--csv") {
            args.csv = true;
        }
        else if (a == "--help" || a == "-h") {
            std::cout
                << "用法：msb_delta_half_fail_trace [选项]\n"
                << "  --level 1|2       输入密文等级：1=TLWELvl1，2=TLWELvl2（仍输出 lvl1 bit）\n"
                << "  --trials N        每个 plain_bits 的试验次数（默认 1000）\n"
                << "  --min bits        扫描 plain_bits 最小值（默认 6）\n"
                << "  --max bits        扫描 plain_bits 最大值（默认 24）\n"
                << "  --err-bits e      三段格式里的误差区低位 bits（默认 8）\n"
                << "  --case msb0|msb1  worst-case 明文选择（默认 msb1）\n"
                << "     msb1: m=2^t-1 (111..11, 最接近 2^W 回绕边界)\n"
                << "     msb0: m=2^(t-1)-1 (011..11, 更接近 half 边界)\n"
                << "  --K K             每个 plain_bits 最多打印 K 个失败样本 trace（默认 5）\n"
                << "  --max-nonzero-a N BlindRotateTrace 内最多打印 N 个非零 ai_q（默认 10）\n"
                << "  --base-seed S     随机种子基数（默认 0xC0FFEE，可用 0x...）\n"
                << "  --no-control      不跑对照输入 m3=3*2^(t-2)\n"
                << "  --stop            某个 plain_bits 出现失败后停止继续扫描\n";
            std::cout
                << "  --collect-qphase  仅收集 BlindRotate 量化相位统计（q_phase_direct/q_phase_parts），不运行 PBS\n"
                << "  --csv             与 --collect-qphase 搭配输出 CSV\n";
            std::exit(0);
        }
        else {
            throw std::runtime_error("unknown arg: " + std::string(a));
        }
    }
    return args;
}

HEDB::TLWELvl1 encrypt_int_lvl1(uint64_t m, uint32_t plain_bits,
                               const HEDB::TFHESecretKey &sk)
{
    const uint32_t W = std::numeric_limits<HEDB::Lvl1::T>::digits;
    const uint32_t scale_bits = W - plain_bits;
    const double scale = std::pow(2.0, static_cast<double>(scale_bits));
    return TFHEpp::tlweSymInt32Encrypt<HEDB::Lvl1>(
        static_cast<HEDB::Lvl1::T>(m), HEDB::Lvl1::α, scale, sk.key.lvl1);
}

HEDB::TLWELvl2 encrypt_int_lvl2(uint64_t m, uint32_t plain_bits,
                               const HEDB::TFHESecretKey &sk)
{
    const uint32_t W = std::numeric_limits<HEDB::Lvl2::T>::digits;
    const uint32_t scale_bits = W - plain_bits;
    const double scale = std::pow(2.0, static_cast<double>(scale_bits));
    return TFHEpp::tlweSymInt32Encrypt<HEDB::Lvl2>(
        static_cast<HEDB::Lvl2::T>(m), HEDB::Lvl2::α, scale, sk.key.lvl2);
}

bool decrypt_bit_lvl1(const HEDB::TLWELvl1 &ct, const HEDB::TFHESecretKey &sk)
{
    return TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(ct, sk.key.lvl1);
}

uint64_t worst_m(uint32_t plain_bits, Args::CaseMode mode)
{
    if (plain_bits == 0) return 0;
    if (mode == Args::CaseMode::msb1) return (1ULL << plain_bits) - 1;
    // msb0: 0111..11
    return (1ULL << (plain_bits - 1)) - 1;
}

uint32_t floor_log2_u32(uint32_t x)
{
    if (x == 0) return 0;
#if defined(__GNUG__)
    return 31U - static_cast<uint32_t>(__builtin_clz(x));
#else
    uint32_t r = 0;
    while (x >>= 1) ++r;
    return r;
#endif
}

struct VecStats {
    uint32_t min = 0;
    uint32_t p10 = 0;
    uint32_t p50 = 0;
    uint32_t p90 = 0;
    uint32_t max = 0;
    double mean = 0.0;
};

VecStats compute_stats(std::vector<uint32_t> v)
{
    VecStats s{};
    if (v.empty()) return s;
    uint64_t sum = 0;
    uint32_t mn = v[0];
    uint32_t mx = v[0];
    for (uint32_t x : v) {
        sum += x;
        if (x < mn) mn = x;
        if (x > mx) mx = x;
    }
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    auto q = [&](double p) -> uint32_t {
        const size_t idx = static_cast<size_t>(p * static_cast<double>(n - 1));
        return v[idx];
    };
    s.min = mn;
    s.max = mx;
    s.p10 = q(0.10);
    s.p50 = q(0.50);
    s.p90 = q(0.90);
    s.mean = static_cast<double>(sum) / static_cast<double>(n);
    return s;
}

}  // namespace

int main(int argc, char **argv)
{
    try {
        const Args args = parse_args(argc, argv);
        if (args.level != 1 && args.level != 2)
            throw std::runtime_error("--level 只能是 1 或 2");
        if (args.trials == 0) throw std::runtime_error("--trials 必须 >= 1");
        if (args.plain_bits_min == 0 || args.plain_bits_min > args.plain_bits_max)
            throw std::runtime_error("plain_bits 范围非法");

        // Keygen
        HEDB::TFHESecretKey sk;
        HEDB::TFHEEvalKey ek;
        ek.emplaceiksk<TFHEpp::lvl10param>(sk);
        ek.emplacebkfft<TFHEpp::lvl01param>(sk);
        if (args.level == 2) {
            ek.emplaceiksk<TFHEpp::lvl20param>(sk);
            ek.emplacebkfft<TFHEpp::lvl02param>(sk);
        }

        const uint32_t W =
            (args.level == 1) ? std::numeric_limits<HEDB::Lvl1::T>::digits
                              : std::numeric_limits<HEDB::Lvl2::T>::digits;

        if (args.collect_qphase) {
            if (args.csv) {
                std::cout << "level,case,plain_bits,trials,shift_bs,mod,half,"
                             "q_direct_p50,q_parts_p50,dist2N_direct_p50,dist2N_parts_p50,"
                             "msb_direct_rate,msb_parts_rate,msb_diff_rate\n";
            }
            else {
                std::cout << "HE3DB experiments：收集 BlindRotate 量化相位统计（不运行 PBS）\n";
                std::cout << "参数：level=" << args.level << "，trials=" << args.trials
                          << "，base_seed=" << args.base_seed << "\n";
                std::cout << "扫描 plain_bits 范围=[" << args.plain_bits_min << ", "
                          << args.plain_bits_max << "]，W=" << W << "\n";
                std::cout << "case="
                          << ((args.case_mode == Args::CaseMode::msb1) ? "msb1" : "msb0")
                          << "（仅影响 m 的选择）\n";
                std::cout << "说明：这里统计的是进入 BlindRotate 前（IKS 后 lvl0 TLWE）的\n"
                             "      q_phase_direct（精确phase后量化）与 q_phase_parts（分量量化路径）。\n";
            }

            // Quantization domain for lvl01param is fixed (domain=lvl0, target=lvl1).
            for (uint32_t plain_bits = args.plain_bits_min;
                 plain_bits <= args.plain_bits_max; ++plain_bits) {
                const uint64_t m = worst_m(plain_bits, args.case_mode);

                std::vector<uint32_t> q_direct_v;
                std::vector<uint32_t> q_parts_v;
                std::vector<uint32_t> dist_direct_v;
                std::vector<uint32_t> dist_parts_v;
                q_direct_v.reserve(args.trials);
                q_parts_v.reserve(args.trials);
                dist_direct_v.reserve(args.trials);
                dist_parts_v.reserve(args.trials);

                uint32_t msb_direct_ones = 0;
                uint32_t msb_parts_ones = 0;
                uint32_t msb_diff = 0;
                uint32_t shift_bs_out = 0;
                uint32_t mod_out = 0;
                uint32_t half_out = 0;
                bool params_set = false;

                for (uint32_t trial = 0; trial < args.trials; ++trial) {
                    const uint64_t seed =
                        args.base_seed + (static_cast<uint64_t>(args.level) << 48) +
                        (static_cast<uint64_t>(plain_bits) << 16) + trial;
                    TFHEpp::generator.seed(seed);

                    TFHEpp::TLWE<TFHEpp::lvl0param> tlwelvl0{};
                    if (args.level == 1) {
                        using InP = TFHEpp::lvl1param;
                        const auto ct = encrypt_int_lvl1(m, plain_bits, sk);
                        const uint64_t offset =
                            HE3DBExperiment::delta_half_offset_u64<typename InP::T>(plain_bits);
                        TFHEpp::TLWE<InP> ct_off = ct;
                        ct_off[InP::k * InP::n] += static_cast<typename InP::T>(offset);
                        TFHEpp::IdentityKeySwitch<TFHEpp::lvl10param>(
                            tlwelvl0, ct_off, *ek.iksklvl10);
                    }
                    else {
                        using InP = TFHEpp::lvl2param;
                        const auto ct = encrypt_int_lvl2(m, plain_bits, sk);
                        const uint64_t offset =
                            HE3DBExperiment::delta_half_offset_u64<typename InP::T>(plain_bits);
                        TFHEpp::TLWE<InP> ct_off = ct;
                        ct_off[InP::k * InP::n] += static_cast<typename InP::T>(offset);
                        TFHEpp::IdentityKeySwitch<TFHEpp::lvl20param>(
                            tlwelvl0, ct_off, *ek.iksklvl20);
                    }

                    const auto rec =
                        TFHEppInstrumented::QuantizePhaseRecord<TFHEpp::lvl01param>(
                            tlwelvl0, sk.key.lvl0, plain_bits, m, seed);

                    const uint32_t mod = static_cast<uint32_t>(rec.mod);
                    const uint32_t half = static_cast<uint32_t>(rec.half);
                    if (!params_set) {
                        shift_bs_out = rec.shift_bs;
                        mod_out = mod;
                        half_out = half;
                        params_set = true;
                    }
                    const uint32_t qd = static_cast<uint32_t>(rec.q_phase_direct);
                    const uint32_t qp = static_cast<uint32_t>(rec.q_phase_parts);

                    // “距离 2N（从下方）”：值越接近 2N(=mod)，dist 越小；
                    // 若发生回绕（qp 很小甚至 0），则 dist 会很大（接近 mod）。
                    const uint32_t dist_d = (qd == 0) ? mod : (mod - qd);
                    const uint32_t dist_p = (qp == 0) ? mod : (mod - qp);

                    q_direct_v.push_back(qd);
                    q_parts_v.push_back(qp);
                    dist_direct_v.push_back(dist_d);
                    dist_parts_v.push_back(dist_p);

                    const bool msb_d = (qd >= half);
                    const bool msb_p = (qp >= half);
                    if (msb_d) ++msb_direct_ones;
                    if (msb_p) ++msb_parts_ones;
                    if (msb_d != msb_p) ++msb_diff;
                }

                const VecStats qd_s = compute_stats(q_direct_v);
                const VecStats qp_s = compute_stats(q_parts_v);
                const VecStats dd_s = compute_stats(dist_direct_v);
                const VecStats dp_s = compute_stats(dist_parts_v);

                // A compact “order of magnitude” indicator for the median distance.
                const uint32_t dd_k = floor_log2_u32(dd_s.p50);
                const uint32_t dp_k = floor_log2_u32(dp_s.p50);

                const double msb_d_rate =
                    static_cast<double>(msb_direct_ones) / static_cast<double>(args.trials);
                const double msb_p_rate =
                    static_cast<double>(msb_parts_ones) / static_cast<double>(args.trials);
                const double msb_diff_rate =
                    static_cast<double>(msb_diff) / static_cast<double>(args.trials);

                if (args.csv) {
                    std::cout << args.level << ","
                              << ((args.case_mode == Args::CaseMode::msb1) ? "msb1" : "msb0")
                              << "," << plain_bits << "," << args.trials << ","
                              << shift_bs_out << "," << mod_out << "," << half_out << ","
                              << qd_s.p50 << "," << qp_s.p50 << ","
                              << dd_s.p50 << "," << dp_s.p50 << ","
                              << msb_d_rate << "," << msb_p_rate << ","
                              << msb_diff_rate << "\n";
                }
                else {
                    std::cout << "\n[plain_bits=" << plain_bits << "] m=" << m << "\n";
                    std::cout << "  q_direct(p50)=" << qd_s.p50
                              << "，dist_to_2N(p50)=" << dd_s.p50 << "≈2^" << dd_k
                              << "（mean=" << dd_s.mean << "）\n";
                    std::cout << "  q_parts (p50)=" << qp_s.p50
                              << "，dist_to_2N(p50)=" << dp_s.p50 << "≈2^" << dp_k
                              << "（mean=" << dp_s.mean << "）\n";
                    std::cout << "  MSB率：direct=" << msb_d_rate
                              << "，parts=" << msb_p_rate
                              << "，MSB不一致率=" << msb_diff_rate << "\n";
                }

                if (args.stop_on_first_bits_with_failure && msb_diff > 0) break;
            }

            return 0;
        }

        std::cout << "HE3DB experiments：验证 offset=Δ/2 的通用 MSB 门失败与 trace\n";
        std::cout << "参数：level=" << args.level << "，trials=" << args.trials
                  << "，每个 plain_bits 最多打印 K=" << args.max_fail_traces_per_bits
                  << "，base_seed=" << args.base_seed << "\n";
        std::cout << "扫描 plain_bits 范围=[" << args.plain_bits_min << ", "
                  << args.plain_bits_max << "]，W=" << W << "，err_bits="
                  << args.err_bits << "\n";
        std::cout << "worst-case 选择："
                  << ((args.case_mode == Args::CaseMode::msb1) ? "msb1" : "msb0")
                  << "\n";

        for (uint32_t plain_bits = args.plain_bits_min;
             plain_bits <= args.plain_bits_max; ++plain_bits) {
            const uint64_t m = worst_m(plain_bits, args.case_mode);
            const bool expected = ((m >> (plain_bits - 1)) & 1ULL) != 0;

            uint32_t fails = 0;
            uint32_t printed = 0;

            for (uint32_t trial = 0; trial < args.trials; ++trial) {
                const uint64_t seed =
                    args.base_seed + (static_cast<uint64_t>(args.level) << 48) +
                    (static_cast<uint64_t>(plain_bits) << 16) + trial;
                TFHEpp::generator.seed(seed);

                HEDB::TLWELvl1 out{};
                if (args.level == 1) {
                    const auto ct = encrypt_int_lvl1(m, plain_bits, sk);
                    HE3DBExperiment::MSBGateBootstrappingDeltaHalfTrace(
                        out, ct, ek, LOGIC, plain_bits, args.err_bits, m, seed,
                        /*sk_dbg=*/nullptr, /*trace_enable=*/false,
                        /*max_nonzero_a_print=*/0);
                }
                else {
                    const auto ct = encrypt_int_lvl2(m, plain_bits, sk);
                    HE3DBExperiment::MSBGateBootstrappingDeltaHalfTrace(
                        out, ct, ek, LOGIC, plain_bits, args.err_bits, m, seed,
                        /*sk_dbg=*/nullptr, /*trace_enable=*/false,
                        /*max_nonzero_a_print=*/0);
                }

                const bool got = decrypt_bit_lvl1(out, sk);
                if (got != expected) {
                    ++fails;
                    if (printed < args.max_fail_traces_per_bits) {
                        ++printed;
                        std::cout << "\n========== 失败样本 TRACE（plain_bits=" << plain_bits
                                  << "） ==========\n";
                        std::cout << "trial=" << trial << "，seed=" << seed << "，m=" << m
                                  << "，expected=" << (expected ? 1 : 0)
                                  << "，got=" << (got ? 1 : 0) << "\n";

                        // 失败样本复跑：同 seed 重加密 + 开启 trace
                        TFHEpp::generator.seed(seed);
                        if (args.level == 1) {
                            const auto ct = encrypt_int_lvl1(m, plain_bits, sk);
                            HE3DBExperiment::MSBGateBootstrappingDeltaHalfTrace(
                                out, ct, ek, LOGIC, plain_bits, args.err_bits, m,
                                seed, &sk, /*trace_enable=*/true,
                                args.max_nonzero_a_print);
                        }
                        else {
                            const auto ct = encrypt_int_lvl2(m, plain_bits, sk);
                            HE3DBExperiment::MSBGateBootstrappingDeltaHalfTrace(
                                out, ct, ek, LOGIC, plain_bits, args.err_bits, m,
                                seed, &sk, /*trace_enable=*/true,
                                args.max_nonzero_a_print);
                        }

                        std::cout
                            << "\n提示：重点看\n"
                               "  1) 加offset后phase 是否靠近/越过 2^W 回绕点\n"
                               "  2) IKS后(lvl0) phase（进入BlindRotate前）是否已经“贴边”\n"
                               "  3) BlindRotateTrace 里 q_phase_direct vs q_phase_parts 是否出现\n"
                               "     msb 分歧且 q_phase_parts 跳到接近 0（量化层回绕/进位）\n";
                    }
                }
            }

            const double rate =
                static_cast<double>(fails) / static_cast<double>(args.trials);
            std::cout << "\n[plain_bits=" << plain_bits << "] worst m=" << m
                      << "（expected=" << (expected ? 1 : 0) << "）\n";
            std::cout << "统计：试验次数=" << args.trials << "，失败数=" << fails
                      << "，失败率=" << rate << "\n";

            if (args.run_control && plain_bits >= 3) {
                const uint64_t m3 = 3ULL << (plain_bits - 2);  // 11..00
                const bool expected3 = ((m3 >> (plain_bits - 1)) & 1ULL) != 0;
                uint32_t fails3 = 0;
                for (uint32_t trial = 0; trial < args.trials; ++trial) {
                    const uint64_t seed =
                        args.base_seed + (static_cast<uint64_t>(args.level) << 48) +
                        (static_cast<uint64_t>(plain_bits) << 16) + trial +
                        0x10000000ULL;
                    TFHEpp::generator.seed(seed);

                    HEDB::TLWELvl1 out{};
                    if (args.level == 1) {
                        const auto ct = encrypt_int_lvl1(m3, plain_bits, sk);
                        HE3DBExperiment::MSBGateBootstrappingDeltaHalfTrace(
                            out, ct, ek, LOGIC, plain_bits, args.err_bits, m3, seed,
                            /*sk_dbg=*/nullptr, /*trace_enable=*/false,
                            /*max_nonzero_a_print=*/0);
                    }
                    else {
                        const auto ct = encrypt_int_lvl2(m3, plain_bits, sk);
                        HE3DBExperiment::MSBGateBootstrappingDeltaHalfTrace(
                            out, ct, ek, LOGIC, plain_bits, args.err_bits, m3, seed,
                            /*sk_dbg=*/nullptr, /*trace_enable=*/false,
                            /*max_nonzero_a_print=*/0);
                    }
                    if (decrypt_bit_lvl1(out, sk) != expected3) ++fails3;
                }
                const double rate3 =
                    static_cast<double>(fails3) / static_cast<double>(args.trials);
                std::cout << "对照输入 m3=3*2^(t-2)=" << m3 << "（expected=" << (expected3 ? 1 : 0) << "）\n";
                std::cout << "  统计：试验次数=" << args.trials << "，失败数=" << fails3
                          << "，失败率=" << rate3 << "\n";
            }

            if (args.stop_on_first_bits_with_failure && fails > 0) break;
        }

        return 0;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
