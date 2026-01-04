// Experiment: show that naive ExtractMSB5 on large-precision integer ciphertexts
// can fail due to wrap-around (carry-chain) when m=2^plain_bits-1.
//
// Hard constraints:
// - Do NOT modify any existing HE3DB/TFHEpp sources.
// - Instrumentation via copied functions in experiments/he3db_msb_trace.hpp.
// - No unconditional printing inside trial loops: only on failure replay with
//   same seed, limited to K samples per plain_bits.

#include "experiments/he3db_msb_trace.hpp"

#include "HEDB/comparison/extract_msb.h"
#include "HEDB/comparison/tfhepp_utils.h"
#include "HEDB/utils/types.h"

#include <tfhe++.hpp>

#include <cmath>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

#ifndef USE_RANDEN
#error "This experiment requires -DUSE_RANDEN so TFHEpp::generator is seedable."
#endif

struct Args {
    uint32_t trials = 500;
    uint32_t plain_bits_min = 6;
    uint32_t plain_bits_max = 10;  // HomMSB(TLWELvl1) only supports up to 10
    uint32_t err_bits = 8;
    uint32_t max_fail_traces_per_bits = 5;
    uint64_t base_seed = 0xC0FFEEULL;
    bool run_control = true;
    bool stop_on_first_bits_with_failure = false;
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
        if (a == "--trials") {
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
        else if (a == "--help" || a == "-h") {
            std::cout
                << "用法：msb5_fail_trace [选项]\n"
                << "  --trials N        试验次数（默认 500）\n"
                << "  --min bits        plain_bits 最小值（默认 6）\n"
                << "  --max bits        plain_bits 最大值（默认 10）\n"
                << "  --err-bits e      误差区低位 bits（默认 8，可设 3）\n"
                << "  --case msb0|msb1  worst-case 明文选择（默认 msb1）\n"
                << "                   msb0: m=2^(t-1)-1  (0 后接长串 1)\n"
                << "                   msb1: m=2^t-1      (全 1，接近回绕边界)\n"
                << "  --K K             每个 plain_bits 最多打印 K 个失败 trace（默认 5）\n"
                << "  --base-seed S     随机种子基数（十进制或 0x..，默认 0xC0FFEE）\n"
                << "  --no-control      不跑对照输入 m3=3*2^(t-2)\n"
                << "  --stop            某个 plain_bits 出现失败后停止继续扫描\n";
            std::exit(0);
        }
        else {
            throw std::runtime_error("unknown arg: " + std::string(a));
        }
    }
    return args;
}

struct Stats {
    uint32_t fails_msb5 = 0;
    uint32_t fails_correct = 0;
    uint32_t correct_trials = 0;
};

HEDB::TLWELvl1 encrypt_int_lvl1(uint64_t m, uint32_t plain_bits,
                               const HEDB::TFHESecretKey &sk)
{
    const uint32_t W = std::numeric_limits<HEDB::Lvl1::T>::digits;
    const uint32_t scale_bits = W - plain_bits;
    const double scale = std::pow(2.0, static_cast<double>(scale_bits));
    return TFHEpp::tlweSymInt32Encrypt<HEDB::Lvl1>(
        static_cast<HEDB::Lvl1::T>(m), HEDB::Lvl1::α, scale, sk.key.lvl1);
}

bool decrypt_bit(const HEDB::TLWELvl1 &ct, const HEDB::TFHESecretKey &sk)
{
    return TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(ct, sk.key.lvl1);
}

Stats run_one_plain_bits(const Args &args, uint32_t plain_bits,
                         const HEDB::TFHESecretKey &sk,
                         const HEDB::TFHEEvalKey &ek, uint64_t m,
                         std::string_view label)
{
    Stats stats{};
    uint32_t printed = 0;
    const bool correct_available = (plain_bits <= 10);
    const uint32_t W = std::numeric_limits<HEDB::Lvl1::T>::digits;
    const uint32_t scale_bits = W - plain_bits;

    for (uint32_t trial = 0; trial < args.trials; ++trial) {
        const uint64_t seed =
            args.base_seed + (static_cast<uint64_t>(plain_bits) << 16) + trial;
        TFHEpp::generator.seed(seed);

        const HEDB::TLWELvl1 ct = encrypt_int_lvl1(m, plain_bits, sk);

        // Path 1 (naive / intentionally fragile): ExtractMSB5 on a multi-bit
        // integer ciphertext.
        HEDB::TLWELvl1 out_msb5;
        HE3DBTrace::ExtractMSB5Trace(out_msb5, ct, ek, LOGIC, plain_bits,
                                     args.err_bits, m, &sk,
                                     /*trace_enable=*/false, /*max_print=*/0);
        const bool got_msb5 = decrypt_bit(out_msb5, sk);

        // Path 2 (HE3DB "correct" for TLWELvl1): HomMSB dispatches to
        // ExtractMSB9/10.
        std::optional<bool> got_correct;
        if (correct_available) {
            HEDB::TLWELvl1 out_correct;
            HEDB::HomMSB(out_correct, ct, plain_bits, ek, LOGIC);
            got_correct = decrypt_bit(out_correct, sk);
            ++stats.correct_trials;
        }

        const bool expected =
            (plain_bits == 0) ? false : (((m >> (plain_bits - 1)) & 1ULL) != 0);

        if (got_msb5 != expected) {
            ++stats.fails_msb5;
            if (printed < args.max_fail_traces_per_bits) {
                ++printed;
                std::cout << "\n========== 失败样本 TRACE（" << label << "） ==========\n";
                std::cout << "plain_bits=" << plain_bits << "，trial=" << trial
                          << "，seed=" << seed << "，m=" << m
                          << "，expected_msb=" << (expected ? 1 : 0)
                          << "，got_msb5=" << (got_msb5 ? 1 : 0);
                if (got_correct.has_value())
                    std::cout << "，got_correct=" << (*got_correct ? 1 : 0);
                else
                    std::cout << "，got_correct=N/A";
                std::cout << "\n";

                // 复跑：同 seed 重加密 + 打开 trace 打印（仅在失败时打印，避免刷屏）
                TFHEpp::generator.seed(seed);
                const HEDB::TLWELvl1 ct_re = encrypt_int_lvl1(m, plain_bits, sk);

                // (0) 复跑失败路径本身：naive ExtractMSB5（打开trace）
                //     这会在 MSBGateBootstrappingTrace 内打印：
                //     - MSB门：输入phase
                //     - MSB门：加offset后phase
                {
                    HEDB::TLWELvl1 out_msb5_trace;
                    HE3DBTrace::ExtractMSB5Trace(
                        out_msb5_trace, ct_re, ek, LOGIC, plain_bits, args.err_bits,
                        m, &sk, /*trace_enable=*/true, /*max_print=*/0);
                }

                // (1) identity refresh：输出仍承载同一个整数编码
                HEDB::TLWELvl1 ct_bs;
                HE3DBTrace::IdeGateBootstrappingTrace(
                    ct_bs, ct_re, scale_bits, ek, plain_bits, args.err_bits, m,
                    &sk, /*trace_enable=*/true, /*max_print=*/0);

                // (2) 用同样三段格式做一行对比（更方便肉眼看 MSB 是否翻转）
                const auto phase_in =
                    HE3DBTrace::tlwe_phase<HEDB::Lvl1>(ct_re, sk.key.lvl1);
                const auto phase_bs =
                    HE3DBTrace::tlwe_phase<HEDB::Lvl1>(ct_bs, sk.key.lvl1);
                const uint64_t phase_in_u = static_cast<uint64_t>(phase_in);
                const uint64_t phase_bs_u = static_cast<uint64_t>(phase_bs);
                const std::string trip_in = HE3DBTrace::format_phase_triplet(
                    phase_in_u, W, plain_bits, args.err_bits);
                const std::string trip_bs = HE3DBTrace::format_phase_triplet(
                    phase_bs_u, W, plain_bits, args.err_bits);
                std::cout << "\n三段对比：\n  " << trip_in << "  ->  " << trip_bs
                          << "\n";

                const uint64_t m_top_in = phase_in_u >> (W - plain_bits);
                const uint64_t m_top_bs = phase_bs_u >> (W - plain_bits);
                const uint64_t msb_in = (m_top_in >> (plain_bits - 1)) & 1ULL;
                const uint64_t msb_bs = (m_top_bs >> (plain_bits - 1)) & 1ULL;
                std::cout
                    << "对比结论：加密后phase 的 m区MSB=" << msb_in
                    << "，刷新后phase 的 m区MSB=" << msb_bs << "。";
                if (msb_in != msb_bs)
                    std::cout << "m区最高位发生变化 -> 属于 MSB 反转现象。\n";
                else
                    std::cout << "m区最高位未变化。\n";
            }
        }

        if (got_correct.has_value() && (*got_correct != expected))
            ++stats.fails_correct;
    }

    return stats;
}

}  // namespace

int main(int argc, char **argv)
{
    try {
        const Args args = parse_args(argc, argv);

        if (args.plain_bits_min < 1 || args.plain_bits_min > args.plain_bits_max)
            throw std::runtime_error("invalid plain_bits range");

        const uint32_t W = std::numeric_limits<HEDB::Lvl1::T>::digits;
        if (args.err_bits == 0) throw std::runtime_error("--err-bits 必须 >= 1");
        if (args.plain_bits_max >= W)
            throw std::runtime_error("plain_bits 不能 >= W");
        if (args.plain_bits_max + args.err_bits >= W)
            throw std::runtime_error(
                "err_bits 太大：需要满足 plain_bits + err_bits < W");

        std::cout << "HE3DB：验证 naive ExtractMSB5 的失败\n";
        std::cout << "参数：trials=" << args.trials
                  << "，每个 plain_bits 最多打印 K=" << args.max_fail_traces_per_bits
                  << "，base_seed=" << args.base_seed << "\n";
        std::cout << "扫描 plain_bits 范围=[" << args.plain_bits_min << ", "
                  << args.plain_bits_max << "]，W=" << W
                  << "，err_bits=" << args.err_bits << "\n";
        std::cout << "worst-case 选择："
                  << ((args.case_mode == Args::CaseMode::msb0) ? "msb0" : "msb1")
                  << "\n";

        // Keygen (HE3DB aliases TFHEpp types).
        HEDB::TFHESecretKey sk;
        HEDB::TFHEEvalKey ek;
        ek.emplaceiksk<TFHEpp::lvl10param>(sk);
        ek.emplacebkfft<TFHEpp::lvl01param>(sk);

        // Sanity: check output polarity for ExtractMSB5 on a bit-like input.
        {
            constexpr uint32_t bits = 5;
            const uint64_t m0 = 0;
            const uint64_t m1 = (1ULL << (bits - 1));
            TFHEpp::generator.seed(0x1111ULL);
            const auto c0 = encrypt_int_lvl1(m0, bits, sk);
            TFHEpp::generator.seed(0x2222ULL);
            const auto c1 = encrypt_int_lvl1(m1, bits, sk);
            HEDB::TLWELvl1 o0, o1;
            HE3DBTrace::ExtractMSB5Trace(o0, c0, ek, LOGIC, bits, args.err_bits,
                                         m0, &sk, false, 0);
            HE3DBTrace::ExtractMSB5Trace(o1, c1, ek, LOGIC, bits, args.err_bits,
                                         m1, &sk, false, 0);
            std::cout << "sanity(bits=5)：m=0 -> " << (decrypt_bit(o0, sk) ? 1 : 0)
                      << "，m=2^(bits-1) -> " << (decrypt_bit(o1, sk) ? 1 : 0)
                      << "（期望 0,1）\n";
        }

        for (uint32_t plain_bits = args.plain_bits_min;
             plain_bits <= args.plain_bits_max; ++plain_bits) {
            if (plain_bits >= W) break;

            uint64_t m_worst = 0;
            if (args.case_mode == Args::CaseMode::msb1) {
                m_worst = (1ULL << plain_bits) - 1;  // 111..11
            }
            else {
                m_worst = (1ULL << (plain_bits - 1)) - 1;  // 011..11
            }
            const uint64_t expected_msb =
                (m_worst >> (plain_bits - 1)) & 1ULL;

            std::cout << "\n[plain_bits=" << plain_bits << "]\n";
            std::cout << "worst m=" << m_worst << "（二进制="
                      << HE3DBTrace::bits_fixed_u64(m_worst, plain_bits)
                      << "，expected_msb=" << expected_msb << "）\n";

            const Stats worst =
                run_one_plain_bits(args, plain_bits, sk, ek, m_worst, "worst");
            std::cout << "统计：试验次数=" << args.trials
                      << "，ExtractMSB5 失败数=" << worst.fails_msb5
                      << "，失败率="
                      << (static_cast<double>(worst.fails_msb5) /
                          static_cast<double>(args.trials))
                      << "\n";
            if (worst.correct_trials > 0) {
                std::cout << "对照（HE3DB 正确方法）失败数=" << worst.fails_correct
                          << " / " << worst.correct_trials << "，失败率="
                          << (static_cast<double>(worst.fails_correct) /
                              static_cast<double>(worst.correct_trials))
                          << "\n";
            }
            else {
                std::cout << "对照（HE3DB 正确方法）=N/A\n";
            }

            if (args.run_control && plain_bits >= 2) {
                const uint64_t m3 = 3ULL << (plain_bits - 2);  // 11..00
                std::cout << "对照输入 m3=3*2^(t-2)=" << m3 << "（二进制="
                          << HE3DBTrace::bits_fixed_u64(m3, plain_bits) << "）\n";
                const Stats ctrl =
                    run_one_plain_bits(args, plain_bits, sk, ek, m3, "control");
                std::cout << "  统计：试验次数=" << args.trials
                          << "，ExtractMSB5 失败数=" << ctrl.fails_msb5
                          << "，失败率="
                          << (static_cast<double>(ctrl.fails_msb5) /
                              static_cast<double>(args.trials))
                          << "\n";
                if (ctrl.correct_trials > 0) {
                    std::cout << "  对照（HE3DB 正确方法）失败数="
                              << ctrl.fails_correct << " / " << ctrl.correct_trials
                              << "，失败率="
                              << (static_cast<double>(ctrl.fails_correct) /
                                  static_cast<double>(ctrl.correct_trials))
                              << "\n";
                }
                else {
                    std::cout << "  对照（HE3DB 正确方法）=N/A\n";
                }
            }

            if (args.stop_on_first_bits_with_failure && worst.fails_msb5 > 0) {
                break;
            }
        }

        return 0;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
