// Benchmark: compare HomMSB runtime across different plain_bits (legal ranges).
//
// Hard constraints:
// - Do NOT modify any existing HE3DB/TFHEpp sources.
// - Only add/modify files under experiments/.

#include "HEDB/comparison/extract_msb.h"
#include "HEDB/utils/types.h"

// Optional: allow benchmarking the experimental MSB6 implementation.
#include "experiments/he3db_msb6.hpp"

#include <tfhe++.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

#ifndef USE_RANDEN
#error "This benchmark requires -DUSE_RANDEN so TFHEpp::generator is seedable."
#endif

using Clock = std::chrono::steady_clock;

struct Args {
    uint32_t level = 1;  // 1: TLWELvl1 input (plain_bits<=10); 2: TLWELvl2 input (plain_bits<=33)
    uint32_t plain_bits_min = 5;
    uint32_t plain_bits_max = 10;
    uint32_t iters = 20;
    uint32_t warmup = 2;
    uint64_t base_seed = 0xC0FFEEULL;
    bool csv = false;
    bool also_msb6 = false;  // only meaningful when level==1 and plain_bits==6
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
        else if (a == "--iters") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --iters");
            args.iters = static_cast<uint32_t>(*v);
        }
        else if (a == "--warmup") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --warmup");
            args.warmup = static_cast<uint32_t>(*v);
        }
        else if (a == "--base-seed") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --base-seed");
            args.base_seed = *v;
        }
        else if (a == "--csv") {
            args.csv = true;
        }
        else if (a == "--also-msb6") {
            args.also_msb6 = true;
        }
        else if (a == "--help" || a == "-h") {
            std::cout
                << "用法：msb_timing_bench [选项]\n"
                << "  --level 1|2       输入密文等级：1=TLWELvl1(plain_bits<=10)，2=TLWELvl2(plain_bits<=33)\n"
                << "  --min bits        plain_bits 最小值（默认：level=1 用 5；level=2 用 5）\n"
                << "  --max bits        plain_bits 最大值（默认：level=1 用 10；level=2 用 33）\n"
                << "  --iters N         每个 plain_bits 的计时迭代次数（默认 20）\n"
                << "  --warmup N        每个 plain_bits 的预热次数（默认 2，不计时）\n"
                << "  --base-seed S     生成不同输入密文的随机种子基数（默认 0xC0FFEE）\n"
                << "  --also-msb6       额外测一次实验版 MSB6(offset=Δ/2)，只对 plain_bits=6/level=1 生效\n"
                << "  --csv             输出 CSV（便于画图）\n";
            std::exit(0);
        }
        else {
            throw std::runtime_error("unknown arg: " + std::string(a));
        }
    }
    return args;
}

std::string impl_name(uint32_t level, uint32_t plain_bits)
{
    if (level == 1) {
        if (plain_bits <= 5) return "ExtractMSB5";
        if (plain_bits <= 9) return "ExtractMSB9";
        if (plain_bits == 10) return "ExtractMSB10";
        return "PlainBitsOut";
    }
    if (level == 2) {
        if (plain_bits <= 5) return "ImExtractMSB5";
        if (plain_bits <= 9) return "ImExtractMSB9";
        if (plain_bits <= 14) return "ImExtractMSB14";
        if (plain_bits <= 19) return "ImExtractMSB19";
        if (plain_bits <= 24) return "ImExtractMSB24";
        if (plain_bits <= 29) return "ImExtractMSB29";
        if (plain_bits <= 33) return "ImExtractMSB33";
        return "PlainBitsOut";
    }
    return "BadLevel";
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

struct BenchStats {
    double mean_ms = 0.0;
    double median_ms = 0.0;
    double min_ms = 0.0;
    double max_ms = 0.0;
    bool got = false;
};

BenchStats compute_stats(std::vector<double> ms, bool got)
{
    if (ms.empty()) return BenchStats{0, 0, 0, 0, got};
    double sum = 0.0;
    double mn = ms[0];
    double mx = ms[0];
    for (double v : ms) {
        sum += v;
        mn = std::min(mn, v);
        mx = std::max(mx, v);
    }
    std::sort(ms.begin(), ms.end());
    const double median = ms[ms.size() / 2];
    return BenchStats{sum / static_cast<double>(ms.size()), median, mn, mx, got};
}

template <class CT>
BenchStats bench_hommsb(const std::vector<CT> &cts, uint32_t warmup,
                        uint32_t plain_bits, const HEDB::TFHEEvalKey &ek,
                        const HEDB::TFHESecretKey &sk)
{
    HEDB::TLWELvl1 out{};
    for (uint32_t i = 0; i < warmup; ++i) {
        HEDB::HomMSB(out, cts[i], plain_bits, ek, LOGIC);
    }

    std::vector<double> ms;
    ms.reserve(cts.size() - warmup);
    for (size_t i = warmup; i < cts.size(); ++i) {
        const auto t0 = Clock::now();
        HEDB::HomMSB(out, cts[i], plain_bits, ek, LOGIC);
        const auto t1 = Clock::now();
        const double dt =
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        ms.push_back(dt);
    }
    const bool got = decrypt_bit_lvl1(out, sk);
    return compute_stats(std::move(ms), got);
}

BenchStats bench_msb6(const std::vector<HEDB::TLWELvl1> &cts, uint32_t warmup,
                      const HEDB::TFHEEvalKey &ek, const HEDB::TFHESecretKey &sk)
{
    HEDB::TLWELvl1 out{};
    for (uint32_t i = 0; i < warmup; ++i) {
        HE3DBExperiment::ExtractMSB6(out, cts[i], ek, LOGIC);
    }

    std::vector<double> ms;
    ms.reserve(cts.size() - warmup);
    for (size_t i = warmup; i < cts.size(); ++i) {
        const auto t0 = Clock::now();
        HE3DBExperiment::ExtractMSB6(out, cts[i], ek, LOGIC);
        const auto t1 = Clock::now();
        const double dt =
            std::chrono::duration<double, std::milli>(t1 - t0).count();
        ms.push_back(dt);
    }
    const bool got = decrypt_bit_lvl1(out, sk);
    return compute_stats(std::move(ms), got);
}

}  // namespace

int main(int argc, char **argv)
{
    try {
        Args args = parse_args(argc, argv);
        if (args.level != 1 && args.level != 2)
            throw std::runtime_error("--level 只能是 1 或 2");
        if (args.iters == 0) throw std::runtime_error("--iters 必须 >= 1");

        // Fill defaults based on level if user kept the old defaults.
        if (args.level == 2 && args.plain_bits_max == 10) {
            args.plain_bits_max = 33;
        }

        const uint32_t legal_max = (args.level == 1) ? 10 : 33;
        if (args.plain_bits_min < 1 || args.plain_bits_min > args.plain_bits_max)
            throw std::runtime_error("plain_bits 范围非法");
        if (args.plain_bits_max > legal_max)
            throw std::runtime_error("plain_bits 超出该 level 的合法范围");

        // Keygen
        HEDB::TFHESecretKey sk;
        HEDB::TFHEEvalKey ek;
        ek.emplaceiksk<TFHEpp::lvl10param>(sk);
        ek.emplacebkfft<TFHEpp::lvl01param>(sk);
        if (args.level == 2) {
            ek.emplaceiksk<TFHEpp::lvl20param>(sk);
            ek.emplaceiksk<TFHEpp::lvl21param>(sk);
            ek.emplacebkfft<TFHEpp::lvl02param>(sk);
        }

        if (args.csv) {
            std::cout << "level,plain_bits,impl,iters,mean_ms,median_ms,min_ms,max_ms,got,expected\n";
        }
        else {
            std::cout << "HE3DB experiments：MSB 提取耗时对比（按 plain_bits）\n";
            std::cout << "level=" << args.level << "，扫描 plain_bits=["
                      << args.plain_bits_min << ", " << args.plain_bits_max
                      << "]，iters=" << args.iters << "，warmup=" << args.warmup
                      << "，base_seed=" << args.base_seed << "\n";
            if (args.also_msb6)
                std::cout << "附加：plain_bits=6 时额外测 MSB6(offset=Δ/2)\n";
        }

        for (uint32_t plain_bits = args.plain_bits_min;
             plain_bits <= args.plain_bits_max; ++plain_bits) {
            const uint64_t m = (1ULL << plain_bits) - 1;  // worst-case: 111..11
            const bool expected = ((m >> (plain_bits - 1)) & 1ULL) != 0;
            const std::string impl = impl_name(args.level, plain_bits);

            const uint32_t total = args.warmup + args.iters;
            const uint64_t seed_base =
                args.base_seed + (static_cast<uint64_t>(args.level) << 40) +
                (static_cast<uint64_t>(plain_bits) << 20);

            if (args.level == 1) {
                std::vector<HEDB::TLWELvl1> cts;
                cts.reserve(total);
                for (uint32_t i = 0; i < total; ++i) {
                    TFHEpp::generator.seed(seed_base + i);
                    cts.push_back(encrypt_int_lvl1(m, plain_bits, sk));
                }

                const BenchStats s = bench_hommsb(
                    cts, args.warmup, plain_bits, ek, sk);
                if (args.csv) {
                    std::cout << args.level << "," << plain_bits << ",HomMSB("
                              << impl << ")," << args.iters << ","
                              << s.mean_ms << "," << s.median_ms << ","
                              << s.min_ms << "," << s.max_ms << ","
                              << (s.got ? 1 : 0) << "," << (expected ? 1 : 0)
                              << "\n";
                }
                else {
                    std::cout << "\n[plain_bits=" << plain_bits << "] 方法=HomMSB("
                              << impl << ")\n";
                    std::cout << "  平均=" << s.mean_ms << " ms，中位数="
                              << s.median_ms << " ms，min=" << s.min_ms
                              << " ms，max=" << s.max_ms << " ms\n";
                    std::cout << "  解密输出 got=" << (s.got ? 1 : 0)
                              << "，expected=" << (expected ? 1 : 0) << "\n";
                }

                if (args.also_msb6 && plain_bits == 6) {
                    const BenchStats s6 = bench_msb6(cts, args.warmup, ek, sk);
                    if (args.csv) {
                        std::cout << args.level << "," << plain_bits
                                  << ",MSB6(offset=Δ/2)," << args.iters << ","
                                  << s6.mean_ms << "," << s6.median_ms << ","
                                  << s6.min_ms << "," << s6.max_ms << ","
                                  << (s6.got ? 1 : 0) << ","
                                  << (expected ? 1 : 0) << "\n";
                    }
                    else {
                        std::cout << "  [附加] MSB6(offset=Δ/2)\n";
                        std::cout << "    平均=" << s6.mean_ms
                                  << " ms，中位数=" << s6.median_ms
                                  << " ms，min=" << s6.min_ms
                                  << " ms，max=" << s6.max_ms << " ms\n";
                        std::cout << "    解密输出 got=" << (s6.got ? 1 : 0)
                                  << "，expected=" << (expected ? 1 : 0) << "\n";
                    }
                }
            }
            else {
                std::vector<HEDB::TLWELvl2> cts;
                cts.reserve(total);
                for (uint32_t i = 0; i < total; ++i) {
                    TFHEpp::generator.seed(seed_base + i);
                    cts.push_back(encrypt_int_lvl2(m, plain_bits, sk));
                }

                const BenchStats s = bench_hommsb(
                    cts, args.warmup, plain_bits, ek, sk);
                if (args.csv) {
                    std::cout << args.level << "," << plain_bits << ",HomMSB("
                              << impl << ")," << args.iters << ","
                              << s.mean_ms << "," << s.median_ms << ","
                              << s.min_ms << "," << s.max_ms << ","
                              << (s.got ? 1 : 0) << "," << (expected ? 1 : 0)
                              << "\n";
                }
                else {
                    std::cout << "\n[plain_bits=" << plain_bits
                              << "] 方法=HomMSB(" << impl << ")\n";
                    std::cout << "  平均=" << s.mean_ms << " ms，中位数="
                              << s.median_ms << " ms，min=" << s.min_ms
                              << " ms，max=" << s.max_ms << " ms\n";
                    std::cout << "  解密输出 got=" << (s.got ? 1 : 0)
                              << "，expected=" << (expected ? 1 : 0) << "\n";
                }
            }
        }

        return 0;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

