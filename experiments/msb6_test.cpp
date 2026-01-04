// Experiment: MSB6 (offset=Δ/2) vs original ExtractMSB5 on 6-bit inputs.
//
// Goal:
// - Check whether 6-bit inputs (especially worst-case m=2^6-1) still fail.
//
// Hard constraints:
// - Do NOT modify any existing HE3DB/TFHEpp sources.
// - Only add/modify files under experiments/.

#include "experiments/he3db_msb6.hpp"

#include "HEDB/comparison/extract_msb.h"
#include "HEDB/utils/types.h"

#include <tfhe++.hpp>

#include <bitset>
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
    uint64_t base_seed = 0xC0FFEEULL;
    enum class CaseMode { msb1, msb0 } case_mode = CaseMode::msb1;
    std::optional<uint64_t> m_override;
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
            if (v > (std::numeric_limits<uint64_t>::max() >> 4))
                return std::nullopt;
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
        else if (a == "--base-seed") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --base-seed");
            args.base_seed = *v;
        }
        else if (a == "--case") {
            const std::string_view v = require(i);
            if (v == "msb1")
                args.case_mode = Args::CaseMode::msb1;
            else if (v == "msb0")
                args.case_mode = Args::CaseMode::msb0;
            else
                throw std::runtime_error("invalid --case (use msb0/msb1)");
        }
        else if (a == "--m") {
            const auto v = parse_u64(require(i));
            if (!v) throw std::runtime_error("invalid --m");
            args.m_override = *v;
        }
        else if (a == "--help" || a == "-h") {
            std::cout
                << "用法：msb6_test [选项]\n"
                << "  --trials N        试验次数（默认 500）\n"
                << "  --base-seed S     随机种子基值（默认 0xC0FFEE）\n"
                << "  --case msb1|msb0  默认 msb1：m=63(111111)；msb0：m=31(011111)\n"
                << "  --m M             指定明文 m（覆盖 --case）\n";
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

bool decrypt_bit(const HEDB::TLWELvl1 &ct, const HEDB::TFHESecretKey &sk)
{
    return TFHEpp::tlweSymDecrypt<HEDB::Lvl1>(ct, sk.key.lvl1);
}

}  // namespace

int main(int argc, char **argv)
{
    try {
        const Args args = parse_args(argc, argv);
        constexpr uint32_t plain_bits = 6;
        constexpr uint32_t W = std::numeric_limits<HEDB::Lvl1::T>::digits;
        static_assert(W == 32, "This experiment assumes lvl1 torus is 32-bit.");

        uint64_t m = 0;
        if (args.m_override.has_value())
            m = *args.m_override;
        else if (args.case_mode == Args::CaseMode::msb1)
            m = (1ULL << plain_bits) - 1;  // 111111
        else
            m = (1ULL << (plain_bits - 1)) - 1;  // 011111

        const bool expected =
            (plain_bits == 0) ? false : (((m >> (plain_bits - 1)) & 1ULL) != 0);

        const uint32_t offset_bits = W - plain_bits - 1;
        const uint64_t offset = 1ULL << offset_bits;
        const uint64_t delta = 1ULL << (W - plain_bits);

        std::cout << "HE3DB experiments：MSB6(offset=Δ/2) 对比 ExtractMSB5\n";
        std::cout << "plain_bits=6，trials=" << args.trials << "，base_seed="
                  << args.base_seed << "\n";
        std::cout << "m=" << m << "（6位二进制=" << std::bitset<6>(m)
                  << "），expected_msb=" << (expected ? 1 : 0) << "\n";
        std::cout << "W=" << W << "，Δ=2^" << (W - plain_bits) << "=" << delta
                  << "，offset=Δ/2=2^" << offset_bits << "=" << offset << "\n";

        // Keygen (HE3DB aliases TFHEpp types).
        HEDB::TFHESecretKey sk;
        HEDB::TFHEEvalKey ek;
        ek.emplaceiksk<TFHEpp::lvl10param>(sk);
        ek.emplacebkfft<TFHEpp::lvl01param>(sk);

        // Sanity for polarity.
        {
            const uint64_t m0 = 0;
            const uint64_t m1 = 1ULL << (plain_bits - 1);  // 100000
            TFHEpp::generator.seed(0x1111ULL);
            const auto c0 = encrypt_int_lvl1(m0, plain_bits, sk);
            TFHEpp::generator.seed(0x2222ULL);
            const auto c1 = encrypt_int_lvl1(m1, plain_bits, sk);

            HEDB::TLWELvl1 o0_5, o1_5, o0_6, o1_6;
            HEDB::ExtractMSB5(o0_5, c0, ek, LOGIC);
            HEDB::ExtractMSB5(o1_5, c1, ek, LOGIC);
            HE3DBExperiment::ExtractMSB6(o0_6, c0, ek, LOGIC);
            HE3DBExperiment::ExtractMSB6(o1_6, c1, ek, LOGIC);

            std::cout << "sanity：m=0 -> ExtractMSB5="
                      << (decrypt_bit(o0_5, sk) ? 1 : 0)
                      << "，MSB6=" << (decrypt_bit(o0_6, sk) ? 1 : 0) << "\n";
            std::cout << "sanity：m=32 -> ExtractMSB5="
                      << (decrypt_bit(o1_5, sk) ? 1 : 0)
                      << "，MSB6=" << (decrypt_bit(o1_6, sk) ? 1 : 0)
                      << "（期望 1）\n";
        }

        uint32_t fails_msb5 = 0;
        uint32_t fails_msb6 = 0;
        uint64_t first_fail_seed_msb5 = 0;
        uint64_t first_fail_seed_msb6 = 0;
        bool recorded5 = false;
        bool recorded6 = false;

        for (uint32_t trial = 0; trial < args.trials; ++trial) {
            const uint64_t seed = args.base_seed + trial;
            TFHEpp::generator.seed(seed);
            const HEDB::TLWELvl1 ct = encrypt_int_lvl1(m, plain_bits, sk);

            HEDB::TLWELvl1 out5, out6;
            HEDB::ExtractMSB5(out5, ct, ek, LOGIC);
            HE3DBExperiment::ExtractMSB6(out6, ct, ek, LOGIC);
            const bool got5 = decrypt_bit(out5, sk);
            const bool got6 = decrypt_bit(out6, sk);

            if (got5 != expected) {
                ++fails_msb5;
                if (!recorded5) {
                    recorded5 = true;
                    first_fail_seed_msb5 = seed;
                }
            }
            if (got6 != expected) {
                ++fails_msb6;
                if (!recorded6) {
                    recorded6 = true;
                    first_fail_seed_msb6 = seed;
                }
            }
        }

        const double rate5 =
            static_cast<double>(fails_msb5) / static_cast<double>(args.trials);
        const double rate6 =
            static_cast<double>(fails_msb6) / static_cast<double>(args.trials);

        std::cout << "\n统计结果：\n";
        std::cout << "ExtractMSB5 失败数=" << fails_msb5 << " / " << args.trials
                  << "，失败率=" << rate5;
        if (recorded5) std::cout << "，首个失败seed=" << first_fail_seed_msb5;
        std::cout << "\n";

        std::cout << "MSB6(offset=Δ/2) 失败数=" << fails_msb6 << " / "
                  << args.trials << "，失败率=" << rate6;
        if (recorded6) std::cout << "，首个失败seed=" << first_fail_seed_msb6;
        std::cout << "\n";

        return 0;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

