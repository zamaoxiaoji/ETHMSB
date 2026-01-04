// Experiment: MSB wrap-around (carry-chain) in TFHEpp bootstrapping.
// Constraints: new standalone file; no modification to TFHEpp sources.

#include <tfhe++.hpp>

#include "experiments/tfhepp_instrumented_pbs.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace {

#ifndef USE_RANDEN
#error "This experiment requires -DUSE_RANDEN so TFHEpp::generator is seedable."
#endif

template <class T>
constexpr int torus_bits_v = std::numeric_limits<T>::digits;

template <class T>
using WideT =
    std::conditional_t<(torus_bits_v<T> <= 32), uint64_t, unsigned __int128>;

template <class P>
typename P::T tlwe_phase(const TFHEpp::TLWE<P> &ct, const TFHEpp::Key<P> &sk)
{
    using T = typename P::T;
    using Wide = WideT<T>;
    constexpr int W = torus_bits_v<T>;
    const Wide mask = (Wide{1} << W) - 1;

    Wide acc = static_cast<Wide>(ct[P::k * P::n]) & mask;
    for (size_t i = 0; i < P::k * P::n; ++i) {
        acc = (acc - (static_cast<Wide>(ct[i]) * static_cast<Wide>(sk[i]))) &
              mask;
    }
    return static_cast<T>(acc);
}

template <class T>
std::string torus_bits_with_sep(T x, uint32_t t)
{
    constexpr int W = torus_bits_v<T>;
    std::string s;
    s.reserve(static_cast<size_t>(W + 1));
    for (int i = W - 1; i >= 0; --i) {
        const bool bit = (static_cast<std::make_unsigned_t<T>>(x) >> i) & 1U;
        s.push_back(bit ? '1' : '0');
        if (t > 0 && t < static_cast<uint32_t>(W) &&
            i == static_cast<int>(W - t)) {
            s.push_back('|');
        }
    }
    return s;
}

template <class T>
uint64_t torus_to_u64(T x)
{
    return static_cast<uint64_t>(static_cast<std::make_unsigned_t<T>>(x));
}

template <class P>
TFHEpp::TLWE<P> tlwe_encrypt_int_scaled(uint64_t m, uint32_t shift, double alpha,
                                       const TFHEpp::Key<P> &sk)
{
    using T = typename P::T;
    constexpr int W = torus_bits_v<T>;
    if (shift >= static_cast<uint32_t>(W)) {
        throw std::runtime_error("shift must be < torus bit width");
    }

    const uint64_t scale = 1ULL << shift;
    const T center = static_cast<T>(m * scale);  // mod 2^W

    std::uniform_int_distribution<T> torus_dist(
        0, std::numeric_limits<T>::max());
    TFHEpp::TLWE<P> ct = {};
    ct[P::k * P::n] = TFHEpp::ModularGaussian<P>(center, alpha);
    for (size_t i = 0; i < P::k * P::n; ++i) {
        ct[i] = torus_dist(TFHEpp::generator);
        ct[P::k * P::n] += ct[i] * sk[i];
    }
    return ct;
}

struct Args {
    uint32_t trials = 2000;
    uint32_t t_min = 8;
    uint32_t t_max = 32;
    uint32_t fail_print_limit = 5;
    double alpha_in = 0.0;
    bool stop_on_first_t_with_failure = true;
    bool run_control = true;  // run m2=2^(t-1) and m3=3*2^(t-2)
};

std::optional<uint64_t> parse_u64(std::string_view s)
{
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

std::optional<double> parse_double(std::string_view s)
{
    try {
        size_t idx = 0;
        const double v = std::stod(std::string{s}, &idx);
        if (idx != s.size()) return std::nullopt;
        return v;
    }
    catch (...) {
        return std::nullopt;
    }
}

Args parse_args(int argc, char **argv)
{
    Args args;
    auto require_value = [&](int &i) -> std::string_view {
        if (i + 1 >= argc) {
            throw std::runtime_error("missing value for argument");
        }
        ++i;
        return argv[i];
    };

    for (int i = 1; i < argc; ++i) {
        const std::string_view a = argv[i];
        if (a == "--trials") {
            const auto v = parse_u64(require_value(i));
            if (!v) throw std::runtime_error("invalid --trials");
            args.trials = static_cast<uint32_t>(*v);
        }
        else if (a == "--t-min") {
            const auto v = parse_u64(require_value(i));
            if (!v) throw std::runtime_error("invalid --t-min");
            args.t_min = static_cast<uint32_t>(*v);
        }
        else if (a == "--t-max") {
            const auto v = parse_u64(require_value(i));
            if (!v) throw std::runtime_error("invalid --t-max");
            args.t_max = static_cast<uint32_t>(*v);
        }
        else if (a == "--fail-print-limit") {
            const auto v = parse_u64(require_value(i));
            if (!v) throw std::runtime_error("invalid --fail-print-limit");
            args.fail_print_limit = static_cast<uint32_t>(*v);
        }
        else if (a == "--alpha-in") {
            const auto v = parse_double(require_value(i));
            if (!v) throw std::runtime_error("invalid --alpha-in");
            args.alpha_in = *v;
        }
        else if (a == "--no-stop") {
            args.stop_on_first_t_with_failure = false;
        }
        else if (a == "--no-control") {
            args.run_control = false;
        }
        else if (a == "--help" || a == "-h") {
            std::cout
                << "Usage: msb_wraparound_experiment [options]\n"
                << "  --trials N            (default 2000)\n"
                << "  --t-min N             (default 8)\n"
                << "  --t-max N             (default 32)\n"
                << "  --alpha-in X          (default 0.0)\n"
                << "  --fail-print-limit K  (default 5)\n"
                << "  --no-stop             don't stop at first failing t\n"
                << "  --no-control          skip control m2=2^(t-1)\n";
            std::exit(0);
        }
        else {
            throw std::runtime_error("unknown argument: " + std::string(a));
        }
    }
    return args;
}

struct CaseResult {
    uint32_t failures = 0;
    uint32_t printed_failures = 0;
};

template <class DomainP, class TargetP, class BkP>
CaseResult run_case(const Args &args, const std::string_view label, uint32_t t,
                    uint64_t m, bool expected_msb, uint32_t shift,
                    const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek,
                    const TFHEpp::Polynomial<TargetP> &testvector)
{
    using Torus = typename DomainP::T;
    constexpr int W = torus_bits_v<Torus>;
    const uint64_t delta_u64 = (shift < 64) ? (1ULL << shift) : 0ULL;
    const Torus half = static_cast<Torus>(Torus{1} << (W - 1));

    const TFHEpp::Key<DomainP> *sk_debug = nullptr;
    if constexpr (std::is_same_v<DomainP, TFHEpp::lvl0param>)
        sk_debug = &sk.key.lvl0;
    else if constexpr (std::is_same_v<DomainP, TFHEpp::lvl1param>)
        sk_debug = &sk.key.lvl1;
    else if constexpr (std::is_same_v<DomainP, TFHEpp::lvl2param>)
        sk_debug = &sk.key.lvl2;

    CaseResult result{};
    for (uint32_t trial = 0; trial < args.trials; ++trial) {
        const uint64_t seed =
            0xC0FFEEULL + (static_cast<uint64_t>(t) << 16) + trial;
        TFHEpp::generator.seed(seed);

        const TFHEpp::TLWE<DomainP> ct_in =
            tlwe_encrypt_int_scaled<DomainP>(m, shift, args.alpha_in,
                                             sk.key.get<DomainP>());
        const Torus phase_in = tlwe_phase<DomainP>(ct_in, sk.key.get<DomainP>());
        const uint64_t m_in_decoded =
            (shift == 0) ? torus_to_u64(phase_in)
                         : (torus_to_u64(phase_in) >> shift);
        const bool msb_in_decoded =
            (t == 0) ? false : ((m_in_decoded >> (t - 1)) & 1ULL);

        TFHEpp::TLWE<TargetP> ct_out;
        TFHEppInstrumented::GateBootstrappingTLWE2TLWEFFTTrace<BkP>(
            ct_out, ct_in, *ek.bkfftlvl01, testvector);
        const bool got_msb =
            TFHEpp::tlweSymDecrypt<TargetP>(ct_out, sk.key.get<TargetP>());

        if (got_msb == expected_msb) continue;
        ++result.failures;

        if (result.printed_failures >= args.fail_print_limit) continue;
        ++result.printed_failures;

        const Torus ideal = static_cast<Torus>(m << shift);
        const Torus gap_to_wrap = static_cast<Torus>(Torus{0} - phase_in);
        const Torus dist_to_half =
            (phase_in >= half) ? (phase_in - half) : (half - phase_in);
        const double gap_over_delta =
            (delta_u64 == 0) ? 0.0
                             : (static_cast<double>(torus_to_u64(gap_to_wrap)) /
                                static_cast<double>(delta_u64));

        std::cout << "\n[FAIL] case=" << label << "\n";
        std::cout << "  trial=" << trial << " seed=" << seed << "\n";
        std::cout << "  W=" << W << " t=" << t << " m=" << m
                  << " expected_msb=" << (expected_msb ? 1 : 0)
                  << " got_msb=" << (got_msb ? 1 : 0) << "\n";
        std::cout << "  shift=" << shift << " Δ=2^(W-t)=" << delta_u64 << "\n";
        std::cout << "  Δm(ideal)= " << torus_bits_with_sep(ideal, t) << "\n";
        std::cout << "  phase_in = " << torus_bits_with_sep(phase_in, t) << "\n";
        std::cout << "  gap_to_wrap=(0-phase_in) = " << torus_to_u64(gap_to_wrap)
                  << " (gap/Δ=" << std::setprecision(6) << gap_over_delta
                  << ", gap_to_wrap<=Δ ? "
                  << ((torus_to_u64(gap_to_wrap) <= delta_u64) ? "yes" : "no")
                  << ")\n";
        std::cout << "  dist_to_half=|phase_in-2^(W-1)| = "
                  << torus_to_u64(dist_to_half) << "\n";
        std::cout << "  m_in_decoded=(phase_in>>shift) = " << m_in_decoded
                  << " msb(m_in_decoded)=" << (msb_in_decoded ? 1 : 0) << "\n";
        const uint64_t gap_u64 = torus_to_u64(gap_to_wrap);
        const uint64_t dist_u64 = torus_to_u64(dist_to_half);
        const uint64_t half_u64 = 1ULL << (W - 1);
        if (gap_u64 <= delta_u64 && dist_u64 >= (half_u64 / 4)) {
            std::cout
                << "  Note: dist_to_half is huge while gap_to_wrap≈Δ, so MSB\n"
                << "        flip is more consistent with wrap-around past 2^W\n"
                << "        (carry-chain 111..11 + small + -> 000..00).\n";
        }
        else if (dist_u64 <= delta_u64) {
            std::cout
                << "  Note: phase_in is very close to 2^(W-1), so MSB flip can\n"
                << "        also be explained by ordinary half-boundary crossing.\n";
        }
        else {
            std::cout << "  Note: failure observed; see gap_to_wrap/dist_to_half.\n";
        }

        // Failure replay (same seed) with trace enabled.
        {
            TFHEpp::generator.seed(seed);
            const TFHEpp::TLWE<DomainP> ct_in_re =
                tlwe_encrypt_int_scaled<DomainP>(m, shift, args.alpha_in,
                                                 sk.key.get<DomainP>());

            TFHEppInstrumented::TraceConfig cfg;
            cfg.enable = true;
            cfg.max_nonzero_a_print = 10;
            cfg.os = &std::cout;

            TFHEppInstrumented::TraceRecord rec;
            rec.t = t;
            rec.m = m;
            rec.seed = seed;

            TFHEpp::TLWE<TargetP> ct_out_re;
            TFHEppInstrumented::GateBootstrappingTLWE2TLWEFFTTrace<BkP>(
                ct_out_re, ct_in_re, *ek.bkfftlvl01, testvector, sk_debug, &cfg,
                &rec);
            const bool got_msb_re =
                TFHEpp::tlweSymDecrypt<TargetP>(ct_out_re, sk.key.get<TargetP>());
            if (got_msb_re != got_msb) {
                std::cout << "  [WARN] replay got_msb changed: got_msb_re="
                          << (got_msb_re ? 1 : 0) << "\n";
            }
        }
    }
    return result;
}

}  // namespace

int main(int argc, char **argv)
{
    try {
        const Args args = parse_args(argc, argv);

        using DomainP = TFHEpp::lvl0param;
        using TargetP = TFHEpp::lvl1param;
        using BkP = TFHEpp::lvl01param;
        using Torus = typename DomainP::T;
        constexpr int W = torus_bits_v<Torus>;

        if (args.t_min < 1 || args.t_max > static_cast<uint32_t>(W) ||
            args.t_min > args.t_max) {
            throw std::runtime_error("invalid t range");
        }

        std::cout << "TFHEpp MSB wrap-around experiment\n";
        std::cout << "  Torus W=" << W << " (Domain=lvl0param uint32)\n";
        std::cout << "  trials=" << args.trials << " alpha_in=" << args.alpha_in
                  << " stop_on_first_t_with_failure="
                  << (args.stop_on_first_t_with_failure ? "yes" : "no") << "\n";

        TFHEpp::SecretKey sk;
        TFHEpp::EvalKey ek;
        ek.emplacebkfft<BkP>(sk);

        auto make_testvector = [](bool use_pos_mu) {
            TFHEpp::Polynomial<TargetP> tv;
            const typename TargetP::T mu = TargetP::μ;
            const typename TargetP::T v =
                use_pos_mu ? mu : static_cast<typename TargetP::T>(0) - mu;
            tv.fill(v);
            return tv;
        };

        auto run_msb_once = [&](uint32_t t, uint64_t m,
                                const TFHEpp::Polynomial<TargetP> &tv,
                                uint64_t seed) -> bool {
            const uint32_t shift = static_cast<uint32_t>(W) - t;
            TFHEpp::generator.seed(seed);
            const TFHEpp::TLWE<DomainP> ct_in =
                tlwe_encrypt_int_scaled<DomainP>(m, shift, 0.0,
                                                 sk.key.get<DomainP>());
            TFHEpp::TLWE<TargetP> ct_out;
            TFHEppInstrumented::GateBootstrappingTLWE2TLWEFFTTrace<BkP>(
                ct_out, ct_in, *ek.bkfftlvl01, tv);
            return TFHEpp::tlweSymDecrypt<TargetP>(ct_out, sk.key.get<TargetP>());
        };

        // Quick sanity check for PBS output polarity (use alpha=0.0).
        bool use_pos_mu = false;
        {
            constexpr uint32_t t_check = 8;
            const uint64_t seed0 = 0x1234ULL;
            const uint64_t seed1 = 0x5678ULL;
            const uint64_t m0 = 0;
            const uint64_t m_hi = 3ULL << (t_check - 2);  // far from half boundary

            const auto tv_neg = make_testvector(false);
            const bool out0_neg = run_msb_once(t_check, m0, tv_neg, seed0);
            const bool out_hi_neg = run_msb_once(t_check, m_hi, tv_neg, seed1);

            if (!out0_neg && out_hi_neg) {
                use_pos_mu = false;
            }
            else {
                const auto tv_pos = make_testvector(true);
                const bool out0_pos = run_msb_once(t_check, m0, tv_pos, seed0);
                const bool out_hi_pos = run_msb_once(t_check, m_hi, tv_pos, seed1);
                if (!out0_pos && out_hi_pos)
                    use_pos_mu = true;
                else
                    std::cout << "  [WARN] sanity check ambiguous; keep -mu\n";
            }

            std::cout << "  sanity(t=8,alpha=0): use_pos_mu="
                      << (use_pos_mu ? "yes" : "no") << "\n";
            const auto tv = make_testvector(use_pos_mu);
            const uint64_t m2 = 1ULL << (t_check - 1);
            const uint64_t mmax = (1ULL << t_check) - 1;
            std::cout << "    m=0 got=" << (run_msb_once(t_check, 0, tv, 0x1ULL) ? 1 : 0)
                      << " expected=0\n";
            std::cout << "    m=2^(t-1) got="
                      << (run_msb_once(t_check, m2, tv, 0x2ULL) ? 1 : 0)
                      << " expected=1 (boundary)\n";
            std::cout << "    m=2^t-1 got="
                      << (run_msb_once(t_check, mmax, tv, 0x3ULL) ? 1 : 0)
                      << " expected=1\n";
        }

        TFHEpp::Polynomial<TargetP> msb_testvector = make_testvector(use_pos_mu);

        for (uint32_t t = args.t_min; t <= args.t_max; ++t) {
            const uint32_t shift = static_cast<uint32_t>(W) - t;
            const uint64_t m_worst =
                (t == 64) ? std::numeric_limits<uint64_t>::max()
                          : ((1ULL << t) - 1);
            const bool expected_msb = true;

            std::cout << "\n[t=" << t << "]"
                      << " worst m=2^t-1=" << m_worst << " (expected msb=1)\n";

            const CaseResult worst = run_case<DomainP, TargetP, BkP>(
                args, "worst(m=2^t-1)", t, m_worst, expected_msb, shift, sk, ek,
                msb_testvector);
            std::cout << "  worst failure_rate=" << worst.failures << "/"
                      << args.trials << " = "
                      << std::setprecision(6)
                      << (static_cast<double>(worst.failures) /
                          static_cast<double>(args.trials))
                      << "\n";

            if (args.run_control) {
                const uint64_t m2 = (t == 0) ? 0ULL : (1ULL << (t - 1));
                std::cout << "  control m2=2^(t-1)=" << m2
                          << " (expected msb=1)\n";
                const CaseResult control = run_case<DomainP, TargetP, BkP>(
                    args, "control(m2=2^(t-1))", t, m2, true, shift, sk, ek,
                    msb_testvector);
                std::cout << "  control failure_rate=" << control.failures
                          << "/" << args.trials << " = " << std::setprecision(6)
                          << (static_cast<double>(control.failures) /
                              static_cast<double>(args.trials))
                          << "\n";

                if (t >= 2) {
                    const uint64_t m3 = 3ULL << (t - 2);  // 11..00 (far from wrap)
                    std::cout << "  control2 m3=3*2^(t-2)=" << m3
                              << " (expected msb=1)\n";
                    const CaseResult control2 = run_case<DomainP, TargetP, BkP>(
                        args, "control2(m3=3*2^(t-2))", t, m3, true, shift, sk,
                        ek, msb_testvector);
                    std::cout << "  control2 failure_rate=" << control2.failures
                              << "/" << args.trials << " = "
                              << std::setprecision(6)
                              << (static_cast<double>(control2.failures) /
                                  static_cast<double>(args.trials))
                              << "\n";
                }
            }

            if (args.stop_on_first_t_with_failure && worst.failures > 0) break;
        }

        return 0;
    }
    catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
