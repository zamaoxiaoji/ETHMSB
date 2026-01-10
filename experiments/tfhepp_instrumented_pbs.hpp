#pragma once

// Instrumented copies of TFHEpp bootstrapping internals for debugging.
// Hard constraint: we do NOT modify TFHEpp original sources; we copy functions
// into a new namespace and only add trace/record logic.

#include <tfhe++.hpp>

#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <ostream>
#include <type_traits>
#include <utility>
#include <vector>

namespace TFHEppInstrumented {

struct TraceConfig {
    bool enable = false;
    uint32_t max_nonzero_a_print = 10;
    std::ostream *os = &std::cout;
};

struct TraceRecord {
    uint32_t t = 0;
    uint64_t m = 0;
    uint64_t seed = 0;

    uint32_t bitwidth = 0;
    uint32_t shift_bs = 0;
    uint64_t mod = 0;
    uint64_t half = 0;
    uint64_t roundoffset = 0;

    uint32_t qb = 0;
    uint32_t bbar = 0;

    uint32_t nonzero_ai_q_total = 0;

    uint64_t sum_ai_s = 0;
    uint64_t q_phase_parts = 0;
    uint64_t q_phase_direct = 0;
    int64_t q_phase_diff_centered = 0;
    uint64_t rot_exponent = 0;

    struct NonzeroAEntry {
        uint32_t i = 0;
        uint32_t ai_q = 0;
        uint32_t s_i = 0;
    };
    std::vector<NonzeroAEntry> nonzero_a_samples;
};

namespace detail {

inline uint64_t mod_u64(int64_t x, uint64_t mod)
{
    const int64_t m = static_cast<int64_t>(mod);
    int64_t r = x % m;
    if (r < 0) r += m;
    return static_cast<uint64_t>(r);
}

inline int64_t centered_diff_u64(uint64_t a, uint64_t b, uint64_t mod)
{
    const uint64_t raw = (a >= b) ? (a - b) : (a + mod - b);
    const uint64_t half = mod / 2;
    if (raw <= half) return static_cast<int64_t>(raw);
    return -static_cast<int64_t>(mod - raw);
}

template <class DomainP>
typename DomainP::T tlwe_phase_exact(const TFHEpp::TLWE<DomainP> &tlwe,
                                     const TFHEpp::Key<DomainP> &sk)
{
    typename DomainP::T phase = tlwe[DomainP::k * DomainP::n];
    for (int i = 0; i < DomainP::k * DomainP::n; ++i)
        phase -= tlwe[i] * sk[i];
    return phase;
}

}  // namespace detail

// Compute the same quantized indices used by BlindRotate (qb/ai_q path),
// WITHOUT running any CMUX/FFT operations. This is useful for collecting
// statistics of modulus switching / carry(wrap-around) behavior efficiently.
//
// The definitions are consistent with TFHEpp::BlindRotate:
//   shift_bs = digits(domainT) - 1 - targetP::nbit + bitwidth
//   roundoffset = 1 << (digits(domainT) - 2 - targetP::nbit + bitwidth)
//   qb = ((b >> shift_bs) << bitwidth)
//   ai_q = (((a_i + roundoffset) >> shift_bs) << bitwidth)
//   sum_ai_s = Σ(ai_q * s_i) mod MOD, MOD=(2N)<<bitwidth
//   q_phase_parts  = (qb - sum_ai_s) mod MOD
//   q_phase_direct = ((phase_exact >> shift_bs) << bitwidth) mod MOD
template <class P, uint32_t num_out = 1>
TraceRecord QuantizePhaseRecord(const TFHEpp::TLWE<typename P::domainP> &tlwe,
                                const TFHEpp::Key<typename P::domainP> &sk,
                                uint32_t t, uint64_t m, uint64_t seed)
{
    TraceRecord rec;
    rec.t = t;
    rec.m = m;
    rec.seed = seed;

    constexpr uint32_t bitwidth = TFHEpp::bits_needed<num_out - 1>();
    constexpr uint32_t shift_bs =
        std::numeric_limits<typename P::domainP::T>::digits - 1 -
        P::targetP::nbit + bitwidth;

    const uint64_t mod = (2ULL * static_cast<uint64_t>(P::targetP::n))
                         << bitwidth;
    const uint64_t half = (static_cast<uint64_t>(P::targetP::n) << bitwidth);

    rec.bitwidth = bitwidth;
    rec.shift_bs = shift_bs;
    rec.mod = mod;
    rec.half = half;
    rec.roundoffset = static_cast<uint64_t>(
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                 P::targetP::nbit + bitwidth));

    const uint32_t qb = static_cast<uint32_t>(
        (tlwe[P::domainP::k * P::domainP::n] >> shift_bs) << bitwidth);
    const uint32_t bbar = 2 * P::targetP::n - qb;
    rec.qb = qb;
    rec.bbar = bbar;

    uint32_t nonzero_cnt = 0;
    uint64_t sum_ai_s = 0;
    for (int i = 0; i < P::domainP::k * P::domainP::n; ++i) {
        const uint32_t ai_q =
            static_cast<uint32_t>(((tlwe[i] + rec.roundoffset) >> shift_bs)
                                  << bitwidth);
        if (ai_q == 0) continue;
        ++nonzero_cnt;
        const uint32_t s_i = static_cast<uint32_t>(sk[i]);
        sum_ai_s = (sum_ai_s +
                    (static_cast<uint64_t>(ai_q) * static_cast<uint64_t>(s_i))) %
                   mod;
    }
    rec.nonzero_ai_q_total = nonzero_cnt;
    rec.sum_ai_s = sum_ai_s;

    const typename P::domainP::T phase_exact =
        detail::tlwe_phase_exact<typename P::domainP>(tlwe, sk);
    const uint64_t q_phase_direct = static_cast<uint64_t>(
        static_cast<uint32_t>((phase_exact >> shift_bs) << bitwidth)) % mod;
    const uint64_t q_phase_parts =
        detail::mod_u64(static_cast<int64_t>(qb) - static_cast<int64_t>(sum_ai_s),
                        mod);
    const uint64_t rot_exponent =
        (static_cast<uint64_t>(bbar) + sum_ai_s) % mod;
    const int64_t diff_centered =
        detail::centered_diff_u64(q_phase_parts, q_phase_direct, mod);

    rec.q_phase_direct = q_phase_direct;
    rec.q_phase_parts = q_phase_parts;
    rec.q_phase_diff_centered = diff_centered;
    rec.rot_exponent = rot_exponent;
    return rec;
}

// Copied from: thirdparty/TFHEpp/include/gatebootstrapping.hpp
// Original name: TFHEpp::BlindRotate (FFT + Polynomial testvector)
// We only add optional tracing/recording controlled by (cfg, rec).
template <class P, uint32_t num_out = 1>
void BlindRotateTrace(
    TFHEpp::TRLWE<typename P::targetP> &res,
    const TFHEpp::TLWE<typename P::domainP> &tlwe,
    const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
    const TFHEpp::Polynomial<typename P::targetP> &testvector,
    const TFHEpp::Key<typename P::domainP> *sk_debug = nullptr,
    TraceConfig *cfg = nullptr, TraceRecord *rec = nullptr)
{
    constexpr uint32_t bitwidth = TFHEpp::bits_needed<num_out - 1>();
    constexpr uint32_t shift_bs =
        std::numeric_limits<typename P::domainP::T>::digits - 1 -
        P::targetP::nbit + bitwidth;

    const uint32_t qb = static_cast<uint32_t>(
        (tlwe[P::domainP::k * P::domainP::n] >> shift_bs) << bitwidth);
    const uint32_t bbar = 2 * P::targetP::n - qb;

    res = {};
    TFHEpp::PolynomialMulByXai<typename P::targetP>(res[P::targetP::k],
                                                   testvector, bbar);

    const bool trace_enabled = (cfg != nullptr) && cfg->enable &&
                               (cfg->os != nullptr) && (sk_debug != nullptr) &&
                               (rec != nullptr);
    const uint32_t max_keep =
        (cfg != nullptr) ? cfg->max_nonzero_a_print : 0;
    const bool want_compute = (sk_debug != nullptr) && (trace_enabled || (rec != nullptr));

    // MOD/half live in the "rotation domain" of BlindRotate.
    const uint64_t mod = (2ULL * static_cast<uint64_t>(P::targetP::n))
                         << bitwidth;
    const uint64_t half = (static_cast<uint64_t>(P::targetP::n) << bitwidth);

    uint32_t nonzero_cnt = 0;
    uint64_t sum_ai_s = 0;
    std::vector<TraceRecord::NonzeroAEntry> samples;
    if (trace_enabled) samples.reserve(std::min<uint32_t>(max_keep, 64U));

    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);

        const uint32_t ai_q =
            static_cast<uint32_t>(((tlwe[i] + roundoffset) >> shift_bs)
                                  << bitwidth);

        if (ai_q != 0) {
            ++nonzero_cnt;
            if (want_compute) {
                const uint32_t s_i = static_cast<uint32_t>((*sk_debug)[i]);
                sum_ai_s = (sum_ai_s +
                            (static_cast<uint64_t>(ai_q) *
                             static_cast<uint64_t>(s_i))) %
                           mod;
                if (trace_enabled && samples.size() < max_keep)
                    samples.push_back({static_cast<uint32_t>(i), ai_q, s_i});
            }
        }

        if (ai_q == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
        TFHEpp::CMUXFFTwithPolynomialMulByXaiMinusOne<typename P::targetP>(
            res, bkfft[i], ai_q);
    }

    if (rec == nullptr) return;

    // Compute direct quantized phase (after exact phase computation).
    uint64_t q_phase_direct = 0;
    if (sk_debug != nullptr) {
        const typename P::domainP::T phase_exact =
            detail::tlwe_phase_exact<typename P::domainP>(tlwe, *sk_debug);
        q_phase_direct = static_cast<uint64_t>(
            static_cast<uint32_t>((phase_exact >> shift_bs) << bitwidth));
        q_phase_direct %= mod;
    }

    const uint64_t q_phase_parts =
        detail::mod_u64(static_cast<int64_t>(qb) - static_cast<int64_t>(sum_ai_s),
                        mod);
    const uint64_t rot_exponent =
        (static_cast<uint64_t>(bbar) + sum_ai_s) % mod;
    const int64_t diff_centered =
        detail::centered_diff_u64(q_phase_parts, q_phase_direct, mod);

    rec->bitwidth = bitwidth;
    rec->shift_bs = shift_bs;
    rec->mod = mod;
    rec->half = half;
    rec->roundoffset = static_cast<uint64_t>(
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                 P::targetP::nbit + bitwidth));
    rec->qb = qb;
    rec->bbar = bbar;
    rec->nonzero_ai_q_total = nonzero_cnt;
    rec->sum_ai_s = sum_ai_s;
    rec->q_phase_parts = q_phase_parts;
    rec->q_phase_direct = q_phase_direct;
    rec->q_phase_diff_centered = diff_centered;
    rec->rot_exponent = rot_exponent;
    if (trace_enabled) rec->nonzero_a_samples = samples;

    if (!trace_enabled) return;

    std::ostream &os = *cfg->os;
    os << "\n[BlindRotateTrace]\n";
    os << "  t=" << rec->t << " m=" << rec->m << " seed=" << rec->seed << "\n";
    os << "  domainW=" << std::numeric_limits<typename P::domainP::T>::digits
       << " targetN=" << P::targetP::n << " targetNbit=" << P::targetP::nbit
       << " bitwidth=" << bitwidth << "\n";
    os << "  shift_bs=" << shift_bs << " roundoffset=2^("
       << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
           P::targetP::nbit + bitwidth)
       << ") mod=" << mod << " half=" << half << "\n";
    os << "  qb=" << qb << " (msb=" << ((qb >= half) ? 1 : 0) << ")"
       << " bbar=" << bbar << " (msb=" << ((bbar >= half) ? 1 : 0) << ")\n";
    os << "  nonzero_ai_q_total=" << nonzero_cnt
       << " printed_nonzero_ai_q=" << samples.size() << "\n";
    os << "  sum_ai_s=" << sum_ai_s << "\n";
    os << "  q_phase_direct=" << q_phase_direct
       << " (msb=" << ((q_phase_direct >= half) ? 1 : 0) << ")\n";
    os << "  q_phase_parts =" << q_phase_parts
       << " (msb=" << ((q_phase_parts >= half) ? 1 : 0) << ")\n";
    os << "  diff_centered(q_parts-q_direct)=" << diff_centered << "\n";
    os << "  rot_exponent=(bbar+sum_ai_s) mod mod = " << rot_exponent << "\n";

    for (const auto &e : samples) {
        os << "    i=" << e.i << " ai_q=" << e.ai_q << " s_i=" << e.s_i
           << "\n";
    }

    const bool direct_msb = (q_phase_direct >= half);
    const bool parts_msb = (q_phase_parts >= half);
    if (direct_msb && !parts_msb && (q_phase_parts <= (half / 8))) {
        os << "  Conclusion: q_phase_direct in upper half (msb=1) but\n"
              "              q_phase_parts wrapped into low region (msb=0,\n"
              "              near 0) => quantized carry/wrap-around in\n"
              "              modulus switching (bbar/ai_q path).\n";
    }
    else {
        os << "  Conclusion: see msb(q_phase_direct) vs msb(q_phase_parts).\n";
    }
}

// Copied from: thirdparty/TFHEpp/src/gatebootstrapping.cpp
// Original name: TFHEpp::GateBootstrappingTLWE2TLWEFFT
template <class P>
void GateBootstrappingTLWE2TLWEFFTTrace(
    TFHEpp::TLWE<typename P::targetP> &res,
    const TFHEpp::TLWE<typename P::domainP> &tlwe,
    const TFHEpp::BootstrappingKeyFFT<P> &bkfft,
    const TFHEpp::Polynomial<typename P::targetP> &testvector,
    const TFHEpp::Key<typename P::domainP> *sk_debug = nullptr,
    TraceConfig *cfg = nullptr, TraceRecord *rec = nullptr)
{
    TFHEpp::TRLWE<typename P::targetP> acc;
    BlindRotateTrace<P>(acc, tlwe, bkfft, testvector, sk_debug, cfg, rec);
    TFHEpp::SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

}  // namespace TFHEppInstrumented
