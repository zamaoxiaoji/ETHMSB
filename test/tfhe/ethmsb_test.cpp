#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <tfhe++.hpp>
#include <tfhe/ethmsb.hpp>

#if defined(USE_80BIT_SECURITY) || defined(USE_CONCRETE) || \
    defined(USE_CGGI19) || defined(USE_TFHE_RS) || defined(USE_TERNARY) || \
    defined(USE_COMPRESS)
#error "ETHMSB artifact tests must use params/128bit.hpp"
#endif

namespace {

template <class P>
constexpr std::uint32_t torus_bits()
{
    return std::numeric_limits<typename P::T>::digits;
}

std::uint64_t mask_bits(std::uint32_t bits)
{
    assert(bits > 0 && bits < 64);
    return (std::uint64_t{1} << bits) - 1;
}

template <class P>
TFHEpp::TLWE<P> encrypt_mod(std::uint64_t value, std::uint32_t full_bits,
                            const TFHEpp::SecretKey &sk)
{
    assert(full_bits > 0 && full_bits < torus_bits<P>());
    using T = typename P::T;
    const T phase = static_cast<T>(value) << (torus_bits<P>() - full_bits);
    TFHEpp::TLWE<P> ct{};
    TFHEpp::tlweSymEncrypt<P>(ct, phase, P::α, sk.key.get<P>());
    return ct;
}

template <class P>
std::uint64_t decrypt_mod(const TFHEpp::TLWE<P> &ct, std::uint32_t full_bits,
                          const TFHEpp::SecretKey &sk)
{
    assert(full_bits > 0 && full_bits < torus_bits<P>());
    using T = typename P::T;
    const T phase = TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
    const std::uint32_t shift = torus_bits<P>() - full_bits;
    const T rounded = phase + (T{1} << (shift - 1));
    return static_cast<std::uint64_t>(rounded >> shift) & mask_bits(full_bits);
}

bool decrypt_logic(const TFHEpp::TLWE<TFHEpp::lvl1param> &ct,
                   const TFHEpp::SecretKey &sk)
{
    return TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(
        ct, sk.key.get<TFHEpp::lvl1param>());
}

bool decrypt_arithmetic_bool(const TFHEpp::TLWE<TFHEpp::lvl1param> &ct,
                             const TFHEpp::SecretKey &sk)
{
    using P = TFHEpp::lvl1param;
    using T = P::T;
    const T phase = TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
    constexpr T half_delta = static_cast<T>(P::μ) << 1;
    constexpr std::uint32_t shift = std::numeric_limits<T>::digits - 1;
    return ((phase + half_delta) >> shift) & T{1};
}

template <class P>
double phase_error_units(const TFHEpp::TLWE<P> &ct, std::uint64_t expected,
                         std::uint32_t full_bits,
                         const TFHEpp::SecretKey &sk)
{
    using T = typename P::T;
    const T phase = TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
    const T expected_phase =
        static_cast<T>(expected) << (torus_bits<P>() - full_bits);
    const T diff = phase - expected_phase;
    using Signed = std::make_signed_t<T>;
    const auto signed_diff = static_cast<Signed>(diff);
    const long double delta =
        static_cast<long double>(T{1} << (torus_bits<P>() - full_bits));
    return static_cast<double>(static_cast<long double>(signed_diff) / delta);
}

[[noreturn]] void fail(const std::string &message)
{
    std::cerr << "ETHMSB test failure: " << message << std::endl;
    std::exit(1);
}

void test_plaintext_schedule()
{
    for (std::uint32_t full_bits = 1; full_bits <= 16; ++full_bits) {
        const std::uint64_t limit = std::uint64_t{1} << full_bits;
        for (std::uint64_t value = 0; value < limit; ++value) {
            bool has_previous_guard = false;
            std::uint32_t previous_guard_pos = 0;
            std::uint64_t previous_guard_bit = 0;
            const bool original_msb = (value >> (full_bits - 1)) & 1U;
            for (std::uint32_t active_bits = TFHEpp::ethmsb::kappa;
                 active_bits < full_bits; active_bits += TFHEpp::ethmsb::kappa) {
                const std::uint32_t guard_pos = active_bits - 1;
                std::uint64_t local = value & mask_bits(active_bits);
                if (has_previous_guard) {
                    local -= previous_guard_bit << previous_guard_pos;
                    if (((local >> previous_guard_pos) & 1U) != 0)
                        fail("rolling plaintext did not clear previous guard");
                }
                const std::uint64_t current_guard_bit =
                    (local >> guard_pos) & 1U;
                if (current_guard_bit != ((value >> guard_pos) & 1U))
                    fail("rolling plaintext changed current guard bit");
                if (((value >> (full_bits - 1)) & 1U) != original_msb)
                    fail("rolling plaintext changed global MSB");
                previous_guard_bit = current_guard_bit;
                previous_guard_pos = guard_pos;
                has_previous_guard = true;
            }
            std::uint64_t final_value = value;
            if (has_previous_guard)
                final_value -= previous_guard_bit << previous_guard_pos;
            const bool result = (final_value >> (full_bits - 1)) & 1U;
            if (result != original_msb)
                fail("rolling plaintext final MSB mismatch");
        }
    }

    std::mt19937_64 rng(0x4554484d5342554cULL);
    for (std::uint32_t full_bits = 17; full_bits <= 33; ++full_bits) {
        const std::uint64_t modulus = std::uint64_t{1} << full_bits;
        for (int i = 0; i < 200; ++i) {
            const std::uint64_t value = rng() & (modulus - 1);
            bool has_previous_guard = false;
            std::uint32_t previous_guard_pos = 0;
            std::uint64_t previous_guard_bit = 0;
            const bool original_msb = (value >> (full_bits - 1)) & 1U;
            for (std::uint32_t active_bits = TFHEpp::ethmsb::kappa;
                 active_bits < full_bits; active_bits += TFHEpp::ethmsb::kappa) {
                const std::uint32_t guard_pos = active_bits - 1;
                std::uint64_t local = value & mask_bits(active_bits);
                if (has_previous_guard)
                    local -= previous_guard_bit << previous_guard_pos;
                previous_guard_bit = (local >> guard_pos) & 1U;
                previous_guard_pos = guard_pos;
                has_previous_guard = true;
            }
            std::uint64_t final_value = value;
            if (has_previous_guard)
                final_value -= previous_guard_bit << previous_guard_pos;
            if (((final_value >> (full_bits - 1)) & 1U) != original_msb)
                fail("random rolling plaintext final MSB mismatch");
        }
    }
}

void test_offset_formulas()
{
    using L1 = TFHEpp::lvl1param;
    using L2 = TFHEpp::lvl2param;
    const auto round10 =
        TFHEpp::ethmsb::make_round_parameters<L1>(10, 5, false, 0);
    if (round10.decision_offset != (L1::T{1} << 26))
        fail("full_bits=10 first decision offset is not Q/64");
    const L1::T final10 =
        TFHEpp::ethmsb::rolling_final_decision_offset<L1>(10, true, 4);
    if (final10 != (static_cast<L1::T>(17) << 21))
        fail("full_bits=10 final offset is not 17Q/2048");
    const L1::T final15 =
        TFHEpp::ethmsb::rolling_final_decision_offset<L1>(15, true, 9);
    if (final15 != (static_cast<L1::T>(513) << 16))
        fail("rolling full_bits=15 final offset is not 513Q/2^16");

    const auto round33_5 =
        TFHEpp::ethmsb::make_round_parameters<L2>(33, 5, false, 0);
    if (round33_5.target_weight != (L2::T{1} << 58) ||
        round33_5.decision_offset != (L2::T{1} << 58))
        fail("full_bits=33 active=5 formula mismatch");

    const auto round33_10 =
        TFHEpp::ethmsb::make_round_parameters<L2>(33, 10, true, 4);
    if (round33_10.decision_offset !=
            (static_cast<L2::T>(17) << 53) ||
        round33_10.target_weight != (L2::T{1} << 58))
        fail("full_bits=33 active=10 guarded formula mismatch");

    const auto round33_30 =
        TFHEpp::ethmsb::make_round_parameters<L2>(33, 30, true, 24);
    if (round33_30.decision_offset !=
            ((L2::T{1} << 57) + (L2::T{1} << 33)) ||
        round33_30.target_weight != (L2::T{1} << 60))
        fail("full_bits=33 active=30 guarded formula mismatch");

    const L2::T final33 =
        TFHEpp::ethmsb::rolling_final_decision_offset<L2>(33, true, 29);
    if (final33 != ((L2::T{1} << 59) + (L2::T{1} << 30)))
        fail("full_bits=33 final rolling offset mismatch");
}

void make_eval_key(TFHEpp::EvalKey &ek, const TFHEpp::SecretKey &sk)
{
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(sk);
    ek.emplaceiksk<TFHEpp::lvl21param>(sk);
    ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    ek.emplacebkfft<TFHEpp::lvl02param>(sk);
}

template <class P>
void rolling_guard_sequence_case(const TFHEpp::SecretKey &sk,
                                 const TFHEpp::EvalKey &ek,
                                 std::uint32_t full_bits,
                                 std::uint64_t value)
{
    const auto input = encrypt_mod<P>(value, full_bits, sk);
    TFHEpp::TLWE<P> previous_guard{};
    bool has_previous_guard = false;
    std::uint32_t previous_guard_pos = 0;
    for (std::uint32_t active_bits = TFHEpp::ethmsb::kappa;
         active_bits < full_bits; active_bits += TFHEpp::ethmsb::kappa) {
        TFHEpp::TLWE<P> local_input{};
        TFHEpp::TLWE<P> next_guard{};
        TFHEpp::ethmsb::RoundParameters<P> params{};
        TFHEpp::ethmsb::ExtractGuardForNextRound<P>(
            next_guard, input, previous_guard, has_previous_guard,
            previous_guard_pos, full_bits, active_bits, ek, &params);

        TFHEpp::ethmsb::left_shift_copy<P>(local_input, input,
                                           params.left_shift);
        if (has_previous_guard)
            TFHEpp::ethmsb::subtract_in_place<P>(local_input, previous_guard);

        std::uint64_t expected_local = value & mask_bits(params.active_bits);
        if (has_previous_guard) {
            const std::uint64_t previous_bit =
                (value >> previous_guard_pos) & 1U;
            expected_local -= previous_bit << previous_guard_pos;
        }
        const std::uint64_t decoded_local =
            decrypt_mod<P>(local_input, params.active_bits, sk);
        const double local_error = phase_error_units<P>(
            local_input, expected_local, params.active_bits, sk);
        const long double local_half_gap =
            has_previous_guard
                ? (static_cast<long double>(std::uint64_t{1}
                                            << previous_guard_pos) +
                   1.0L) /
                      2.0L
                : 0.5L;
        if (!(std::abs(local_error) < local_half_gap)) {
            std::ostringstream os;
            os << "rolling local input mismatch full_bits=" << full_bits
               << " active_bits=" << active_bits
               << " next_active_bits=" << params.next_active_bits
               << " previous_guard_pos=" << previous_guard_pos
               << " has_previous_guard=" << has_previous_guard
               << " value=" << value << " expected=" << expected_local
               << " actual=" << decoded_local
               << " phase_error_delta=" << local_error
               << " required_half_gap=" << local_half_gap;
            fail(os.str());
        }

        const std::uint64_t expected_guard =
            ((value >> params.guard_pos) & 1U) << params.guard_pos;
        const std::uint64_t decoded_guard =
            decrypt_mod<P>(next_guard, params.next_active_bits, sk);
        const double guard_error = phase_error_units<P>(
            next_guard, expected_guard, params.next_active_bits, sk);
        const long double guard_half_gap =
            (static_cast<long double>(std::uint64_t{1}
                                      << params.guard_pos) +
             1.0L) /
            2.0L;
        if (!(std::abs(guard_error) < guard_half_gap)) {
            std::ostringstream os;
            os << "rolling guard mismatch full_bits=" << full_bits
               << " active_bits=" << active_bits
               << " next_active_bits=" << params.next_active_bits
               << " guard_pos=" << params.guard_pos
               << " previous_guard_pos=" << previous_guard_pos
               << " has_previous_guard=" << has_previous_guard
               << " decision_offset=" << params.decision_offset
               << " target_weight=" << params.target_weight
               << " value=" << value << " expected=" << expected_guard
               << " actual=" << decoded_guard
               << " phase_error_delta=" << guard_error
               << " required_half_gap=" << guard_half_gap;
            fail(os.str());
        }

        previous_guard = next_guard;
        previous_guard_pos = params.guard_pos;
        has_previous_guard = true;
    }

    if (has_previous_guard) {
        TFHEpp::TLWE<P> final_input = input;
        TFHEpp::ethmsb::subtract_in_place<P>(final_input, previous_guard);
        const std::uint64_t previous_bit =
            (value >> previous_guard_pos) & 1U;
        const std::uint64_t expected_final =
            value - (previous_bit << previous_guard_pos);
        const std::uint64_t decoded_final =
            decrypt_mod<P>(final_input, full_bits, sk);
        const double final_error =
            phase_error_units<P>(final_input, expected_final, full_bits, sk);
        const long double final_half_gap =
            (static_cast<long double>(std::uint64_t{1}
                                      << previous_guard_pos) +
             1.0L) /
            2.0L;
        if (!(std::abs(final_error) < final_half_gap)) {
            std::ostringstream os;
            os << "rolling final input MSB mismatch full_bits=" << full_bits
               << " previous_guard_pos=" << previous_guard_pos
               << " value=" << value << " expected=" << expected_final
               << " actual=" << decoded_final
               << " phase_error_delta=" << final_error
               << " required_half_gap=" << final_half_gap;
            fail(os.str());
        }
    }
}

void test_real_rolling_guards(const TFHEpp::SecretKey &sk,
                              const TFHEpp::EvalKey &ek)
{
    for (std::uint64_t value : {0ULL, 1ULL, 15ULL, 16ULL, 17ULL, 31ULL,
                               495ULL, 512ULL, 513ULL, 1023ULL}) {
        rolling_guard_sequence_case<TFHEpp::lvl1param>(sk, ek, 10, value);
    }

    const std::array<std::uint64_t, 8> lvl2_truth_values{
        0ULL, 15ULL, 495ULL, 512ULL, 513ULL, 1023ULL,
        (std::uint64_t{1} << 16) - 1, std::uint64_t{1} << 16};
    for (std::uint64_t value : lvl2_truth_values)
        rolling_guard_sequence_case<TFHEpp::lvl2param>(sk, ek, 17, value);

    for (std::uint64_t value : {(std::uint64_t{1} << 29) - 1,
                               std::uint64_t{1} << 29,
                               (std::uint64_t{1} << 29) + 1})
        rolling_guard_sequence_case<TFHEpp::lvl2param>(sk, ek, 30, value);

    const std::uint64_t all_guard_bits =
        (std::uint64_t{1} << 4) | (std::uint64_t{1} << 9) |
        (std::uint64_t{1} << 14) | (std::uint64_t{1} << 19) |
        (std::uint64_t{1} << 24) | (std::uint64_t{1} << 29);
    const std::array<std::uint64_t, 9> full33_values{
        0ULL,
        1ULL,
        (std::uint64_t{1} << 32) - 1,
        std::uint64_t{1} << 32,
        (std::uint64_t{1} << 32) + 1,
        (std::uint64_t{1} << 33) - 1,
        all_guard_bits,
        0x155555555ULL,
        0x0AAAAAAAAULL};
    for (std::uint64_t value : full33_values)
        rolling_guard_sequence_case<TFHEpp::lvl2param>(sk, ek, 33, value);
    std::mt19937_64 rng(0x4554484d53425233ULL);
    for (int i = 0; i < 8; ++i)
        rolling_guard_sequence_case<TFHEpp::lvl2param>(
            sk, ek, 33, rng() & ((std::uint64_t{1} << 33) - 1));
}

template <class P>
void direct_ethmsb_case(const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek,
                        std::uint32_t full_bits, std::uint64_t value,
                        TFHEpp::ethmsb::ResultEncoding result_encoding)
{
    const auto ct = encrypt_mod<P>(value, full_bits, sk);
    TFHEpp::TLWE<TFHEpp::lvl1param> msb{};
    TFHEpp::TLWE<TFHEpp::lvl1param> nonnegative{};
    TFHEpp::ethmsb::HomETHMSB(msb, ct, full_bits, ek, result_encoding);
    TFHEpp::ethmsb::HomETHNonNegative(nonnegative, ct, full_bits, ek,
                                      result_encoding);
    const bool expected_msb = (value >> (full_bits - 1)) & 1U;
    const bool decoded_msb =
        result_encoding == TFHEpp::ethmsb::ResultEncoding::Logic
            ? decrypt_logic(msb, sk)
            : decrypt_arithmetic_bool(msb, sk);
    const bool decoded_nonnegative =
        result_encoding == TFHEpp::ethmsb::ResultEncoding::Logic
            ? decrypt_logic(nonnegative, sk)
            : decrypt_arithmetic_bool(nonnegative, sk);
    if (decoded_msb != expected_msb || decoded_nonnegative == expected_msb) {
        const auto msb_phase = TFHEpp::tlweSymPhase<TFHEpp::lvl1param>(
            msb, sk.key.get<TFHEpp::lvl1param>());
        const auto nonnegative_phase =
            TFHEpp::tlweSymPhase<TFHEpp::lvl1param>(
                nonnegative, sk.key.get<TFHEpp::lvl1param>());
        std::ostringstream os;
        os << "direct ETHMSB polarity mismatch full_bits=" << full_bits
           << " value=" << value
           << " encoding="
           << TFHEpp::ethmsb::result_encoding_name(result_encoding)
           << " expected_msb=" << expected_msb
           << " decoded_msb=" << decoded_msb
           << " decoded_nonnegative=" << decoded_nonnegative
           << " msb_phase=" << msb_phase
           << " nonnegative_phase=" << nonnegative_phase;
        fail(os.str());
    }
}

void test_direct_ethmsb(const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek)
{
    const std::array<std::uint32_t, 15> bits{5,  6,  9,  10, 11,
                                             14, 15, 16, 19, 20,
                                             24, 25, 29, 30, 33};
    for (std::uint32_t full_bits : bits) {
        const std::array<std::uint64_t, 4> values{
            0ULL,
            (std::uint64_t{1} << (full_bits - 1)) - 1,
            std::uint64_t{1} << (full_bits - 1),
            (std::uint64_t{1} << full_bits) - 1};
        for (std::uint64_t value : values) {
            if (full_bits <= 10) {
                direct_ethmsb_case<TFHEpp::lvl1param>(
                    sk, ek, full_bits, value,
                    TFHEpp::ethmsb::ResultEncoding::Logic);
                direct_ethmsb_case<TFHEpp::lvl1param>(
                    sk, ek, full_bits, value,
                    TFHEpp::ethmsb::ResultEncoding::Arithmetic);
            }
            else {
                direct_ethmsb_case<TFHEpp::lvl2param>(
                    sk, ek, full_bits, value,
                    TFHEpp::ethmsb::ResultEncoding::Logic);
                direct_ethmsb_case<TFHEpp::lvl2param>(
                    sk, ek, full_bits, value,
                    TFHEpp::ethmsb::ResultEncoding::Arithmetic);
            }
        }
    }
}

template <class P>
void comparison_case(const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek,
                     std::uint32_t operand_bits, std::uint64_t lhs,
                     std::uint64_t rhs,
                     TFHEpp::ethmsb::ResultEncoding result_encoding)
{
    const std::uint32_t full_bits = operand_bits + 1;
    const auto clhs = encrypt_mod<P>(lhs, full_bits, sk);
    const auto crhs = encrypt_mod<P>(rhs, full_bits, sk);
    TFHEpp::TLWE<TFHEpp::lvl1param> result{};
    TFHEpp::ethmsb::HomGreaterThan<P>(result, clhs, crhs, operand_bits, ek,
                                      result_encoding);
    const bool decoded =
        result_encoding == TFHEpp::ethmsb::ResultEncoding::Logic
            ? decrypt_logic(result, sk)
            : decrypt_arithmetic_bool(result, sk);
    const bool expected = lhs > rhs;
    if (decoded != expected) {
        const auto result_phase = TFHEpp::tlweSymPhase<TFHEpp::lvl1param>(
            result, sk.key.get<TFHEpp::lvl1param>());
        std::ostringstream os;
        os << "comparison mismatch operand_bits=" << operand_bits
           << " full_bits=" << full_bits << " lhs=" << lhs << " rhs=" << rhs
           << " encoding="
           << TFHEpp::ethmsb::result_encoding_name(result_encoding)
           << " expected=" << expected << " actual=" << decoded
           << " result_phase=" << result_phase;
        fail(os.str());
    }
}

void test_real_comparisons(const TFHEpp::SecretKey &sk,
                           const TFHEpp::EvalKey &ek)
{
    const std::array<std::uint32_t, 7> operand_bits_list{4, 8, 10, 12,
                                                        16, 24, 32};
    for (const std::uint32_t operand_bits : operand_bits_list) {
        const std::uint64_t max_value =
            operand_bits == 64 ? std::numeric_limits<std::uint64_t>::max()
                               : (std::uint64_t{1} << operand_bits) - 1;
        std::vector<std::pair<std::uint64_t, std::uint64_t>> cases{
            {0, 0},
            {1, 0},
            {0, 1},
            {max_value, 0},
            {0, max_value},
            {max_value, max_value},
            {max_value, max_value >> 1},
            {0xAAAAAAAAULL & max_value, 0x55555555ULL & max_value},
        };
        for (auto [lhs, rhs] : cases) {
            if (operand_bits + 1 <= 10) {
                comparison_case<TFHEpp::lvl1param>(
                    sk, ek, operand_bits, lhs, rhs,
                    TFHEpp::ethmsb::ResultEncoding::Logic);
                comparison_case<TFHEpp::lvl1param>(
                    sk, ek, operand_bits, lhs, rhs,
                    TFHEpp::ethmsb::ResultEncoding::Arithmetic);
            }
            else {
                comparison_case<TFHEpp::lvl2param>(
                    sk, ek, operand_bits, lhs, rhs,
                    TFHEpp::ethmsb::ResultEncoding::Logic);
                comparison_case<TFHEpp::lvl2param>(
                    sk, ek, operand_bits, lhs, rhs,
                    TFHEpp::ethmsb::ResultEncoding::Arithmetic);
            }
        }
    }
}

}  // namespace

int main()
{
    std::cout << "ETHMSB parameter header: "
              << TFHEpp::ethmsb::parameter_header_name() << std::endl;
    std::cout << "Lvl0 n=" << TFHEpp::lvl0param::n
              << " T=" << torus_bits<TFHEpp::lvl0param>() << std::endl;
    std::cout << "Lvl1 n=" << TFHEpp::lvl1param::n
              << " T=" << torus_bits<TFHEpp::lvl1param>()
              << " Bgbit=" << TFHEpp::lvl1param::Bgbit
              << " l=" << TFHEpp::lvl1param::l << std::endl;
    std::cout << "Lvl2 n=" << TFHEpp::lvl2param::n
              << " T=" << torus_bits<TFHEpp::lvl2param>()
              << " Bgbit=" << TFHEpp::lvl2param::Bgbit
              << " l=" << TFHEpp::lvl2param::l << std::endl;
    std::cout << "KS lvl10 t=" << TFHEpp::lvl10param::t
              << " basebit=" << TFHEpp::lvl10param::basebit
              << " lvl20 t=" << TFHEpp::lvl20param::t
              << " basebit=" << TFHEpp::lvl20param::basebit
              << " lvl21 t=" << TFHEpp::lvl21param::t
              << " basebit=" << TFHEpp::lvl21param::basebit << std::endl;

    test_plaintext_schedule();
    test_offset_formulas();

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    make_eval_key(ek, sk);

    test_real_rolling_guards(sk, ek);
    test_direct_ethmsb(sk, ek);
    test_real_comparisons(sk, ek);

    std::cout << "ETHMSB tests passed" << std::endl;
}
