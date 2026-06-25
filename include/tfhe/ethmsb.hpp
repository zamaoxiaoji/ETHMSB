#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "cloudkey.hpp"
#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "tlwe.hpp"

namespace TFHEpp::ethmsb {

inline constexpr std::uint32_t kappa = 5;

enum class ResultEncoding { Logic, Arithmetic };

enum class OutputPolarity {
    UpperHalfIsOne,
    LowerHalfIsOne,
};

template <class P>
inline constexpr std::uint32_t torus_bits_v =
    std::numeric_limits<typename P::T>::digits;

template <class P>
inline constexpr std::size_t tlwe_size_v = P::k * P::n + 1;

template <class P>
inline constexpr std::size_t tlwe_b_index_v = P::k * P::n;

template <class P>
struct RoundParameters {
    std::uint32_t full_bits{};
    std::uint32_t active_bits{};
    std::uint32_t next_active_bits{};
    std::uint32_t guard_pos{};
    std::uint32_t previous_guard_pos{};
    std::uint32_t left_shift{};
    bool has_previous_guard{};
    typename P::T decision_offset{};
    typename P::T target_weight{};
    typename P::T lut_half_amplitude{};
};

struct ScheduleCounts {
    std::uint32_t guard_rounds{};
    std::uint32_t pbs_count{};
    std::uint32_t iks_count{};
};

template <class P>
Polynomial<P> constant_polynomial(typename P::T value)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    Polynomial<P> poly{};
    poly.fill(value);
    return poly;
}

template <class P>
Polynomial<P> signed_amplitude_lut(typename P::T lut_half_amplitude,
                                   OutputPolarity output_polarity)
{
    if (output_polarity == OutputPolarity::UpperHalfIsOne)
        return constant_polynomial<P>(typename P::T{} - lut_half_amplitude);
    return constant_polynomial<P>(lut_half_amplitude);
}

template <class P>
void validate_full_bits(std::uint32_t full_bits,
                        std::uint32_t max_supported_bits)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    constexpr std::uint32_t T = torus_bits_v<P>;
    if (full_bits == 0)
        throw std::invalid_argument("ETHMSB full_bits must be positive");
    if (full_bits > max_supported_bits)
        throw std::invalid_argument("ETHMSB full_bits exceeds supported range");
    if (full_bits >= T)
        throw std::invalid_argument("ETHMSB full_bits leaves no offset bit");
}

template <class P>
RoundParameters<P> make_round_parameters(std::uint32_t full_bits,
                                          std::uint32_t active_bits,
                                          bool has_previous_guard,
                                          std::uint32_t previous_guard_pos)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    constexpr std::uint32_t T = torus_bits_v<P>;
    if (full_bits == 0 || full_bits >= T)
        throw std::invalid_argument("ETHMSB full_bits is invalid");
    if (active_bits < kappa)
        throw std::invalid_argument("ETHMSB active_bits is below kappa");
    if (active_bits >= full_bits)
        throw std::invalid_argument(
            "ETHMSB active_bits must be smaller than full_bits");
    if (active_bits >= T)
        throw std::invalid_argument("ETHMSB active_bits leaves no offset bit");
    const std::uint32_t guard_pos = active_bits - 1;
    if (guard_pos >= 64)
        throw std::invalid_argument("ETHMSB guard_pos exceeds uint64 range");
    const std::uint32_t next_active_bits =
        std::min(active_bits + kappa, full_bits);
    const std::uint32_t step = next_active_bits - active_bits;
    if (step == 0 || step > kappa)
        throw std::invalid_argument("ETHMSB rolling step is invalid");
    if (has_previous_guard &&
        previous_guard_pos + kappa + 1 != active_bits)
        throw std::invalid_argument("ETHMSB rolling guard position is invalid");
    const std::uint32_t left_shift = full_bits - active_bits;
    if (left_shift >= T)
        throw std::invalid_argument("ETHMSB left_shift exceeds Torus width");

    const std::uint32_t base_shift = T - active_bits - 1;
    const std::uint32_t target_shift = T - step - 1;
    if (target_shift == 0 || target_shift >= T)
        throw std::invalid_argument("ETHMSB target weight shift is invalid");

    RoundParameters<P> params{};
    params.full_bits = full_bits;
    params.active_bits = active_bits;
    params.next_active_bits = next_active_bits;
    params.guard_pos = guard_pos;
    params.previous_guard_pos = previous_guard_pos;
    params.left_shift = left_shift;
    params.has_previous_guard = has_previous_guard;
    params.decision_offset = typename P::T{1} << base_shift;
    if (has_previous_guard) {
        const std::uint32_t guard_offset_shift =
            base_shift + previous_guard_pos;
        if (guard_offset_shift >= T)
            throw std::invalid_argument("ETHMSB guard offset overflows Torus");
        params.decision_offset += typename P::T{1} << guard_offset_shift;
    }
    params.target_weight = typename P::T{1} << target_shift;
    params.lut_half_amplitude = params.target_weight >> 1;
    return params;
}

template <class P>
typename P::T rolling_final_decision_offset(std::uint32_t full_bits,
                                            bool has_previous_guard,
                                            std::uint32_t previous_guard_pos)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    constexpr std::uint32_t T = torus_bits_v<P>;
    if (full_bits >= T)
        throw std::invalid_argument("ETHMSB final offset shift is invalid");
    const std::uint32_t base_shift = T - full_bits - 1;
    typename P::T offset = typename P::T{1} << base_shift;
    if (has_previous_guard) {
        const std::uint32_t guard_offset_shift =
            base_shift + previous_guard_pos;
        if (guard_offset_shift >= T)
            throw std::invalid_argument(
                "ETHMSB final guard offset overflows Torus");
        offset += typename P::T{1} << guard_offset_shift;
    }
    return offset;
}

template <class P>
void left_shift_copy(TLWE<P> &shifted, const TLWE<P> &state,
                     std::uint32_t left_shift)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    for (std::size_t i = 0; i < tlwe_size_v<P>; ++i)
        shifted[i] = state[i] << left_shift;
}

template <class P>
void subtract_in_place(TLWE<P> &state, const TLWE<P> &weighted_bit)
{
    for (std::size_t i = 0; i < tlwe_size_v<P>; ++i)
        state[i] -= weighted_bit[i];
}

inline void MSBPBSWithOffset(TLWE<lvl1param> &res,
                             const TLWE<lvl1param> &input,
                             lvl1param::T decision_offset,
                             lvl1param::T lut_half_amplitude,
                             bool add_output_recenter, const EvalKey &ek,
                             OutputPolarity output_polarity =
                                 OutputPolarity::UpperHalfIsOne)
{
    // decision_offset is the one input-side phase shift used for the MSB
    // decision.  output_recenter below only maps the signed LUT output to
    // 0/target_weight arithmetic form.
    TLWE<lvl1param> input_with_offset = input;
    input_with_offset[tlwe_b_index_v<lvl1param>] += decision_offset;
    TLWE<lvl0param> switched{};
    IdentityKeySwitch<lvl10param>(switched, input_with_offset,
                                  ek.getiksk<lvl10param>());
    GateBootstrappingTLWE2TLWE<lvl01param>(
        res, switched, ek.getbkfft<lvl01param>(),
        signed_amplitude_lut<lvl1param>(lut_half_amplitude, output_polarity));
    if (add_output_recenter)
        res[tlwe_b_index_v<lvl1param>] += lut_half_amplitude;
}

inline void MSBPBSWithOffset(TLWE<lvl2param> &res,
                             const TLWE<lvl2param> &input,
                             lvl2param::T decision_offset,
                             lvl2param::T lut_half_amplitude,
                             bool add_output_recenter, const EvalKey &ek,
                             OutputPolarity output_polarity =
                                 OutputPolarity::UpperHalfIsOne)
{
    // Keep the offset at the input level, then use the HE3DB-aligned
    // lvl2->lvl0->lvl2 PBS route.  No historical Torus offset is accumulated.
    TLWE<lvl2param> input_with_offset = input;
    input_with_offset[tlwe_b_index_v<lvl2param>] += decision_offset;
    TLWE<lvl0param> switched{};
    IdentityKeySwitch<lvl20param>(switched, input_with_offset,
                                  ek.getiksk<lvl20param>());
    GateBootstrappingTLWE2TLWE<lvl02param>(
        res, switched, ek.getbkfft<lvl02param>(),
        signed_amplitude_lut<lvl2param>(lut_half_amplitude, output_polarity));
    if (add_output_recenter)
        res[tlwe_b_index_v<lvl2param>] += lut_half_amplitude;
}

inline void MSBPBSWithOffset(TLWE<lvl1param> &res,
                             const TLWE<lvl2param> &input,
                             lvl2param::T decision_offset,
                             lvl1param::T lut_half_amplitude,
                             bool add_output_recenter, const EvalKey &ek,
                             OutputPolarity output_polarity =
                                 OutputPolarity::UpperHalfIsOne)
{
    // Final Lvl2->Lvl1 extraction uses the same single decision_offset rule;
    // subtract/update stages never add or remove offsets.
    TLWE<lvl2param> input_with_offset = input;
    input_with_offset[tlwe_b_index_v<lvl2param>] += decision_offset;
    TLWE<lvl0param> switched{};
    IdentityKeySwitch<lvl20param>(switched, input_with_offset,
                                  ek.getiksk<lvl20param>());
    GateBootstrappingTLWE2TLWE<lvl01param>(
        res, switched, ek.getbkfft<lvl01param>(),
        signed_amplitude_lut<lvl1param>(lut_half_amplitude, output_polarity));
    if (add_output_recenter)
        res[tlwe_b_index_v<lvl1param>] += lut_half_amplitude;
}

template <class P>
void ExtractGuardForNextRound(TLWE<P> &next_guard,
                              const TLWE<P> &original_input,
                              const TLWE<P> &previous_guard,
                              bool has_previous_guard,
                              std::uint32_t previous_guard_pos,
                              std::uint32_t full_bits,
                              std::uint32_t active_bits, const EvalKey &ek,
                              RoundParameters<P> *round = nullptr)
{
    const RoundParameters<P> params =
        make_round_parameters<P>(full_bits, active_bits, has_previous_guard,
                                 previous_guard_pos);
    TLWE<P> local_input{};
    left_shift_copy<P>(local_input, original_input, params.left_shift);
    if (has_previous_guard) subtract_in_place<P>(local_input, previous_guard);
    MSBPBSWithOffset(next_guard, local_input, params.decision_offset,
                     params.lut_half_amplitude, true, ek,
                     OutputPolarity::UpperHalfIsOne);
    if (round != nullptr) *round = params;
}

inline lvl1param::T final_output_amplitude(ResultEncoding result_encoding)
{
    if (result_encoding == ResultEncoding::Logic)
        return static_cast<lvl1param::T>(lvl1param::μ);
    return static_cast<lvl1param::T>(lvl1param::μ) << 1;
}

inline bool final_output_recenter(ResultEncoding result_encoding)
{
    return result_encoding == ResultEncoding::Arithmetic;
}

inline void FinalMSBWithOffset(TLWE<lvl1param> &res,
                               const TLWE<lvl1param> &state,
                               std::uint32_t full_bits,
                               bool has_previous_guard,
                               std::uint32_t previous_guard_pos,
                               const EvalKey &ek,
                               ResultEncoding result_encoding,
                               OutputPolarity output_polarity =
                                   OutputPolarity::UpperHalfIsOne)
{
    const lvl1param::T offset =
        rolling_final_decision_offset<lvl1param>(
            full_bits, has_previous_guard, previous_guard_pos);
    MSBPBSWithOffset(res, state, offset, final_output_amplitude(result_encoding),
                     final_output_recenter(result_encoding), ek,
                     output_polarity);
}

inline void FinalMSBWithOffset(TLWE<lvl1param> &res,
                               const TLWE<lvl2param> &state,
                               std::uint32_t full_bits,
                               bool has_previous_guard,
                               std::uint32_t previous_guard_pos,
                               const EvalKey &ek,
                               ResultEncoding result_encoding,
                               OutputPolarity output_polarity =
                                   OutputPolarity::UpperHalfIsOne)
{
    const lvl2param::T offset =
        rolling_final_decision_offset<lvl2param>(
            full_bits, has_previous_guard, previous_guard_pos);
    MSBPBSWithOffset(res, state, offset, final_output_amplitude(result_encoding),
                     final_output_recenter(result_encoding), ek,
                     output_polarity);
}

inline void ClassicMSBLikeHE3DB(TLWE<lvl1param> &res,
                                const TLWE<lvl1param> &input,
                                std::uint32_t full_bits, const EvalKey &ek,
                                ResultEncoding result_encoding,
                                OutputPolarity output_polarity)
{
    FinalMSBWithOffset(res, input, full_bits, false, 0, ek, result_encoding,
                       output_polarity);
}

inline void ClassicMSBLikeHE3DB(TLWE<lvl1param> &res,
                                const TLWE<lvl2param> &input,
                                std::uint32_t full_bits, const EvalKey &ek,
                                ResultEncoding result_encoding,
                                OutputPolarity output_polarity)
{
    TLWE<lvl1param> input_lvl1{};
    IdentityKeySwitch<lvl21param>(input_lvl1, input, ek.getiksk<lvl21param>());
    ClassicMSBLikeHE3DB(res, input_lvl1, full_bits, ek, result_encoding,
                        output_polarity);
}

inline void HomETHMSBWithPolarity(TLWE<lvl1param> &res,
                                  const TLWE<lvl1param> &input,
                                  std::uint32_t full_bits, const EvalKey &ek,
                                  ResultEncoding result_encoding,
                                  OutputPolarity output_polarity)
{
    validate_full_bits<lvl1param>(full_bits, 10);
    if (full_bits <= kappa) {
        ClassicMSBLikeHE3DB(res, input, full_bits, ek, result_encoding,
                            output_polarity);
        return;
    }
    TLWE<lvl1param> previous_guard{};
    bool has_previous_guard = false;
    std::uint32_t previous_guard_pos = 0;
    for (std::uint32_t active_bits = kappa; active_bits < full_bits;
         active_bits += kappa) {
        TLWE<lvl1param> next_guard{};
        ExtractGuardForNextRound<lvl1param>(
            next_guard, input, previous_guard, has_previous_guard,
            previous_guard_pos, full_bits, active_bits, ek);
        previous_guard = next_guard;
        previous_guard_pos = active_bits - 1;
        has_previous_guard = true;
    }
    TLWE<lvl1param> final_input = input;
    subtract_in_place<lvl1param>(final_input, previous_guard);
    FinalMSBWithOffset(res, final_input, full_bits, has_previous_guard,
                       previous_guard_pos, ek, result_encoding,
                       output_polarity);
}

inline void HomETHMSBWithPolarity(TLWE<lvl1param> &res,
                                  const TLWE<lvl2param> &input,
                                  std::uint32_t full_bits, const EvalKey &ek,
                                  ResultEncoding result_encoding,
                                  OutputPolarity output_polarity)
{
    validate_full_bits<lvl2param>(full_bits, 33);
    if (full_bits <= kappa) {
        ClassicMSBLikeHE3DB(res, input, full_bits, ek, result_encoding,
                            output_polarity);
        return;
    }
    TLWE<lvl2param> previous_guard{};
    bool has_previous_guard = false;
    std::uint32_t previous_guard_pos = 0;
    for (std::uint32_t active_bits = kappa; active_bits < full_bits;
         active_bits += kappa) {
        TLWE<lvl2param> next_guard{};
        ExtractGuardForNextRound<lvl2param>(
            next_guard, input, previous_guard, has_previous_guard,
            previous_guard_pos, full_bits, active_bits, ek);
        previous_guard = next_guard;
        previous_guard_pos = active_bits - 1;
        has_previous_guard = true;
    }
    TLWE<lvl2param> final_input = input;
    subtract_in_place<lvl2param>(final_input, previous_guard);
    FinalMSBWithOffset(res, final_input, full_bits, has_previous_guard,
                       previous_guard_pos, ek, result_encoding,
                       output_polarity);
}

// Returns the mathematical unsigned MSB:
// lower half -> false/0, upper half -> true/arithmetic-one.
// This polarity is opposite to HE3DB's original MSBGateBootstrapping wrapper.
inline void HomETHMSB(TLWE<lvl1param> &res, const TLWE<lvl1param> &input,
                      std::uint32_t full_bits, const EvalKey &ek,
                      ResultEncoding result_encoding)
{
    HomETHMSBWithPolarity(res, input, full_bits, ek, result_encoding,
                          OutputPolarity::UpperHalfIsOne);
}

inline void HomETHMSB(TLWE<lvl1param> &res, const TLWE<lvl2param> &input,
                      std::uint32_t full_bits, const EvalKey &ek,
                      ResultEncoding result_encoding)
{
    HomETHMSBWithPolarity(res, input, full_bits, ek, result_encoding,
                          OutputPolarity::UpperHalfIsOne);
}

// HE3DB-compatible polarity: lower half -> true/arithmetic-one,
// upper half -> false/0.
inline void HomETHNonNegative(TLWE<lvl1param> &res,
                              const TLWE<lvl1param> &input,
                              std::uint32_t full_bits, const EvalKey &ek,
                              ResultEncoding result_encoding)
{
    HomETHMSBWithPolarity(res, input, full_bits, ek, result_encoding,
                          OutputPolarity::LowerHalfIsOne);
}

inline void HomETHNonNegative(TLWE<lvl1param> &res,
                              const TLWE<lvl2param> &input,
                              std::uint32_t full_bits, const EvalKey &ek,
                              ResultEncoding result_encoding)
{
    HomETHMSBWithPolarity(res, input, full_bits, ek, result_encoding,
                          OutputPolarity::LowerHalfIsOne);
}

template <class P>
void HomGreaterThan(TLWE<lvl1param> &res, const TLWE<P> &lhs,
                    const TLWE<P> &rhs, std::uint32_t operand_bits,
                    const EvalKey &ek, ResultEncoding result_encoding)
{
    // Precondition: lhs and rhs encode unsigned operand_bits-bit values using
    // Delta = Q / 2^(operand_bits + 1).  The extra bit is reserved for the
    // two's-complement subtraction result, matching HE3DB comparison inputs.
    const std::uint32_t full_bits = operand_bits + 1;
    TLWE<P> diff{};
    for (std::size_t i = 0; i < tlwe_size_v<P>; ++i) diff[i] = rhs[i] - lhs[i];
    HomETHMSB(res, diff, full_bits, ek, result_encoding);
}

template <class P>
void HomLessThan(TLWE<lvl1param> &res, const TLWE<P> &lhs,
                 const TLWE<P> &rhs, std::uint32_t operand_bits,
                 const EvalKey &ek, ResultEncoding result_encoding)
{
    // Same encoding precondition as HomGreaterThan.
    const std::uint32_t full_bits = operand_bits + 1;
    TLWE<P> diff{};
    for (std::size_t i = 0; i < tlwe_size_v<P>; ++i) diff[i] = lhs[i] - rhs[i];
    HomETHMSB(res, diff, full_bits, ek, result_encoding);
}

inline ScheduleCounts schedule_counts(std::uint32_t full_bits)
{
    if (full_bits == 0)
        throw std::invalid_argument("ETHMSB full_bits must be positive");
    ScheduleCounts counts{};
    counts.guard_rounds = (full_bits - 1) / kappa;
    counts.pbs_count = counts.guard_rounds + 1;
    counts.iks_count = counts.pbs_count;
    return counts;
}

template <class P>
constexpr const char *level_name()
{
    if constexpr (std::is_same_v<P, lvl1param>)
        return "Lvl1";
    else if constexpr (std::is_same_v<P, lvl2param>)
        return "Lvl2";
    else if constexpr (std::is_same_v<P, lvlhalfparam>)
        return "LvlHalf";
    else if constexpr (std::is_same_v<P, lvl0param>)
        return "Lvl0";
    else
        return "Unknown";
}

inline const char *result_encoding_name(ResultEncoding result_encoding)
{
    return result_encoding == ResultEncoding::Logic ? "LOGIC" : "ARITHMETIC";
}

inline std::string parameter_header_name()
{
#if defined(USE_80BIT_SECURITY)
    return "params/CGGI16.hpp";
#elif defined(USE_COMPRESS)
    return "params/compress.hpp";
#elif defined(USE_CGGI19)
    return "params/CGGI19.hpp";
#elif defined(USE_CONCRETE)
    return "params/concrete.hpp";
#elif defined(USE_TFHE_RS)
    return "params/tfhe-rs.hpp";
#elif defined(USE_TERNARY)
    return "params/ternary.hpp";
#else
    return "params/128bit.hpp";
#endif
}

}  // namespace TFHEpp::ethmsb
