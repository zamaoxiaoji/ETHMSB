#pragma once

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include "HEDB/comparison/comparison.h"

namespace HEDB
{
namespace ethmsb_detail
{
inline constexpr std::uint32_t kappa = 5;

inline std::uint32_t next_segment_active_bits(std::uint32_t active_bits,
                                              std::uint32_t full_bits)
{
    return std::min(active_bits + kappa, full_bits);
}

enum class ResultEncoding { Logic, Arithmetic };

enum class OutputPolarity {
    UpperHalfIsOne,
    LowerHalfIsOne,
};

inline ResultEncoding result_encoding_from_he3db(bool result_type)
{
    return IS_LOGIC(result_type) ? ResultEncoding::Logic
                                 : ResultEncoding::Arithmetic;
}

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

template <class P>
TFHEpp::Polynomial<P> constant_polynomial(typename P::T value)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    TFHEpp::Polynomial<P> poly{};
    poly.fill(value);
    return poly;
}

template <class P>
TFHEpp::Polynomial<P> signed_amplitude_lut(
    typename P::T lut_half_amplitude, OutputPolarity output_polarity)
{
    if (output_polarity == OutputPolarity::UpperHalfIsOne)
        return constant_polynomial<P>(
            static_cast<typename P::T>(0) - lut_half_amplitude);
    return constant_polynomial<P>(lut_half_amplitude);
}

template <class P>
void validate_full_bits(std::uint32_t full_bits,
                        std::uint32_t max_supported_bits)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    constexpr std::uint32_t torus_bits = torus_bits_v<P>;
    if (full_bits == 0)
        throw std::invalid_argument("ETHMSB full_bits must be positive");
    if (full_bits > max_supported_bits)
        throw std::invalid_argument("ETHMSB full_bits exceeds supported range");
    if (full_bits >= torus_bits)
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
    constexpr std::uint32_t torus_bits = torus_bits_v<P>;
    if (full_bits == 0 || full_bits >= torus_bits)
        throw std::invalid_argument("ETHMSB full_bits is invalid");
    if (active_bits < kappa)
        throw std::invalid_argument("ETHMSB active_bits is below kappa");
    if (active_bits >= full_bits)
        throw std::invalid_argument(
            "ETHMSB active_bits must be smaller than full_bits");
    if (active_bits >= torus_bits)
        throw std::invalid_argument("ETHMSB active_bits leaves no offset bit");

    const std::uint32_t guard_pos = active_bits - 1;
    const std::uint32_t next_active_bits =
        next_segment_active_bits(active_bits, full_bits);
    const std::uint32_t step = next_active_bits - active_bits;
    if (step == 0 || step > kappa)
        throw std::invalid_argument("ETHMSB rolling step is invalid");
    if (has_previous_guard &&
        next_segment_active_bits(previous_guard_pos + 1, full_bits) !=
            active_bits)
        throw std::invalid_argument("ETHMSB rolling guard position is invalid");

    const std::uint32_t left_shift = full_bits - active_bits;
    const std::uint32_t base_shift = torus_bits - active_bits - 1;
    const std::uint32_t target_shift = torus_bits - step - 1;
    if (left_shift >= torus_bits || target_shift == 0 ||
        target_shift >= torus_bits)
        throw std::invalid_argument("ETHMSB shift is invalid");

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
        if (guard_offset_shift >= torus_bits)
            throw std::invalid_argument("ETHMSB guard offset overflows Torus");
        params.decision_offset += typename P::T{1} << guard_offset_shift;
    }
    params.target_weight = typename P::T{1} << target_shift;
    params.lut_half_amplitude = params.target_weight >> 1;
    return params;
}

template <class P>
typename P::T rolling_final_decision_offset(
    std::uint32_t full_bits, bool has_previous_guard,
    std::uint32_t previous_guard_pos)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    constexpr std::uint32_t torus_bits = torus_bits_v<P>;
    if (full_bits >= torus_bits)
        throw std::invalid_argument("ETHMSB final offset shift is invalid");
    const std::uint32_t base_shift = torus_bits - full_bits - 1;
    typename P::T offset = typename P::T{1} << base_shift;
    if (has_previous_guard) {
        const std::uint32_t guard_offset_shift =
            base_shift + previous_guard_pos;
        if (guard_offset_shift >= torus_bits)
            throw std::invalid_argument(
                "ETHMSB final guard offset overflows Torus");
        offset += typename P::T{1} << guard_offset_shift;
    }
    return offset;
}

template <class P>
void left_shift_copy(TFHEpp::TLWE<P> &shifted, const TFHEpp::TLWE<P> &state,
                     std::uint32_t left_shift)
{
    static_assert(std::is_unsigned_v<typename P::T>,
                  "ETHMSB requires an unsigned Torus type");
    for (std::size_t i = 0; i < tlwe_size_v<P>; ++i)
        shifted[i] = state[i] << left_shift;
}

template <class P>
void subtract_in_place(TFHEpp::TLWE<P> &state,
                       const TFHEpp::TLWE<P> &weighted_bit)
{
    for (std::size_t i = 0; i < tlwe_size_v<P>; ++i)
        state[i] -= weighted_bit[i];
}

inline void MSBPBSWithOffset(TLWELvl1 &res, const TLWELvl1 &input,
                             Lvl1::T decision_offset,
                             Lvl1::T lut_half_amplitude,
                             bool add_output_recenter,
                             const TFHEEvalKey &ek,
                             OutputPolarity output_polarity =
                                 OutputPolarity::UpperHalfIsOne)
{
    TLWELvl1 input_with_offset = input;
    input_with_offset[tlwe_b_index_v<Lvl1>] += decision_offset;
    TLWELvl0 switched{};
    TFHEpp::IdentityKeySwitch<Lvl10>(switched, input_with_offset,
                                     ek.getiksk<Lvl10>());
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
        res, switched, ek.getbkfft<Lvl01>(),
        signed_amplitude_lut<Lvl1>(lut_half_amplitude, output_polarity));
    if (add_output_recenter)
        res[tlwe_b_index_v<Lvl1>] += lut_half_amplitude;
}

inline void MSBPBSWithOffset(TLWELvl2 &res, const TLWELvl2 &input,
                             Lvl2::T decision_offset,
                             Lvl2::T lut_half_amplitude,
                             bool add_output_recenter,
                             const TFHEEvalKey &ek,
                             OutputPolarity output_polarity =
                                 OutputPolarity::UpperHalfIsOne)
{
    TLWELvl2 input_with_offset = input;
    input_with_offset[tlwe_b_index_v<Lvl2>] += decision_offset;
    TLWELvl0 switched{};
    TFHEpp::IdentityKeySwitch<Lvl20>(switched, input_with_offset,
                                     ek.getiksk<Lvl20>());
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl02>(
        res, switched, ek.getbkfft<Lvl02>(),
        signed_amplitude_lut<Lvl2>(lut_half_amplitude, output_polarity));
    if (add_output_recenter)
        res[tlwe_b_index_v<Lvl2>] += lut_half_amplitude;
}

inline void MSBPBSWithOffset(TLWELvl1 &res, const TLWELvl2 &input,
                             Lvl2::T decision_offset,
                             Lvl1::T lut_half_amplitude,
                             bool add_output_recenter,
                             const TFHEEvalKey &ek,
                             OutputPolarity output_polarity =
                                 OutputPolarity::UpperHalfIsOne)
{
    TLWELvl2 input_with_offset = input;
    input_with_offset[tlwe_b_index_v<Lvl2>] += decision_offset;
    TLWELvl0 switched{};
    TFHEpp::IdentityKeySwitch<Lvl20>(switched, input_with_offset,
                                     ek.getiksk<Lvl20>());
    TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(
        res, switched, ek.getbkfft<Lvl01>(),
        signed_amplitude_lut<Lvl1>(lut_half_amplitude, output_polarity));
    if (add_output_recenter)
        res[tlwe_b_index_v<Lvl1>] += lut_half_amplitude;
}

template <class P>
void ExtractGuardForNextRound(TFHEpp::TLWE<P> &next_guard,
                              const TFHEpp::TLWE<P> &original_input,
                              const TFHEpp::TLWE<P> &previous_guard,
                              bool has_previous_guard,
                              std::uint32_t previous_guard_pos,
                              std::uint32_t full_bits,
                              std::uint32_t active_bits,
                              const TFHEEvalKey &ek)
{
    const RoundParameters<P> params =
        make_round_parameters<P>(full_bits, active_bits, has_previous_guard,
                                 previous_guard_pos);
    TFHEpp::TLWE<P> local_input{};
    left_shift_copy<P>(local_input, original_input, params.left_shift);
    if (has_previous_guard)
        subtract_in_place<P>(local_input, previous_guard);
    MSBPBSWithOffset(next_guard, local_input, params.decision_offset,
                     params.lut_half_amplitude, true, ek,
                     OutputPolarity::UpperHalfIsOne);
}

inline Lvl1::T final_output_amplitude(ResultEncoding result_encoding)
{
    if (result_encoding == ResultEncoding::Logic)
        return static_cast<Lvl1::T>(Lvl1::μ);
    return static_cast<Lvl1::T>(Lvl1::μ) << 1;
}

inline bool final_output_recenter(ResultEncoding result_encoding)
{
    return result_encoding == ResultEncoding::Arithmetic;
}

inline bool he3db_result_type(ResultEncoding result_encoding)
{
    return result_encoding == ResultEncoding::Logic ? LOGIC : ARITHMETIC;
}

inline void FinalMSBWithOffset(TLWELvl1 &res, const TLWELvl1 &state,
                               std::uint32_t full_bits,
                               bool has_previous_guard,
                               std::uint32_t previous_guard_pos,
                               const TFHEEvalKey &ek,
                               ResultEncoding result_encoding,
                               OutputPolarity output_polarity =
                                   OutputPolarity::UpperHalfIsOne)
{
    const Lvl1::T offset = rolling_final_decision_offset<Lvl1>(
        full_bits, has_previous_guard, previous_guard_pos);
    MSBPBSWithOffset(res, state, offset, final_output_amplitude(result_encoding),
                     final_output_recenter(result_encoding), ek,
                     output_polarity);
}

inline void FinalMSBWithOffset(TLWELvl1 &res, const TLWELvl2 &state,
                               std::uint32_t full_bits,
                               bool has_previous_guard,
                               std::uint32_t previous_guard_pos,
                               const TFHEEvalKey &ek,
                               ResultEncoding result_encoding,
                               OutputPolarity output_polarity =
                                   OutputPolarity::UpperHalfIsOne)
{
    const Lvl2::T offset = rolling_final_decision_offset<Lvl2>(
        full_bits, has_previous_guard, previous_guard_pos);
    MSBPBSWithOffset(res, state, offset, final_output_amplitude(result_encoding),
                     final_output_recenter(result_encoding), ek,
                     output_polarity);
}

inline void BaseMSB(TLWELvl1 &res, const TLWELvl1 &input,
                    std::uint32_t full_bits, const TFHEEvalKey &ek,
                    ResultEncoding result_encoding,
                    OutputPolarity output_polarity)
{
    FinalMSBWithOffset(res, input, full_bits, false, 0, ek, result_encoding,
                       output_polarity);
}

inline void BaseMSB(TLWELvl1 &res, const TLWELvl2 &input,
                    std::uint32_t full_bits, const TFHEEvalKey &ek,
                    ResultEncoding result_encoding,
                    OutputPolarity output_polarity)
{
    TLWELvl1 input_lvl1{};
    TFHEpp::IdentityKeySwitch<Lvl21>(input_lvl1, input, ek.getiksk<Lvl21>());
    BaseMSB(res, input_lvl1, full_bits, ek, result_encoding, output_polarity);
}

inline void HomETHMSBWithPolarity(TLWELvl1 &res, const TLWELvl1 &input,
                                  std::uint32_t full_bits,
                                  const TFHEEvalKey &ek,
                                  ResultEncoding result_encoding,
                                  OutputPolarity output_polarity)
{
    validate_full_bits<Lvl1>(full_bits, 10);
    if (full_bits <= kappa) {
        BaseMSB(res, input, full_bits, ek, result_encoding, output_polarity);
        return;
    }
    TLWELvl1 previous_guard{};
    bool has_previous_guard = false;
    std::uint32_t previous_guard_pos = 0;
    for (std::uint32_t active_bits = kappa; active_bits < full_bits;
         active_bits = next_segment_active_bits(active_bits, full_bits)) {
        TLWELvl1 next_guard{};
        ExtractGuardForNextRound<Lvl1>(
            next_guard, input, previous_guard, has_previous_guard,
            previous_guard_pos, full_bits, active_bits, ek);
        previous_guard = next_guard;
        previous_guard_pos = active_bits - 1;
        has_previous_guard = true;
    }

    TLWELvl1 final_input = input;
    subtract_in_place<Lvl1>(final_input, previous_guard);
    FinalMSBWithOffset(res, final_input, full_bits, has_previous_guard,
                       previous_guard_pos, ek, result_encoding,
                       output_polarity);
}

inline void HomETHMSBWithPolarity(TLWELvl1 &res, const TLWELvl2 &input,
                                  std::uint32_t full_bits,
                                  const TFHEEvalKey &ek,
                                  ResultEncoding result_encoding,
                                  OutputPolarity output_polarity)
{
    validate_full_bits<Lvl2>(full_bits, 33);
    if (full_bits <= kappa) {
        BaseMSB(res, input, full_bits, ek, result_encoding, output_polarity);
        return;
    }

    TLWELvl2 previous_guard{};
    bool has_previous_guard = false;
    std::uint32_t previous_guard_pos = 0;
    for (std::uint32_t active_bits = kappa; active_bits < full_bits;
         active_bits = next_segment_active_bits(active_bits, full_bits)) {
        TLWELvl2 next_guard{};
        ExtractGuardForNextRound<Lvl2>(
            next_guard, input, previous_guard, has_previous_guard,
            previous_guard_pos, full_bits, active_bits, ek);
        previous_guard = next_guard;
        previous_guard_pos = active_bits - 1;
        has_previous_guard = true;
    }

    TLWELvl2 final_input = input;
    subtract_in_place<Lvl2>(final_input, previous_guard);
    FinalMSBWithOffset(res, final_input, full_bits, has_previous_guard,
                       previous_guard_pos, ek, result_encoding,
                       output_polarity);
}

inline void HomETHMSB(TLWELvl1 &res, const TLWELvl1 &input,
                      std::uint32_t full_bits, const TFHEEvalKey &ek,
                      ResultEncoding result_encoding)
{
    HomETHMSBWithPolarity(res, input, full_bits, ek, result_encoding,
                          OutputPolarity::UpperHalfIsOne);
}

inline void HomETHMSB(TLWELvl1 &res, const TLWELvl2 &input,
                      std::uint32_t full_bits, const TFHEEvalKey &ek,
                      ResultEncoding result_encoding)
{
    HomETHMSBWithPolarity(res, input, full_bits, ek, result_encoding,
                          OutputPolarity::UpperHalfIsOne);
}

}  // namespace ethmsb_detail

template <typename P>
void ethmsb_greater_than(TFHEpp::TLWE<P> &cipher1,
                         TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                         uint32_t plain_bits, TFHEEvalKey &ek,
                         bool result_type)
{
    TFHEpp::TLWE<P> sub_tlwe{};
    for (std::size_t i = 0; i < ethmsb_detail::tlwe_size_v<P>; ++i)
        sub_tlwe[i] = cipher2[i] - cipher1[i];
    ethmsb_detail::HomETHMSB(
        res, sub_tlwe, plain_bits + 1, ek,
        ethmsb_detail::result_encoding_from_he3db(result_type));
}

template <typename P>
void ethmsb_greater_than_equal(TFHEpp::TLWE<P> &cipher1,
                               TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                               uint32_t plain_bits, TFHEEvalKey &ek,
                               bool result_type)
{
    TFHEpp::TLWE<P> sub_tlwe{};
    for (std::size_t i = 0; i < ethmsb_detail::tlwe_size_v<P>; ++i)
        sub_tlwe[i] = cipher1[i] - cipher2[i];
    ethmsb_detail::HomETHMSB(
        res, sub_tlwe, plain_bits + 1, ek,
        ethmsb_detail::ResultEncoding::Logic);
    HomNOT<Lvl1>(res, res);
    if (IS_ARITHMETIC(result_type))
        TFHEpp::LOG_to_ARI(res, res, ek);
}

template <typename P>
void ethmsb_less_than(TFHEpp::TLWE<P> &cipher1,
                      TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                      uint32_t plain_bits, TFHEEvalKey &ek, bool result_type)
{
    TFHEpp::TLWE<P> sub_tlwe{};
    for (std::size_t i = 0; i < ethmsb_detail::tlwe_size_v<P>; ++i)
        sub_tlwe[i] = cipher1[i] - cipher2[i];
    ethmsb_detail::HomETHMSB(
        res, sub_tlwe, plain_bits + 1, ek,
        ethmsb_detail::result_encoding_from_he3db(result_type));
}

template <typename P>
void ethmsb_less_than_equal(TFHEpp::TLWE<P> &cipher1,
                            TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,
                            uint32_t plain_bits, TFHEEvalKey &ek,
                            bool result_type)
{
    TFHEpp::TLWE<P> sub_tlwe{};
    for (std::size_t i = 0; i < ethmsb_detail::tlwe_size_v<P>; ++i)
        sub_tlwe[i] = cipher2[i] - cipher1[i];
    ethmsb_detail::HomETHMSB(
        res, sub_tlwe, plain_bits + 1, ek,
        ethmsb_detail::ResultEncoding::Logic);
    HomNOT<Lvl1>(res, res);
    if (IS_ARITHMETIC(result_type))
        TFHEpp::LOG_to_ARI(res, res, ek);
}

}  // namespace HEDB
