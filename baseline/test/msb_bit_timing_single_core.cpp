#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"
#include "ethmsb_compare_adapter.h"

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace HEDB;

namespace
{
using Clock = std::chrono::steady_clock;

enum class Operation { GreaterThan, GreaterThanEqual, Equal, NotEqual };

struct Options {
    Operation operation = Operation::GreaterThan;
    std::uint32_t from_bit = 0;
    std::uint32_t to_bit = 32;
    std::uint32_t trials = 10;
    std::uint64_t seed = 20260625;
    std::string output = "results/msb_bit_timing_single_core_1_33.csv";
};

void configure_single_thread()
{
#ifdef _OPENMP
    omp_set_num_threads(1);
#endif
}

void make_eval_key(TFHEEvalKey &ek, TFHESecretKey &sk)
{
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);
}

std::uint64_t mask_bits(std::uint32_t bits)
{
    if (bits >= 64) return std::numeric_limits<std::uint64_t>::max();
    return (std::uint64_t{1} << bits) - 1;
}

template <class P>
TFHEpp::TLWE<P> encrypt_full_bits(std::uint64_t value, std::uint32_t full_bits,
                                  const TFHEpp::Key<P> &key)
{
    using T = typename P::T;
    constexpr std::uint32_t torus_bits = std::numeric_limits<T>::digits;
    if (full_bits == 0 || full_bits >= torus_bits)
        throw std::invalid_argument("invalid full_bits");

    const T phase = static_cast<T>(value) << (torus_bits - full_bits);
    return TFHEpp::tlweSymEncrypt<P>(phase, P::α, key);
}

bool decrypt_logic(const TLWELvl1 &cipher, const TFHESecretKey &sk)
{
    return TFHEpp::tlweSymDecrypt<Lvl1>(cipher, sk.key.lvl1);
}

template <class P>
struct TimingResult {
    double he3db_ms = 0.0;
    double ethmsb_ms = 0.0;
    std::uint64_t he3db_errors = 0;
    std::uint64_t ethmsb_errors = 0;
};

template <class P>
void he3db_compare(TFHEpp::TLWE<P> &lhs_cipher, TFHEpp::TLWE<P> &rhs_cipher,
                   TLWELvl1 &res, std::uint32_t bit, TFHEEvalKey &ek,
                   Operation operation)
{
    if (operation == Operation::GreaterThan)
    {
        greater_than<P>(lhs_cipher, rhs_cipher, res, bit, ek, LOGIC);
        return;
    }

    if (operation == Operation::GreaterThanEqual)
    {
        greater_than_equal<P>(lhs_cipher, rhs_cipher, res, bit, ek, LOGIC);
        return;
    }

    if (operation == Operation::Equal)
    {
        equal<P>(lhs_cipher, rhs_cipher, res, bit, ek, LOGIC);
        return;
    }

    TLWELvl1 greater_res{}, less_res{};
    greater_than<P>(lhs_cipher, rhs_cipher, greater_res, bit, ek, LOGIC);
    less_than<P>(lhs_cipher, rhs_cipher, less_res, bit, ek, LOGIC);
    HomOR(res, greater_res, less_res, ek, LOGIC);
}

template <class P>
void ethmsb_equal(TFHEpp::TLWE<P> &lhs_cipher, TFHEpp::TLWE<P> &rhs_cipher,
                  TLWELvl1 &res, std::uint32_t bit, TFHEEvalKey &ek)
{
    TLWELvl1 greater_equal_res{}, less_equal_res{};
    ethmsb_greater_than_equal<P>(lhs_cipher, rhs_cipher, greater_equal_res, bit,
                                 ek, LOGIC);
    ethmsb_less_than_equal<P>(lhs_cipher, rhs_cipher, less_equal_res, bit, ek,
                              LOGIC);
    HomAND(res, greater_equal_res, less_equal_res, ek, LOGIC);
}

template <class P>
void ethmsb_compare(TFHEpp::TLWE<P> &lhs_cipher, TFHEpp::TLWE<P> &rhs_cipher,
                    TLWELvl1 &res, std::uint32_t bit, TFHEEvalKey &ek,
                    Operation operation)
{
    if (operation == Operation::GreaterThan)
    {
        ethmsb_greater_than<P>(lhs_cipher, rhs_cipher, res, bit, ek, LOGIC);
        return;
    }

    if (operation == Operation::GreaterThanEqual)
    {
        ethmsb_greater_than_equal<P>(lhs_cipher, rhs_cipher, res, bit, ek,
                                     LOGIC);
        return;
    }

    if (operation == Operation::Equal)
    {
        ethmsb_equal<P>(lhs_cipher, rhs_cipher, res, bit, ek);
        return;
    }

    TLWELvl1 greater_res{}, less_res{};
    ethmsb_greater_than<P>(lhs_cipher, rhs_cipher, greater_res, bit, ek,
                           LOGIC);
    ethmsb_less_than<P>(lhs_cipher, rhs_cipher, less_res, bit, ek, LOGIC);
    HomOR(res, greater_res, less_res, ek, LOGIC);
}

bool expected_result(std::uint64_t lhs, std::uint64_t rhs,
                     Operation operation)
{
    if (operation == Operation::GreaterThan) return lhs > rhs;
    if (operation == Operation::GreaterThanEqual) return lhs >= rhs;
    if (operation == Operation::Equal) return lhs == rhs;
    return lhs != rhs;
}

template <class P>
TimingResult<P> run_comparison_bit(std::uint32_t bit, std::uint32_t trials,
                                   std::uint64_t seed, TFHESecretKey &sk,
                                   TFHEEvalKey &ek, Operation operation)
{
    std::mt19937_64 rng(seed ^ (std::uint64_t{0x9e3779b97f4a7c15ULL} * bit));
    std::uniform_int_distribution<std::uint64_t> dist(0, mask_bits(bit));
    const auto key = sk.key.get<P>();

    TLWELvl1 he3db_res{}, ethmsb_res{};
    TimingResult<P> result{};

    for (std::uint32_t trial = 0; trial < trials; ++trial)
    {
        const std::uint64_t lhs = dist(rng);
        const std::uint64_t rhs = dist(rng);
        const bool expected = expected_result(lhs, rhs, operation);
        auto lhs_cipher = encrypt_full_bits<P>(lhs, bit + 1, key);
        auto rhs_cipher = encrypt_full_bits<P>(rhs, bit + 1, key);

        if ((trial & 1U) == 0)
        {
            auto start = Clock::now();
            he3db_compare<P>(lhs_cipher, rhs_cipher, he3db_res, bit, ek,
                             operation);
            auto end = Clock::now();
            result.he3db_ms +=
                std::chrono::duration<double, std::milli>(end - start).count();

            start = Clock::now();
            ethmsb_compare<P>(lhs_cipher, rhs_cipher, ethmsb_res, bit, ek,
                              operation);
            end = Clock::now();
            result.ethmsb_ms +=
                std::chrono::duration<double, std::milli>(end - start).count();
        }
        else
        {
            auto start = Clock::now();
            ethmsb_compare<P>(lhs_cipher, rhs_cipher, ethmsb_res, bit, ek,
                              operation);
            auto end = Clock::now();
            result.ethmsb_ms +=
                std::chrono::duration<double, std::milli>(end - start).count();

            start = Clock::now();
            he3db_compare<P>(lhs_cipher, rhs_cipher, he3db_res, bit, ek,
                             operation);
            end = Clock::now();
            result.he3db_ms +=
                std::chrono::duration<double, std::milli>(end - start).count();
        }

        result.he3db_errors += decrypt_logic(he3db_res, sk) == expected ? 0 : 1;
        result.ethmsb_errors +=
            decrypt_logic(ethmsb_res, sk) == expected ? 0 : 1;
    }

    result.he3db_ms /= trials;
    result.ethmsb_ms /= trials;
    return result;
}

Operation parse_operation(const std::string &name)
{
    if (name == "gt" || name == "greater-than") return Operation::GreaterThan;
    if (name == "ge" || name == "greater-equal" ||
        name == "greater-than-equal")
        return Operation::GreaterThanEqual;
    if (name == "eq" || name == "equal") return Operation::Equal;
    if (name == "ne" || name == "neq" || name == "not-equal")
        return Operation::NotEqual;
    throw std::invalid_argument("unknown operation: " + name);
}

Options parse_args(int argc, char **argv)
{
    Options options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];
        if (arg == "--op" && i + 1 < argc)
            options.operation = parse_operation(argv[++i]);
        else if ((arg == "--from-bit" || arg == "--from") && i + 1 < argc)
            options.from_bit = static_cast<std::uint32_t>(
                std::stoul(argv[++i]));
        else if ((arg == "--to-bit" || arg == "--to") && i + 1 < argc)
            options.to_bit =
                static_cast<std::uint32_t>(std::stoul(argv[++i]));
        else if ((arg == "--trials" || arg == "-n") && i + 1 < argc)
            options.trials =
                static_cast<std::uint32_t>(std::stoul(argv[++i]));
        else if (arg == "--seed" && i + 1 < argc)
            options.seed = std::stoull(argv[++i]);
        else if ((arg == "--out" || arg == "-o") && i + 1 < argc)
            options.output = argv[++i];
        else if (arg == "--help" || arg == "-h")
        {
            std::cout
                << "Usage: msb_bit_timing_single_core [--op gt|ge|eq|ne] "
                   "[--from-bit N] [--to-bit N] [--trials N] [--out FILE]\n";
            std::exit(0);
        }
        else
            throw std::invalid_argument("unknown or incomplete argument: " +
                                        arg);
    }

    if (options.to_bit < options.from_bit || options.to_bit > 32 ||
        options.trials == 0)
        throw std::invalid_argument("invalid benchmark range or trial count");
    return options;
}
}  // namespace

int main(int argc, char **argv)
{
    try
    {
        const Options options = parse_args(argc, argv);
        TFHESecretKey sk;
        TFHEEvalKey ek;
        make_eval_key(ek, sk);

        std::filesystem::path out_path(options.output);
        if (out_path.has_parent_path())
            std::filesystem::create_directories(out_path.parent_path());
        std::ofstream out(out_path);
        if (!out) throw std::runtime_error("failed to open output csv");
        out << "bit,he3db_ms,ethmsb_ms,speedup\n";

        for (std::uint32_t bit = options.from_bit; bit <= options.to_bit; ++bit)
        {
            if (bit <= 9)
            {
                const auto r =
                    run_comparison_bit<Lvl1>(bit, options.trials, options.seed,
                                             sk, ek, options.operation);
                if (r.he3db_errors != 0 || r.ethmsb_errors != 0)
                    throw std::runtime_error("correctness failure at bit " +
                                             std::to_string(bit));
                out << bit << ',' << r.he3db_ms << ',' << r.ethmsb_ms << ','
                    << (r.he3db_ms / r.ethmsb_ms) << '\n';
                std::cout << "bit=" << bit << " he3db_ms=" << r.he3db_ms
                          << " ethmsb_ms=" << r.ethmsb_ms
                          << " speedup=" << (r.he3db_ms / r.ethmsb_ms)
                          << '\n';
            }
            else
            {
                const auto r =
                    run_comparison_bit<Lvl2>(bit, options.trials, options.seed,
                                             sk, ek, options.operation);
                if (r.he3db_errors != 0 || r.ethmsb_errors != 0)
                    throw std::runtime_error("correctness failure at bit " +
                                             std::to_string(bit));
                out << bit << ',' << r.he3db_ms << ',' << r.ethmsb_ms << ','
                    << (r.he3db_ms / r.ethmsb_ms) << '\n';
                std::cout << "bit=" << bit << " he3db_ms=" << r.he3db_ms
                          << " ethmsb_ms=" << r.ethmsb_ms
                          << " speedup=" << (r.he3db_ms / r.ethmsb_ms)
                          << '\n';
            }
        }
        std::cout << "wrote " << options.output << '\n';
    }
    catch (const std::exception &e)
    {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}
