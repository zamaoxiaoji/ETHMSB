#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../include/tfhe++.hpp"

#if defined(USE_80BIT_SECURITY) || defined(USE_CONCRETE) || \
    defined(USE_CGGI19) || defined(USE_TFHE_RS) || defined(USE_TERNARY) || \
    defined(USE_COMPRESS)
#error "ETHMSB artifact benchmark must use params/128bit.hpp"
#endif

namespace {

using Clock = std::chrono::steady_clock;

template <class P>
constexpr std::uint32_t torus_bits()
{
    return std::numeric_limits<typename P::T>::digits;
}

template <class P>
TFHEpp::TLWE<P> encrypt_mod(std::uint64_t value, std::uint32_t full_bits,
                            const TFHEpp::SecretKey &sk)
{
    using T = typename P::T;
    const T phase = static_cast<T>(value) << (torus_bits<P>() - full_bits);
    TFHEpp::TLWE<P> ct{};
    TFHEpp::tlweSymEncrypt<P>(ct, phase, P::α, sk.key.get<P>());
    return ct;
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
    using T = TFHEpp::lvl1param::T;
    const T phase = TFHEpp::tlweSymPhase<TFHEpp::lvl1param>(
        ct, sk.key.get<TFHEpp::lvl1param>());
    constexpr T half_delta = static_cast<T>(TFHEpp::lvl1param::μ) << 1;
    return ((phase + half_delta) >>
            (std::numeric_limits<T>::digits - 1)) &
           T{1};
}

void make_eval_key(TFHEpp::EvalKey &ek, const TFHEpp::SecretKey &sk)
{
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(sk);
    ek.emplaceiksk<TFHEpp::lvl21param>(sk);
    ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    ek.emplacebkfft<TFHEpp::lvl02param>(sk);
}

std::string shell_one_line(const char *cmd)
{
    std::array<char, 256> buffer{};
    std::string result;
    FILE *pipe = popen(cmd, "r");
    if (pipe == nullptr) return "unknown";
    while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe) !=
           nullptr)
        result += buffer.data();
    pclose(pipe);
    while (!result.empty() &&
           (result.back() == '\n' || result.back() == '\r'))
        result.pop_back();
    return result.empty() ? "unknown" : result;
}

std::string thread_count_string()
{
    const char *omp = std::getenv("OMP_NUM_THREADS");
    return omp == nullptr ? "environment-default" : omp;
}

struct Stats {
    double mean_ms{};
    double median_ms{};
    double stddev_ms{};
    double min_ms{};
    double max_ms{};
};

struct CountBreakdown {
    std::uint32_t pbs_l1{};
    std::uint32_t pbs_l2{};
    std::uint32_t iks_lvl1_to_lvl0{};
    std::uint32_t iks_lvl2_to_lvl0{};
    std::uint32_t iks_lvl2_to_lvl1{};
    const char *final_path{};
};

template <class P>
CountBreakdown count_breakdown(std::uint32_t full_bits)
{
    const auto counts = TFHEpp::ethmsb::schedule_counts(full_bits);
    CountBreakdown out{};
    if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
        out.pbs_l1 = counts.pbs_count;
        out.iks_lvl1_to_lvl0 = counts.iks_count;
        out.final_path = "Lvl1->Lvl0->Lvl1";
    }
    else {
        if (full_bits <= TFHEpp::ethmsb::kappa) {
            out.pbs_l1 = 1;
            out.iks_lvl1_to_lvl0 = 1;
            out.iks_lvl2_to_lvl1 = 1;
            out.final_path = "Lvl2->Lvl1->Lvl0->Lvl1";
        }
        else {
            out.pbs_l1 = 1;
            out.pbs_l2 = counts.guard_rounds;
            out.iks_lvl2_to_lvl0 = counts.iks_count;
            out.final_path = "Lvl2->Lvl0->Lvl1";
        }
    }
    return out;
}

Stats summarize(std::vector<double> values)
{
    Stats s{};
    if (values.empty()) return s;
    std::sort(values.begin(), values.end());
    s.min_ms = values.front();
    s.max_ms = values.back();
    s.median_ms = values[values.size() / 2];
    s.mean_ms =
        std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double variance = 0.0;
    for (double v : values) variance += (v - s.mean_ms) * (v - s.mean_ms);
    s.stddev_ms = std::sqrt(variance / values.size());
    return s;
}

template <class P>
struct EncryptedPair {
    TFHEpp::TLWE<P> lhs;
    TFHEpp::TLWE<P> rhs;
    bool expected{};
};

template <class P>
std::vector<EncryptedPair<P>> prepare_pairs(std::uint32_t operand_bits,
                                            std::uint32_t trials,
                                            std::uint64_t seed,
                                            const TFHEpp::SecretKey &sk)
{
    const std::uint32_t full_bits = operand_bits + 1;
    std::mt19937_64 rng(seed ^ (std::uint64_t{0x9e3779b97f4a7c15} *
                                static_cast<std::uint64_t>(operand_bits)));
    const std::uint64_t max_value =
        (std::uint64_t{1} << operand_bits) - 1;
    std::uniform_int_distribution<std::uint64_t> dist(0, max_value);
    std::vector<EncryptedPair<P>> pairs;
    pairs.reserve(trials);
    for (std::uint32_t i = 0; i < trials; ++i) {
        const std::uint64_t lhs = dist(rng);
        const std::uint64_t rhs = dist(rng);
        pairs.push_back({encrypt_mod<P>(lhs, full_bits, sk),
                         encrypt_mod<P>(rhs, full_bits, sk), lhs > rhs});
    }
    return pairs;
}

template <class P>
void run_one(std::ofstream &csv, std::uint32_t operand_bits,
             TFHEpp::ethmsb::ResultEncoding result_encoding,
             std::uint32_t trials, std::uint64_t seed,
             const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek,
             const std::string &git_commit)
{
    const std::uint32_t full_bits = operand_bits + 1;
    auto pairs = prepare_pairs<P>(operand_bits, trials, seed, sk);
    std::vector<double> timings;
    timings.reserve(pairs.size());
    std::uint32_t failures = 0;

    for (const auto &pair : pairs) {
        TFHEpp::TLWE<TFHEpp::lvl1param> result{};
        const auto start = Clock::now();
        TFHEpp::ethmsb::HomGreaterThan<P>(
            result, pair.lhs, pair.rhs, operand_bits, ek, result_encoding);
        const auto end = Clock::now();
        timings.push_back(
            std::chrono::duration<double, std::milli>(end - start).count());
        const bool decoded =
            result_encoding == TFHEpp::ethmsb::ResultEncoding::Logic
                ? decrypt_logic(result, sk)
                : decrypt_arithmetic_bool(result, sk);
        if (decoded != pair.expected) ++failures;
    }

    const auto counts = TFHEpp::ethmsb::schedule_counts(full_bits);
    const auto detailed_counts = count_breakdown<P>(full_bits);
    const Stats stats = summarize(timings);
    csv << "ETHMSB-current-branch," << operand_bits << ',' << full_bits << ','
        << TFHEpp::ethmsb::level_name<P>() << ','
        << TFHEpp::ethmsb::level_name<P>() << ",Lvl1,"
        << TFHEpp::ethmsb::kappa << ','
        << TFHEpp::ethmsb::result_encoding_name(result_encoding) << ','
        << trials << ',' << seed << ',' << failures << ','
        << counts.pbs_count << ',' << counts.iks_count << ','
        << detailed_counts.pbs_l1 << ',' << detailed_counts.pbs_l2 << ','
        << detailed_counts.iks_lvl1_to_lvl0 << ','
        << detailed_counts.iks_lvl2_to_lvl0 << ','
        << detailed_counts.iks_lvl2_to_lvl1 << ','
        << detailed_counts.final_path << ','
        << stats.mean_ms << ',' << stats.median_ms << ',' << stats.stddev_ms
        << ',' << stats.min_ms << ',' << stats.max_ms << ',' << git_commit
        << ',' << TFHEpp::ethmsb::parameter_header_name() << ",\""
        << __VERSION__ << "\",\""
#ifdef NDEBUG
        << "Release"
#else
        << "Debug"
#endif
        << "\"," << thread_count_string() << '\n';
}

struct Cli {
    std::uint32_t trials = 5;
    std::uint64_t seed = 0x4554484d5342ULL;
    std::string output = "results/ethmsb_vs_he3db.csv";
};

Cli parse_cli(int argc, char **argv)
{
    Cli cli;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--trials" && i + 1 < argc)
            cli.trials = static_cast<std::uint32_t>(std::stoul(argv[++i]));
        else if (arg == "--seed" && i + 1 < argc)
            cli.seed = std::stoull(argv[++i]);
        else if (arg == "--out" && i + 1 < argc)
            cli.output = argv[++i];
        else {
            std::cerr << "Usage: ethmsb_bench [--trials N] [--seed S]"
                         " [--out results/ethmsb_vs_he3db.csv]\n";
            std::exit(2);
        }
    }
    return cli;
}

}  // namespace

int main(int argc, char **argv)
{
    const Cli cli = parse_cli(argc, argv);
    std::filesystem::create_directories(
        std::filesystem::path(cli.output).parent_path());
    const std::string git_commit = shell_one_line("git rev-parse HEAD");

    std::cout << "ETHMSB benchmark parameter header: "
              << TFHEpp::ethmsb::parameter_header_name() << '\n'
              << "TFHEpp commit: " << git_commit << '\n'
              << "thread_count: " << thread_count_string() << '\n'
              << "trials: " << cli.trials << " seed: " << cli.seed << '\n';

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    make_eval_key(ek, sk);

    std::ofstream csv(cli.output);
    csv << "scheme,operand_bits,full_bits,input_level,intermediate_level,"
           "output_level,kappa,result_type,trials,seed,failures,pbs_count,"
           "iks_count,pbs_l1_count,pbs_l2_count,iks_lvl1_to_lvl0_count,"
           "iks_lvl2_to_lvl0_count,iks_lvl2_to_lvl1_count,final_path,"
           "mean_ms,median_ms,stddev_ms,min_ms,max_ms,"
           "TFHEpp_commit,parameter_header,compiler,compile_flags,"
           "thread_count\n";

    for (std::uint32_t operand_bits : {4U, 8U}) {
        run_one<TFHEpp::lvl1param>(csv, operand_bits,
                                   TFHEpp::ethmsb::ResultEncoding::Logic,
                                   cli.trials, cli.seed, sk, ek, git_commit);
        run_one<TFHEpp::lvl1param>(
            csv, operand_bits, TFHEpp::ethmsb::ResultEncoding::Arithmetic,
            cli.trials, cli.seed, sk, ek, git_commit);
    }
    for (std::uint32_t operand_bits : {10U, 12U, 16U, 24U, 32U}) {
        run_one<TFHEpp::lvl2param>(csv, operand_bits,
                                   TFHEpp::ethmsb::ResultEncoding::Logic,
                                   cli.trials, cli.seed, sk, ek, git_commit);
        run_one<TFHEpp::lvl2param>(
            csv, operand_bits, TFHEpp::ethmsb::ResultEncoding::Arithmetic,
            cli.trials, cli.seed, sk, ek, git_commit);
    }

    std::cout << "wrote " << cli.output << std::endl;
}
