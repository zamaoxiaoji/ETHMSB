#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <tfhe/ethmsb.hpp>

#if defined(USE_80BIT_SECURITY) || defined(USE_CONCRETE) || \
    defined(USE_CGGI19) || defined(USE_TFHE_RS) || defined(USE_TERNARY) || \
    defined(USE_COMPRESS)
#error "ETHMSB correctness stress runner must use params/128bit.hpp"
#endif

namespace {

using Clock = std::chrono::steady_clock;

enum class Mode { Msb, GreaterThan };
enum class EncodingSelection { Logic, Arithmetic, Both };

struct Cli {
    Mode mode = Mode::Msb;
    EncodingSelection encodings = EncodingSelection::Both;
    std::uint32_t from_bit = 6;
    std::uint32_t to_bit = 33;
    std::uint64_t trials = std::uint64_t{1} << 20;
    std::uint64_t seed = 123;
    std::uint32_t threads = 0;
    std::uint32_t chunk_size = 16;
    std::uint64_t progress_interval = std::uint64_t{1} << 16;
    std::string output = "results/ethmsb_correctness_6_33.csv";
};

struct FirstFailure {
    bool present = false;
    std::uint64_t index = 0;
    std::uint64_t lhs = 0;
    std::uint64_t rhs = 0;
    std::uint64_t value = 0;
    bool expected = false;
    bool actual = false;
};

struct Stats {
    std::uint64_t failures = 0;
    FirstFailure first_failure{};
    double elapsed_seconds = 0.0;
};

struct TrialResult {
    bool expected = false;
    bool actual = false;
    std::uint64_t lhs = 0;
    std::uint64_t rhs = 0;
    std::uint64_t value = 0;
};

template <class P>
constexpr std::uint32_t torus_bits()
{
    return std::numeric_limits<typename P::T>::digits;
}

std::uint64_t mask_bits(std::uint32_t bits)
{
    if (bits >= 64) return std::numeric_limits<std::uint64_t>::max();
    return (std::uint64_t{1} << bits) - 1;
}

std::uint64_t splitmix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

template <class P>
TFHEpp::TLWE<P> encrypt_mod(std::uint64_t value, std::uint32_t full_bits,
                            const TFHEpp::Key<P> &key)
{
    using T = typename P::T;
    const T phase = static_cast<T>(value) << (torus_bits<P>() - full_bits);
    TFHEpp::TLWE<P> ct{};
    TFHEpp::tlweSymEncrypt<P>(ct, phase, P::α, key);
    return ct;
}

bool decrypt_logic(const TFHEpp::TLWE<TFHEpp::lvl1param> &ct,
                   const TFHEpp::Key<TFHEpp::lvl1param> &key)
{
    return TFHEpp::tlweSymDecrypt<TFHEpp::lvl1param>(ct, key);
}

bool decrypt_arithmetic_bool(const TFHEpp::TLWE<TFHEpp::lvl1param> &ct,
                             const TFHEpp::Key<TFHEpp::lvl1param> &key)
{
    using P = TFHEpp::lvl1param;
    using T = P::T;
    const T phase = TFHEpp::tlweSymPhase<P>(ct, key);
    constexpr T half_delta = static_cast<T>(P::μ) << 1;
    constexpr std::uint32_t shift = std::numeric_limits<T>::digits - 1;
    return ((phase + half_delta) >> shift) & T{1};
}

bool decode_result(const TFHEpp::TLWE<TFHEpp::lvl1param> &ct,
                   const TFHEpp::Key<TFHEpp::lvl1param> &key,
                   TFHEpp::ethmsb::ResultEncoding encoding)
{
    return encoding == TFHEpp::ethmsb::ResultEncoding::Logic
               ? decrypt_logic(ct, key)
               : decrypt_arithmetic_bool(ct, key);
}

void make_eval_key(TFHEpp::EvalKey &ek, const TFHEpp::SecretKey &sk)
{
    ek.emplaceiksk<TFHEpp::lvl10param>(sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(sk);
    ek.emplaceiksk<TFHEpp::lvl21param>(sk);
    ek.emplacebkfft<TFHEpp::lvl01param>(sk);
    ek.emplacebkfft<TFHEpp::lvl02param>(sk);
}

const char *mode_name(Mode mode)
{
    return mode == Mode::Msb ? "MSB" : "GREATER_THAN";
}

std::string compile_flags_name()
{
    std::string flags =
#ifdef NDEBUG
        "Release";
#else
        "Debug";
#endif
#ifdef USE_AVX512
    flags += " USE_AVX512=ON";
#else
    flags += " USE_AVX512=OFF";
#endif
#ifdef USE_SPQLIOS_INTL
    flags += " USE_SPQLIOS_INTL=ON";
#endif
    return flags;
}

std::uint32_t default_threads()
{
#ifdef _OPENMP
    return static_cast<std::uint32_t>(std::max(1, omp_get_num_procs()));
#else
    return 1;
#endif
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

template <class P>
TrialResult run_msb_trial(
    std::uint64_t trial_index, std::uint32_t full_bits, std::uint64_t seed,
    const TFHEpp::Key<P> &input_key,
    const TFHEpp::Key<TFHEpp::lvl1param> &result_key, const TFHEpp::EvalKey &ek,
    TFHEpp::ethmsb::ResultEncoding encoding)
{
    const std::uint64_t value =
        splitmix64(seed ^ (std::uint64_t{0x4d5342} << 32) ^
                   (static_cast<std::uint64_t>(full_bits) << 48) ^
                   trial_index) &
        mask_bits(full_bits);
    const auto ct = encrypt_mod<P>(value, full_bits, input_key);
    TFHEpp::TLWE<TFHEpp::lvl1param> result{};
    TFHEpp::ethmsb::HomETHMSB(result, ct, full_bits, ek, encoding);
    TrialResult out{};
    out.value = value;
    out.expected = ((value >> (full_bits - 1)) & 1U) != 0;
    out.actual = decode_result(result, result_key, encoding);
    return out;
}

template <class P>
TrialResult run_gt_trial(
    std::uint64_t trial_index, std::uint32_t operand_bits, std::uint64_t seed,
    const TFHEpp::Key<P> &input_key,
    const TFHEpp::Key<TFHEpp::lvl1param> &result_key, const TFHEpp::EvalKey &ek,
    TFHEpp::ethmsb::ResultEncoding encoding)
{
    const std::uint32_t full_bits = operand_bits + 1;
    const std::uint64_t mask = mask_bits(operand_bits);
    const std::uint64_t lhs =
        splitmix64(seed ^ (std::uint64_t{0x47544c} << 32) ^
                   (static_cast<std::uint64_t>(operand_bits) << 48) ^
                   (trial_index << 1)) &
        mask;
    const std::uint64_t rhs =
        splitmix64(seed ^ (std::uint64_t{0x475452} << 32) ^
                   (static_cast<std::uint64_t>(operand_bits) << 48) ^
                   ((trial_index << 1) | 1ULL)) &
        mask;
    const auto clhs = encrypt_mod<P>(lhs, full_bits, input_key);
    const auto crhs = encrypt_mod<P>(rhs, full_bits, input_key);
    TFHEpp::TLWE<TFHEpp::lvl1param> result{};
    TFHEpp::ethmsb::HomGreaterThan<P>(result, clhs, crhs, operand_bits, ek,
                                      encoding);
    TrialResult out{};
    out.lhs = lhs;
    out.rhs = rhs;
    out.expected = lhs > rhs;
    out.actual = decode_result(result, result_key, encoding);
    return out;
}

void merge_first_failure(FirstFailure &dst, const FirstFailure &src)
{
    if (!src.present) return;
    if (!dst.present || src.index < dst.index) dst = src;
}

template <class P>
Stats run_case(Mode mode, std::uint32_t bit, const Cli &cli,
               const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek,
               TFHEpp::ethmsb::ResultEncoding encoding)
{
    const auto input_key = sk.key.get<P>();
    const auto result_key = sk.key.get<TFHEpp::lvl1param>();
    Stats stats{};
    const auto started = Clock::now();
    std::atomic<std::uint64_t> completed{0};
    std::atomic<std::uint64_t> failures_seen{0};

#ifdef _OPENMP
    omp_set_num_threads(static_cast<int>(cli.threads));
#endif

#pragma omp parallel
    {
        std::uint64_t local_failures = 0;
        FirstFailure local_first{};

#pragma omp for schedule(dynamic, cli.chunk_size)
        for (std::int64_t trial = 0;
             trial < static_cast<std::int64_t>(cli.trials); ++trial) {
            const TrialResult result =
                mode == Mode::Msb
                    ? run_msb_trial<P>(static_cast<std::uint64_t>(trial), bit,
                                       cli.seed, input_key, result_key, ek,
                                       encoding)
                    : run_gt_trial<P>(static_cast<std::uint64_t>(trial), bit,
                                      cli.seed, input_key, result_key, ek,
                                      encoding);
            if (result.actual != result.expected) {
                ++local_failures;
                failures_seen.fetch_add(1);
                if (!local_first.present) {
                    local_first.present = true;
                    local_first.index = static_cast<std::uint64_t>(trial);
                    local_first.lhs = result.lhs;
                    local_first.rhs = result.rhs;
                    local_first.value = result.value;
                    local_first.expected = result.expected;
                    local_first.actual = result.actual;
                }
            }

            if (cli.progress_interval != 0) {
                const std::uint64_t done = completed.fetch_add(1) + 1;
                if (done % cli.progress_interval == 0 || done == cli.trials) {
#pragma omp critical(ethmsb_progress)
                    {
                        std::cerr << mode_name(mode) << " bit=" << bit
                                  << " "
                                  << TFHEpp::ethmsb::result_encoding_name(
                                         encoding)
                                  << " progress=" << done << '/'
                                  << cli.trials << " failures_so_far="
                                  << failures_seen.load() << '\n';
                    }
                }
            }
        }

#pragma omp critical(ethmsb_merge)
        {
            stats.failures += local_failures;
            merge_first_failure(stats.first_failure, local_first);
        }
    }

    const auto finished = Clock::now();
    stats.elapsed_seconds =
        std::chrono::duration<double>(finished - started).count();
    return stats;
}

std::vector<TFHEpp::ethmsb::ResultEncoding> selected_encodings(
    EncodingSelection selection)
{
    if (selection == EncodingSelection::Logic)
        return {TFHEpp::ethmsb::ResultEncoding::Logic};
    if (selection == EncodingSelection::Arithmetic)
        return {TFHEpp::ethmsb::ResultEncoding::Arithmetic};
    return {TFHEpp::ethmsb::ResultEncoding::Logic,
            TFHEpp::ethmsb::ResultEncoding::Arithmetic};
}

std::string usage()
{
    return "Usage: ethmsb_correctness [--mode msb|gt] [--from-bit N] "
           "[--to-bit N] [--trials N] [--seed N] [--threads N] "
           "[--chunk-size N] [--encoding logic|arithmetic|both] "
           "[--progress-interval N] [--out FILE]\n";
}

Cli parse_cli(int argc, char **argv)
{
    Cli cli;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto need_value = [&](const char *name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << '\n'
                          << usage();
                std::exit(2);
            }
            return argv[++i];
        };
        if (arg == "--mode") {
            const std::string value = need_value("--mode");
            if (value == "msb")
                cli.mode = Mode::Msb;
            else if (value == "gt" || value == "greater-than")
                cli.mode = Mode::GreaterThan;
            else {
                std::cerr << "Unknown mode: " << value << '\n' << usage();
                std::exit(2);
            }
        }
        else if (arg == "--from-bit")
            cli.from_bit = static_cast<std::uint32_t>(
                std::stoul(need_value("--from-bit")));
        else if (arg == "--to-bit")
            cli.to_bit =
                static_cast<std::uint32_t>(std::stoul(need_value("--to-bit")));
        else if (arg == "--trials")
            cli.trials = std::stoull(need_value("--trials"));
        else if (arg == "--seed")
            cli.seed = std::stoull(need_value("--seed"));
        else if (arg == "--threads")
            cli.threads = static_cast<std::uint32_t>(
                std::stoul(need_value("--threads")));
        else if (arg == "--chunk-size")
            cli.chunk_size = std::max<std::uint32_t>(
                1, static_cast<std::uint32_t>(
                       std::stoul(need_value("--chunk-size"))));
        else if (arg == "--progress-interval")
            cli.progress_interval =
                std::stoull(need_value("--progress-interval"));
        else if (arg == "--out")
            cli.output = need_value("--out");
        else if (arg == "--encoding") {
            const std::string value = need_value("--encoding");
            if (value == "logic")
                cli.encodings = EncodingSelection::Logic;
            else if (value == "arithmetic")
                cli.encodings = EncodingSelection::Arithmetic;
            else if (value == "both")
                cli.encodings = EncodingSelection::Both;
            else {
                std::cerr << "Unknown encoding: " << value << '\n'
                          << usage();
                std::exit(2);
            }
        }
        else if (arg == "--help" || arg == "-h") {
            std::cout << usage();
            std::exit(0);
        }
        else {
            std::cerr << "Unknown argument: " << arg << '\n' << usage();
            std::exit(2);
        }
    }

    if (cli.threads == 0) cli.threads = default_threads();
    if (cli.from_bit > cli.to_bit) {
        std::cerr << "--from-bit must be <= --to-bit\n" << usage();
        std::exit(2);
    }
    if (cli.mode == Mode::GreaterThan && cli.to_bit > 32) {
        std::cerr << "greater-than mode supports operand bits only up to 32 "
                     "(full_bits=operand_bits+1 <= 33)\n";
        std::exit(2);
    }
    return cli;
}

void write_header(std::ofstream &csv)
{
    csv << "mode,bit,operand_bits,full_bits,input_level,result_type,trials,"
           "failures,correct_rate,first_failure_index,first_lhs,first_rhs,"
           "first_value,first_expected,first_actual,seed,threads,chunk_size,"
           "elapsed_seconds,trials_per_second,parameter_header,tfhepp_commit,"
           "compiler,compile_flags\n";
}

template <class P>
void write_row(std::ofstream &csv, Mode mode, std::uint32_t bit,
               TFHEpp::ethmsb::ResultEncoding encoding, const Cli &cli,
               const Stats &stats, const std::string &git_commit)
{
    const std::uint32_t operand_bits =
        mode == Mode::Msb ? bit - 1 : bit;
    const std::uint32_t full_bits = mode == Mode::Msb ? bit : bit + 1;
    const long double correct_rate =
        static_cast<long double>(cli.trials - stats.failures) /
        static_cast<long double>(cli.trials);
    const double trials_per_second =
        stats.elapsed_seconds == 0.0
            ? 0.0
            : static_cast<double>(cli.trials) / stats.elapsed_seconds;
    csv << mode_name(mode) << ',' << bit << ',' << operand_bits << ','
        << full_bits << ',' << TFHEpp::ethmsb::level_name<P>() << ','
        << TFHEpp::ethmsb::result_encoding_name(encoding) << ','
        << cli.trials << ',' << stats.failures << ',' << std::setprecision(17)
        << static_cast<double>(correct_rate) << ',';
    if (stats.first_failure.present) {
        csv << stats.first_failure.index << ',' << stats.first_failure.lhs
            << ',' << stats.first_failure.rhs << ','
            << stats.first_failure.value << ','
            << stats.first_failure.expected << ','
            << stats.first_failure.actual << ',';
    }
    else {
        csv << "NA,NA,NA,NA,NA,NA,";
    }
    csv << cli.seed << ',' << cli.threads << ',' << cli.chunk_size << ','
        << std::setprecision(6) << stats.elapsed_seconds << ','
        << trials_per_second << ','
        << TFHEpp::ethmsb::parameter_header_name() << ',' << git_commit
        << ",\"" << __VERSION__ << "\",\"" << compile_flags_name() << "\"\n";
    csv.flush();
}

template <class P>
void run_and_record(std::ofstream &csv, Mode mode, std::uint32_t bit,
                    TFHEpp::ethmsb::ResultEncoding encoding, const Cli &cli,
                    const TFHEpp::SecretKey &sk, const TFHEpp::EvalKey &ek,
                    const std::string &git_commit)
{
    std::cerr << "Starting " << mode_name(mode) << " bit=" << bit << " "
              << TFHEpp::ethmsb::result_encoding_name(encoding)
              << " trials=" << cli.trials << " threads=" << cli.threads
              << '\n';
    const Stats stats = run_case<P>(mode, bit, cli, sk, ek, encoding);
    write_row<P>(csv, mode, bit, encoding, cli, stats, git_commit);
    const long double correct_rate =
        static_cast<long double>(cli.trials - stats.failures) /
        static_cast<long double>(cli.trials);
    std::cout << mode_name(mode) << " bit=" << bit << " "
              << TFHEpp::ethmsb::result_encoding_name(encoding)
              << " failures=" << stats.failures << '/' << cli.trials
              << " correct_rate=" << std::setprecision(17)
              << static_cast<double>(correct_rate)
              << " elapsed_s=" << std::setprecision(6)
              << stats.elapsed_seconds << '\n';
}

}  // namespace

int main(int argc, char **argv)
{
    const Cli cli = parse_cli(argc, argv);
    const std::filesystem::path out_path(cli.output);
    if (out_path.has_parent_path())
        std::filesystem::create_directories(out_path.parent_path());

    std::cout << "ETHMSB correctness stress runner\n"
              << "mode=" << mode_name(cli.mode)
              << " bits=" << cli.from_bit << ".." << cli.to_bit
              << " trials=" << cli.trials << " seed=" << cli.seed
              << " threads=" << cli.threads << '\n'
              << "parameter_header="
              << TFHEpp::ethmsb::parameter_header_name() << '\n'
              << "compile_flags=" << compile_flags_name() << '\n';

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    make_eval_key(ek, sk);

    const std::string git_commit = shell_one_line("git rev-parse HEAD");
    std::ofstream csv(cli.output);
    write_header(csv);

    const auto encodings = selected_encodings(cli.encodings);
    for (std::uint32_t bit = cli.from_bit; bit <= cli.to_bit; ++bit) {
        for (const auto encoding : encodings) {
            const std::uint32_t full_bits =
                cli.mode == Mode::Msb ? bit : bit + 1;
            if (full_bits <= 10) {
                run_and_record<TFHEpp::lvl1param>(csv, cli.mode, bit,
                                                  encoding, cli, sk, ek,
                                                  git_commit);
            }
            else {
                run_and_record<TFHEpp::lvl2param>(csv, cli.mode, bit,
                                                  encoding, cli, sk, ek,
                                                  git_commit);
            }
        }
    }

    std::cout << "wrote " << cli.output << std::endl;
}
