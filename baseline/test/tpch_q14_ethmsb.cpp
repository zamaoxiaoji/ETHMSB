#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"
#include "HEDB/conversion/repack.h"
#include "ethmsb_compare_adapter.h"

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace HEDB;
using namespace std;

namespace
{
int default_threads()
{
#ifdef _OPENMP
    return omp_get_num_procs();
#else
    return 1;
#endif
}

bool is_power_of_two(size_t value)
{
    return value != 0 && (value & (value - 1)) == 0;
}

void configure_threads(int threads)
{
#ifdef _OPENMP
    omp_set_num_threads(threads);
#else
    (void)threads;
#endif
}

class ScopedOmpThreads {
public:
    explicit ScopedOmpThreads(int threads)
    {
#ifdef _OPENMP
        previous_ = omp_get_max_threads();
        omp_set_num_threads(threads > 0 ? threads : 1);
#else
        (void)threads;
#endif
    }

    ~ScopedOmpThreads()
    {
#ifdef _OPENMP
        omp_set_num_threads(previous_);
#endif
    }

private:
    int previous_ = 1;
};
}  // namespace

/*
    TPC-H Query 14
    select
        100.00 * sum(case
            when p_type like 'PROMO%'
                then l_extendedprice * (1 - l_discount)
            else 0
        end) / sum(l_extendedprice * (1 - l_discount)) as promo_revenue
    from
        lineitem,
        part
    where
        l_partkey = p_partkey
        and l_shipdate >= date ':1'
        and l_shipdate < date ':1' + interval '1' month;
    Consider the joined table
*/
bool relational_query14_ethmsb(size_t num, int threads, int agg_threads,
                               int pair_threads)
{
    if (!is_power_of_two(num))
        throw std::invalid_argument("rows must be a power of two for the rotate-sum aggregation");
    configure_threads(threads);
    std::cout << "Relational SQL Query14 ETHMSB Test: "<< std::endl;
    std::cout << "--------------------------------------------------------"<< std::endl;
    std::cout << "Records: " << num << std::endl;
    std::cout << "Threads: " << threads << std::endl;
    std::cout << "Aggregation Threads: " << (agg_threads > 0 ? agg_threads : 1) << std::endl;
    std::cout << "Pair Threads: " << (pair_threads > 0 ? pair_threads : 1) << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    TFHESecretKey sk;
    TFHEEvalKey ek;
    std::uniform_int_distribution<uint32_t> shipdate_message(10000, 20000);
    std::uniform_int_distribution<uint32_t> revenue_message(0, 100);
    std::uniform_int_distribution<uint32_t> ptype_message(0, 100);
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    // Filtering
    std::vector<uint64_t> ship_date(num);
    std::vector<uint32_t> ptype(num);
    std::vector<TLWELvl2> shipdate_ciphers(num);
    std::vector<TLWELvl1> ptype_ciphers(num);

    uint32_t ship_bits = 16, ptype_bits = 7, ship_scale_bits, ptype_scale_bits;
    ship_scale_bits = std::numeric_limits<Lvl2::T>::digits - ship_bits - 1;
    ptype_scale_bits = std::numeric_limits<Lvl1::T>::digits - ptype_bits - 1;

    TLWELvl2 predicate1_cipher, predicate2_cipher;
    TLWELvl1 predicate3_cipher, predicate4_cipher;
    uint64_t predicate1_value = 10592, predicate2_value = 10957;
    uint32_t predicate3_value = 30, predicate4_value = 70;
    predicate1_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl2>(predicate1_value, Lvl2::α, pow(2., ship_scale_bits), sk.key.get<Lvl2>());
    predicate2_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl2>(predicate2_value, Lvl2::α, pow(2., ship_scale_bits), sk.key.get<Lvl2>());
    predicate3_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(predicate3_value, Lvl1::α, pow(2., ptype_scale_bits), sk.key.get<Lvl1>());
    predicate4_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(predicate4_value, Lvl1::α, pow(2., ptype_scale_bits), sk.key.get<Lvl1>());


    // Start sql evaluation
    std::vector<TLWELvl1> filter_res(num), filter_case_res(num);

    std::vector<double> revenue(num);

    for (size_t i = 0; i < num; i++)
    {
        revenue[i] = revenue_message(engine);
    }

    for (size_t i = 0; i < num; i++)
    {
        // Generate data
        ship_date[i] = shipdate_message(engine);
        ptype[i] = ptype_message(engine);
    }
    if (num > 0)
    {
        ship_date[0] = 10600;
        ptype[0] = 50;
        revenue[0] = 100;
    }

#pragma omp parallel for num_threads(threads) schedule(static)
    for (long i = 0; i < static_cast<long>(num); i++)
    {
        shipdate_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(ship_date[i], Lvl2::α, pow(2., ship_scale_bits), sk.key.get<Lvl2>());
        ptype_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(ptype[i], Lvl1::α, pow(2., ptype_scale_bits), sk.key.get<Lvl1>());
    }

    std::chrono::system_clock::time_point start, end;
    double filtering_time = 0, aggregation_time;
    start = std::chrono::system_clock::now();

#pragma omp parallel for num_threads(threads) schedule(static)
    for (long i = 0; i < static_cast<long>(num); i++)
    {

        TLWELvl1 pre_res;
        ethmsb_greater_than<Lvl2>(shipdate_ciphers[i], predicate1_cipher, filter_res[i], ship_bits, ek, LOGIC);
        ethmsb_less_than<Lvl2>(shipdate_ciphers[i], predicate2_cipher, pre_res, ship_bits, ek, LOGIC);
        HomAND(filter_res[i], pre_res, filter_res[i], ek, LOGIC);
        ethmsb_greater_than<Lvl1>(ptype_ciphers[i], predicate3_cipher, pre_res, ptype_bits, ek, LOGIC);
        HomAND(filter_case_res[i], pre_res, filter_res[i], ek, LOGIC);
        ethmsb_less_than<Lvl1>(ptype_ciphers[i], predicate4_cipher, pre_res, ptype_bits, ek, LOGIC);
        HomAND(filter_case_res[i], pre_res, filter_case_res[i], ek, ARITHMETIC);
        HomAND(filter_res[i], filter_res[i], filter_res[i], ek, ARITHMETIC);

    }
    end = std::chrono::system_clock::now();

    filtering_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "Filtering Time: " << filtering_time << " ms" << std::endl;

    std::vector<uint64_t> plain_filter_res(num), plain_filter_case_res(num);
    double plain_agg_res = 0, plain_agg_case_res = 0;
#pragma omp parallel for num_threads(threads) schedule(static) reduction(+:plain_agg_res, plain_agg_case_res)
    for (long i = 0; i < static_cast<long>(num); i++)
    {
        if (ship_date[i] > predicate1_value && ship_date[i] < predicate2_value)
        {
            plain_filter_res[i] = 1;
            plain_agg_res += revenue[i];
            if (ptype[i] > predicate3_value && ptype[i] < predicate4_value)
            {
                plain_filter_case_res[i] = 1;
                plain_agg_case_res += revenue[i];
            }
            else
            {
                plain_filter_case_res[i] = 0;
            }

        }
        else
        {
            plain_filter_res[i] = 0;
            plain_filter_case_res[i] = 0;
        }

    }

    uint32_t rlwe_scale_bits = 29;
#pragma omp parallel for num_threads(threads) schedule(static)
    for (long i = 0; i < static_cast<long>(num); i++)
    {
        TFHEpp::ari_rescale(filter_res[i], filter_res[i], rlwe_scale_bits, ek);
        TFHEpp::ari_rescale(filter_case_res[i], filter_case_res[i], rlwe_scale_bits, ek);
    }

    std::cout << "Filtering finish" << std::endl;

    std::cout << "Aggregation :" << std::endl;
    uint64_t scale_bits = 29;
    uint64_t modq_bits = 32;
    uint64_t modulus_bits = 45;
    uint64_t repack_scale_bits = modulus_bits + scale_bits - modq_bits;
    uint64_t slots_count = filter_res.size();
    ScopedOmpThreads scoped_threads(agg_threads);
    const int effective_agg_threads = agg_threads > 0 ? agg_threads : 1;
    std::cout << "Generating Parameters..." << std::endl;
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    size_t poly_modulus_degree = 65536;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_modulus_degree, {59, 42, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59}));
    double scale = std::pow(2.0, scale_bits);

    //context instance
    seal::SEALContext context(parms, true, seal::sec_level_type::none);

    //key generation
    seal::KeyGenerator keygen(context);
    seal::SecretKey seal_secret_key = keygen.secret_key();
    seal::RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    seal::GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);


    //utils
    seal::Encryptor encryptor(context, seal_secret_key);
    seal::Evaluator evaluator(context);
    seal::Evaluator evaluator_case(context);
    seal::Decryptor decryptor(context, seal_secret_key);
    seal::Decryptor decryptor_case(context, seal_secret_key);

    //encoder
    seal::CKKSEncoder ckks_encoder(context);
    seal::CKKSEncoder ckks_encoder_case(context);



    //generate evaluation key
    std::cout << "Generating Conversion Key..." << std::endl;
    LTPreKey pre_key;
    LWEsToRLWEKeyGen(pre_key, std::pow(2., modulus_bits), seal_secret_key, sk, Lvl1::n, ckks_encoder, encryptor, context);


    // conversion
    std::cout << "Starting Conversion..." << std::endl;
    seal::Ciphertext result, result_case;
    const int paired_threads = std::min(2, std::max(1, pair_threads));
    start = std::chrono::system_clock::now();
#pragma omp parallel sections num_threads(paired_threads)
    {
#pragma omp section
        {
            LWEsToRLWE(result, filter_res, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder, galois_keys, relin_keys, evaluator, context);
            HomRound(result, result.scale(), ckks_encoder, relin_keys, evaluator, decryptor, context);
        }
#pragma omp section
        {
            LWEsToRLWE(result_case, filter_case_res, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits), ckks_encoder_case, galois_keys, relin_keys, evaluator_case, context);
            HomRound(result_case, result_case.scale(), ckks_encoder_case, relin_keys, evaluator_case, decryptor_case, context);
        }
    }
    end = std::chrono::system_clock::now();
    aggregation_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    seal::Plaintext plain, plain_case;
    std::vector<double> computed(slots_count), computed_case(slots_count);
#pragma omp parallel sections num_threads(paired_threads)
    {
#pragma omp section
        {
            decryptor.decrypt(result, plain);
            seal::pack_decode(computed, plain, ckks_encoder);
        }
#pragma omp section
        {
            decryptor_case.decrypt(result_case, plain_case);
            seal::pack_decode(computed_case, plain_case, ckks_encoder_case);
        }
    }

    double err1 = 0., err2 = 0.;

#pragma omp parallel for num_threads(threads) schedule(static) reduction(+:err1, err2)
    for (long i = 0; i < static_cast<long>(slots_count); ++i)
    {
        err1 += std::abs(computed[i] - plain_filter_res[i]);
        err2 += std::abs(computed_case[i] - plain_filter_case_res[i]);
    }

    const double repack_avg_error = err1 / slots_count;
    const double repack_case_avg_error = err2 / slots_count;
    printf("Repack average error = %f ~ 2^%.1f\n", repack_avg_error, std::log2(repack_avg_error));
    printf("Repack average error = %f ~ 2^%.1f\n", repack_case_avg_error, std::log2(repack_case_avg_error));


    // Filter result * data
    seal::Ciphertext revenue_cipher;
    double qd = parms.coeff_modulus()[result.coeff_modulus_size() - 1].value();
    seal::pack_encode(revenue, qd, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, revenue_cipher);

    std::cout << "Aggregating price and discount .." << std::endl;
    start = std::chrono::system_clock::now();
#pragma omp parallel sections num_threads(paired_threads)
    {
#pragma omp section
        {
            seal::Ciphertext local_revenue_cipher = revenue_cipher;
            seal::multiply_and_relinearize(result, local_revenue_cipher, result, evaluator, relin_keys);
            evaluator.rescale_to_next_inplace(result);
        }
#pragma omp section
        {
            seal::Ciphertext local_revenue_cipher = revenue_cipher;
            seal::multiply_and_relinearize(result_case, local_revenue_cipher, result_case, evaluator_case, relin_keys);
            evaluator_case.rescale_to_next_inplace(result_case);
        }
    }
    std::cout << "Remian modulus: " << result.coeff_modulus_size() << std::endl;
    int logrow = log2(num);

    auto rotate_sum = [&](seal::Ciphertext &cipher, seal::Evaluator &local_evaluator)
    {
        seal::Ciphertext temp;
        for (size_t i = 0; i < static_cast<size_t>(logrow); i++)
        {
            temp = cipher;
            size_t step = 1 << (logrow - i - 1);
            local_evaluator.rotate_vector_inplace(temp, step, galois_keys);
            local_evaluator.add_inplace(cipher, temp);
        }
    };
#pragma omp parallel sections num_threads(paired_threads)
    {
#pragma omp section
        {
            rotate_sum(result, evaluator);
        }
#pragma omp section
        {
            rotate_sum(result_case, evaluator_case);
        }
    }
    end = std::chrono::system_clock::now();
    aggregation_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::vector<double> agg_result(slots_count), agg_case_result(slots_count);
#pragma omp parallel sections num_threads(paired_threads)
    {
#pragma omp section
        {
            decryptor.decrypt(result, plain);
            seal::pack_decode(agg_result, plain, ckks_encoder);
        }
#pragma omp section
        {
            decryptor_case.decrypt(result_case, plain_case);
            seal::pack_decode(agg_case_result, plain_case, ckks_encoder_case);
        }
    }

    std::cout << "Aggregation Time: " << aggregation_time << " ms" << std::endl;
    std::cout << "Query Evaluation Time: " << filtering_time + aggregation_time << " ms" << std::endl;

    const double encrypted_ratio = agg_case_result[0] / agg_result[0];
    const double plain_ratio = plain_agg_case_res / plain_agg_res;
    const double final_abs_error = std::abs(encrypted_ratio - plain_ratio);
    std::cout << "Encrypted query result: " << std::endl;
    std::cout << std::setw(12) <<"promo_revenue" << std::endl;
    std::cout << std::setw(12) << encrypted_ratio << std::endl;
    std::cout << "Plain query result: " << std::endl;
    std::cout << std::setw(12) <<"promo_revenue" << std::endl;
    std::cout << std::setw(12) << plain_ratio << std::endl;
    std::cout << "Final absolute error: " << final_abs_error << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    const bool correct = std::isfinite(encrypted_ratio) && std::isfinite(plain_ratio) &&
                         repack_avg_error <= 0.25 && repack_case_avg_error <= 0.25 &&
                         final_abs_error <= 0.01;
    std::cout << "BENCH_RESULT,Q14," << num << "," << threads << ","
              << filtering_time << "," << aggregation_time << "," << (filtering_time + aggregation_time) << ","
              << repack_avg_error << "," << repack_case_avg_error << ","
              << plain_ratio << "," << encrypted_ratio << "," << final_abs_error << ","
              << (correct ? 1 : 0) << std::endl;
    return correct;

}

int main(int argc, char **argv)
{
    size_t rows = 1024;
    int threads = default_threads();
    int agg_threads = threads;
    int pair_threads = -1;
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if ((arg == "--rows" || arg == "-n") && i + 1 < argc)
            rows = std::stoull(argv[++i]);
        else if ((arg == "--threads" || arg == "-t") && i + 1 < argc)
            threads = std::stoi(argv[++i]);
        else if (arg == "--agg-threads" && i + 1 < argc)
            agg_threads = std::stoi(argv[++i]);
        else if (arg == "--pair-threads" && i + 1 < argc)
            pair_threads = std::stoi(argv[++i]);
        else if (arg == "--help" || arg == "-h")
        {
            cout << "Usage: " << argv[0] << " [--rows power_of_two] [--threads n] [--agg-threads n] [--pair-threads n]" << endl;
            return 0;
        }
        else
        {
            cerr << "Unknown argument: " << arg << endl;
            return 1;
        }
    }
    if (pair_threads <= 0)
        pair_threads = std::min(2, std::max(1, agg_threads));
    return relational_query14_ethmsb(rows, threads, agg_threads, pair_threads) ? 0 : 2;
}
