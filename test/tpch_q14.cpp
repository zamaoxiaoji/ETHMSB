#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"
#include "HEDB/conversion/repack.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace HEDB;
using namespace std;

/***
 * TPC-H Query 14
 * select
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
*/


uint64_t generate_date(uint64_t down, uint64_t up)
{
    uint64_t dyear, dmonth, dday, uyear, umonth, uday;
    dyear = down / 10000;
    dmonth = (down / 100) % 100;
    dday = down % 100;
    uyear = up / 10000;
    umonth = (up / 100) % 100;
    uday = up % 100;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    uniform_int_distribution<Lvl1::T> day_message(dday, uday);
    uniform_int_distribution<Lvl1::T> month_message(dmonth, umonth);
    uniform_int_distribution<Lvl1::T> year_message(dyear, uyear);
    return day_message(engine) + 100 * month_message(engine) + 10000 * year_message(engine);

}

void predicate_evaluation(std::vector<TLWELvl1> &pred_cres_total, std::vector<uint32_t> &pred_res_total,
    std::vector<TLWELvl1> &pred_cres_part, std::vector<uint32_t> &pred_res_part, size_t rows,
    std::vector<std::vector<uint32_t>> p_type_data, std::vector<uint64_t> ship_data,
    std::vector<uint64_t> l_partkey_data, std::vector<uint64_t> p_partkey_data, TFHESecretKey &sk, double &filter_time)
{

    std::cout<< "Predicate evaluation: " << std::endl;
    using P = Lvl2;
    TFHEEvalKey ek;
    std::cout<< "Generating evaluation key..." << std::endl;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    uint32_t  ship_bits = 16, p_type_bits = 8, l_partkey_bits = 16, p_partkey_bits = 16,
              ship_scale_bits, p_type_scale_bits, l_partkey_scale_bits, p_partkey_scale_bits;

    ship_scale_bits = std::numeric_limits<Lvl2::T>::digits - ship_bits - 1;
    p_type_scale_bits = std::numeric_limits<Lvl1::T>::digits - p_type_bits - 1;
    l_partkey_scale_bits = std::numeric_limits<Lvl2::T>::digits - l_partkey_bits - 1;
    p_partkey_scale_bits = std::numeric_limits<Lvl2::T>::digits - p_partkey_bits - 1;

    //Encrypt database
    std::cout<< "Encrypting Database..." << std::endl;
    std::vector<std::vector<TLWELvl1>> p_type_ciphers(rows , std::vector<TLWELvl1>(7));
    std::vector<TLWELvl2> ship_ciphers(rows), l_partkey_ciphers(rows), p_partkey_ciphers(rows);
#ifdef _OPENMP
    {
        int omp_threads = omp_get_max_threads();
        bool use_parallel_enc = (rows >= static_cast<size_t>(omp_threads * 4));
#pragma omp parallel for schedule(static) if(use_parallel_enc) num_threads(omp_threads)
        for (size_t i = 0; i < rows; i++)
        {
            for(int j = 0; j<7; j++){
                p_type_ciphers[i][j] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(p_type_data[i][j], Lvl1::α, pow(2., p_type_scale_bits), sk.key.get<Lvl1>());
            }
            ship_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(ship_data[i], Lvl2::α, pow(2., ship_scale_bits), sk.key.get<Lvl2>());
            l_partkey_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(l_partkey_data[i], Lvl2::α, pow(2., l_partkey_scale_bits), sk.key.get<Lvl2>());
            p_partkey_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(p_partkey_data[i], Lvl2::α, pow(2., p_partkey_scale_bits), sk.key.get<Lvl2>());
        }
    }
#else
    for (size_t i = 0; i < rows; i++)
    {
        for(int j = 0; j<7; j++){
            p_type_ciphers[i][j] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(p_type_data[i][j], Lvl1::α, pow(2., p_type_scale_bits), sk.key.get<Lvl1>());
        }
        ship_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(ship_data[i], Lvl2::α, pow(2., ship_scale_bits), sk.key.get<Lvl2>());
        l_partkey_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(l_partkey_data[i], Lvl2::α, pow(2., l_partkey_scale_bits), sk.key.get<Lvl2>());
        p_partkey_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(p_partkey_data[i], Lvl2::α, pow(2., p_partkey_scale_bits), sk.key.get<Lvl2>());
    }
#endif

    //Encrypt Predicate values
    std::cout<< "Encrypting Predicate Values..." << std::endl;

    Lvl2::T pred1 = 20101, pred2 = 20201;
    std::vector<uint32_t> pred3 = {80,82,79,77,79}; //PROMO 的 ascii码


    std::vector<Lvl1::T> pred_res1(rows, 0), pred_res2(rows, 0), pred_res3(rows, 1), pred_res4(rows, 0);
    for (size_t i = 0; i < rows; i++)
    {
        pred_res1[i] = (ship_data[i] >= pred1) ? 1 : 0;
        pred_res2[i] = (ship_data[i] < pred2) ? 1 : 0;
        for(int j =0; j< 5;j++){
            if( p_type_data[i][j] ==  pred3[j]){
                pred_res3[i] *= 1;
            }
            else{
                pred_res3[i] *= 0;
            }
        }
        pred_res4[i] = (l_partkey_data[i] == p_partkey_data[i]) ? 1 : 0;

        pred_res_part[i] = pred_res1[i] * pred_res2[i] * pred_res3[i] * pred_res4[i];
        pred_res_total[i] = pred_res1[i] * pred_res2[i] * pred_res4[i];
    }

    TLWELvl2 pred_cipher1, pred_cipher2;
    std::vector<TLWELvl1> pred_cipher3(5);

    pred_cipher1 = TFHEpp::tlweSymInt32Encrypt<Lvl2>(pred1, Lvl2::α, pow(2., ship_scale_bits), sk.key.get<Lvl2>());
    pred_cipher2 = TFHEpp::tlweSymInt32Encrypt<Lvl2>(pred2, Lvl2::α, pow(2., ship_scale_bits), sk.key.get<Lvl2>());
    for(int i=0; i< 5 ;i++){
        pred_cipher3[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(pred3[i], Lvl1::α, pow(2., p_type_scale_bits), sk.key.get<Lvl1>());
    }

    // Predicate Evaluation
    std::cout<< "Start Predicate Evaluation..." << std::endl;
    std::vector<TLWELvl1> pred_cres1(rows), pred_cres2(rows), pred_cres3(rows), pred_cres4(rows), pred_cres_total_logic(rows);
    TLWELvl1 pred_cres_str;

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
#ifdef _OPENMP
    {
        int omp_threads = omp_get_max_threads();
        bool use_parallel_filter = (rows >= static_cast<size_t>(omp_threads * 2));
#pragma omp parallel for schedule(static) if(use_parallel_filter) num_threads(omp_threads)
        for (size_t i = 0; i < rows; i++)
        {
            greater_than_equal<Lvl2>(ship_ciphers[i], pred_cipher1, pred_cres1[i], ship_bits, ek, LOGIC);
            less_than<Lvl2>(ship_ciphers[i], pred_cipher2, pred_cres2[i], ship_bits, ek, LOGIC);
            HomAND(pred_cres_total[i], pred_cres1[i], pred_cres2[i], ek, LOGIC);

            equal<Lvl2>(l_partkey_ciphers[i] ,p_partkey_ciphers[i], pred_cres4[i], l_partkey_bits, ek, LOGIC);
            HomAND(pred_cres_total_logic[i], pred_cres_total[i], pred_cres4[i], ek, LOGIC);

            HomAND(pred_cres_total[i], pred_cres_total[i], pred_cres4[i], ek, ARITHMETIC);

            equal<Lvl1>(p_type_ciphers[i][0] ,pred_cipher3[0] , pred_cres3[i], p_type_bits, ek, LOGIC);
            for(int j = 1; j < 5; j++){
                equal<Lvl1>(p_type_ciphers[i][j] ,pred_cipher3[j] , pred_cres_str, p_type_bits, ek, LOGIC);
                HomAND(pred_cres3[i], pred_cres_str, pred_cres3[i], ek, LOGIC);
            }
            HomAND(pred_cres_part[i], pred_cres3[i], pred_cres_total_logic[i], ek, ARITHMETIC);
        }
    }
#else
    for (size_t i = 0; i < rows; i++)
    {
        greater_than_equal<Lvl2>(ship_ciphers[i], pred_cipher1, pred_cres1[i], ship_bits, ek, LOGIC);
        less_than<Lvl2>(ship_ciphers[i], pred_cipher2, pred_cres2[i], ship_bits, ek, LOGIC);
        HomAND(pred_cres_total[i], pred_cres1[i], pred_cres2[i], ek, LOGIC);

        equal<Lvl2>(l_partkey_ciphers[i] ,p_partkey_ciphers[i], pred_cres4[i], l_partkey_bits, ek, LOGIC);
        HomAND(pred_cres_total_logic[i], pred_cres_total[i], pred_cres4[i], ek, LOGIC);

        HomAND(pred_cres_total[i], pred_cres_total[i], pred_cres4[i], ek, ARITHMETIC);

        //std::cout<< "11111111111111111111111111111111111" << std::endl;
        equal<Lvl1>(p_type_ciphers[i][0] ,pred_cipher3[0] , pred_cres3[i], p_type_bits, ek, LOGIC);
        for(int j = 1; j < 5; j++){
            equal<Lvl1>(p_type_ciphers[i][j] ,pred_cipher3[j] , pred_cres_str, p_type_bits, ek, LOGIC);
            HomAND(pred_cres3[i], pred_cres_str, pred_cres3[i], ek, LOGIC);
        }
        HomAND(pred_cres_part[i], pred_cres3[i], pred_cres_total_logic[i], ek, ARITHMETIC);
    }
#endif
    end = std::chrono::system_clock::now();
    std::vector<uint32_t> pred_cres_total_de(rows), pred_cres_part_de(rows), pred_cres1_de(rows), pred_cres2_de(rows), pred_cres3_de(rows), pred_cres4_de(rows), pred_cres5_de(rows);
    // for (size_t i = 0; i < rows; i++)
    // {
    //     pred_cres_de[i] = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres[i], pow(2., 31), sk.key.get<Lvl1>());
    //     pred_cres1_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres1[i], sk.key.lvl1);
    //     pred_cres2_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres2[i], sk.key.lvl1);
    //     pred_cres3_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres3[i], sk.key.lvl1);
    //     pred_cres4_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres4[i], sk.key.lvl1);
    //     pred_cres5_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres5[i], sk.key.lvl1);
    //}

    size_t error_time_total = 0;
    size_t error_time_part = 0;

    uint32_t rlwe_scale_bits = 29;
    for (size_t i = 0; i < rows; i++)
    {
        TFHEpp::ari_rescale(pred_cres_total[i], pred_cres_total[i], rlwe_scale_bits, ek);
        TFHEpp::ari_rescale(pred_cres_part[i], pred_cres_part[i], rlwe_scale_bits, ek);
    }
    for (size_t i = 0; i < rows; i++)
    {
        pred_cres_total_de[i] = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres_total[i], pow(2., 29), sk.key.get<Lvl1>());
        pred_cres_part_de[i] = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres_part[i], pow(2., 29), sk.key.get<Lvl1>());
    }
    for (size_t i = 0; i < rows; i++)
    {
        error_time_total += (pred_cres_total_de[i] == pred_res_total[i])? 0 : 1;
        error_time_part += (pred_cres_part_de[i] == pred_res_part[i])? 0 : 1;
    }
    cout << "Predicate Evaluaton Time (s): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000 << std::endl;
    cout << "Predicate Total Error: " << error_time_total << std::endl;
    cout << "Predicate Part Error: " << error_time_part << std::endl;
    filter_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

}


void aggregation(std::vector<TLWELvl1> &pred_cres, std::vector<uint32_t> &pred_res, size_t tfhe_n,
            std::vector<double> &extendedprice_data, std::vector<double> &discount_data_double,
             size_t rows, TFHESecretKey &sk, double &aggregation_time,
             double &plain_aggregation_result , double &de_aggregation_result)
{

    std::cout << "Aggregation :" << std::endl;
    uint64_t scale_bits = 29;
    uint64_t modq_bits = 32;
    uint64_t modulus_bits = 45;
    uint64_t repack_scale_bits = modulus_bits + scale_bits - modq_bits;
    uint64_t slots_count = pred_cres.size();
    std::cout << "Generating Parameters..." << std::endl;
    seal::EncryptionParameters parms(seal::scheme_type::ckks);
    size_t poly_modulus_degree = 65536;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_modulus_degree, {59, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 59}));
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
    seal::Decryptor decryptor(context, seal_secret_key);

    //encoder
    seal::CKKSEncoder ckks_encoder(context);



    //generate evaluation key
    std::cout << "Generating Conversion Key..." << std::endl;
    LTPreKey pre_key;
    LWEsToRLWEKeyGen(pre_key, std::pow(2., modulus_bits), seal_secret_key, sk, tfhe_n, ckks_encoder, encryptor, context);


    // conversion (P0: chunked parallel LWEsToRLWE + HomRound)
    std::cout << "Starting Conversion..." << std::endl;
    struct ChunkAggRes { seal::Ciphertext cipher; double time_ms = 0.0; };
    auto convert_chunk = [&](size_t begin, size_t count) -> ChunkAggRes {
        ChunkAggRes r;
        seal::CKKSEncoder enc_local(context);
        seal::Evaluator eval_local(context);
        seal::Decryptor dec_local(context, seal_secret_key);
        std::vector<TLWELvl1> slice(pred_cres.begin() + begin, pred_cres.begin() + begin + count);
        auto t0 = std::chrono::system_clock::now();
        LWEsToRLWE(r.cipher, slice, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits),
                   enc_local, galois_keys, relin_keys, eval_local, context);
        HomRound(r.cipher, r.cipher.scale(), enc_local, relin_keys, eval_local, dec_local, context);
        auto t1 = std::chrono::system_clock::now();
        r.time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
        return r;
    };
    const size_t chunk_size = 16384;
    const size_t chunk_count = (rows + chunk_size - 1) / chunk_size;
    std::vector<ChunkAggRes> chunks(chunk_count);
#ifdef _OPENMP
    {
        int omp_threads = omp_get_max_threads();
        bool use_parallel_chunks = (!omp_in_parallel()) && (chunk_count > 1);
#pragma omp parallel for schedule(static) if(use_parallel_chunks) num_threads(std::min<int>(omp_threads, (int)chunk_count))
        for (size_t c = 0; c < chunk_count; ++c)
        {
            size_t begin = c * chunk_size;
            size_t cnt = std::min(chunk_size, rows - begin);
            chunks[c] = convert_chunk(begin, cnt);
        }
    }
#else
    for (size_t c = 0; c < chunk_count; ++c)
    {
        size_t begin = c * chunk_size;
        size_t cnt = std::min(chunk_size, rows - begin);
        chunks[c] = convert_chunk(begin, cnt);
    }
#endif
    seal::Ciphertext result = chunks.front().cipher;
    for (size_t c = 1; c < chunk_count; ++c) evaluator.add_inplace(result, chunks[c].cipher);
    aggregation_time = 0.0;
    for (auto &c : chunks) aggregation_time += c.time_ms;
    seal::Plaintext plain;
    std::vector<double> computed(slots_count);
    decryptor.decrypt(result, plain);
    seal::pack_decode(computed, plain, ckks_encoder);
    double err = 0.;
    for (size_t i = 0; i < slots_count; ++i) err += std::abs(computed[i] - pred_res[i]);
    printf("Repack average error = %f ~ 2^%.1f\n", err / slots_count, std::log2(err / slots_count));


    // Filter result * data
    std::vector<double> price_discount(extendedprice_data.size());
    seal::Ciphertext price_discount_cipher;
#ifdef _OPENMP
    {
        int omp_threads = omp_get_max_threads();
        bool use_parallel_pd = (rows >= static_cast<size_t>(omp_threads * 16));
#pragma omp parallel for schedule(static) if(use_parallel_pd) num_threads(omp_threads)
        for (size_t i = 0; i < rows; i++)
        {
            price_discount[i] = extendedprice_data[i] * (100-discount_data_double[i]);
        }
    }
#else
    for (size_t i = 0; i < rows; i++)
    {
        price_discount[i] = extendedprice_data[i] * (100-discount_data_double[i]);
    }
#endif
    double qd = parms.coeff_modulus()[result.coeff_modulus_size() - 1].value();
    seal::pack_encode(price_discount, qd, plain, ckks_encoder);
    encryptor.encrypt_symmetric(plain, price_discount_cipher);

    std::cout << "Aggregating price and discount .." << std::endl;
    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    seal::multiply_and_relinearize(result, price_discount_cipher, result, evaluator, relin_keys);
    evaluator.rescale_to_next_inplace(result);
    int logrow = log2(rows);

    seal::Ciphertext temp;
    for (size_t i = 0; i < logrow; i++)
    {
        temp = result;
        size_t step = 1 << (logrow - i - 1);
        evaluator.rotate_vector_inplace(temp, step, galois_keys);
        evaluator.add_inplace(result, temp);
    }
    end = std::chrono::system_clock::now();
    aggregation_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::vector<double> agg_result(slots_count);
    decryptor.decrypt(result, plain);
    seal::pack_decode(agg_result, plain, ckks_encoder);
    double plain_result = 0;
#ifdef _OPENMP
    {
        int omp_threads = omp_get_max_threads();
        bool use_parallel_plain = (rows >= static_cast<size_t>(omp_threads * 16));
#pragma omp parallel for schedule(static) reduction(+:plain_result) if(use_parallel_plain) num_threads(omp_threads)
        for (size_t i = 0; i < rows; i++)
        {
            plain_result += extendedprice_data[i] * (100 - discount_data_double[i]) * pred_res[i];
        }
    }
#else
    for (size_t i = 0; i < rows; i++)
    {
        plain_result += extendedprice_data[i] * (100 - discount_data_double[i]) * pred_res[i];
    }
#endif
    plain_aggregation_result = plain_result;
    de_aggregation_result = std::round(agg_result[0]);

    cout << "aggregation_time: " << aggregation_time / 1000 << " s" << endl;
    cout << "Plain_result: " << plain_result << endl;
    cout << "Encrypted query result: " << std::round(agg_result[0]) <<endl;

}

void query_evaluation(size_t rows)
{
    cout << "-------tpch_q14-------"  <<endl;
    cout << "data rows: " << rows <<endl;
    TFHESecretKey sk;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
     // Generate database
    std::vector<uint32_t> discount_data(rows);
    std::vector<uint64_t> ship_data(rows);
    std::vector<uint64_t> l_partkey_data(rows);
    std::vector<uint64_t> p_partkey_data(rows);
    std::vector<std::vector<uint32_t>> p_type(rows , std::vector<uint32_t>(7));

    uint32_t  ship_bits = 16, p_type_bits = 7, discount_bits = 6, l_partkey_bits = 16 , p_partkey_bits = 16 ,
              ship_scale_bits, p_type_scale_bits, discount_scale_bits, l_partkey_scale_bits, p_partkey_scale_bits;
    uniform_int_distribution<Lvl1::T> p_type_message(0, 127);
    uniform_int_distribution<Lvl1::T> discount_message(0, (1 << discount_bits) - 1);
    uniform_int_distribution<Lvl1::T> possibility_message(0, 9);
    uniform_int_distribution<Lvl1::T> partkey_message(0, (1 << l_partkey_bits) - 1);

    ship_scale_bits = std::numeric_limits<Lvl2::T>::digits - ship_bits - 1;
    p_type_scale_bits = std::numeric_limits<Lvl1::T>::digits - p_type_bits - 1;
    discount_scale_bits = std::numeric_limits<Lvl1::T>::digits - discount_bits - 1;

    uint32_t possibility1 , possibility2;
    for (size_t i = 0; i < rows; i++)
    {


        possibility1 = possibility_message(engine);
        if(possibility1 < 5){
            for(int j =0; j< 7; j++){
                p_type[i][j] = p_type_message(engine);
            }
        }
        else{
            p_type[i] = {80,82,79,77,79,60,60};
        }

        ship_data[i] = generate_date(10101, 21230);
        discount_data[i] = discount_message(engine);

        possibility2 = possibility_message(engine);
        if(possibility2 < 5){
            l_partkey_data[i] = partkey_message(engine);
            p_partkey_data[i] = partkey_message(engine);
        }
        else{
            l_partkey_data[i] = partkey_message(engine);
            p_partkey_data[i] = l_partkey_data[i];
        }
    }

    std::vector<double> extendedprice_data(rows), discount_data_double(rows);
    std::uniform_int_distribution<uint64_t> extendedprice_message(1, 5);
    for (size_t i = 0; i < rows; i++)
    {
        extendedprice_data[i] = (double)extendedprice_message(engine);
        discount_data_double[i] = (double)discount_data[i];
    }

    //quanlity_data[0] = 1;
    discount_data[0] = 9;
    ship_data[0] = 21215;
    double filter_time, aggregation_time_total, aggregation_time_part;
    std::vector<TLWELvl1> pred_cres_total(rows);
    std::vector<uint32_t> pred_res_total(rows, 0);
    std::vector<TLWELvl1> pred_cres_part(rows);
    std::vector<uint32_t> pred_res_part(rows, 0);

    double plain_aggregation_total_result, de_aggregation_total_result, plain_aggregation_part_result, de_aggregation_part_result ;

    predicate_evaluation(pred_cres_total, pred_res_total, pred_cres_part, pred_res_part, rows, p_type, ship_data, l_partkey_data,  p_partkey_data, sk, filter_time);
    aggregation(pred_cres_total, pred_res_total, Lvl1::n, extendedprice_data, discount_data_double, rows, sk, aggregation_time_total, plain_aggregation_total_result, de_aggregation_total_result);
    aggregation(pred_cres_part, pred_res_part, Lvl1::n, extendedprice_data, discount_data_double, rows, sk, aggregation_time_part, plain_aggregation_part_result, de_aggregation_part_result);


    double plain_result = 100 * plain_aggregation_part_result / plain_aggregation_total_result ;
    double de_result = 100 * de_aggregation_part_result / de_aggregation_total_result ;

    cout << "plain_aggregation_part_result: " << plain_aggregation_part_result << endl;
    cout << "plain_aggregation_total_result: " << plain_aggregation_total_result << endl;

    cout << "plain_result: " << plain_result << endl;

    cout << "de_aggregation_part_result: " << de_aggregation_part_result << endl;
    cout << "de_aggregation_total_result: " << de_aggregation_total_result << endl;

    cout << "de_result: " << de_result  << endl;
    cout << "End-to-End Time: " << (filter_time + aggregation_time_part + aggregation_time_total) / 1000 << " s" << endl;

}

int main()
{
    // query_evaluation(256);
    // query_evaluation(1024);
    // query_evaluation(4096);
    // query_evaluation(8192);
    query_evaluation(32768);
    // query_evaluation(32400);
}
