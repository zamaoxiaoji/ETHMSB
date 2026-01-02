#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"
#include "HEDB/conversion/repack.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace HEDB;
using namespace std;

/***
 * TPC-H Query 19
select
    sum(l_extendedprice * (1 - l_discount) ) as revenue
from
    lineitem, part
where
(
    p_partkey = l_partkey
    and p_brand = ‘[BRAND1]’ /*特定品牌。BRAND1、BRAND2、BRAND3＝‘Brand＃MN’，M和N是两个字母，代表两个数值，相互独立，取值在1到5之间
    and p_container in ( ‘SM CASE’, ‘SM BOX’, ‘SM PACK’, ‘SM PKG’) //包装范围
    and l_quantity >= [QUANTITY1] and l_quantity <= [QUANTITY1] + 10 /* QUANTITY1 是1到10之间的任意取值
    and p_size between 1 and 5 //尺寸范围
    and l_shipmode in (‘AIR’, ‘AIR REG’) //运输模式，可以是(AIR, AIR REG, SHIP, SHIP REG, TRUCK, TRUCK REG, RAIL, RAIL REG)
    and l_shipinstruct = ‘DELIVER IN PERSON’ //运输指令，可以是（DELIVER IN PERSON，COLLECT COD，SHIP BY AIR，SHIP BY RAIL，SHIP BY TRUCK，PICKUP）
)
or
(
    p_partkey = l_partkey
    and p_brand = ‘[BRAND2]’
    and p_container in (‘MED BAG’, ‘MED BOX’, ‘MED PKG’, ‘MED PACK’)
    and l_quantity >= [QUANTITY2] and l_quantity <= [QUANTITY2] + 10 /* QUANTITY2 是10到20之间的任意取值
    and p_size between 1 and 10
    and l_shipmode in (‘AIR’, ‘AIR REG’)
    and l_shipinstruct = ‘DELIVER IN PERSON’
)
or
(
    p_partkey = l_partkey
    and p_brand = ‘[BRAND3]’
    and p_container in ( ‘LG CASE’, ‘LG BOX’, ‘LG PACK’, ‘LG PKG’)
    and l_quantity >= [QUANTITY3] and l_quantity <= [QUANTITY3] + 10 /* QUANTITY3 是20到30之间的任意取值
    and p_size between 1 and 15
    and l_shipmode in (‘AIR’, ‘AIR REG’)
    and l_shipinstruct = ‘DELIVER IN PERSON’
);
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


void predicate_evaluation(std::vector<TLWELvl1>& pred_cres, std::vector<uint32_t>& pred_res, size_t rows,
                          std::vector<uint32_t>& p_brand_data, std::vector<uint32_t>& p_container_data,
                          std::vector<uint32_t>& l_quantity_data,
                          std::vector<uint32_t>& p_size_data, std::vector<uint32_t>& l_shipmode_data,
                          std::vector<uint32_t>& l_shipinstruct_data,
                          std::vector<uint64_t> l_partkey_data, std::vector<uint64_t> p_partkey_data, TFHESecretKey& sk,
                          double& filter_time)

{
    std::cout << "Predicate evaluation: " << std::endl;
    //using P = Lvl2;
    TFHEEvalKey ek;
    std::cout << "Generating evaluation key..." << std::endl;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    uint32_t p_brand_bits = 6, p_container_bits = 4, l_quantity_bits = 6, p_size_bits = 6, l_shipmode_bits = 4,
             l_shipinstruct_bits = 4, l_partkey_bits = 16, p_partkey_bits = 16;
    uint32_t p_brand_scale_bits, p_container_scale_bits, l_quantity_scale_bits, p_size_scale_bits, l_shipmode_scale_bits
             , l_shipinstruct_scale_bits, l_partkey_scale_bits, p_partkey_scale_bits;
    p_brand_scale_bits = std::numeric_limits<Lvl1::T>::digits - p_brand_bits - 1;
    p_container_scale_bits = std::numeric_limits<Lvl1::T>::digits - p_container_bits - 1;
    l_quantity_scale_bits = std::numeric_limits<Lvl1::T>::digits - l_quantity_bits - 1;
    p_size_scale_bits = std::numeric_limits<Lvl1::T>::digits - p_size_bits - 1;
    l_shipmode_scale_bits = std::numeric_limits<Lvl1::T>::digits - l_shipmode_bits - 1;
    l_shipinstruct_scale_bits = std::numeric_limits<Lvl1::T>::digits - l_shipinstruct_bits - 1;
    l_partkey_scale_bits = std::numeric_limits<Lvl2::T>::digits - l_partkey_bits - 1;
    p_partkey_scale_bits = std::numeric_limits<Lvl2::T>::digits - p_partkey_bits - 1;

    //Encrypt database
    std::cout << "Encrypting Database..." << std::endl;

    std::vector<TLWELvl1> p_brand_ciphers(rows);
    std::vector<TLWELvl1> p_container_ciphers(rows);
    std::vector<TLWELvl1> l_quantity_ciphers(rows);
    std::vector<TLWELvl1> p_size_ciphers(rows);
    std::vector<TLWELvl1> l_shipmode_ciphers(rows);
    std::vector<TLWELvl1> l_shipinstruct_ciphers(rows);
    std::vector<TLWELvl2> l_partkey_ciphers(rows), p_partkey_ciphers(rows);

#ifdef _OPENMP
    int omp_threads = omp_get_max_threads();
    bool use_parallel_enc = (rows >= static_cast<size_t>(omp_threads * 4));
#pragma omp parallel for schedule(static) if(use_parallel_enc) num_threads(omp_threads)
#endif
    for (size_t i = 0; i < rows; i++)
    {
        p_brand_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(p_brand_data[i], Lvl1::α, pow(2., p_brand_scale_bits),
                                                               sk.key.get<Lvl1>());
        p_container_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(p_container_data[i], Lvl1::α,
                                                                   pow(2., p_container_scale_bits), sk.key.get<Lvl1>());
        l_quantity_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(l_quantity_data[i], Lvl1::α,
                                                                  pow(2., l_quantity_scale_bits), sk.key.get<Lvl1>());
        p_size_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(p_size_data[i], Lvl1::α, pow(2., p_size_scale_bits),
                                                              sk.key.get<Lvl1>());
        l_shipmode_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(l_shipmode_data[i], Lvl1::α,
                                                                  pow(2., l_shipmode_scale_bits), sk.key.get<Lvl1>());
        l_shipinstruct_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(l_shipinstruct_data[i], Lvl1::α,
                                                                      pow(2., l_shipinstruct_scale_bits),
                                                                      sk.key.get<Lvl1>());
        l_partkey_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(l_partkey_data[i], Lvl2::α,
                                                                 pow(2., l_partkey_scale_bits), sk.key.get<Lvl2>());
        p_partkey_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<Lvl2>(p_partkey_data[i], Lvl2::α,
                                                                 pow(2., p_partkey_scale_bits), sk.key.get<Lvl2>());
    }

    //Encrypt Predicate values
    std::cout << "Encrypting Predicate Values..." << std::endl;

    Lvl1::T brand1 = 0;
    std::vector<uint32_t> container1 = {0, 1, 2, 3};
    Lvl1::T quantity1 = 10;
    std::vector<uint32_t> shipmode = {0, 1};
    Lvl1::T shipinstruct = 0;

    Lvl1::T brand2 = 1;
    std::vector<uint32_t> container2 = {4, 5, 6, 7};
    Lvl1::T quantity2 = 15;

    Lvl1::T brand3 = 2;
    std::vector<uint32_t> container3 = {8, 9, 10, 11};
    Lvl1::T quantity3 = 20;

    std::vector<Lvl1::T> brand_res1(rows, 0), brand_res2(rows, 0), brand_res3(rows, 1);
    std::vector<Lvl1::T> container_res1(rows, 0), container_res2(rows, 0), container_res3(rows, 1);
    std::vector<Lvl1::T> quantity_res1(rows, 0), quantity_res2(rows, 0), quantity_res3(rows, 1);
    std::vector<Lvl1::T> size_res1(rows, 0), size_res2(rows, 0), size_res3(rows, 1);
    std::vector<Lvl1::T> shipmode_res1(rows, 0);
    std::vector<Lvl1::T> shipinstruct_res1(rows, 0);
    std::vector<Lvl1::T> partkey_res1(rows, 0);
    std::vector<Lvl1::T> pred_res1(rows, 0), pred_res2(rows, 0), pred_res3(rows, 0);

    for (size_t i = 0; i < rows; i++)
    {
        shipmode_res1[i] = (l_shipmode_data[i] == shipmode[0] || l_shipmode_data[i] == shipmode[1]) ? 1 : 0;
        shipinstruct_res1[i] = (l_shipinstruct_data[i] == shipinstruct) ? 1 : 0;
        partkey_res1[i] = (l_partkey_data[i] == p_partkey_data[i]) ? 1 : 0;

        brand_res1[i] = (p_brand_data[i] == brand1) ? 1 : 0;
        container_res1[i] = (p_container_data[i] == container1[0] || p_container_data[i] == container1[1] ||
                                p_container_data[i] == container1[2] || p_container_data[i] == container1[3])
                                ? 1
                                : 0;
        quantity_res1[i] = (l_quantity_data[i] >= quantity1 && l_quantity_data[i] <= quantity1 + 10) ? 1 : 0;
        size_res1[i] = (p_size_data[i] >= 1 && p_size_data[i] <= 5) ? 1 : 0;


        pred_res1[i] = (partkey_res1[i] == 1 && shipmode_res1[i] == 1 && shipinstruct_res1[i] == 1 && brand_res1[i] == 1
                           && container_res1[i] == 1 && quantity_res1[i] == 1 && size_res1[i])
                           ? 1
                           : 0;

        brand_res2[i] = (p_brand_data[i] == brand2) ? 1 : 0;
        container_res2[i] = (p_container_data[i] == container2[0] || p_container_data[i] == container2[1] ||
                                p_container_data[i] == container2[2] || p_container_data[i] == container2[3])
                                ? 1
                                : 0;
        quantity_res2[i] = (l_quantity_data[i] >= quantity2 && l_quantity_data[i] <= quantity2 + 10) ? 1 : 0;
        size_res2[i] = (p_size_data[i] >= 1 && p_size_data[i] <= 10) ? 1 : 0;

        pred_res2[i] = (partkey_res1[i] == 1 && shipmode_res1[i] + shipinstruct_res1[i] + brand_res2[i] + container_res2
                           [i] + quantity_res2[i] + size_res2[i] > 5)
                           ? 1
                           : 0;

        brand_res3[i] = (p_brand_data[i] == brand3) ? 1 : 0;
        container_res3[i] = (p_container_data[i] == container3[0] || p_container_data[i] == container3[1] ||
                                p_container_data[i] == container3[2] || p_container_data[i] == container3[3])
                                ? 1
                                : 0;
        quantity_res3[i] = (l_quantity_data[i] >= quantity3 && l_quantity_data[i] <= quantity3 + 10) ? 1 : 0;
        size_res3[i] = (p_size_data[i] >= 1 && p_size_data[i] <= 15) ? 1 : 0;

        pred_res3[i] = (partkey_res1[i] == 1 && shipmode_res1[i] + shipinstruct_res1[i] + brand_res3[i] + container_res3
                           [i] + quantity_res3[i] + size_res3[i] > 5)
                           ? 1
                           : 0;

        pred_res[i] = (pred_res1[i] == 1 || pred_res2[i] == 1 || pred_res3[i] == 1) ? 1 : 0;

        // cout << "pred_res1[" << i <<"]: " <<  pred_res1[i] << std::endl;
        // cout << "pred_res[" << i <<"]: " <<  pred_res[i] << std::endl;
    }


    TLWELvl1 pred_brand1_cipher, pred_brand2_cipher, pred_brand3_cipher;
    std::vector<TLWELvl1> pred_container1_cipher(4), pred_container2_cipher(4), pred_container3_cipher(4);
    TLWELvl1 pred_quantity1_left_cipher, pred_quantity1_right_cipher, pred_quantity2_left_cipher,
             pred_quantity2_right_cipher, pred_quantity3_left_cipher, pred_quantity3_right_cipher;
    TLWELvl1 pred_size1_left_cipher, pred_size1_right_cipher, pred_size2_left_cipher, pred_size2_right_cipher,
             pred_size3_left_cipher, pred_size3_right_cipher;
    std::vector<TLWELvl1> pred_shipmode_cipher(2);
    TLWELvl1 pred_shipinstruct_cipher;


    pred_brand1_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(brand1, Lvl1::α, pow(2., p_brand_scale_bits),
                                                           sk.key.get<Lvl1>());
    for (int i = 0; i < 4; i++)
    {
        pred_container1_cipher[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(container1[i], Lvl1::α,
                                                                      pow(2., p_container_scale_bits),
                                                                      sk.key.get<Lvl1>());
    }
    pred_quantity1_left_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(quantity1, Lvl1::α, pow(2., l_quantity_scale_bits),
                                                                   sk.key.get<Lvl1>());
    pred_quantity1_right_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(quantity1 + 10, Lvl1::α,
                                                                    pow(2., l_quantity_scale_bits), sk.key.get<Lvl1>());
    pred_size1_left_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, Lvl1::α, pow(2., p_size_scale_bits),
                                                               sk.key.get<Lvl1>());
    pred_size1_right_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(5, Lvl1::α, pow(2., p_size_scale_bits),
                                                                sk.key.get<Lvl1>());
    pred_shipmode_cipher[0] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(shipmode[0], Lvl1::α, pow(2., l_shipmode_scale_bits),
                                                                sk.key.get<Lvl1>());
    pred_shipmode_cipher[1] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(shipmode[1], Lvl1::α, pow(2., l_shipmode_scale_bits),
                                                                sk.key.get<Lvl1>());
    pred_shipinstruct_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(shipinstruct, Lvl1::α,
                                                                 pow(2., l_shipinstruct_scale_bits),
                                                                 sk.key.get<Lvl1>());

    pred_brand2_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(brand2, Lvl1::α, pow(2., p_brand_scale_bits),
                                                           sk.key.get<Lvl1>());
    for (int i = 0; i < 4; i++)
    {
        pred_container2_cipher[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(container2[i], Lvl1::α,
                                                                      pow(2., p_container_scale_bits),
                                                                      sk.key.get<Lvl1>());
    }
    pred_quantity2_left_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(quantity2, Lvl1::α, pow(2., l_quantity_scale_bits),
                                                                   sk.key.get<Lvl1>());
    pred_quantity2_right_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(quantity2 + 10, Lvl1::α,
                                                                    pow(2., l_quantity_scale_bits), sk.key.get<Lvl1>());
    pred_size2_left_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, Lvl1::α, pow(2., p_size_scale_bits),
                                                               sk.key.get<Lvl1>());
    pred_size2_right_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(10, Lvl1::α, pow(2., p_size_scale_bits),
                                                                sk.key.get<Lvl1>());

    pred_brand3_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(brand3, Lvl1::α, pow(2., p_brand_scale_bits),
                                                           sk.key.get<Lvl1>());
    for (int i = 0; i < 4; i++)
    {
        pred_container3_cipher[i] = TFHEpp::tlweSymInt32Encrypt<Lvl1>(container3[i], Lvl1::α,
                                                                      pow(2., p_container_scale_bits),
                                                                      sk.key.get<Lvl1>());
    }
    pred_quantity3_left_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(quantity3, Lvl1::α, pow(2., l_quantity_scale_bits),
                                                                   sk.key.get<Lvl1>());
    pred_quantity3_right_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(quantity3 + 10, Lvl1::α,
                                                                    pow(2., l_quantity_scale_bits), sk.key.get<Lvl1>());
    pred_size3_left_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, Lvl1::α, pow(2., p_size_scale_bits),
                                                               sk.key.get<Lvl1>());
    pred_size3_right_cipher = TFHEpp::tlweSymInt32Encrypt<Lvl1>(15, Lvl1::α, pow(2., p_size_scale_bits),
                                                                sk.key.get<Lvl1>());

    // Predicate Evaluation
    std::cout << "Start Predicate Evaluation..." << std::endl;
    std::vector<TLWELvl1> pred_cres1(rows), pred_cres2(rows), pred_cres3(rows);
    std::vector<TLWELvl1> partkey_cres(rows), brand_cres(rows), container_cres(rows), quantity_cres(rows),
                          size_cres(rows), shipmode_cres(rows), shipinstruct_cres(rows);
    std::vector<TLWELvl1> temp_cres(rows);

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

#ifdef _OPENMP
    {
        int omp_threads = omp_get_max_threads();
        bool use_parallel_filter = (rows >= static_cast<size_t>(omp_threads * 2));
#pragma omp parallel for schedule(static) if(use_parallel_filter) num_threads(omp_threads)
        for (size_t i = 0; i < rows; i++)
        {
            equal<Lvl2>(l_partkey_ciphers[i], p_partkey_ciphers[i], partkey_cres[i], l_partkey_bits, ek, LOGIC);
            equal<Lvl1>(p_brand_ciphers[i], pred_brand1_cipher, brand_cres[i], p_brand_bits, ek, LOGIC);
            HomAND(pred_cres1[i], partkey_cres[i], brand_cres[i], ek, LOGIC);
            equal<Lvl1>(p_container_ciphers[i], pred_container1_cipher[0], container_cres[i], p_container_bits, ek, LOGIC);
            for (uint32_t j = 1; j < 4; j++)
            {
                equal<Lvl1>(p_container_ciphers[i], pred_container1_cipher[j], temp_cres[i], p_container_bits, ek, LOGIC);
                HomOR(container_cres[i], container_cres[i], temp_cres[i], ek, LOGIC);
            }
            HomAND(pred_cres1[i], pred_cres1[i], container_cres[i], ek, LOGIC);
            greater_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity1_left_cipher, quantity_cres[i], l_quantity_bits,
                                     ek, LOGIC);
            HomAND(pred_cres1[i], pred_cres1[i], quantity_cres[i], ek, LOGIC);
            less_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity1_right_cipher, quantity_cres[i], l_quantity_bits, ek,
                                  LOGIC);
            HomAND(pred_cres1[i], pred_cres1[i], quantity_cres[i], ek, LOGIC);
            greater_than_equal<Lvl1>(p_size_ciphers[i], pred_size1_left_cipher, size_cres[i], p_size_bits, ek, LOGIC);
            HomAND(pred_cres1[i], pred_cres1[i], size_cres[i], ek, LOGIC);
            less_than_equal<Lvl1>(p_size_ciphers[i], pred_size1_right_cipher, size_cres[i], p_size_bits, ek, LOGIC);
            HomAND(pred_cres1[i], pred_cres1[i], size_cres[i], ek, LOGIC);
            equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[0], shipmode_cres[i], l_shipmode_bits, ek, LOGIC);
            equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[1], temp_cres[i], l_shipmode_bits, ek, LOGIC);
            HomOR(shipmode_cres[i], shipmode_cres[i], temp_cres[i], ek, LOGIC);
            HomAND(pred_cres1[i], pred_cres1[i], shipmode_cres[i], ek, LOGIC);
            equal<Lvl1>(l_shipinstruct_ciphers[i], pred_shipinstruct_cipher, shipinstruct_cres[i], l_shipinstruct_bits,
                        ek, LOGIC);
            HomAND(pred_cres1[i], pred_cres1[i], shipinstruct_cres[i], ek, LOGIC);
            equal<Lvl2>(l_partkey_ciphers[i], p_partkey_ciphers[i], partkey_cres[i], l_partkey_bits, ek, LOGIC);
            equal<Lvl1>(p_brand_ciphers[i], pred_brand2_cipher, brand_cres[i], p_brand_bits, ek, LOGIC);
            HomAND(pred_cres2[i], partkey_cres[i], brand_cres[i], ek, LOGIC);
            equal<Lvl1>(p_container_ciphers[i], pred_container2_cipher[0], container_cres[i], p_container_bits, ek, LOGIC);
            for (uint32_t j = 1; j < 4; j++)
            {
                equal<Lvl1>(p_container_ciphers[i], pred_container2_cipher[j], temp_cres[i], p_container_bits, ek, LOGIC);
                HomOR(container_cres[i], container_cres[i], temp_cres[i], ek, LOGIC);
            }
            HomAND(pred_cres2[i], pred_cres2[i], container_cres[i], ek, LOGIC);
            greater_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity2_left_cipher, quantity_cres[i], l_quantity_bits,
                                     ek, LOGIC);
            HomAND(pred_cres2[i], pred_cres2[i], quantity_cres[i], ek, LOGIC);
            less_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity2_right_cipher, quantity_cres[i], l_quantity_bits, ek,
                                  LOGIC);
            HomAND(pred_cres2[i], pred_cres2[i], quantity_cres[i], ek, LOGIC);
            greater_than_equal<Lvl1>(p_size_ciphers[i], pred_size2_left_cipher, size_cres[i], p_size_bits, ek, LOGIC);
            HomAND(pred_cres2[i], pred_cres2[i], size_cres[i], ek, LOGIC);
            less_than_equal<Lvl1>(p_size_ciphers[i], pred_size2_right_cipher, size_cres[i], p_size_bits, ek, LOGIC);
            HomAND(pred_cres2[i], pred_cres2[i], size_cres[i], ek, LOGIC);
            equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[0], shipmode_cres[i], l_shipmode_bits, ek, LOGIC);
            equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[1], temp_cres[i], l_shipmode_bits, ek, LOGIC);
            HomOR(shipmode_cres[i], shipmode_cres[i], temp_cres[i], ek, LOGIC);
            HomAND(pred_cres2[i], pred_cres2[i], shipmode_cres[i], ek, LOGIC);
            equal<Lvl1>(l_shipinstruct_ciphers[i], pred_shipinstruct_cipher, shipinstruct_cres[i], l_shipinstruct_bits,
                        ek, LOGIC);
            HomAND(pred_cres2[i], pred_cres2[i], shipinstruct_cres[i], ek, LOGIC);
            equal<Lvl2>(l_partkey_ciphers[i], p_partkey_ciphers[i], partkey_cres[i], l_partkey_bits, ek, LOGIC);
            equal<Lvl1>(p_brand_ciphers[i], pred_brand3_cipher, brand_cres[i], p_brand_bits, ek, LOGIC);
            HomAND(pred_cres3[i], partkey_cres[i], brand_cres[i], ek, LOGIC);
            equal<Lvl1>(p_container_ciphers[i], pred_container3_cipher[0], container_cres[i], p_container_bits, ek, LOGIC);
            for (uint32_t j = 1; j < 4; j++)
            {
                equal<Lvl1>(p_container_ciphers[i], pred_container3_cipher[j], temp_cres[i], p_container_bits, ek, LOGIC);
                HomOR(container_cres[i], container_cres[i], temp_cres[i], ek, LOGIC);
            }
            HomAND(pred_cres3[i], pred_cres3[i], container_cres[i], ek, LOGIC);
            greater_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity3_left_cipher, quantity_cres[i], l_quantity_bits,
                                     ek, LOGIC);
            HomAND(pred_cres3[i], pred_cres3[i], quantity_cres[i], ek, LOGIC);
            less_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity3_right_cipher, quantity_cres[i], l_quantity_bits, ek,
                                  LOGIC);
            HomAND(pred_cres3[i], pred_cres3[i], quantity_cres[i], ek, LOGIC);
            greater_than_equal<Lvl1>(p_size_ciphers[i], pred_size3_left_cipher, size_cres[i], p_size_bits, ek, LOGIC);
            HomAND(pred_cres3[i], pred_cres3[i], size_cres[i], ek, LOGIC);
            less_than_equal<Lvl1>(p_size_ciphers[i], pred_size3_right_cipher, size_cres[i], p_size_bits, ek, LOGIC);
            HomAND(pred_cres3[i], pred_cres3[i], size_cres[i], ek, LOGIC);
            equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[0], shipmode_cres[i], l_shipmode_bits, ek, LOGIC);
            equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[1], temp_cres[i], l_shipmode_bits, ek, LOGIC);
            HomOR(shipmode_cres[i], shipmode_cres[i], temp_cres[i], ek, LOGIC);
            HomAND(pred_cres3[i], pred_cres3[i], shipmode_cres[i], ek, LOGIC);
            equal<Lvl1>(l_shipinstruct_ciphers[i], pred_shipinstruct_cipher, shipinstruct_cres[i], l_shipinstruct_bits,
                        ek, LOGIC);
            HomAND(pred_cres3[i], pred_cres3[i], shipinstruct_cres[i], ek, LOGIC);
    }
    }
#else
    for (size_t i = 0; i < rows; i++)
    {
        equal<Lvl2>(l_partkey_ciphers[i], p_partkey_ciphers[i], partkey_cres[i], l_partkey_bits, ek, LOGIC);
        equal<Lvl1>(p_brand_ciphers[i], pred_brand1_cipher, brand_cres[i], p_brand_bits, ek, LOGIC);
        HomAND(pred_cres1[i], partkey_cres[i], brand_cres[i], ek, LOGIC);
        // cout << "----------------------------brand----------------------------"  << std::endl;
        // cout << "p_brand_data[" << i <<"]: " <<  p_brand_data[i] << std::endl;
        // cout << "brand_res1[" << i <<"]: " <<  brand_res1[i] << std::endl;
        // uint32_t brand_cres_de = TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(brand_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "brand_cres_de[" << i <<"]: " <<  brand_cres_de << std::endl;

        equal<Lvl1>(p_container_ciphers[i], pred_container1_cipher[0], container_cres[i], p_container_bits, ek, LOGIC);
        for (uint32_t j = 1; j < 4; j++)
        {
            equal<Lvl1>(p_container_ciphers[i], pred_container1_cipher[j], temp_cres[i], p_container_bits, ek, LOGIC);
            HomOR(container_cres[i], container_cres[i], temp_cres[i], ek, LOGIC);
        }
        // cout << "----------------------------container----------------------------"  << std::endl;
        // cout << "p_container_data[" << i <<"]: " <<  p_container_data[i] << std::endl;
        // cout << "container_res1[" << i <<"]: " <<  container_res1[i] << std::endl;
        // uint32_t container_cres_de = TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(container_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "container_cres_de[" << i <<"]: " <<  container_cres_de << std::endl;

        HomAND(pred_cres1[i], pred_cres1[i], container_cres[i], ek, LOGIC);
        // cout << "---------------------------brand-and-container---------------------------"  << std::endl;
        // uint32_t pred_cres1_de1 = TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de_1 : " <<  pred_cres1_de1 << std::endl;

        greater_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity1_left_cipher, quantity_cres[i], l_quantity_bits,
                                 ek, LOGIC);
        // cout << "----------------------------quantity_left----------------------------"  << std::endl;
        // cout << "l_quantity_data[" << i <<"]: " <<  l_quantity_data[i] << std::endl;
        // cout << "quantity_res1[" << i <<"]: " <<  quantity_res1[i] << std::endl;
        // uint32_t quantity_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(quantity_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "quantity_>=_cres_de : " <<  quantity_cres_de << std::endl;

        HomAND(pred_cres1[i], pred_cres1[i], quantity_cres[i], ek, LOGIC);
        // cout << "---------------------------pred_cres1-andd-quantity-left---------------------------"  << std::endl;
        // uint32_t pred_cres1_de2 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de_2 : " <<  pred_cres1_de2 << std::endl;

        less_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity1_right_cipher, quantity_cres[i], l_quantity_bits, ek,
                              LOGIC);
        // cout << "----------------------------quantity_right----------------------------"  << std::endl;
        // quantity_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(quantity_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "quantity_<=_cres_de : " <<  quantity_cres_de << std::endl;

        HomAND(pred_cres1[i], pred_cres1[i], quantity_cres[i], ek, LOGIC);
        // cout << "---------------------------pred_cres1-and-quantity-left---------------------------"  << std::endl;
        // uint32_t pred_cres1_de3 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de_3 : " <<  pred_cres1_de3 << std::endl;

        greater_than_equal<Lvl1>(p_size_ciphers[i], pred_size1_left_cipher, size_cres[i], p_size_bits, ek, LOGIC);
        // cout << "----------------------------size_left----------------------------"  << std::endl;
        // cout << "p_size_data[" << i <<"]: " <<  p_size_data[i] << std::endl;
        // cout << "size_res1[" << i <<"]: " <<  size_res1[i] << std::endl;
        // uint32_t size_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(size_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "size_>=_cres_de : " <<  size_cres_de << std::endl;

        HomAND(pred_cres1[i], pred_cres1[i], size_cres[i], ek, LOGIC);
        // cout << "------------------------pred_cres1-and-size_left----------------------------"  << std::endl;
        // uint32_t pred_cres1_de4 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de_4 : " <<  pred_cres1_de4 << std::endl;

        less_than_equal<Lvl1>(p_size_ciphers[i], pred_size1_right_cipher, size_cres[i], p_size_bits, ek, LOGIC);
        // cout << "----------------------------size_right----------------------------"  << std::endl;
        // size_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(size_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "size_<=_cres_de : " <<  size_cres_de << std::endl;

        HomAND(pred_cres1[i], pred_cres1[i], size_cres[i], ek, LOGIC);
        // cout << "------------------------pred_cres1-and-size_right----------------------------"  << std::endl;
        // uint32_t pred_cres1_de5 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de_5 : " <<  pred_cres1_de5 << std::endl;

        equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[0], shipmode_cres[i], l_shipmode_bits, ek, LOGIC);
        // cout << "----------------------------shipmode_1----------------------------"  << std::endl;
        // cout << "l_shipmode_data[" << i <<"]: " <<  l_shipmode_data[i] << std::endl;
        // cout << "shipmode_res1[" << i <<"]: " <<  shipmode_res1[i] << std::endl;
        // uint32_t shipmode_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(shipmode_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "shipmode_cres_de_1 : " <<  shipmode_cres_de << std::endl;

        equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[1], temp_cres[i], l_shipmode_bits, ek, LOGIC);
        // cout << "----------------------------shipmode_2----------------------------"  << std::endl;
        // shipmode_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(temp_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "shipmode_cres_de_2 : " <<  shipmode_cres_de << std::endl;

        HomOR(shipmode_cres[i], shipmode_cres[i], temp_cres[i], ek, LOGIC);
        // cout << "----------------------------shipmode1_OR_shipmode2----------------------------"  << std::endl;
        // shipmode_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(shipmode_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "shipmode_cres_de : " <<  shipmode_cres_de << std::endl;

        HomAND(pred_cres1[i], pred_cres1[i], shipmode_cres[i], ek, LOGIC);
        // cout << "------------------------pred_cres1-and-shipmode----------------------------"  << std::endl;
        // uint32_t pred_cres1_de6 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de_6 : " <<  pred_cres1_de6 << std::endl;

        equal<Lvl1>(l_shipinstruct_ciphers[i], pred_shipinstruct_cipher, shipinstruct_cres[i], l_shipinstruct_bits, ek,
                    LOGIC);
        // cout << "------------------------shipinstruct----------------------------"  << std::endl;
        // cout << "l_shipinstruct_data[" << i <<"]: " <<  l_shipinstruct_data[i] << std::endl;
        // cout << "shipinstruct_res1[" << i <<"]: " <<  shipinstruct_res1[i] << std::endl;
        // uint32_t shipinstruct_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(shipinstruct_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "shipinstruct_cres_de : " <<  shipinstruct_cres_de << std::endl;

        HomAND(pred_cres1[i], pred_cres1[i], shipinstruct_cres[i], ek, LOGIC);
        // cout << "------------------------pred_cres1-and-shipinstruct----------------------------"  << std::endl;
        // uint32_t pred_cres1_de7 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de_7 : " <<  pred_cres1_de7 << std::endl;


        /////////////////////////////////////////////////////////////
        equal<Lvl2>(l_partkey_ciphers[i], p_partkey_ciphers[i], partkey_cres[i], l_partkey_bits, ek, LOGIC);
        equal<Lvl1>(p_brand_ciphers[i], pred_brand2_cipher, brand_cres[i], p_brand_bits, ek, LOGIC);
        HomAND(pred_cres2[i], partkey_cres[i], brand_cres[i], ek, LOGIC);

        equal<Lvl1>(p_container_ciphers[i], pred_container2_cipher[0], container_cres[i], p_container_bits, ek, LOGIC);
        for (uint32_t j = 1; j < 4; j++)
        {
            equal<Lvl1>(p_container_ciphers[i], pred_container2_cipher[j], temp_cres[i], p_container_bits, ek, LOGIC);
            HomOR(container_cres[i], container_cres[i], temp_cres[i], ek, LOGIC);
        }
        HomAND(pred_cres2[i], pred_cres2[i], container_cres[i], ek, LOGIC);

        greater_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity2_left_cipher, quantity_cres[i], l_quantity_bits,
                                 ek, LOGIC);
        HomAND(pred_cres2[i], pred_cres2[i], quantity_cres[i], ek, LOGIC);
        less_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity2_right_cipher, quantity_cres[i], l_quantity_bits, ek,
                              LOGIC);
        HomAND(pred_cres2[i], pred_cres2[i], quantity_cres[i], ek, LOGIC);

        greater_than_equal<Lvl1>(p_size_ciphers[i], pred_size2_left_cipher, size_cres[i], p_size_bits, ek, LOGIC);
        HomAND(pred_cres2[i], pred_cres2[i], size_cres[i], ek, LOGIC);
        less_than_equal<Lvl1>(p_size_ciphers[i], pred_size2_right_cipher, size_cres[i], p_size_bits, ek, LOGIC);
        HomAND(pred_cres2[i], pred_cres2[i], size_cres[i], ek, LOGIC);

        equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[0], shipmode_cres[i], l_shipmode_bits, ek, LOGIC);
        equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[1], temp_cres[i], l_shipmode_bits, ek, LOGIC);
        HomOR(shipmode_cres[i], shipmode_cres[i], temp_cres[i], ek, LOGIC);

        HomAND(pred_cres2[i], pred_cres2[i], shipmode_cres[i], ek, LOGIC);

        equal<Lvl1>(l_shipinstruct_ciphers[i], pred_shipinstruct_cipher, shipinstruct_cres[i], l_shipinstruct_bits, ek,
                    LOGIC);
        HomAND(pred_cres2[i], pred_cres2[i], shipinstruct_cres[i], ek, LOGIC);


        // //////////////////////////////////////////////////////////////////////
        equal<Lvl2>(l_partkey_ciphers[i], p_partkey_ciphers[i], partkey_cres[i], l_partkey_bits, ek, LOGIC);
        equal<Lvl1>(p_brand_ciphers[i], pred_brand3_cipher, brand_cres[i], p_brand_bits, ek, LOGIC);
        HomAND(pred_cres3[i], partkey_cres[i], brand_cres[i], ek, LOGIC);

        equal<Lvl1>(p_container_ciphers[i], pred_container3_cipher[0], container_cres[i], p_container_bits, ek, LOGIC);
        for (uint32_t j = 1; j < 4; j++)
        {
            equal<Lvl1>(p_container_ciphers[i], pred_container3_cipher[j], temp_cres[i], p_container_bits, ek, LOGIC);
            HomOR(container_cres[i], container_cres[i], temp_cres[i], ek, LOGIC);
        }
        HomAND(pred_cres3[i], pred_cres3[i], container_cres[i], ek, LOGIC);

        greater_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity3_left_cipher, quantity_cres[i], l_quantity_bits,
                                 ek, LOGIC);
        HomAND(pred_cres3[i], pred_cres3[i], quantity_cres[i], ek, LOGIC);
        less_than_equal<Lvl1>(l_quantity_ciphers[i], pred_quantity3_right_cipher, quantity_cres[i], l_quantity_bits, ek,
                              LOGIC);
        HomAND(pred_cres3[i], pred_cres3[i], quantity_cres[i], ek, LOGIC);

        greater_than_equal<Lvl1>(p_size_ciphers[i], pred_size3_left_cipher, size_cres[i], p_size_bits, ek, LOGIC);
        HomAND(pred_cres3[i], pred_cres3[i], size_cres[i], ek, LOGIC);
        less_than_equal<Lvl1>(p_size_ciphers[i], pred_size3_right_cipher, size_cres[i], p_size_bits, ek, LOGIC);
        HomAND(pred_cres3[i], pred_cres3[i], size_cres[i], ek, LOGIC);

        equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[0], shipmode_cres[i], l_shipmode_bits, ek, LOGIC);
        equal<Lvl1>(l_shipmode_ciphers[i], pred_shipmode_cipher[1], temp_cres[i], l_shipmode_bits, ek, LOGIC);
        HomOR(shipmode_cres[i], shipmode_cres[i], temp_cres[i], ek, LOGIC);

        HomAND(pred_cres3[i], pred_cres3[i], shipmode_cres[i], ek, LOGIC);

        equal<Lvl1>(l_shipinstruct_ciphers[i], pred_shipinstruct_cipher, shipinstruct_cres[i], l_shipinstruct_bits, ek,
                    LOGIC);
        HomAND(pred_cres3[i], pred_cres3[i], shipinstruct_cres[i], ek, LOGIC);


        // ///////////////////////////////////////
        HomOR(pred_cres[i], pred_cres1[i], pred_cres2[i], ek, LOGIC);
        HomOR(pred_cres[i], pred_cres[i], pred_cres3[i], ek, ARITHMETIC);

        // pred_cres[i] = pred_cres1[i];
        // uint32_t pred_cres1_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres1[i], pow(2., 29), sk.key.get<Lvl1>());
        // uint32_t pred_cres_de = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        // cout << "pred_cres1_de[" << i <<"]: " <<  pred_cres1_de << std::endl;
        // cout << "pred_cres_de[" << i <<"]: " <<  pred_cres_de << std::endl;
    }
#endif
    end = std::chrono::system_clock::now();
    std::vector<uint32_t> pred_cres_total_de(rows), pred_cres_part_de(rows), pred_cres1_de(rows), pred_cres2_de(rows),
                          pred_cres3_de(rows), pred_cres4_de(rows), pred_cres5_de(rows);
    // for (size_t i = 0; i < rows; i++)
    // {
    //     pred_cres_de[i] = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres[i], pow(2., 31), sk.key.get<Lvl1>());
    //     pred_cres1_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres1[i], sk.key.lvl1);
    //     pred_cres2_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres2[i], sk.key.lvl1);
    //     pred_cres3_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres3[i], sk.key.lvl1);
    //     pred_cres4_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres4[i], sk.key.lvl1);
    //     pred_cres5_de[i] =  TFHEpp::tlweSymDecrypt<Lvl1>(pred_cres5[i], sk.key.lvl1);
    //}

    size_t error_time = 0;
    std::vector<uint32_t> pred_cres_de(rows);

    uint32_t rlwe_scale_bits = 29;
    for (size_t i = 0; i < rows; i++)
    {
        TFHEpp::ari_rescale(pred_cres[i], pred_cres[i], rlwe_scale_bits, ek);
    }
    for (size_t i = 0; i < rows; i++)
    {
        pred_cres_de[i] = TFHEpp::tlweSymInt32Decrypt<Lvl1>(pred_cres[i], pow(2., 29), sk.key.get<Lvl1>());
        //cout << "pred_cres_de[" << i <<"]: " <<  pred_cres_de[i] << std::endl;
    }
    for (size_t i = 0; i < rows; i++)
    {
        error_time += (pred_cres_de[i] == pred_res[i]) ? 0 : 1;
    }
    cout << "Predicate Evaluaton Time (s): " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).
        count() / 1000 << std::endl;
    cout << "Predicate Error: " << error_time << std::endl;
    filter_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}


void aggregation(std::vector<TLWELvl1>& pred_cres, std::vector<uint32_t>& pred_res, size_t tfhe_n,
                 std::vector<double>& extendedprice_data, std::vector<double>& discount_data,
                 size_t rows, TFHESecretKey& sk, double& aggregation_time)
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
    parms.set_coeff_modulus(seal::CoeffModulus::Create(poly_modulus_degree, {
                                                           59, 42, 42, 42, 42, 42, 42, 42, 45, 45, 45, 45, 45, 45, 45,
                                                           45, 45, 45, 45, 59
                                                       }));
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
    LWEsToRLWEKeyGen(pre_key, std::pow(2., modulus_bits), seal_secret_key, sk, tfhe_n, ckks_encoder, encryptor,
                     context);


    // P0: Chunked parallel LWEsToRLWE + HomRound
    std::cout << "Starting Conversion..." << std::endl;

    struct ChunkAggRes {
        seal::Ciphertext cipher;
        double time_ms = 0.0;
    };

    auto convert_chunk = [&](size_t begin, size_t count) -> ChunkAggRes {
        ChunkAggRes r;
        seal::CKKSEncoder enc_local(context);
        seal::Evaluator eval_local(context);
        seal::Decryptor dec_local(context, seal_secret_key);
        // Share pre_key, relin_keys, galois_keys as read-only
        std::vector<TLWELvl1> slice(pred_cres.begin() + begin, pred_cres.begin() + begin + count);
        auto s0 = std::chrono::system_clock::now();
        LWEsToRLWE(r.cipher, slice, pre_key, scale, std::pow(2., modq_bits), std::pow(2., modulus_bits - modq_bits),
                   enc_local, galois_keys, relin_keys, eval_local, context);
        HomRound(r.cipher, r.cipher.scale(), enc_local, relin_keys, eval_local, dec_local, context);
        auto s1 = std::chrono::system_clock::now();
        r.time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(s1 - s0).count();
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

    // Merge chunk ciphertexts
    seal::Ciphertext result = chunks.front().cipher;
    for (size_t c = 1; c < chunk_count; ++c)
    {
        evaluator.add_inplace(result, chunks[c].cipher);
    }
    aggregation_time = 0.0;
    for (auto &c : chunks) aggregation_time += c.time_ms;

    // Optional correctness check (avg error)
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
            price_discount[i] = extendedprice_data[i] * (100 - discount_data[i]);
        }
    }
#else
    for (size_t i = 0; i < rows; i++)
    {
        price_discount[i] = extendedprice_data[i] * (100 - discount_data[i]);
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
            plain_result += extendedprice_data[i] * (100 - discount_data[i]) * pred_res[i];
        }
    }
#else
    for (size_t i = 0; i < rows; i++)
    {
        plain_result += extendedprice_data[i] * (100 - discount_data[i]) * pred_res[i];
    }
#endif

    cout << "aggregation_time: " << aggregation_time / 1000 << " s" << endl;
    cout << "Plain_result: " << plain_result << endl;
    cout << "Encrypted query result: " << std::round(agg_result[0]) << endl;
}

void query_evaluation(size_t rows)
{
    cout << "-------tpch_q19-------" << endl;
    cout << "data rows: " << rows << endl;
    TFHESecretKey sk;

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    // Generate database
    std::vector<uint32_t> p_brand_data(rows);
    std::vector<uint32_t> p_container_data(rows);
    std::vector<uint32_t> l_quantity_data(rows);
    std::vector<uint32_t> p_size_data(rows);
    std::vector<uint32_t> l_shipmode_data(rows);
    std::vector<uint32_t> l_shipinstruct_data(rows);
    std::vector<uint64_t> l_partkey_data(rows);
    std::vector<uint64_t> p_partkey_data(rows);

    uint32_t p_brand_bits = 6, p_container_bits = 4, l_quantity_bits = 6, p_size_bits = 6, l_shipmode_bits = 4,
             l_shipinstruct = 4, discount_bits = 4, l_partkey_bits = 16, p_partkey_bits = 16;
    uint32_t p_brand_scale_bits, p_container_scale_bits, l_quantity_scale_bits, p_size_scale_bits, l_shipmode_scale_bits
             , l_shipinstruct_scale_bits, discount_scale_bits, l_partkey_scale_bits, p_partkey_scale_bits;
    uniform_int_distribution<Lvl1::T> brand_message(0, 55);
    uniform_int_distribution<Lvl1::T> container_message(0, 11);
    uniform_int_distribution<Lvl1::T> quantity_message(10, 20);
    uniform_int_distribution<Lvl1::T> size_message(1, 15);
    uniform_int_distribution<Lvl1::T> shipmode_message(0, 11);
    uniform_int_distribution<Lvl1::T> shipinstruct_message(0, 7);
    uniform_int_distribution<Lvl1::T> discount_message(0, (1 << discount_bits) - 1);
    uniform_int_distribution<Lvl1::T> partkey_message(0, (1 << l_partkey_bits) - 1);
    uniform_int_distribution<Lvl1::T> possibility_message(0, 9);


    p_brand_scale_bits = std::numeric_limits<Lvl2::T>::digits - p_brand_bits - 1;
    p_container_scale_bits = std::numeric_limits<Lvl1::T>::digits - p_container_bits - 1;
    l_quantity_scale_bits = std::numeric_limits<Lvl1::T>::digits - l_quantity_bits - 1;
    p_size_scale_bits = std::numeric_limits<Lvl1::T>::digits - p_size_bits - 1;
    l_shipmode_scale_bits = std::numeric_limits<Lvl1::T>::digits - l_shipmode_bits - 1;
    l_shipinstruct_scale_bits = std::numeric_limits<Lvl1::T>::digits - l_shipinstruct - 1;

    uint32_t possibility = 0;
    for (size_t i = 0; i < rows; i++)
    {
        p_brand_data[i] = brand_message(engine);
        //p_brand_data[i] = 1;
        p_container_data[i] = container_message(engine);
        l_quantity_data[i] = quantity_message(engine);
        p_size_data[i] = size_message(engine);
        l_shipmode_data[i] = shipmode_message(engine);
        //l_shipmode_data[i] = 0;
        l_shipinstruct_data[i] = shipinstruct_message(engine);

        possibility = possibility_message(engine);
        if (possibility < 2)
        {
            l_partkey_data[i] = partkey_message(engine);
            p_partkey_data[i] = partkey_message(engine);
        }
        else
        {
            l_partkey_data[i] = partkey_message(engine);
            p_partkey_data[i] = l_partkey_data[i];
        }
    }

    std::vector<double> extendedprice_data(rows), discount_data(rows);
    std::uniform_int_distribution<uint64_t> extendedprice_message(1, 5);
    for (size_t i = 0; i < rows; i++)
    {
        extendedprice_data[i] = (double)extendedprice_message(engine);
        discount_data[i] = (double)discount_message(engine);
    }

    //quanlity_data[0] = 1;
    // discount_data[0] = 9;
    // ship_data[0] = 21215;
    double filter_time, aggregation_time;
    std::vector<TLWELvl1> pred_cres(rows);
    std::vector<uint32_t> pred_res(rows, 0);

    predicate_evaluation(pred_cres, pred_res, rows, p_brand_data, p_container_data, l_quantity_data, p_size_data,
                         l_shipmode_data, l_shipinstruct_data, l_partkey_data, p_partkey_data, sk, filter_time);
    aggregation(pred_cres, pred_res, Lvl1::n, extendedprice_data, discount_data, rows, sk, aggregation_time);


    cout << "filter_time: " << filter_time / 1000 << " s" << endl;
    cout << "aggregation_time: " << aggregation_time / 1000 << " s" << endl;
    cout << "End-to-End Time: " << (filter_time + aggregation_time) / 1000 << " s" << endl;
}

int main()
{
    // query_evaluation(1000);
    // cout << " " << endl;
    // query_evaluation(6000);
    // cout << " " << endl;

    // cout << " " << endl;
    //query_evaluation(32768);

    // query_evaluation(256);
    // cout << " " << endl;
    // query_evaluation(512);
    // cout << " " << endl;

    // cout << " " << endl;
    // query_evaluation(2048);
    // cout << " " << endl;
    query_evaluation(1024);
    query_evaluation(4096);
    query_evaluation(8192);
    query_evaluation(16384);
}
