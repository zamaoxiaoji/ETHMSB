#include <iostream>
#include <chrono>
#include <random>
#include<bitset>
#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"

#include "HEDB/comparison/tfhepp_utils.h"
#include <gatebootstrapping.hpp>
#include "detwfa.hpp"

#include <cassert>
#include <tfhe++.hpp>


#include <bits/stdint-uintn.h>
#include <iostream>
#include "circuitbootstrapping.hpp"


#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"
#include "HEDB/conversion/repack.h"

using namespace std;
using namespace HEDB;
using namespace TFHEpp;

int test(uint32_t num_test)
{
    //constexpr uint32_t num_test = 10;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl21param;

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    std::vector<std::array<uint8_t, privksP::targetP::n>> pa(num_test);
    std::vector<std::array<typename privksP::targetP::T, privksP::targetP::n>>
        pmu(num_test);
    std::vector<uint8_t> pones(num_test);
    std::array<bool, privksP::targetP::n> pres;

    //std::array<typename privksP::targetP::T, privksP::targetP::n> pres_uint;
    Polynomial<privksP::targetP> pres_uint;


    for (std::array<uint8_t, privksP::targetP::n> &i : pa)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < privksP::targetP::n; j++)
            pmu[i][j] = pa[i][j] ? privksP::targetP::μ : -privksP::targetP::μ;
    for (int i = 0; i < num_test; i++) pones[i] = true;
    std::vector<TFHEpp::TRLWE<typename privksP::targetP>> ca(num_test);
    std::vector<TFHEpp::TLWE<typename iksP::domainP>> cones(num_test);
    std::vector<TFHEpp::TRGSWFFT<typename privksP::targetP>> bootedTGSW(
        num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = TFHEpp::trlweSymEncrypt<typename privksP::targetP>(
            pmu[i], privksP::targetP::α,
            sk->key.get<typename privksP::targetP>());
    cones = TFHEpp::bootsSymEncrypt(pones, *sk);

    std::chrono::system_clock::time_point start, end;
    #ifdef USE_PERF
    ProfilerStart("cb.prof");
    #endif
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW[test],
                                                            cones[test], ek);
    }
    end = std::chrono::system_clock::now();
    #ifdef USE_PERF
    ProfilerStop();
    #endif
    for (int test = 0; test < num_test; test++) {
        TFHEpp::trgswfftExternalProduct<typename privksP::targetP>(
            ca[test], ca[test], bootedTGSW[test]);
        pres = TFHEpp::trlweSymDecrypt<typename privksP::targetP>(ca[test],
                                                                  sk->key.lvl1);
        
        pres_uint = TFHEpp::trlweSymInt32Decrypt<typename privksP::targetP>(ca[test], 
                            pow(2., 27), sk->key.lvl1);

        std::cout << std::endl;
        std::cout << "pres: " << std::endl;
        for(int i = 0 ; i< privksP::targetP::n ; i++){
            std::cout << pres_uint[i]*pow(2., 27) ;
        }

        std::cout << std::endl;
        std::cout << "pa: " << std::endl;
        for(int i = 0 ; i< privksP::targetP::n ; i++){
            std::cout << (int)pa[test][i] ;
        }

        std::cout << std::endl;
        std::cout << "pmu: " << std::endl;
        for(int i = 0 ; i< privksP::targetP::n ; i++){
            std::cout << pmu[test][i] ;
        }

        for (int i = 0 ; i < privksP::targetP::n; i++)
            assert(pres[i] == pa[test][i]);
    }
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;
}

void lwe_to_gsw_test( uint32_t k, int num_test){

    using P = Lvl1;
    // using iksP = Lvl20;
    // using bkP = Lvl02;

    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl21param;

    uint32_t scale_bits = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<P::T>::digits - k - 1;
    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);


    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, 1);
    

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    TLWELvl1 res;
    int choose = 1;
    int d_res = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    std::chrono::system_clock::time_point start, end;
    
    TFHEpp::Polynomial<P> p0, p1, pres;

    for(int time = 0; time < num_test; time++){
        for (typename P::T &i : p0) i = message(engine);
        for (typename P::T &i : p1) i = message(engine);

        TFHEpp::TRLWE<P> c0 = TFHEpp::trlweSymIntEncrypt<P>(p0, P::α, sk->key.get<P>());
        TFHEpp::TRLWE<P> c1 = TFHEpp::trlweSymIntEncrypt<P>(p1, P::α, sk->key.get<P>());
        // 0 或者 1 scale_bits = 29 ，进行加密；
        TLWE<typename iksP::domainP> k_1 = TFHEpp::tlweSymInt32Encrypt<P>( choose, P::α, pow(2., scale_bits_2), sk->key.get<P>());
        //TLWE<typename iksP::domainP> k_1 = TFHEpp::tlweSymEncrypt<Lvl1>(Lvl1::μ, Lvl1::α, sk->key.get<P>());
        
        
        TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;
        TFHEpp::TLWE<typename iksP::domainP> encaddress;
        std::chrono::system_clock::time_point start, end;

        ///////////////////////////
        TFHEpp::TRLWE<TFHEpp::lvl1param> r_k_1;
        TRGSW<typename privksP::targetP> gsw;
        
        TFHEpp::CircuitBootstrapping<iksP, bkP, privksP>(gsw, k_1, ek);
        //TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl11param>(r_k_1, k_1, *iksk);
        bootedTGSW = ApplyFFT2trgsw<P>(gsw);
        for (TRLWE<P> &trlwe0 : gsw){
            Polynomial<P> de_rlwe =  TFHEpp::trlweSymInt32Decrypt<P>(trlwe0, pow(2., 10), sk->key.get<P>());
            std::cout << "----------de_boot_gsw: " << std::endl;
            for(int i = 0 ; i<P::n ; i++){
                std::cout << de_rlwe[i]<< ' ';
            }
            std::cout << std::endl;
        } 

        // 将TLWE转化为TGSW类型的密文
        start = std::chrono::system_clock::now();
        //TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
        end = std::chrono::system_clock::now();
        lwe_to_gsw_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        // 直接将首项为0 或者 1的向量，加密为TGSW类型的密文
        const Polynomial<TFHEpp::lvl1param> plainpoly = {static_cast<typename lvl1param::T>(choose)};
        ////////////////////////////////
        std::cout << "plainpoly: " << std::endl;
            for(int i = 0 ; i<P::n ; i++){
                std::cout << plainpoly[i]<< ' ';
            }
            std::cout << std::endl;

        TRGSWFFT<lvl1param> trgswfft = trgswfftSymEncrypt<lvl1param>(plainpoly, lvl1param::α, sk->key.get<P>());
        
        //////////////////////////////////////
        TRLWE<P> rlwe = TFHEpp::trlweSymInt32Encrypt<P>(plainpoly,lvl1param::α, pow(2., 29), sk->key.get<P>());
        Polynomial<P> de_rlwe0 =  TFHEpp::trlweSymInt32Decrypt<P>(rlwe, pow(2., 29), sk->key.get<P>());
        std::cout << "de_rlwe0: " << std::endl;
        for(int i = 0 ; i<P::n ; i++){
            std::cout << de_rlwe0[i]<< ' ';
        }
        std::cout << std::endl;

        ////////////////////////////////////////////
        TRGSW<lvl1param> trgsw = trgswSymEncrypt<lvl1param>(plainpoly, lvl1param::α, sk->key.get<P>());
        trgswfft =  ApplyFFT2trgsw<P>(trgsw);
        
        for (TRLWE<P> &trlwe1 : trgsw){
            Polynomial<P> de_rlwe =  TFHEpp::trlweSymInt32Decrypt<P>(trlwe1, pow(2., 10), sk->key.get<P>());
            std::cout << "de_enc_gsw: " << std::endl;
            for(int i = 0 ; i<P::n ; i++){
                std::cout << de_rlwe[i]<< ' ';
            }
            std::cout << std::endl;
        } 


        //uint32_t plain_bits = 2;
        uint32_t plain_bits = k;
        std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
        
        Lvl1::T a = message_2(engine);
        Lvl1::T b = message_2(engine);
        
    
        // std::cout << "a : " << a << std::endl;
        // std::cout << "b : " << b << std::endl;
        
        TLWELvl1 ca = TFHEpp::tlweSymInt32Encrypt<P>(a, P::α, pow(2., scale_bits_2), sk->key.get<P>());
        TLWELvl1 cb = TFHEpp::tlweSymInt32Encrypt<P>(b, P::α, pow(2., scale_bits_2), sk->key.get<P>());
        
        TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl11param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl11param>();
        TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl11param>(*iksk, *sk);
        
        TFHEpp::TRLWE<TFHEpp::lvl1param> r_ca, r_cb;
        
        TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl11param>(r_ca, ca, *iksk);
        start = std::chrono::system_clock::now();
        TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl11param>(r_cb, cb, *iksk);
        end = std::chrono::system_clock::now();
        lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        
        c0 = r_ca;
        c1 = r_cb;

        // 做cmux选择
        TRLWELvl1 cres;
        
        start = std::chrono::system_clock::now();
        CMUXFFT<lvl1param>(cres, bootedTGSW, c0, c1);
        //CMUXFFT<lvl1param>(cres, trgswfft, c0, c1);
        end = std::chrono::system_clock::now();
        cmux_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        Polynomial<P> de_cres = TFHEpp::trlweSymInt32Decrypt<P>(cres, pow(2., scale_bits_2), sk->key.get<P>());
        Polynomial<P> de_c0 = TFHEpp::trlweSymInt32Decrypt<P>(c0, pow(2., scale_bits_2), sk->key.get<P>());
        Polynomial<P> de_c1 = TFHEpp::trlweSymInt32Decrypt<P>(c1, pow(2., scale_bits_2), sk->key.get<P>());

        TLWE<Lvl1> extract_cres;
        SampleExtractIndex<lvl1param>(extract_cres, cres, 0);
        uint32_t de_extract_cres = TFHEpp::tlweSymInt32Decrypt<P>(extract_cres, pow(2., scale_bits_2), sk->key.get<P>());

        bool is_c0 = true;
        bool is_c1 = true;

        std::cout << "de_cres: " << std::endl;
        for(int i = 0 ; i<P::n ; i++){
            std::cout << de_cres[i]<< ' ';
        }
        std::cout << std::endl;
        std::cout << "de_c0: " << std::endl;
        for(int i = 0 ; i<P::n ; i++){
            std::cout << de_c0[i]<< ' ';
        }
        std::cout << std::endl;
        std::cout << "de_c1: " << std::endl;
        for(int i = 0 ; i<P::n ; i++){
            std::cout << de_c1[i]<< ' ';
        }
        std::cout << std::endl;

        for(int i = 0 ; i<P::n ; i++){
            if(de_cres[i] != de_c0[i]){
                is_c0 = false;
            }
            if(de_cres[i] != de_c1[i]){
                is_c1 = false;
            }
        }

        // if(de_extract_cres == a){
        //     is_c0 = true;
        //     is_c1 = false;
        // }
        // else if(de_extract_cres == b){
        //     is_c0 = false;
        //     is_c1 = true;
        // }
        // else{
        //     is_c0 = false;
        //     is_c1 = false;
        // }

        // std::cout << "is_c0: " << is_c0 <<  std::endl;
        // std::cout << "is_c1: " << is_c1 <<  std::endl;
        if( !((choose == 1 && is_c0) || (choose == 0 && is_c1)) ){
            error_times ++;
            std::cout << "a : " << a << std::endl;
            std::cout << "b : " << b << std::endl;
            //std::cout << "de_extract_cres : " << de_extract_cres << std::endl;

            // std::cout << "de_cres: " << std::endl;
            // for(int i = 0 ; i<P::n ; i++){
            //     std::cout << de_cres[i]<< ' ';
            // }
            // std::cout << std::endl;
            // std::cout << "de_c0: " << std::endl;
            // for(int i = 0 ; i<P::n ; i++){
            //     std::cout << de_c0[i]<< ' ' ;
            // }
            // std::cout << std::endl;
            // std::cout << "de_c1: " << std::endl;
            // for(int i = 0 ; i<P::n ; i++){
            //     std::cout << de_c1[i]<< ' ' ;
            // }
            // std::cout << std::endl;

            std::cout << "is_c0: " << is_c0 <<  std::endl;
            std::cout << "is_c1: " << is_c1 <<  std::endl;
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "lwe_to_gsw_time: " << lwe_to_gsw_time / num_test << "ms" <<  std::endl;
    std::cout << "lwe_to_rlwe_time: " << lwe_to_rlwe_time / num_test << "ms" <<  std::endl;
    std::cout << "my_and_time: " << cmux_time / num_test << "ms" <<  std::endl;
}

void lwe_to_gsw_test2( uint32_t k, int num_test){

    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t plain_bits = k;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - k - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 2));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    std::chrono::system_clock::time_point start, end;
    
    Lvl1::T p0 , p1;
    TLWELvl2 c0, c1 ;
    TFHEpp::TRLWE<P> r_c0, r_c1, r_cres;
    Polynomial<P> d_r_cres, d_r_c0, d_r_c1;
    TLWE<P> extract_cres;
    P::T d_extract_cres;

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;
    TRGSWFFT<P> trgswfft;
    
    int choose = 1;
    const Polynomial<P> plainpoly = {static_cast<typename P::T>(choose)};

    bool is_c0,is_c1;
    
    for(int time = 0; time < num_test; time++){
        
        // p0 = 901;
        // p1 = 900;
        p0 = message_2(engine);
        // p1 = message_2(engine);
        p1 = p0 +2;
        std::cout << "p0 : " << p0 << std::endl;
        std::cout << "p1 : " << p1 << std::endl;
        
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits_2), sk->key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits_2), sk->key.get<P>());

        // 0 或者 1 scale_bits = 29 ，进行加密；
        //k_1 = TFHEpp::tlweSymInt32Encrypt<iksP::domainP>( choose, iksP::domainP::α, pow(2., 29), sk->key.get<iksP::domainP>());
        
        greater_than<P>(c0, c1, k_1, k, ek, LOGIC);
        //my_greater_than<Lvl2>(c0, c1, k_1, k, ek, LOGIC, 29);

        // 将TLWE转化为TGSW类型的密文
        start = std::chrono::system_clock::now();
        TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
        end = std::chrono::system_clock::now();
        lwe_to_gsw_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        // 直接将首项为0 或者 1的向量，加密为TGSW类型的密文
        //trgswfft = trgswfftSymEncrypt<P>(plainpoly, P::α, sk->key.get<P>());
        

        // TFHEpp::TLWEToTRLWE<P>(c0 , r_c0);
        // Polynomial<P> de_rlwe0 =  TFHEpp::trlweSymInt32Decrypt<P>(r_c0, pow(2., scale_bits_2), sk->key.get<P>());
        // std::cout << "TLWEToTRLWE: " << std::endl;
        // for(int i = 0 ; i<P::n ; i++){
        //     std::cout << de_rlwe0[i]<< ' ';
        // }
        // std::cout << std::endl;

        // LWE密文转RLWE
        
        TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl22param>(r_c0, c0, *iksk);
        start = std::chrono::system_clock::now();
        TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl22param>(r_c1, c1, *iksk);
        end = std::chrono::system_clock::now();
        lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        // de_rlwe0 =  TFHEpp::trlweSymInt32Decrypt<P>(r_c0, pow(2., scale_bits_2), sk->key.get<P>());
        // std::cout << "TLWE2TRLWEIKS: " << std::endl;
        // for(int i = 0 ; i<P::n ; i++){
        //     std::cout << de_rlwe0[i]<< ' ';
        // }
        // std::cout << std::endl;
        
        // 做cmux选择
        start = std::chrono::system_clock::now();
        CMUXFFT<P>(r_cres, bootedTGSW, r_c0, r_c1);
        //CMUXFFT<lvl1param>(cres, trgswfft, r_c0, r_c1);
        end = std::chrono::system_clock::now();
        cmux_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_r_cres = TFHEpp::trlweSymInt32Decrypt<P>(r_cres, pow(2., scale_bits_2), sk->key.get<P>());
        d_r_c0 = TFHEpp::trlweSymInt32Decrypt<P>(r_c0, pow(2., scale_bits_2), sk->key.get<P>());
        d_r_c1 = TFHEpp::trlweSymInt32Decrypt<P>(r_c1, pow(2., scale_bits_2), sk->key.get<P>());

        
        // SampleExtractIndex<P>(extract_cres, r_cres, 0);
        // d_extract_cres = TFHEpp::tlweSymInt32Decrypt<P>(extract_cres, pow(2., scale_bits_2), sk->key.get<P>());

        is_c0 = true;
        is_c1 = true;

        // std::cout << "d_r_cres: " << std::endl;
        // for(int i = 0 ; i<P::n ; i++){
        //     std::cout << d_r_cres[i]<< ' ';
        // }
        // std::cout << std::endl;
        // std::cout << "d_r_c0: " << std::endl;
        // for(int i = 0 ; i<P::n ; i++){
        //     std::cout << d_r_c0[i]<< ' ';
        // }
        // std::cout << std::endl;
        // std::cout << "d_r_c1: " << std::endl;
        // for(int i = 0 ; i<P::n ; i++){
        //     std::cout << d_r_c1[i]<< ' ';
        // }
        // std::cout << std::endl;

        for(int i = 0 ; i<P::n ; i++){
            if(d_r_cres[i] != d_r_c0[i]){
                is_c0 = false;
            }
            if(d_r_cres[i] != d_r_c1[i]){
                is_c1 = false;
            }
        }

        // if(de_extract_cres == p0){
        //     is_c0 = true;
        //     is_c1 = false;
        // }
        // else if(d_extract_cres == p1){
        //     is_c0 = false;
        //     is_c1 = true;
        // }
        // else{
        //     is_c0 = false;
        //     is_c1 = false;
        // }

        // std::cout << "is_c0: " << is_c0 <<  std::endl;
        // std::cout << "is_c1: " << is_c1 <<  std::endl;
        if( !((choose == 1 && is_c0) || (choose == 0 && is_c1)) ){
            error_times ++;
            std::cout << "p0 : " << p0 << std::endl;
            std::cout << "p1 : " << p1 << std::endl;
            //std::cout << "d_extract_cres : " << d_extract_cres << std::endl;

            std::cout << "d_r_cres: " << std::endl;
            for(int i = 0 ; i<P::n ; i++){
                std::cout << d_r_cres[i]<< ' ';
            }
            std::cout << std::endl;
            std::cout << "d_r_c0: " << std::endl;
            for(int i = 0 ; i<P::n ; i++){
                std::cout << d_r_c0[i]<< ' ' ;
            }
            std::cout << std::endl;
            std::cout << "d_r_c1: " << std::endl;
            for(int i = 0 ; i<P::n ; i++){
                std::cout << d_r_c1[i]<< ' ' ;
            }
            std::cout << std::endl;

            // std::cout << "is_c0: " << is_c0 <<  std::endl;
            // std::cout << "is_c1: " << is_c1 <<  std::endl;
        }
        
    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "lwe_to_gsw_time: " << lwe_to_gsw_time / num_test << "ms" <<  std::endl;
    std::cout << "lwe_to_rlwe_time: " << lwe_to_rlwe_time / num_test << "ms" <<  std::endl;
    std::cout << "cmux_time: " << cmux_time / num_test << "ms" <<  std::endl;
    

}

void lwe_to_gsw_test3( int num_test){
    std::cout << "-----lwe_to_gsw_test3(cmux test)--------" << std::endl;
    std::cout << "  num_test: " << num_test <<std::endl;
    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t k =8;
    uint32_t plain_bits = k;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - k - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    // TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    // TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 2));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    std::chrono::system_clock::time_point start, end;
    
    Lvl1::T p0 , p1;
    TLWELvl2 c0, c1 ;
    TFHEpp::TRLWE<P> r_c0, r_c1, r_cres;
    Polynomial<P> d_r_cres, d_r_c0, d_r_c1;
    TLWE<P> extract_cres;
    P::T d_extract_cres;

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;
    TRGSWFFT<P> trgswfft;
    
    int choose = 1;
    const Polynomial<P> plainpoly = {static_cast<typename P::T>(choose)};

    bool is_c0,is_c1;
    
    for(int time = 0; time < num_test; time++){
        
        p0 = 15;
        p1 = 20;
        // p0 = message_2(engine);
        // p1 = message_2(engine);
        
        // std::cout << "p0 : " << p0 << std::endl;
        // std::cout << "p1 : " << p1 << std::endl;
        
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits_2), sk->key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits_2), sk->key.get<P>());

        // 0 或者 1 scale_bits = 29 ，进行加密；
        //k_1 = TFHEpp::tlweSymInt32Encrypt<iksP::domainP>( choose, iksP::domainP::α, pow(2., 29), sk->key.get<iksP::domainP>());
        
        greater_than<Lvl2>(c0, c1, k_1, k, ek, LOGIC);
        //my_greater_than<Lvl2>(c0, c1, k_1, k, ek, LOGIC, 29);

        // 将TLWE转化为TGSW类型的密文
        start = std::chrono::system_clock::now();
        TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
        end = std::chrono::system_clock::now();
        lwe_to_gsw_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        
        TFHEpp::TLWEToTRLWE<P>(c0 , r_c0);
        start = std::chrono::system_clock::now();
        //TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl22param>(r_c1, c1, *iksk);
        TFHEpp::TLWEToTRLWE<P>(c1 , r_c1);
        end = std::chrono::system_clock::now();
        lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        
        // // 做cmux选择
        start = std::chrono::system_clock::now();
        CMUXFFT<P>(r_cres, bootedTGSW, r_c0, r_c1);
        end = std::chrono::system_clock::now();
        cmux_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        SampleExtractIndex<P>(extract_cres, r_cres, 0);
        d_extract_cres = TFHEpp::tlweSymInt32Decrypt<P>(extract_cres, pow(2., scale_bits_2), sk->key.get<P>());
        
        
        if( !((p0 > p1 && d_extract_cres == p0) || (p0 <= p1 && d_extract_cres == p1)) ){
            error_times ++;
            std::cout << "p0 : " << p0 << std::endl;
            std::cout << "p1 : " << p1 << std::endl;
            std::cout << "d_extract_cres : " << d_extract_cres << std::endl;
        }
        

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "lwe_to_gsw_time: " << lwe_to_gsw_time / num_test << "ms" <<  std::endl;
    std::cout << "lwe_to_rlwe_time: " << lwe_to_rlwe_time / num_test << "ms" <<  std::endl;
    std::cout << "cmux_time: " << cmux_time / num_test << "ms" <<  std::endl;
    

}

void n_num_find_max_test( int num_test, int arr_num){

    std::cout << "-----n_num_find_max_test(Lvl2)--------" << std::endl;
    std::cout << "  num_test: " << num_test <<std::endl;
    std::cout << "  arr_num : " << arr_num <<std::endl;
    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t k= 8;
    uint32_t plain_bits = k;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - k - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    double extract_time = 0;
    double end2end_time = 0;
    std::chrono::system_clock::time_point start, end, start2, end2;
    
    Lvl1::T p_result = 0;
    Lvl1::T p_arr[arr_num];
    std::vector<TLWELvl2> c_arr(arr_num);
    std::vector<TFHEpp::TRLWE<P>> r_c_arr(arr_num);

    TLWE<P> cres;
    TFHEpp::TRLWE<P> r_cres;
    P::T d_cres;

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;

    
    for(int time = 0; time < num_test; time++){
        p_result = 0;
        for (size_t i = 0; i < arr_num; ++i) {
            p_arr[i] = message_2(engine);
            if(p_arr[i] > p_result){
                p_result = p_arr[i];
            }
        }

        
        for (size_t i = 0; i < arr_num; ++i) {
            c_arr[i] = TFHEpp::tlweSymInt32Encrypt<P>(p_arr[i], P::α, pow(2., scale_bits_2), sk->key.get<P>());

            // LWE密文转RLWE
            start = std::chrono::system_clock::now();
            //TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl22param>(r_c_arr[i], c_arr[i], *iksk);
            TFHEpp::TLWEToTRLWE<P>(c_arr[i] , r_c_arr[i]);
            end = std::chrono::system_clock::now();
            lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }
        
        
        cres = c_arr[0];
        r_cres = r_c_arr[0];
        start2 = std::chrono::system_clock::now();
        for(int i = 1; i< arr_num; i++){
            
            //greater_than<Lvl2>(cres, c_arr[i], k_1, k, ek, LOGIC);
            my_greater_than<Lvl2>(cres, c_arr[i], k_1,k, ek, LOGIC, 29);
            //将TLWE转化为TGSW类型的密文
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);

            // 做cmux选择
            CMUXFFT<P>(r_cres, bootedTGSW, r_cres, r_c_arr[i]);

            // 密文提取
            SampleExtractIndex<P>(cres, r_cres, 0);
           
        }
        end2 = std::chrono::system_clock::now();
        end2end_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        

        d_cres = TFHEpp::tlweSymInt32Decrypt<P>(cres, pow(2., scale_bits_2), sk->key.get<P>());
        
        if(p_result != d_cres){
            error_times ++;
            std::cout << "p_arr: "  ;
            for(int i =0 ;i< arr_num; i++){
                std::cout << p_arr[i] << "  ";
            }
            std::cout << "  " <<  std::endl;

            std::cout << "p_result: " << p_result <<  std::endl;
            std::cout << "d_cres: " << d_cres <<  std::endl;
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "end2end_time: " << end2end_time/(1000 *num_test)  + lwe_to_rlwe_time / (1000000 * num_test)<< "s" <<  std::endl;
}

void n_num_find_min_test( int num_test, int arr_num){

    std::cout << "-----n_num_find_min_test(Lvl2)--------" << std::endl;
    std::cout << "  num_test: " << num_test <<std::endl;
    std::cout << "  arr_num : " << arr_num <<std::endl;
    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t k = 8;
    uint32_t plain_bits = k;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - k - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    double extract_time = 0;
    double end2end_time = 0;
    std::chrono::system_clock::time_point start, end, start2, end2;
    
    Lvl1::T p_result = (1 << plain_bits);
    Lvl1::T p_arr[arr_num];
    std::vector<TLWELvl2> c_arr(arr_num);
    std::vector<TFHEpp::TRLWE<P>> r_c_arr(arr_num);

    TLWE<P> cres;
    TFHEpp::TRLWE<P> r_cres;
    P::T d_cres;

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;

    
    for(int time = 0; time < num_test; time++){
        p_result = (1 << plain_bits);
        for (size_t i = 0; i < arr_num; ++i) {
            p_arr[i] = message_2(engine);
            if(p_arr[i] < p_result){
                p_result = p_arr[i];
            }
        }

        
        for (size_t i = 0; i < arr_num; ++i) {
            c_arr[i] = TFHEpp::tlweSymInt32Encrypt<P>(p_arr[i], P::α, pow(2., scale_bits_2), sk->key.get<P>());

            // LWE密文转RLWE
            start = std::chrono::system_clock::now();
            //TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl22param>(r_c_arr[i], c_arr[i], *iksk);
            TFHEpp::TLWEToTRLWE<P>(c_arr[i] , r_c_arr[i]);
            end = std::chrono::system_clock::now();
            lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }
        
        
        cres = c_arr[0];
        r_cres = r_c_arr[0];
        start2 = std::chrono::system_clock::now();
        for(int i = 1; i< arr_num; i++){
            
            //greater_than<Lvl2>(cres, c_arr[i], k_1, k, ek, LOGIC);
            my_less_than<Lvl2>(cres, c_arr[i], k_1,k, ek, LOGIC, 29);
            //将TLWE转化为TGSW类型的密文
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);

            // 做cmux选择
            CMUXFFT<P>(r_cres, bootedTGSW, r_cres, r_c_arr[i]);

            // 密文提取
            SampleExtractIndex<P>(cres, r_cres, 0);
           
        }
        end2 = std::chrono::system_clock::now();
        end2end_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        

        d_cres = TFHEpp::tlweSymInt32Decrypt<P>(cres, pow(2., scale_bits_2), sk->key.get<P>());
        
        if(p_result != d_cres){
            error_times ++;
            std::cout << "p_arr: "  ;
            for(int i =0 ;i< arr_num; i++){
                std::cout << p_arr[i] << "  ";
            }
            std::cout << "  " <<  std::endl;

            std::cout << "p_result: " << p_result <<  std::endl;
            std::cout << "d_cres: " << d_cres <<  std::endl;
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "end2end_time: " << end2end_time/(1000 *num_test)  + lwe_to_rlwe_time / (1000000 * num_test)<< "s" <<  std::endl;
}

void do_swap_test( int num_test, int arr_num){
    // 两个数做交换
    arr_num = 2;

    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t k = 8;
    uint32_t plain_bits = k;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - k - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    double extract_time = 0;
    double swap_time = 0;
    double end2end_time = 0;
    std::chrono::system_clock::time_point start, end, start2, end2;
    
    Lvl1::T p_max, p_min;
    Lvl1::T p_arr[arr_num];
    std::vector<TLWELvl2> c_arr(arr_num);
    std::vector<TFHEpp::TRLWE<P>> r_c_arr(arr_num);

    TLWE<P> cres , c_max;
    TFHEpp::TRLWE<P> r_cres, r_c_max, r_sub0, r_sub1;
    P::T d_cres1 , d_cres2;

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;

    
    for(int time = 0; time < num_test; time++){
        p_max = 0;
        p_min = (1 << (plain_bits));
        for (size_t i = 0; i < arr_num; ++i) {
            p_arr[i] = message_2(engine);
            if(p_arr[i] > p_max){
                p_max = p_arr[i];
            }
            if(p_arr[i] < p_min){
                p_min = p_arr[i];
            }
        }

        
        for (size_t i = 0; i < arr_num; ++i) {
            c_arr[i] = TFHEpp::tlweSymInt32Encrypt<P>(p_arr[i], P::α, pow(2., scale_bits_2), sk->key.get<P>());

            // LWE密文转RLWE
            start = std::chrono::system_clock::now();
            TFHEpp::TLWEToTRLWE<P>(c_arr[i] , r_c_arr[i]);
            end = std::chrono::system_clock::now();
            lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }
        
        // 找两个数中的max
        cres = c_arr[0];
        r_cres = r_c_arr[0];
        start2 = std::chrono::system_clock::now();
        for(int i = 1; i< arr_num; i++){
            
            greater_than<Lvl2>(cres, c_arr[i], k_1, k, ek, LOGIC);
            TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
            CMUXFFT<P>(r_cres, bootedTGSW, r_cres, r_c_arr[i]);
            SampleExtractIndex<P>(cres, r_cres, 0);
        }
        
        c_max = cres;
        r_c_max = r_cres;

        start = std::chrono::system_clock::now();
        for (int j = 0; j <= P :: k; j++)
        {
            for (int i = 0; i < P::n; i++){
                r_sub0[j][i] = r_c_max[j][i] - r_c_arr[0][j][i] ;
                r_sub1[j][i] = r_c_max[j][i] - r_c_arr[1][j][i] ;
            } 
        }

        trgswfftExternalProduct<P>(r_sub0, r_sub0, bootedTGSW);
        trgswfftExternalProduct<P>(r_sub1, r_sub1, bootedTGSW);
        
        for (int j = 0; j <= P :: k; j++)
        {
            for (int i = 0; i < P::n; i++){
                r_c_arr[0][j][i] = r_c_arr[0][j][i] - r_sub0[j][i] - r_sub1[j][i] ;
                //r_c_arr[1][j][i] = r_c_arr[1][j][i] - r_sub0[j][i] - r_sub1[j][i] ;
            } 
        }

        r_c_arr[1] = r_c_max;

        SampleExtractIndex<P>(c_arr[0], r_c_arr[0], 0);
        SampleExtractIndex<P>(c_arr[1], r_c_arr[1], 0);

        end = std::chrono::system_clock::now();
        swap_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        end2 = std::chrono::system_clock::now();
        end2end_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        
        d_cres1 = TFHEpp::tlweSymInt32Decrypt<P>(c_arr[0], pow(2., scale_bits_2), sk->key.get<P>());
        //std::cout << "d_cres[0]: " << d_cres1 <<  std::endl;
        d_cres2 = TFHEpp::tlweSymInt32Decrypt<P>(c_arr[1], pow(2., scale_bits_2), sk->key.get<P>());
        //std::cout << "d_cres[1]: " << d_cres2 <<  std::endl;

        std::cout << "max: " << p_max  <<  std::endl;
        std::cout << "min: " << p_min  <<  std::endl;
        std::cout << "p_arr: "  ;
        for(int i =0 ;i< arr_num; i++){
            std::cout << p_arr[i] << "  ";
        }
        std::cout << "  " <<  std::endl;

        std::cout << "swap: " << d_cres1 << "  " << d_cres2 <<  std::endl;
        
        if(p_max != d_cres2 || p_min != d_cres1){
            error_times ++;
            // std::cout << "p_arr: "  ;
            // for(int i =0 ;i< arr_num; i++){
            //     std::cout << p_arr[i] << "  ";
            // }
            // std::cout << "  " <<  std::endl;

            // std::cout << "swap: " << d_cres1 << "  " << d_cres2 <<  std::endl;
            
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "lwe_to_gsw_time(per): " << lwe_to_gsw_time / (num_test * (arr_num-1)) << "ms" <<  std::endl;
    std::cout << "lwe_to_rlwe_time(per): " << lwe_to_rlwe_time / (num_test * arr_num) << "us" <<  std::endl;
    std::cout << "extract_time(per): " << extract_time /(num_test * (arr_num-1)) << "us" <<  std::endl;
    std::cout << "cmux_time(per): " << cmux_time / (num_test * (arr_num-1)) << "us" <<  std::endl;
    std::cout << "swap_time: " << swap_time/num_test << "ms" <<  std::endl;
    
}


void my_order_test( int num_test, int arr_num){
    std::cout << "-------my_order_test-------- " <<  std::endl;
    std::cout << "  num_test: " << num_test<<  std::endl;
    std::cout << " data_szie: " << arr_num<<  std::endl;

    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t plain_bits = 8;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - plain_bits - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    double swap_time = 0;
    double do_order_time = 0;
    double end2end_time = 0;
    std::chrono::system_clock::time_point start, end, start2, end2;
    
    Lvl1::T p_arr[arr_num];
    std::vector<TLWELvl2> c_arr(arr_num);
    std::vector<TFHEpp::TRLWE<P>> r_c_arr(arr_num);

    TLWE<P> cres , c_max;
    TFHEpp::TRLWE<P> r_cres, r_c_max, r_sub0, r_sub1;
    //P::T d_cres1 , d_cres2;
    Lvl1::T d_c_arr[arr_num];

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;

    
    for(int time = 0; time < num_test; time++){

        // 生成数据
        for (size_t i = 0; i < arr_num; ++i) {
            p_arr[i] = message_2(engine);
        }

        
        // 加密
        for (size_t i = 0; i < arr_num; ++i) {
            c_arr[i] = TFHEpp::tlweSymInt32Encrypt<P>(p_arr[i], P::α, pow(2., scale_bits_2), sk->key.get<P>());

            // LWE密文转RLWE
            start = std::chrono::system_clock::now();
            //TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl22param>(r_c_arr[i], c_arr[i], *iksk);
            TFHEpp::TLWEToTRLWE<P>(c_arr[i] , r_c_arr[i]);
            end = std::chrono::system_clock::now();
            lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }

        //明文原始排序
        std::cout << "p_arr: "<< std::endl;
        for(int i = 0 ; i < arr_num ;i++){
            std::cout << p_arr[i] << "  ";
        }
        std::cout << " "<< std::endl;

        //明文排序
        std::sort(p_arr, p_arr + arr_num);
        
        //密文排序
        start2 = std::chrono::system_clock::now();
        for(int j = arr_num; j> 1; j--){
            
            for(int i = 1; i< j; i++){
                
                //greater_than<Lvl2>(c_arr[i-1], c_arr[i], k_1, plain_bits, ek, LOGIC);
                my_greater_than<Lvl2>(c_arr[i-1], c_arr[i], k_1, plain_bits, ek, LOGIC, 29);
                // 将TLWE转化为TGSW类型的密文
                start = std::chrono::system_clock::now();
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
                end = std::chrono::system_clock::now();
                lwe_to_gsw_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

                // 做cmux选择大值
                start = std::chrono::system_clock::now();
                CMUXFFT<P>(r_c_max, bootedTGSW, r_c_arr[i-1], r_c_arr[i]);
                end = std::chrono::system_clock::now();
                cmux_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                //////做交换 - start////////
                start = std::chrono::system_clock::now();
                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        r_sub0[j][k] = r_c_max[j][k] - r_c_arr[i-1][j][k] ;
                        r_sub1[j][k] = r_c_max[j][k] - r_c_arr[i][j][k] ;
                    } 
                }

                trgswfftExternalProduct<P>(r_sub0, r_sub0, bootedTGSW);
                trgswfftExternalProduct<P>(r_sub1, r_sub1, bootedTGSW);
                
                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        r_c_arr[i-1][j][k] = r_c_arr[i-1][j][k] - r_sub0[j][k] - r_sub1[j][k] ;
                       
                    } 
                }

                r_c_arr[i] = r_c_max;

                SampleExtractIndex<P>(c_arr[i-1], r_c_arr[i-1], 0);
                SampleExtractIndex<P>(c_arr[i], r_c_arr[i], 0);
                end = std::chrono::system_clock::now();
                swap_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
                /////做交换 - end/////////
            }
        }
        end2 = std::chrono::system_clock::now();
        do_order_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        
        
        for(int i = 0 ; i < arr_num ;i++){
            d_c_arr[i] = TFHEpp::tlweSymInt32Decrypt<P>(c_arr[i], pow(2., scale_bits_2), sk->key.get<P>());
        }

        std::cout << "p_arr: "<< std::endl;
        for(int i = 0 ; i < arr_num ;i++){
            std::cout << p_arr[i] << "  ";
        }
        std::cout << " "<< std::endl;

        std::cout << "d_c_arr: "<< std::endl;
        for(int i = 0 ; i < arr_num ;i++){
            std::cout << d_c_arr[i] << "  ";
        }
        std::cout << " "<< std::endl;
        
        for(int i = 0 ; i < arr_num ;i++){
            if(d_c_arr[i] != p_arr[i]){
                error_times ++;
                break;
            }
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "do_order_time: " << do_order_time/(1000 * num_test) << "s" <<  std::endl;
    std::cout << "do_order + lwe_to_rlwe time: " << (do_order_time /1000 + lwe_to_rlwe_time/1000000) / num_test << "s" <<  std::endl;
}

void my_order_test_v2( int num_test, int arr_num){
    std::cout << "-------my_order_test_v2-------- " <<  std::endl;
    std::cout << "  num_test: " << num_test<<  std::endl;
    std::cout << " data_szie: " << arr_num<<  std::endl;

    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t plain_bits = 8;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - plain_bits - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    double swap_time = 0;
    double do_order_time = 0;
    double end2end_time = 0;
    std::chrono::system_clock::time_point start, end, start2, end2;
    
    Lvl1::T p_arr[arr_num];
    std::vector<TLWELvl2> c_arr(arr_num);
    std::vector<TFHEpp::TRLWE<P>> r_c_arr(arr_num);

    TLWE<P> cres , c_max;
    TFHEpp::TRLWE<P> r_cres, r_c_max, r_sub0, r_sub1;
    //P::T d_cres1 , d_cres2;
    Lvl1::T d_c_arr[arr_num];

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;

    for(int time = 0; time < num_test; time++){

        // 生成数据
        for (size_t i = 0; i < arr_num; ++i) {
            p_arr[i] = message_2(engine);
        }

        // 加密
        for (size_t i = 0; i < arr_num; ++i) {
            c_arr[i] = TFHEpp::tlweSymInt32Encrypt<P>(p_arr[i], P::α, pow(2., scale_bits_2), sk->key.get<P>());

            start = std::chrono::system_clock::now();
            TFHEpp::TLWEToTRLWE<P>(c_arr[i] , r_c_arr[i]);
            end = std::chrono::system_clock::now();
            lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        }

        //明文未排序
        // std::cout << "p_arr: "<< std::endl;
        // for(int i = 0 ; i < arr_num ;i++){
        //     std::cout << p_arr[i] << "  ";
        // }
        // std::cout << " "<< std::endl;

        //明文排序
        std::sort(p_arr, p_arr + arr_num);
        
        //密文排序（普通冒泡）
        // start2 = std::chrono::system_clock::now();
        // for(int j = arr_num; j> 1; j--){
            
        //     for(int i = 1; i< j; i++){
                
        //         //greater_than<Lvl2>(c_arr[i-1], c_arr[i], k_1, plain_bits, ek, LOGIC);
        //         my_greater_than<Lvl2>(c_arr[i-1], c_arr[i], k_1, plain_bits, ek, LOGIC, 29);
        //         // 将TLWE转化为TGSW类型的密文
        //         start = std::chrono::system_clock::now();
        //         TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
        //         end = std::chrono::system_clock::now();
        //         lwe_to_gsw_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //         // 做cmux选择大值
        //         start = std::chrono::system_clock::now();
        //         CMUXFFT<P>(r_c_max, bootedTGSW, r_c_arr[i-1], r_c_arr[i]);
        //         end = std::chrono::system_clock::now();
        //         cmux_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        //         //////做交换 - start////////
        //         start = std::chrono::system_clock::now();
        //         for (int j = 0; j <= P :: k; j++)
        //         {
        //             for (int k = 0; k < P::n; k++){
        //                 //r_sub0[j][k] = r_c_max[j][k] - r_c_arr[i-1][j][k] ;
        //                 r_sub1[j][k] = r_c_max[j][k] - r_c_arr[i][j][k] ;
        //             } 
        //         }
        //         // trgswfftExternalProduct<P>(r_sub0, r_sub0, bootedTGSW);
        //         // trgswfftExternalProduct<P>(r_sub1, r_sub1, bootedTGSW);
                
        //         for (int j = 0; j <= P :: k; j++)
        //         {
        //             for (int k = 0; k < P::n; k++){
        //                 r_c_arr[i-1][j][k] = r_c_arr[i-1][j][k] - r_sub1[j][k] ;
        //             } 
        //         }

        //         r_c_arr[i] = r_c_max;

        //         SampleExtractIndex<P>(c_arr[i-1], r_c_arr[i-1], 0);
        //         SampleExtractIndex<P>(c_arr[i], r_c_arr[i], 0);
        //         end = std::chrono::system_clock::now();
        //         swap_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        //         /////做交换 - end/////////
        //     }
        // }
        // end2 = std::chrono::system_clock::now();
        // do_order_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();

        //密文排序（奇偶排序）
        start2 = std::chrono::system_clock::now();
        for(int num = 0 ; num < arr_num; num ++){
            
            for(int i = (num%2); i+1< arr_num; i= i+2){
                
                //greater_than<Lvl2>(c_arr[i-1], c_arr[i], k_1, plain_bits, ek, LOGIC);
                my_greater_than<Lvl2>(c_arr[i], c_arr[i+1], k_1, plain_bits, ek, LOGIC, 29);
                // 将TLWE转化为TGSW类型的密文
                start = std::chrono::system_clock::now();
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
                end = std::chrono::system_clock::now();
                lwe_to_gsw_time += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

                // 做cmux选择大值
                start = std::chrono::system_clock::now();
                CMUXFFT<P>(r_c_max, bootedTGSW, r_c_arr[i], r_c_arr[i+1]);
                end = std::chrono::system_clock::now();
                cmux_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

                //////做交换 - start////////
                start = std::chrono::system_clock::now();
                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        //r_sub0[j][k] = r_c_max[j][k] - r_c_arr[i-1][j][k] ;
                        r_sub1[j][k] = r_c_max[j][k] - r_c_arr[i+1][j][k] ;
                    } 
                }
                // trgswfftExternalProduct<P>(r_sub0, r_sub0, bootedTGSW);
                // trgswfftExternalProduct<P>(r_sub1, r_sub1, bootedTGSW);
                
                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        r_c_arr[i][j][k] = r_c_arr[i][j][k] - r_sub1[j][k] ;
                    } 
                }

                r_c_arr[i+1] = r_c_max;

                SampleExtractIndex<P>(c_arr[i], r_c_arr[i], 0);
                SampleExtractIndex<P>(c_arr[i+1], r_c_arr[i+1], 0);
                end = std::chrono::system_clock::now();
                swap_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
                /////做交换 - end/////////
            }
        }
        end2 = std::chrono::system_clock::now();
        do_order_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        
        
        for(int i = 0 ; i < arr_num ;i++){
            d_c_arr[i] = TFHEpp::tlweSymInt32Decrypt<P>(c_arr[i], pow(2., scale_bits_2), sk->key.get<P>());
        }

        // std::cout << "p_arr: "<< std::endl;
        // for(int i = 0 ; i < arr_num ;i++){
        //     std::cout << p_arr[i] << "  ";
        // }
        // std::cout << " "<< std::endl;

        // std::cout << "d_c_arr: "<< std::endl;
        // for(int i = 0 ; i < arr_num ;i++){
        //     std::cout << d_c_arr[i] << "  ";
        // }
        // std::cout << " "<< std::endl;
        
        for(int i = 0 ; i < arr_num ;i++){
            if(d_c_arr[i] != p_arr[i]){
                error_times ++;

                std::cout << "p_arr: "<< std::endl;
                for(int i = 0 ; i < arr_num ;i++){
                    std::cout << p_arr[i] << "  ";
                }
                std::cout << " "<< std::endl;

                std::cout << "d_c_arr: "<< std::endl;
                for(int i = 0 ; i < arr_num ;i++){
                    std::cout << d_c_arr[i] << "  ";
                }
                std::cout << " "<< std::endl;

                break;
            }
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "do_order_time: " << do_order_time/(1000 * num_test) << "s" <<  std::endl;
    std::cout << "do_order + lwe_to_rlwe time: " << (do_order_time /1000 + lwe_to_rlwe_time/1000000) / num_test << "s" <<  std::endl;
}

void two_order_compare_test(int num_test, int arr_num){
    std::cout << "-------two_order_compare_test-------- " <<  std::endl;
    std::cout << "  num_test: " << num_test<<  std::endl;
    std::cout << " data_szie: " << arr_num<<  std::endl;

    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;

    uint32_t plain_bits = 8;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - plain_bits - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    double swap_time = 0;
    double he3db_order_time = 0;
    double my_order_time = 0;
    double end2end_time = 0;
    std::chrono::system_clock::time_point start, end, start2, end2;
    
    Lvl1::T p_arr[arr_num];
    std::vector<TLWELvl2> c_arr(arr_num), c_arr2(arr_num);
    std::vector<TFHEpp::TRLWE<P>> r_c_arr(arr_num), r_c_arr2(arr_num);

    TLWE<P> cres , c_max , c_min;
    TFHEpp::TRLWE<P> r_cres, r_c_max, r_c_min,  r_sub0, r_sub1;
    //P::T d_cres1 , d_cres2;
    Lvl1::T d_c_arr[arr_num], d_c_arr2[arr_num];

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;

    
    for(int time = 0; time < num_test; time++){

        // 生成数据
        for (size_t i = 0; i < arr_num; ++i) {
            p_arr[i] = message_2(engine);
        }

        // 加密
        for (size_t i = 0; i < arr_num; ++i) {
            c_arr[i] = TFHEpp::tlweSymInt32Encrypt<P>(p_arr[i], P::α, pow(2., scale_bits_2), sk->key.get<P>());
            // LWE密文转RLWE
            TFHEpp::TLWEToTRLWE<P>(c_arr[i] , r_c_arr[i]);
            end = std::chrono::system_clock::now();
            lwe_to_rlwe_time += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            c_arr2[i] = c_arr[i];
            r_c_arr2[i] = r_c_arr[i];
        }

        //明文排序
        std::sort(p_arr, p_arr + arr_num);

        //HE3DB密文排序(对c_arr2，r_c_arr2，d_c_arr2操作)
        //std::cout << "-----HE3DB sort--------: "<< std::endl;
        start2 = std::chrono::system_clock::now();
        for(int j = arr_num; j> 1; j--){
            for(int i = 1; i< j; i++){
                greater_than<Lvl2>(c_arr2[i-1], c_arr2[i], k_1, plain_bits, ek, LOGIC);
                
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
               
                CMUXFFT<P>(r_c_max, bootedTGSW, r_c_arr2[i-1], r_c_arr2[i]);
                

                less_than<Lvl2>(c_arr2[i-1], c_arr2[i], k_1, plain_bits, ek, LOGIC);
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
                
                CMUXFFT<P>(r_c_min, bootedTGSW, r_c_arr2[i-1], r_c_arr2[i]);

                r_c_arr2[i] = r_c_max;
                r_c_arr2[i-1] = r_c_min;
                
                SampleExtractIndex<P>(c_arr2[i-1], r_c_arr2[i-1], 0);
                SampleExtractIndex<P>(c_arr2[i], r_c_arr2[i], 0);
            }
        }
        end2 = std::chrono::system_clock::now();
        he3db_order_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        
        
        for(int i = 0 ; i < arr_num ;i++){
            d_c_arr2[i] = TFHEpp::tlweSymInt32Decrypt<P>(c_arr2[i], pow(2., scale_bits_2), sk->key.get<P>());
        }

        // std::cout << "p_arr: "<< std::endl;
        // for(int i = 0 ; i < arr_num ;i++){
        //     std::cout << p_arr[i] << "  ";
        // }
        // std::cout << " "<< std::endl;

        // std::cout << "d_c_arr2: "<< std::endl;
        // for(int i = 0 ; i < arr_num ;i++){
        //     std::cout << d_c_arr2[i] << "  ";
        // }
        // std::cout << " "<< std::endl;
        
        
        //我的密文排序(对c_arr，r_c_arr，d_c_arr操作)
        //std::cout << "-----my sort--------: "<< std::endl;
        start2 = std::chrono::system_clock::now();
        for(int j = arr_num; j> 1; j--){
            
            for(int i = 1; i< j; i++){
                
                my_greater_than<Lvl2>(c_arr[i-1], c_arr[i], k_1, plain_bits, ek, LOGIC, 29);
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);

                CMUXFFT<P>(r_c_max, bootedTGSW, r_c_arr[i-1], r_c_arr[i]);

                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        r_sub0[j][k] = r_c_max[j][k] - r_c_arr[i-1][j][k] ;
                        r_sub1[j][k] = r_c_max[j][k] - r_c_arr[i][j][k] ;
                    } 
                }

                trgswfftExternalProduct<P>(r_sub0, r_sub0, bootedTGSW);
                trgswfftExternalProduct<P>(r_sub1, r_sub1, bootedTGSW);
                
                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        r_c_arr[i-1][j][k] = r_c_arr[i-1][j][k] - r_sub0[j][k] - r_sub1[j][k] ;
                       
                    } 
                }

                r_c_arr[i] = r_c_max;

                SampleExtractIndex<P>(c_arr[i-1], r_c_arr[i-1], 0);
                SampleExtractIndex<P>(c_arr[i], r_c_arr[i], 0);
            }
        }
        end2 = std::chrono::system_clock::now();
        my_order_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        
        
        for(int i = 0 ; i < arr_num ;i++){
            d_c_arr[i] = TFHEpp::tlweSymInt32Decrypt<P>(c_arr[i], pow(2., scale_bits_2), sk->key.get<P>());
        }

        // std::cout << "p_arr: "<< std::endl;
        // for(int i = 0 ; i < arr_num ;i++){
        //     std::cout << p_arr[i] << "  ";
        // }
        // std::cout << " "<< std::endl;

        // std::cout << "d_c_arr: "<< std::endl;
        // for(int i = 0 ; i < arr_num ;i++){
        //     std::cout << d_c_arr[i] << "  ";
        // }
        // std::cout << " "<< std::endl;
        
        for(int i = 0 ; i < arr_num ;i++){
            if(d_c_arr[i] != p_arr[i]){
                error_times ++;

                std::cout << "-----my sort--------: "<< std::endl;
                std::cout << "p_arr: "<< std::endl;
                for(int i = 0 ; i < arr_num ;i++){
                    std::cout << p_arr[i] << "  ";
                }
                std::cout << " "<< std::endl;

                std::cout << "d_c_arr: "<< std::endl;
                for(int i = 0 ; i < arr_num ;i++){
                    std::cout << d_c_arr[i] << "  ";
                }
                std::cout << " "<< std::endl;

                break;
            }
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "he3db_order + lwe_to_rlwe time: " << (he3db_order_time /1000 + lwe_to_rlwe_time/1000000) / num_test << "s" <<  std::endl;
    std::cout << "my_order + lwe_to_rlwe time: " << (my_order_time /1000 + lwe_to_rlwe_time/1000000) / num_test << "s" <<  std::endl;
}


void two_order_compare_test_v2( int num_test, int arr_num){
    std::cout << "-------two_order_compare_test_v2-------- " <<  std::endl;
    std::cout << "  num_test: " << num_test<<  std::endl;
    std::cout << " data_szie: " << arr_num<<  std::endl;

    using P = Lvl2;
    using iksP = TFHEpp::lvl10param;
    using bkP = TFHEpp::lvl02param;
    using privksP = TFHEpp::lvl22param;
    
    uint32_t plain_bits = 8;
    uint32_t scale_bits_1 = 29 ;
    uint32_t scale_bits_2 = std::numeric_limits<privksP::targetP::T>::digits - plain_bits - 1;

    // TFHESecretKey sk;
    // TFHEEvalKey ek;
    // ek.emplacebkfft<Lvl01>(sk);
    // ek.emplacebkfft<Lvl02>(sk);
    // ek.emplaceiksk<Lvl20>(sk);
    // ek.emplaceiksk<Lvl10>(sk);
    // ek.emplaceiksk<Lvl21>(sk);

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey;
    TFHEpp::EvalKey ek;

    ek.emplacebkfft<Lvl01>(*sk);
    ek.emplacebkfft<Lvl02>(*sk);
    ek.emplaceiksk<Lvl20>(*sk);
    ek.emplaceiksk<Lvl10>(*sk);
    ek.emplaceiksk<Lvl21>(*sk);

    ek.emplaceiksk<iksP>(*sk);
    ek.emplacebkfft<bkP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);

    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl22param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl22param>(*iksk, *sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message_2(0, (1 << (plain_bits) - 1));
    
    int error_times = 0;
    double lwe_to_gsw_time = 0;
    double lwe_to_rlwe_time = 0;
    double cmux_time = 0;
    double swap_time = 0;
    double he3db_order_time = 0;
    double my_order_time = 0;
    double end2end_time = 0;
    std::chrono::system_clock::time_point start, end, start2, end2;
    
    Lvl1::T p_arr[arr_num];
    std::vector<TLWELvl2> c_arr(arr_num), c_arr2(arr_num);
    std::vector<TFHEpp::TRLWE<P>> r_c_arr(arr_num), r_c_arr2(arr_num);

    TLWE<P> cres , c_max , c_min;
    TFHEpp::TRLWE<P> r_cres, r_c_max, r_c_min,  r_sub0, r_sub1;
    //P::T d_cres1 , d_cres2;
    Lvl1::T d_c_arr[arr_num], d_c_arr2[arr_num];

    TLWE<typename iksP::domainP> k_1;
    TFHEpp::TRGSWFFT<typename privksP::targetP> bootedTGSW;

    
    for(int time = 0; time < num_test; time++){

        // 生成数据
        for (size_t i = 0; i < arr_num; ++i) {
            p_arr[i] = message_2(engine);
        }

        // 加密
        for (size_t i = 0; i < arr_num; ++i) {
            c_arr[i] = TFHEpp::tlweSymInt32Encrypt<P>(p_arr[i], P::α, pow(2., scale_bits_2), sk->key.get<P>());
            TFHEpp::TLWEToTRLWE<P>(c_arr[i] , r_c_arr[i]);
    
            c_arr2[i] = c_arr[i];
            r_c_arr2[i] = r_c_arr[i];
        }

        //明文排序
        std::sort(p_arr, p_arr + arr_num);

        //HE3DB密文排序(对c_arr2，r_c_arr2，d_c_arr2操作)
        //std::cout << "-----HE3DB sort--------: "<< std::endl;
        start2 = std::chrono::system_clock::now();
        for(int num = 0 ; num < arr_num; num ++){
        
            for(int i = (num%2); i+1< arr_num; i= i+2){
                greater_than<Lvl2>(c_arr2[i], c_arr2[i+1], k_1, plain_bits, ek, LOGIC);
                
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
               
                CMUXFFT<P>(r_c_max, bootedTGSW, r_c_arr2[i], r_c_arr2[i+1]);

                less_than<Lvl2>(c_arr2[i], c_arr2[i+1], k_1, plain_bits, ek, LOGIC);
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);
                
                CMUXFFT<P>(r_c_min, bootedTGSW, r_c_arr2[i], r_c_arr2[i+1]);

                r_c_arr2[i+1] = r_c_max;
                r_c_arr2[i] = r_c_min;
                
                SampleExtractIndex<P>(c_arr2[i], r_c_arr2[i], 0);
                SampleExtractIndex<P>(c_arr2[i+1], r_c_arr2[i+1], 0);
            }
        }
        end2 = std::chrono::system_clock::now();
        he3db_order_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        
        
        for(int i = 0 ; i < arr_num ;i++){
            d_c_arr2[i] = TFHEpp::tlweSymInt32Decrypt<P>(c_arr2[i], pow(2., scale_bits_2), sk->key.get<P>());
        }

        std::cout << "p_arr: "<< std::endl;
        for(int i = 0 ; i < arr_num ;i++){
            std::cout << p_arr[i] << "  ";
        }
        std::cout << " "<< std::endl;

        std::cout << "d_c_arr2: "<< std::endl;
        for(int i = 0 ; i < arr_num ;i++){
            std::cout << d_c_arr2[i] << "  ";
        }
        std::cout << " "<< std::endl;
        
        
        //我的密文排序(对c_arr，r_c_arr，d_c_arr操作)
        //std::cout << "-----my sort--------: "<< std::endl;
        start2 = std::chrono::system_clock::now();
        for(int num = 0 ; num < arr_num; num ++){
        
            for(int i = (num%2); i+1< arr_num; i= i+2){
                
                //greater_than<Lvl2>(c_arr[i-1], c_arr[i], k_1, plain_bits, ek, LOGIC);
                my_greater_than<Lvl2>(c_arr[i], c_arr[i+1], k_1, plain_bits, ek, LOGIC, 29);
                
                TFHEpp::CircuitBootstrappingFFT<iksP, bkP, privksP>(bootedTGSW, k_1, ek);

                
                CMUXFFT<P>(r_c_max, bootedTGSW, r_c_arr[i], r_c_arr[i+1]);
               
                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        
                        r_sub1[j][k] = r_c_max[j][k] - r_c_arr[i+1][j][k] ;
                    } 
                }
                
                for (int j = 0; j <= P :: k; j++)
                {
                    for (int k = 0; k < P::n; k++){
                        r_c_arr[i][j][k] = r_c_arr[i][j][k] - r_sub1[j][k] ;
                    } 
                }

                r_c_arr[i+1] = r_c_max;

                SampleExtractIndex<P>(c_arr[i], r_c_arr[i], 0);
                SampleExtractIndex<P>(c_arr[i+1], r_c_arr[i+1], 0);
               
            }
        }
        end2 = std::chrono::system_clock::now();
        my_order_time += std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();
        
        
        for(int i = 0 ; i < arr_num ;i++){
            d_c_arr[i] = TFHEpp::tlweSymInt32Decrypt<P>(c_arr[i], pow(2., scale_bits_2), sk->key.get<P>());
        }

        std::cout << "p_arr: "<< std::endl;
        for(int i = 0 ; i < arr_num ;i++){
            std::cout << p_arr[i] << "  ";
        }
        std::cout << " "<< std::endl;

        std::cout << "d_c_arr: "<< std::endl;
        for(int i = 0 ; i < arr_num ;i++){
            std::cout << d_c_arr[i] << "  ";
        }
        std::cout << " "<< std::endl;
        

        for(int i = 0 ; i < arr_num ;i++){
            if(d_c_arr[i] != p_arr[i]){
                error_times ++;

                // std::cout << "p_arr: "<< std::endl;
                // for(int i = 0 ; i < arr_num ;i++){
                //     std::cout << p_arr[i] << "  ";
                // }
                // std::cout << " "<< std::endl;

                // std::cout << "d_c_arr: "<< std::endl;
                // for(int i = 0 ; i < arr_num ;i++){
                //     std::cout << d_c_arr[i] << "  ";
                // }
                // std::cout << " "<< std::endl;

                break;
            }
        }

    }
    
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "he3db_order + lwe_to_rlwe time: " << (he3db_order_time /1000 + lwe_to_rlwe_time/1000000) / num_test << "s" <<  std::endl;
    std::cout << "my_order + lwe_to_rlwe time: " << (my_order_time /1000 + lwe_to_rlwe_time/1000000) / num_test << "s" <<  std::endl;
}


int main()
{
    int num_test = 100;
    //do_swap_test(10, 2);
    // my_order_test_v2(100, 8);

    two_order_compare_test_v2(5, 64);

    //lwe_to_gsw_test3(100);
    
    
}