

#include <iostream>
#include <chrono>
#include <random>
#include<bitset>
#include "src/HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"

#include "src/HEDB/comparison/tfhepp_utils.h"
#include <gatebootstrapping.hpp>
#include "detwfa.hpp"

#include <cassert>
#include <tfhe++.hpp>


using namespace std;
using namespace HEDB;
using namespace TFHEpp;


// In comparison, if a and b are x bits, a - b is (x+1) bits
// If we define the input data is p bits, than the scale should be Q/2p
// The database is the same as comparison
// tlwelvl1_comparison_test include 1 - 9 bits, here the plain bits represent the range of input data.


void my_comparison_test(uint32_t plain_bits, int num_test){
    
    std::vector<uint32_t> error_time(4,0) ;
     
    for (int test_time = 0; test_time < num_test; test_time++) 
    {
        std::cout << "------ Test of msb_related tttttt------" << std::endl;
        std::cout << "Test Time : " << num_test << std::endl;
        std::cout << "Plain bits : " << plain_bits << std::endl;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        using P = Lvl1;
        TFHESecretKey sk;
        TFHEEvalKey ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        
        ek.emplaceiksk<Lvl11>(sk);
        
        uint32_t scale_bits;
        
        std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
        std::uniform_int_distribution<typename P::T> type(0, 1);
        scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
        
        std::cout << "scale_bits : " << scale_bits << std::endl;
        
        Lvl1::T a = pow(2., plain_bits) -1;
        Lvl1::T b = message(engine);
        Lvl1::T da, da_offset, dres, dres_offset;
       
        
        std::cout << "a : " << a << std::endl;
        
        
        TLWELvl1 ca = TFHEpp::tlweSymInt32Encrypt<P>(a, P::α, pow(2., scale_bits), sk.key.get<P>());
        
        //uint32_t ca_extra_bit = pair_paf.second;
        uint32_t ca_extra_bit = ca[P::n] & (1UL << 31);
        std::cout << "ca_extra_bit:" << ca_extra_bit << std::endl;
        
        TLWELvl1 res1;
        MSBGateBootstrapping(res1, ca, ek, LOGIC);
        
        std::cout << "-----------原本提取最高有效位------------ : " << std::endl;
        da = TFHEpp::tlweSymInt32Decrypt_print<P>(res1, pow(2., 29), sk.key.lvl1);
    
        TLWELvl1 ca_offset;
        for(int i = 0; i <= P:: n ; i++){
            ca_offset[i] = ca[i] << 1;
        }
        
        TLWELvl1 res2;
        MSBGateBootstrapping(res2, ca_offset, ek, LOGIC);

        std::cout << "-----------左移一位提取最高有效位------------ : " << std::endl;
        da  = TFHEpp::tlweSymInt32Decrypt_print<P>(res2, pow(2., 29), sk.key.lvl1);
        
        TLWELvl1 res3;
        for(int i = 0; i <= P:: n ; i++){
            res3[i] = res1[i] + res2[i];
        }

        std::cout << "-----------新的提取最高有效位------------ : " << std::endl;
        da  = TFHEpp::tlweSymInt32Decrypt_print<P>(res3, pow(2., 29), sk.key.lvl1);
    }
}

void comparison_modify_test(uint32_t plain_bits, int num_test){
    
    std::vector<uint32_t> error_time(4,0) ;
     
    for (int test_time = 0; test_time < num_test; test_time++) 
    {
        std::cout << "------ Test of msb_related tttttt------" << std::endl;
        std::cout << "Test Time : " << num_test << std::endl;
        std::cout << "Plain bits : " << plain_bits << std::endl;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        using P = Lvl1;
        TFHESecretKey sk;
        TFHEEvalKey ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        
        ek.emplaceiksk<Lvl11>(sk);
        
        uint32_t scale_bits;
        
        std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
        std::uniform_int_distribution<typename P::T> type(0, 1);
        scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
        
        std::cout << "scale_bits : " << scale_bits << std::endl;
        
        Lvl1::T a = 15;
        Lvl1::T b = message(engine);
        Lvl1::T da, da_offset, dres, dres_offset;
       
        
        std::cout << "a : " << a << std::endl;
        
        
        pair <TLWELvl1, uint32_t> pair_paf = TFHEpp::new_tlweSymInt32Encrypt<P>(a, P::α, pow(2., scale_bits), sk.key.get<P>());
        TLWELvl1 ca = pair_paf.first;
        TLWELvl1 ca_offset = ca;
        
        //uint32_t ca_extra_bit = pair_paf.second;
        uint32_t ca_extra_bit = ca[P::n] & (1UL << 31);
        std::cout << "ca_extra_bit:" << ca_extra_bit << std::endl;
        
        TLWELvl1 res;
        
        // std::cout << "-----------ca Decrypt------------ : " << std::endl;
        // da = TFHEpp::tlweSymInt32Decrypt_print<P>(ca, pow(2., scale_bits ), sk.key.lvl1);
    
        
        for(int i = 0; i <= P:: n ; i++){
            ca_offset[i] = ca[i] << 1;
        }
        
        std::cout << "-----------ca left offset Decrypt------------ : " << std::endl;
        pair <P::T, uint64_t> pair_paf2  = TFHEpp::new_tlweSymInt32Decrypt<P>(ca_offset, pow(2., scale_bits ), sk.key.lvl1);
        
        TFHEpp::TLWE<TFHEpp::lvl1param> ks_ca;
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl11param>(ks_ca, ca_offset, *ek.iksklvl11);
        
        std::cout << "ca_offset[n]: " << ca_offset[P::n]  << std::endl;
        std::cout << "ks_ca[n]    : " << ks_ca[P::n] << std::endl;
        std::cout << "ca_offset[n]: " << bitset<32>(ca_offset[P::n]) << std::endl;
        std::cout << "ks_ca[n]    : " << bitset<32>(ks_ca[P::n]) << std::endl;
        
        std::cout << "-----------ks_ca Decrypt------------ : " << std::endl;
        pair <P::T, uint64_t> pair_paf3  = TFHEpp::new_tlweSymInt32Decrypt<P>(ks_ca, pow(2., scale_bits ), sk.key.lvl1);
        Lvl1::T ks_da = pair_paf3.first;
        
        std::cout << "-----------ca left switch Decrypt again------------ : " << std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<P>(ca_offset, pow(2., scale_bits ), sk.key.lvl1);
        
        for(int i = 0; i <= P:: n ; i++){
            ca_offset[i] = ks_ca[i] >> 1;
        }
        
        std::cout << "-----------ca right offset Decrypt------------ : " << std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<P>(ca_offset, pow(2., scale_bits ), sk.key.lvl1);
        
        if(ca_extra_bit > 0) ca_offset[P::n] += (1UL << 31);
        
        std::cout << "-----------ca finish offset Decrypt------------ : " << std::endl;
        da_offset = TFHEpp::tlweSymInt32Decrypt_print<P>(ca_offset, pow(2., scale_bits ), sk.key.lvl1);
        
        std::cout << "dada/da_offset:   " << std::round((double)da/da_offset) << std::endl;

        // my_MSBGateBootstrapping(res, ca, ek, ARITHMETIC, 0);
        // std::cout << "-----------MSB(ca) Decrypt------------ : " << std::endl;
        // dres = TFHEpp::tlweSymInt32Decrypt_print<P>(res, pow(2., 31), sk.key.lvl1);
        
        // my_MSBGateBootstrapping(res, ca_offset, ek, ARITHMETIC, 0);
        // std::cout << "----------offset MSB(ca) Decrypt------------ : " << std::endl;
        // dres_offset = TFHEpp::tlweSymInt32Decrypt_print<P>(res, pow(2., 31), sk.key.lvl1);
        
        // std::cout << "MSB : " << dres << std::endl;
        // std::cout << "offset MSB  " << dres_offset << std::endl;
        
        
        
        
        ////////////////my_MSBGateBootstrapping/////////////////////
        // Lvl1::T μ = Lvl1::μ;
        // if (IS_ARITHMETIC(ARITHMETIC)) μ = μ << 1;
        // else μ = μ >> 0;
        // constexpr uint64_t offset = 1ULL << (std::numeric_limits<Lvl1::T>::digits - 6);
        // TFHEpp::TLWE<Lvl1> tlweoffset = ca;
        // tlweoffset[Lvl1::k * Lvl1::n] += offset;
        
        // std::cout << "-----------tlweoffset----Decrpt---------" <<std::endl;  
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(tlweoffset, pow(2., 31), sk.key.lvl1);
        
        // TFHEpp::TLWE<Lvl0> tlwelvl0;
        // TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        
        // std::cout << "-----------KeySwitch-----Decrpt---------" <<std::endl;  
        // TFHEpp::tlweSymDecrypt_print<Lvl0>(tlwelvl0, sk.key.lvl0);
        
        // TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<Lvl1>(μ));
        
        // std::cout << "-----------GateBootstrapping----Decrpt---------" <<std::endl;  
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.lvl1);
        
        //if (IS_ARITHMETIC(ARITHMETIC)) res[Lvl1::k * Lvl1::n] += (μ);
        /////////////////////////////////////////////////////////////////////////////////////

        if( da != a ) error_time[0] ++ ;
        
        if(dres == 0) error_time[1] ++ ;
        
        if(da_offset != 15) error_time[2] ++ ;
        
        if(std::round((double)da_offset/da) != 1) error_time[3] ++ ;
    }
    for(int i =0 ; i< 4 ;i++){
        std::cout << "error_time"<<"[" << i <<"]"<<": " << error_time[i] << std::endl;
    }
    
    
}


void cmux_test(uint32_t plain_bits, int num_tests){
    using P = Lvl1;
    
    constexpr uint32_t num_test = 100;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    uint32_t scale_bits;
    
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    Lvl1::T a = message(engine);
    Lvl1::T b = message(engine);
    Lvl1::T da, db, d_subs, d_res;
   
    
    std::cout << "a : " << a << std::endl;
    // std::cout << "a[0] : " << a[0] << std::endl;
    // std::cout << "a[1] : " << a[1] << std::endl;
    // std::cout << "a[2] : " << a[2] << std::endl;
    // std::cout << "a[3] : " << a[3] << std::endl;
    std::cout << "b : " << b << std::endl;
    
    TLWELvl1 ca = TFHEpp::tlweSymInt32Encrypt<P>(a, P::α, pow(2., scale_bits), sk.key.get<P>());
    TLWELvl1 cb = TFHEpp::tlweSymInt32Encrypt<P>(b, P::α, pow(2., scale_bits), sk.key.get<P>());
    
    TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl11param>* iksk = new TFHEpp::TLWE2TRLWEIKSKey<TFHEpp::lvl11param>();
    TFHEpp::tlwe2trlweikskkgen<TFHEpp::lvl11param>(*iksk, sk);
    
    TFHEpp::TRLWE<TFHEpp::lvl1param> r_ca, r_cb, cres;
    TRGSWFFT<lvl1param> cs;
    int32_t ps = 0;
    array<bool, lvl1param::n> pres;
    
    TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl11param>(r_ca, ca, *iksk);
    TFHEpp::TLWE2TRLWEIKS<TFHEpp::lvl11param>(r_cb, cb, *iksk);
    
    Polynomial<TFHEpp::lvl1param> plainpoly = {};
        plainpoly[0] = ps;
        cs = trgswfftSymEncrypt<lvl1param>(plainpoly, lvl1param::α,
                                              sk.key.lvl1);
    
    chrono::system_clock::time_point start, end;
    
    start = chrono::system_clock::now();
    CMUXFFT<lvl1param>(cres, cs, r_ca, r_cb);
    end = chrono::system_clock::now();
    
    pres = trlweSymDecrypt<lvl1param>(cres, sk.key.lvl1);
    
    for (int i = 0; i < lvl1param::n; i++)
            assert(pres[i] ==
                   (((ps[test] > 0) ? a[i] : b[i]) > 0));
    
    
    // random_device seed_gen;
    // default_random_engine engine(seed_gen());
    // uniform_int_distribution<uint32_t> binary(0, 1);

    // SecretKey *sk = new SecretKey;
    // vector<int32_t> ps(num_test);
    // vector<array<uint8_t, lvl1param::n>> p1(num_test);
    // vector<array<uint8_t, lvl1param::n>> p0(num_test);

    // vector<array<uint32_t, lvl1param::n>> pmu1(num_test);
    // vector<array<uint32_t, lvl1param::n>> pmu0(num_test);
    // array<bool, lvl1param::n> pres;

    // for (int32_t &p : ps) p = binary(engine);
    // for (array<uint8_t, lvl1param::n> &i : p1)
    //     for (uint8_t &p : i) p = binary(engine);
    // for (array<uint8_t, lvl1param::n> &i : p0)
    //     for (uint8_t &p : i) p = binary(engine);

    // for (int i = 0; i < num_test; i++)
    //     for (int j = 0; j < lvl1param::n; j++)
    //         pmu1[i][j] = (p1[i][j] > 0) ? lvl1param::μ : -lvl1param::μ;
    // for (int i = 0; i < num_test; i++)
    //     for (int j = 0; j < lvl1param::n; j++)
    //         pmu0[i][j] = (p0[i][j] > 0) ? lvl1param::μ : -lvl1param::μ;
    // vector<TRGSWFFT<lvl1param>> cs(num_test);
    // vector<TRLWE<lvl1param>> c1(num_test);
    // vector<TRLWE<lvl1param>> c0(num_test);
    // vector<TRLWE<lvl1param>> cres(num_test);

    // for (int i = 0; i < num_test; i++) {
    //     Polynomial<TFHEpp::lvl1param> plainpoly = {};
    //     plainpoly[0] = ps[i];
    //     cs[i] = trgswfftSymEncrypt<lvl1param>(plainpoly, lvl1param::α,
    //                                           sk->key.lvl1);
    // }
    // for (int i = 0; i < num_test; i++)
    //     c1[i] = trlweSymEncrypt<lvl1param>(pmu1[i], lvl1param::α, sk->key.lvl1);
    // for (int i = 0; i < num_test; i++)
    //     c0[i] = trlweSymEncrypt<lvl1param>(pmu0[i], lvl1param::α, sk->key.lvl1);

    // chrono::system_clock::time_point start, end;
    // start = chrono::system_clock::now();
    // for (int test = 0; test < num_test; test++) {
    //     CMUXFFT<lvl1param>(cres[test], cs[test], c1[test], c0[test]);
    // }
    // end = chrono::system_clock::now();

    // for (int test = 0; test < num_test; test++) {
    //     pres = trlweSymDecrypt<lvl1param>(cres[test], sk->key.lvl1);
    //     for (int i = 0; i < lvl1param::n; i++)
    //         assert(pres[i] ==
    //               (((ps[test] > 0) ? p1[test][i] : p0[test][i]) > 0));
    // }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();
    cout << elapsed / num_test << "μs" << endl;
}

// void tfhepp_cmux_test(uint32_t plain_bits, int num_tests){
//     using P = Lvl1;
     
//     constexpr uint32_t num_test = 100;
//     random_device seed_gen;
//     default_random_engine engine(seed_gen());
//     uniform_int_distribution<uint32_t> binary(0, 1);

//     SecretKey *sk = new SecretKey;
//     vector<int32_t> ps(num_test);
//     vector<array<uint8_t, lvl1param::n>> p1(num_test);
//     vector<array<uint8_t, lvl1param::n>> p0(num_test);

//     vector<array<uint32_t, lvl1param::n>> pmu1(num_test);
//     vector<array<uint32_t, lvl1param::n>> pmu0(num_test);
//     array<bool, lvl1param::n> pres;

//     for (int32_t &p : ps) p = binary(engine);
//     for (array<uint8_t, lvl1param::n> &i : p1)
//         for (uint8_t &p : i) p = binary(engine);
//     for (array<uint8_t, lvl1param::n> &i : p0)
//         for (uint8_t &p : i) p = binary(engine);

//     for (int i = 0; i < num_test; i++)
//         for (int j = 0; j < lvl1param::n; j++)
//             pmu1[i][j] = (p1[i][j] > 0) ? lvl1param::μ : -lvl1param::μ;
//     for (int i = 0; i < num_test; i++)
//         for (int j = 0; j < lvl1param::n; j++)
//             pmu0[i][j] = (p0[i][j] > 0) ? lvl1param::μ : -lvl1param::μ;
//     vector<TRGSWFFT<lvl1param>> cs(num_test);
//     vector<TRLWE<lvl1param>> c1(num_test);
//     vector<TRLWE<lvl1param>> c0(num_test);
//     vector<TRLWE<lvl1param>> cres(num_test);

//     for (int i = 0; i < num_test; i++) {
//         Polynomial<TFHEpp::lvl1param> plainpoly = {};
//         plainpoly[0] = ps[i];
//         cs[i] = trgswfftSymEncrypt<lvl1param>(plainpoly, lvl1param::α,
//                                               sk->key.lvl1);
//     }
//     for (int i = 0; i < num_test; i++)
//         c1[i] = trlweSymEncrypt<lvl1param>(pmu1[i], lvl1param::α, sk->key.lvl1);
//     for (int i = 0; i < num_test; i++)
//         c0[i] = trlweSymEncrypt<lvl1param>(pmu0[i], lvl1param::α, sk->key.lvl1);

//     chrono::system_clock::time_point start, end;
//     start = chrono::system_clock::now();
//     for (int test = 0; test < num_test; test++) {
//         CMUXFFT<lvl1param>(cres[test], cs[test], c1[test], c0[test]);
//     }
//     end = chrono::system_clock::now();

//     for (int test = 0; test < num_test; test++) {
//         pres = trlweSymDecrypt<lvl1param>(cres[test], sk->key.lvl1);
//         for (int i = 0; i < lvl1param::n; i++)
//             assert(pres[i] ==
//                   (((ps[test] > 0) ? p1[test][i] : p0[test][i]) > 0));
//     }
//     cout << "Passed" << endl;
//     double elapsed =
//         std::chrono::duration_cast<std::chrono::microseconds>(end - start)
//             .count();
//     cout << elapsed / num_test << "μs" << endl;
    
// }

void add_replace_and_test(uint32_t plain_bits, int num_test)
    {
        std::cout << "------ Test of Add_replace_Andn Function tttttt------" << std::endl;
        std::cout << "Test Time : " << num_test << std::endl;
        std::cout << "Plain bits : " << plain_bits << std::endl;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        using P = Lvl1;
        TFHESecretKey sk;
        TFHEEvalKey ek;
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        uint32_t scale_bits;
        
        std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
        std::uniform_int_distribution<typename P::T> type(0, 1);
        scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
        
        std::cout << "scale_bits : " << scale_bits << std::endl;
        
        int num = 7; 
        std::vector<Lvl1::T> p0_list(num, 0);
        std::vector<Lvl1::T> p1_list(num, 0);
        std::vector<typename std::make_signed<typename P::T>::type> gres_list(num, 0);
        std::vector<Lvl1::T> dgres_list(num, 0);
        std::vector<uint32_t> error_time_list(num +2, 0);
        
        std::vector<TLWELvl1> c0_list(num);
        std::vector<TLWELvl1> c1_list(num);
        std::vector<TLWELvl1> cres_list(num);
        
        uint32_t u_right_move_bit = 2;
        bool result_type = LOGIC;
        
        for (int test_time = 0; test_time < num_test; test_time++) 
        {
            for(int i = 0; i < num; i++)
            {
                p0_list[i] = message(engine);
                p1_list[i] = message(engine);
                
                if(p0_list[i] >= p1_list[i]) gres_list[i] = 1;
                else gres_list[i] = 0;
                
                c0_list[i] = TFHEpp::tlweSymInt32Encrypt<P>(p0_list[i], P::α, pow(2., scale_bits), sk.key.get<P>());
                c1_list[i] = TFHEpp::tlweSymInt32Encrypt<P>(p1_list[i], P::α, pow(2., scale_bits), sk.key.get<P>());
                
                //std::cout<< "gres["<< i << "]: "<< gres_list[i]<<"   p0["<< i << "]:"<< p0_list[i] << "   p1["<< i << "]:"<< p1_list[i]  << std::endl;
                
                my_greater_than_equal<P>(c0_list[i], c1_list[i], cres_list[i], plain_bits, ek, result_type, u_right_move_bit +1);
                
                //std::cout<<"Greater than-----Decrypt--------------------:"<< std::endl;
                if(IS_ARITHMETIC(result_type)) dgres_list[i] = TFHEpp:: tlweSymInt32Decrypt<P>(cres_list[i], pow(2., 31), sk.key.lvl1);
                else dgres_list[i] = TFHEpp::tlweSymDecrypt<P>(cres_list[i], sk.key.lvl1);
                
                //std::cout<<"dgres["<< i << "]: "<< dgres_list[i] << std::endl;            
                if (gres_list[i] != dgres_list[i]){
                    error_time_list[i] += 1;
                }
                
                cres_list[i][Lvl1 :: n] +=  P::μ >> u_right_move_bit;  
                // cres_list[i][Lvl1 :: n] + 2^27, 让最高有效位0 映射到（以2^28权重位为单位1权重位的表示的）0；最高有效位1映射为1
                
            }
           
            int32_t add = 0;
            TLWELvl1 c_add;
            //typename std::make_signed<typename P::T>::type d_add;
            P::T d_add;
            
            for (size_t i = 0; i < num; i++)
            {
                if (gres_list[i] < 1) add -= 0;
                else add +=1;
            }
            
            for (size_t i = 0; i <= Lvl1 :: n; i++)
            {
                c_add[i] = 0;
                for (size_t j = 0; j < num; j++)
                {
                    c_add[i] += cres_list[j][i];
                }
                
            }
            
            //std::cout<<"d_add-----Decrypt----------------------------:"<< std::endl;
            d_add = TFHEpp::tlweSymInt32Decrypt<P>(c_add, pow(2., 29 - u_right_move_bit +1), sk.key.lvl1);
            
            // std::cout<<"d_add:"<< d_add <<std::endl;
            // std::cout<<"add:"<< add <<std::endl;
                    
            if (add != d_add)
            {
                error_time_list[ num ] += 1;
                
                std::cout << "++++++++++++错例打印++++++++++" << std::endl;
                for(size_t i = 0; i < num; i++)
                {
                    std::cout<< "gres["<< i << "]: "<< gres_list[i]<< "dgres["<< i << "]: "<< dgres_list[i]<<"   p0["<< i << "]:"<< p0_list[i] <<  "   p1["<< i << "]:"<< p1_list[i] << std::endl;
                } 
                
                std::cout<<"Greater than-----Decrypt--------------------:"<< std::endl;
                for(size_t i = 0; i < num; i++)
                {
                    TFHEpp::tlweSymDecrypt_print<P>(cres_list[i], sk.key.lvl1);
                    std::cout<<"dgres["<< i << "]:"<< dgres_list[i] <<std::endl;
                }     
                  
                std::cout<<"d_add-----Decrypt----------------------------:"<< std::endl;
                TFHEpp::tlweSymInt32Decrypt_print<P>(c_add, pow(2., 29 - u_right_move_bit +1 ), sk.key.lvl1);
                std::cout<<"d_add:"<< d_add <<std::endl;
                std::cout<<"add  :"<< add <<std::endl;
            } 
            
            //最后再引入一个比较   明文加密，做密态比较；可以省略明文加密，直接做减法，再提取最高有效位
            TLWELvl1 add_cres;
            uint32_t filter = 7;
            uint32_t add_plain_bits = std::numeric_limits<typename P::T>::digits - (29 - u_right_move_bit);
            P::T add_dres;
            
            //uint32_t  threshold = (1ULL << (add_plain_bits-1)) + filter -1;
            uint32_t  threshold = filter-1;
        
            TLWELvl1 c_threshold = TFHEpp::tlweSymInt32Encrypt<P>(threshold, P::α, pow(2., 29 - u_right_move_bit +1), sk.key.get<P>());
            
            //c_add[Lvl1 :: n] += 1ULL << (std::numeric_limits<typename P::T>::digits -1);
            
            // std::cout<<"--------c_add += 2^31:---------- "<<std::endl;
            // TFHEpp::my_tlweSymDecrypt<P>(c_add, sk.key.lvl1);
            
            TLWELvl1 c_subs;
            for (size_t i = 0; i <= Lvl1 :: n; i++)
                {
                    c_subs[i] = c_threshold[i] - c_add[i];
                }
            // std::cout<<"c_threshold[i] - c_add[i]+++++++++ Decrypt +++++++++++++++++:"<<std::endl;
            // TFHEpp::tlweSymInt32Decrypt<P>(c_subs, pow(2., scale_bits), sk.key.lvl1);
            
            my_greater_than<P>(c_add, c_threshold, add_cres, plain_bits, ek, ARITHMETIC,0);
            
            //std::cout<<"Threshold Greater than-----Decrypt--------------------:"<< std::endl;
            add_dres = TFHEpp::tlweSymInt32Decrypt<P>(add_cres, pow(2., 31), sk.key.lvl1);
            
            if(!(d_add == 7 && add_dres == 1 || d_add < 7 && add_dres == 0))
            {
                std::cout<<"--------c_add += 2^31:---------- "<<std::endl;
                TFHEpp::tlweSymDecrypt_print<P>(c_add, sk.key.lvl1);
                std::cout<<"c_threshold[i] - c_add[i]+++++++++ Decrypt +++++++++++++++++:"<<std::endl;
                TFHEpp::tlweSymDecrypt_print<P>(c_subs, sk.key.lvl1);
                //TFHEpp::tlweSymInt32Decrypt<P>(c_subs, pow(2., scale_bits), sk.key.lvl1);
                
                std::cout<<"Threshold Greater than-----Decrypt--------------------:"<< std::endl;
                TFHEpp::tlweSymInt32Decrypt_print<P>(add_cres, pow(2., 31), sk.key.lvl1);
                std::cout<<"d_add: "<< d_add <<std::endl;
                std::cout<<"add_dres: "<< add_dres <<std::endl;
                error_time_list[num +1] += 1;
            }
            
        }
        
        for (size_t i = 0; i < num +2; i++)
        {
            std::cout << "Error time" << i << " : " << error_time_list[i] << std::endl;
        }
        std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    }
    
// 验证 “用源码的算子对 ARITHMETIC 比较的结果用ARITHMETIC And 来做是有问题的。”   
void test (uint32_t plain_bits, int num_test) {
    std::cout << "------ Test of TLWE lvl1param Comparison Function tttttt------" << std::endl;
    std::cout << "Test Time : " << num_test << std::endl;
    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    uint32_t scale_bits;
    
    
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    p0 = 3;
    p1 = 2;
    p2 = 3;
    p3 = 4;
    
    bool result_type = ARITHMETIC;
    
    if (p0 > p1) gres = 1;
    else gres = 0;
    
    if (p2 < p3) lres = 1;
    else lres = 0;
    
    c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
    c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
    c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());
    c3 = TFHEpp::tlweSymInt32Encrypt<P>(p3, P::α, pow(2., scale_bits), sk.key.get<P>());
    
    greater_than<P>(c0, c1, cres, plain_bits, ek, result_type);
    
    std::cout << "--------- dgres ---------" << std::endl;
    if(IS_ARITHMETIC(result_type)) dgres = TFHEpp::tlweSymInt32Decrypt<P>(cres, pow(2., 31), sk.key.lvl1);
    else dgres = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
    
    std::cout << " dgres: " << dgres << std::endl;
    
    
    cres1 = cres;
    cres2 = cres;
    
    less_than<P>(c2, c3, cres, plain_bits, ek, result_type);
    
    std::cout << "--------- dlres ---------" << std::endl;
    if(IS_ARITHMETIC(result_type)) dlres = TFHEpp::tlweSymInt32Decrypt_print<P>(cres, pow(2., 31), sk.key.lvl1);
    else dlres = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
    
    std::cout << " dlres: " << dlres << std::endl;
    
    HomAND(c_ands, cres, cres1, ek, result_type);
    
    std::cout << "--------- and ---------" << std::endl;
    
    if(IS_ARITHMETIC(result_type)) d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c_ands, pow(2., 31), sk.key.lvl1);
    else d_ands = TFHEpp::tlweSymDecrypt_print<P>(c_ands, sk.key.lvl1);
    
    std::cout << " d_ands: " << d_ands << std::endl;
    
    
    
    for (size_t i = 0; i <= Lvl1 :: n; i++)
    {
        c_ands[i] = cres[i] + cres1[i] + cres2[i];
    }
    
    std::cout << "--------- d_add_three ---------" << std::endl;
    d_ands = TFHEpp::tlweSymInt32Decrypt_print<P>(c_ands, pow(2., 29 ), sk.key.lvl1);
    
    
    std::cout << " d_add_three: " << d_ands << std::endl;

}   
    

void tlwelvl1_comparison_test(uint32_t plain_bits, int num_test)
{
    std::cout << "------ Test of TLWE lvl1param Comparison Function tttttt------" << std::endl;
    std::cout << "Test Time : " << num_test << std::endl;
    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    uint32_t scale_bits;
    
    ////////////////////////////////////////////////////////
    //std::vector<uint32_t> error_time(5, 0);
    //std::vector<double> comparison_time(5, 0.);
    std::vector<uint32_t> error_time(7, 0);
    std::vector<double> comparison_time(7, 0.);
    ////////////////////////////////////////////////////////////

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;

    
 ///////////////////////////////////////////////////////////////////////////  
    // uint32_t a = 4294967240;
    // uint32_t d = 3224961240;
    // std::cout<<"a  :"<< bitset<32>(a) << std::endl;
    // std::cout<<"d  :"<< bitset<32>(d) << std::endl;
    // uint32_t ad = a + d;
    // std::cout<<"ad :"<< bitset<32>(ad) << std::endl;
    // //std::cout<<"-a :"<< bitset<32>(-a) << std::endl;
    
    // int32_t  ad1 = static_cast<int32_t>(ad);
    // std::cout<<"ad1:"<< bitset<32>(ad1) << std::endl;
    // int32_t  d1 = static_cast<int32_t>(d);
    // std::cout<<"d1 :"<< bitset<32>(d1) << std::endl;
    // int32_t a_1 = ad1 - d1;
    // std::cout<<"a_1:"<< bitset<32>(a_1) << std::endl;
    // uint32_t a1 = static_cast<uint32_t>(a_1);
    // std::cout<<"a1:"<< bitset<32>(a1) << std::endl;
    // std::cout<<"a1:"<< a1 << std::endl;
 //////////////////////////////////////////////////////////////////////////// 

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine);
        p2 = message(engine);
        p3 = message(engine);
        
        bool result_type = ARITHMETIC;
        if (p0 > p1) gres = 1;
        else gres = 0;
        if (p0 >= p1) geres = 1;
        else geres = 0;
        
        /////////////////////////////////////
        // if (p0 < p1) lres = 1;
        // else lres = 0;
        if (p2 < p3) lres = 1;
        else lres = 0;
        /////////////////////////////////////
        
        if (p0 <= p1) leres = 1;
        else leres = 0;
        if (p0 == p1) eres = 1;
        else eres = 0;
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());
        c3 = TFHEpp::tlweSymInt32Encrypt<P>(p3, P::α, pow(2., scale_bits), sk.key.get<P>());
        
        // std::cout<<"p0:"<< p0 << std::endl;
        // std::cout<<"p1:"<< p1 << std::endl;
        // std::cout<<"p2:"<< p2 << std::endl;
        // std::cout<<"p3:"<< p3 << std::endl;
        
        // std::cout<<"c0[n]:"<< c0[Lvl1 :: n] << std::endl;
        // std::cout<<"c0 +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
        // d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c0, pow(2., scale_bits), sk.key.lvl1);
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     c0[i] = c0[i] << 1;
        // }
        // std::cout<<"c0[n] << 1:"<< c0[Lvl1 :: n] << std::endl;
        //  std::cout<<"c0 << 1 +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
        // d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c0, pow(2., scale_bits-2), sk.key.lvl1);
        
        // std::cout<<"c0 +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
        // d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c0, pow(2., scale_bits), sk.key.lvl1);
        // std::cout<<"c1 +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
        // d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c1, pow(2., scale_bits), sk.key.lvl1);
        
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     c_ands[i] = c0[i] + c1[i];
        // }
        // d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c_ands, pow(2., scale_bits), sk.key.lvl1);
        // std::cout<<"c0+c1解密:"<< d_ands << std::endl;
        // d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c_ands, pow(2., scale_bits+1), sk.key.lvl1);
        // std::cout<<"c0+c1解密+1:"<< d_ands << std::endl;
        // d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c_ands, pow(2., scale_bits-2), sk.key.lvl1);
        // std::cout<<"c0+c1解密+2:"<< d_ands << std::endl;

        //////////////////////////////////
        if (gres == 1 && geres == 1) ands = 1;
        else ands = 0;
        //////////////////////////////////
        
        
        for (size_t i = 0; i <= Lvl1 :: n; i++)
            {
                c_ands[i] = c1[i] - c0[i];
            }
        // std::cout<<"c1 - c0 +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
        // d_ands = TFHEpp::tlweSymInt32Decrypt_print<P>(c_ands, pow(2., scale_bits), sk.key.lvl1);
        
        
        //Greater than
        start = std::chrono::system_clock::now();
        
        //greater_than<P>(c0, c1, cres, plain_bits, ek, result_type);
        result_type = LOGIC;
        my_greater_than<P>(c0, c1, cres, plain_bits, ek, ARITHMETIC, 1);
        
        // greater_than<P>(c0, c1, cres, plain_bits, ek, LOGIC);
        // std::cout<<"LOGIC Greater than-----Decrypt--------------------:"<< std::endl;
        // dgres = TFHEpp:: tlweSymInt32Decrypt<P>(cres, pow(2., 31), sk.key.lvl1);
       


        //cres[Lvl1::k * Lvl1::n] += 1U << 27;             看这里！！！！！这里自己的修改注释掉了哦
         
        end = std::chrono::system_clock::now();
        //ARI_to_LOG(cres, cres, ek);
        //result_type = LOGIC;
        comparison_time[0] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
       
        //std::cout<<"Greater than-----Decrypt--------------------:"<< std::endl;
        if(IS_ARITHMETIC(result_type)) dgres = TFHEpp:: tlweSymInt32Decrypt<P>(cres, pow(2., 31 -1), sk.key.lvl1);
        else dgres = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);

        //std::cout<<"dgres:"<< dgres << std::endl;
        
        if (gres != dgres){
            error_time[0] += 1;
        
            std::cout<<"gres:"<< gres << std::endl;
            std::cout<<"dgres:"<< dgres << std::endl;
            std::cout<<"p0:"<< p0 << std::endl;
            std::cout<<"p1:"<< p1 << std::endl;
            std::cout<<"c1 - c0 +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
            d_ands = TFHEpp::tlweSymInt32Decrypt_print<P>(c_ands, pow(2., scale_bits), sk.key.lvl1);
            
            TLWELvl1 shift_tlwe, sign_tlwe5, res;
            
            for (size_t i = 0; i <= Lvl1 :: n; i++)
            {
                shift_tlwe[i] = c_ands[i] << (plain_bits + 1 - 5);
            }
            
            my_MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC, plain_bits + 1 - 5);
            
            std::cout<<"sign_tlwe5 +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<P>(sign_tlwe5, pow(2., scale_bits), sk.key.lvl1);
            
            for (size_t i = 0; i <= Lvl1 :: n; i++)
            {
                res[i] = c_ands[i] - sign_tlwe5[i];
            }
            
            std::cout<<"c_sub - sign_tlwe5[i]+++++++++ Decrypt +++++++++++++++++:"<<std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<P>(res, pow(2., scale_bits), sk.key.lvl1);
            
            std::cout<<"sign_c_sub +++++++++ Decrypt +++++++++++++++++:"<<std::endl;
            dgres = TFHEpp::tlweSymDecrypt_print<P>(cres, sk.key.lvl1);
        }

        //////////////////////////////////////
        cres1 = cres;
        /////////////////////////////////////

/*
        // greater than equal
        start = std::chrono::system_clock::now();
        greater_than_equal<P>(c0, c1, cres, plain_bits, ek, result_type);
        
        cres[Lvl1::k * Lvl1::n] += 1U << 27;
        
        end = std::chrono::system_clock::now();
        comparison_time[1] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout<<"Greater than or euqal to-------Decrypt-----:"<< std::endl;
        if(IS_ARITHMETIC(result_type)) dgeres = TFHEpp::tlweSymInt32Decrypt<P>(cres, pow(2., 31), sk.key.lvl1);
        else dgeres = TFHEpp::my_tlweSymDecrypt<P>(cres, sk.key.lvl1);
        std::cout<<"dgeres:"<< dgeres << std::endl;

        if (geres != dgeres) error_time[1] += 1;
*/


        // less than  
        start = std::chrono::system_clock::now();
        //less_than<P>(c2, c3, cres, plain_bits, ek, result_type);
        less_than<P>(c2, c3, cres, plain_bits, ek, result_type);
    
        //cres[Lvl1::k * Lvl1::n] += 1U << 27;              看这里！！！！！这里自己的修改注释掉了哦
        
        end = std::chrono::system_clock::now();
        comparison_time[2] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        //std::cout<<"less than-----Decrypt------------------------:"<< std::endl;
        if(IS_ARITHMETIC(result_type)) dlres = TFHEpp::tlweSymInt32Decrypt<P>(cres, pow(2., 31 -1), sk.key.lvl1);
        else dlres = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
        //std::cout<<"dlres:"<< dlres << std::endl;
        
        if (lres != dlres) error_time[2] += 1;

/*        
        //HomAnd
        result_type = LOGIC;
        ARI_to_LOG(cres, cres, ek);
        ARI_to_LOG(cres1, cres1, ek);
        TLWELvl1 res;
        Lvl1::T offset = Lvl1::μ;                                 
        if(IS_ARITHMETIC(result_type)) offset = (offset << 1);
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            res[i] = cres[i] + cres1[i];
        
        std::cout<<"res----------------Decrypt---------------" << std::endl;
        TFHEpp::tlweSymInt32Decrypt<P>(res, pow(2., 31), sk.key.lvl1);
        res[Lvl1::k * Lvl1::n] -= Lvl1::μ >> 1;   // - 1/8
        std::cout<<"res-q/8 ------------Decrypt-------------:" << std::endl;
        TFHEpp::tlweSymInt32Decrypt<P>(res, pow(2., 31), sk.key.lvl1);
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, res, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<Lvl1>(-offset));
        std::cout<<"res_Boot------------Decrypt-------------" << std::endl;
        TFHEpp::tlweSymInt32Decrypt<P>(res, pow(2., 31), sk.key.lvl1);
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += offset;
        
        
        //And
        ////////////////////////////////////
        start = std::chrono::system_clock::now();
        HomAND(c_ands, cres, cres1, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[5] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        std::cout<<"c_ands----------Decrypt-----------------" << std::endl;
        if(IS_ARITHMETIC(result_type)) d_ands = TFHEpp::tlweSymInt32Decrypt<P>(c_ands, pow(2., 31), sk.key.lvl1);
        else d_ands = TFHEpp::tlweSymDecrypt<P>(c_ands, sk.key.lvl1);
        
        std::cout<<"d_ands:"<< d_ands << std::endl;
        if (ands != d_ands) error_time[5] += 1;
        ///////////////////////////////////
*/
        
        
        // typename P::T phase_c1 = c[P::k * P::n];
        // typename P::T phase_c = c[P::k * P::n];
        // for (int k = 0; k < P::k; k++){
        //     for (int i = 0; i < P::n; i++){
        //         phase_c1 -= cres1[k * P::n + i] * sk.key.lvl1[k * P::n + i];
        //         phase_c  -= cres[k * P::n + i] * sk.key.lvl1[k * P::n + i];
        //     }
        // }        
        // std::cout<<"phase_c1:"<< phase_c1 << std::endl;
        // std::cout<<"phase_c:"<< phase_c << std::endl;
        
  
        // Add
        ////////////////////////////////////
        TLWELvl1 c_add;
        //typename P::T d_add;
        typename std::make_signed<typename P::T>::type d_add;
        
        typename std::make_signed<typename P::T>::type gres1 = static_cast<typename std::make_signed<typename P::T>::type>(gres);
        typename std::make_signed<typename P::T>::type lres1 = static_cast<typename std::make_signed<typename P::T>::type>(lres);
        
        start = std::chrono::system_clock::now();
       
        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_add[i] = cres[i] + cres1[i];
        }
        end = std::chrono::system_clock::now();
        comparison_time[6] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout<<"d_add-----Decrypt----------------------------:"<< std::endl;
        
        // if(IS_ARITHMETIC(result_type)) d_add = TFHEpp::my_tlweSymInt32Decrypt<P>(c_add, pow(2., 31), sk.key.lvl1, 4);
        // else d_add = TFHEpp::my_tlweSymDecrypt<P>(c_add, sk.key.lvl1);
        d_add = TFHEpp::my_tlweSymInt32Decrypt_print<P>(c_add, pow(2., 28), sk.key.lvl1);

        std::cout<<"d_add:"<< d_add <<std::endl;
        std::cout << "+++++++++++++++++++++++++++++++" << std::endl;
        if (gres1 < 1) gres1 -= 1;
        if (lres1 < 1) lres1 -= 1;
        if (gres1 + lres1 != d_add){
            error_time[6] += 1;
        } 
        ////////////////////////////////// 

        //less than or equal to
        start = std::chrono::system_clock::now();
        less_than_equal<P>(c0, c1, cres, plain_bits, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[3] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) dleres = TFHEpp::tlweSymInt32Decrypt<P>(cres, pow(2., 31), sk.key.lvl1);
        else dleres = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
        if (leres != dleres) error_time[3] += 1;

        //equal to
        start = std::chrono::system_clock::now();
        equal<P>(c0, c1, cres, plain_bits, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[4] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) deres = TFHEpp::tlweSymInt32Decrypt<P>(cres, pow(2., 31), sk.key.lvl1);
        else deres = TFHEpp::tlweSymDecrypt<P>(cres, sk.key.lvl1);
        if (eres != deres) error_time[4] += 1;

    }
    

    
    std::cout << "Test greater than, greater than or equal, less than, less than or equal , equal to , and ,add" << std::endl;
    for (size_t i = 0; i < 7; i++)
    {
        std::cout << "Error time" << i << " : " << error_time[i] << std::endl;
        std::cout << "Time" << i << " : " << comparison_time[i] / num_test << "ms" << std::endl;
    }
    std::cout << "total_Time" << 6 << " : " << comparison_time[6]  << "ms" << std::endl;
      
}

// Lvl2 include 10 - 32 bits
void tlwelvl2_comparison_test(uint32_t plain_bits, int num_test)
{
    std::cout << "------ Test of TLWE lvl2param Comparison Function------" << std::endl;
    std::cout << "Test Time : " << num_test << std::endl;
    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl2;
    TFHESecretKey sk;
    TFHEEvalKey ek;

    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);
    uint32_t scale_bits;

    ////////////////////////////////////////////////////////
    //std::vector<uint32_t> error_time(5, 0);
    //std::vector<double> comparison_time(5, 0.);
    std::vector<uint32_t> error_time(7, 0);
    std::vector<double> comparison_time(7, 0.);
    ////////////////////////////////////////////////////////////

    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits -1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    typename P::T p0,p1, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0,c1,c;
    TLWELvl1 cres, cres1, c_ands;
    std::chrono::system_clock::time_point start, end;
    for (int test = 0; test < num_test; test++) 
    {
        p0 = message(engine);
        p1 = message(engine);
        bool result_type = LOGIC;
        if (p0 > p1) gres = 1;
        else gres = 0;
        if (p0 >= p1) geres = 1;
        else geres = 0;
        if (p0 < p1) lres = 1;
        else lres = 0;
        if (p0 <= p1) leres = 1;
        else leres = 0;
        if (p0 == p1) eres = 1;
        else eres = 0;
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());

        //////////////////////////////////
        if (gres == 1 && geres == 1) ands = 1;
        else ands = 0;
        //////////////////////////////////


        //Greater than
        start = std::chrono::system_clock::now();
        greater_than<P>(c0, c1, cres, plain_bits, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[0] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) dgres = TFHEpp::tlweSymInt32Decrypt<Lvl1>(cres, pow(2., 31), sk.key.lvl1);
        else dgres = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (gres != dgres) error_time[0] += 1;


        //////////////////////////////////////
        cres1 = cres;
        /////////////////////////////////////


        // Greater than or euqal to
        start = std::chrono::system_clock::now();
        greater_than_equal<P>(c0, c1, cres, plain_bits, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[1] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) dgeres = TFHEpp::tlweSymInt32Decrypt<Lvl1>(cres, pow(2., 31), sk.key.lvl1);
        else dgeres = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (geres != dgeres) error_time[1] += 1;

        // And
        //////////////////////////////////////
        start = std::chrono::system_clock::now();
        HomAND(c_ands, cres, cres1, ek, LOGIC);
        end = std::chrono::system_clock::now();
        comparison_time[5] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) d_ands = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ands, pow(2., 31), sk.key.lvl1);
        else d_ands = TFHEpp::tlweSymDecrypt<Lvl1>(c_ands, sk.key.lvl1);
        if (ands != d_ands) error_time[5] += 1;
        /////////////////////////////////////

        // Add
        ////////////////////////////////////
        TLWELvl2 add;
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            add[i] = c0[i] + c1[i];
        }
        end = std::chrono::system_clock::now();
        comparison_time[6] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        //////////////////////////////////


        // less than 
        start = std::chrono::system_clock::now();
        less_than<P>(c0, c1, cres, plain_bits, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[2] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) dlres = TFHEpp::tlweSymInt32Decrypt<Lvl1>(cres, pow(2., 31), sk.key.lvl1);
        else dlres = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (lres != dlres) error_time[2] += 1;

        //less than or equal to
        start = std::chrono::system_clock::now();
        less_than_equal<P>(c0, c1, cres, plain_bits, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[3] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) dleres = TFHEpp::tlweSymInt32Decrypt<Lvl1>(cres, pow(2., 31), sk.key.lvl1);
        else dleres = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (leres != dleres) error_time[3] += 1;

        //equal to
        start = std::chrono::system_clock::now();
        equal<P>(c0, c1, cres, plain_bits, ek, result_type);
        end = std::chrono::system_clock::now();
        comparison_time[4] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        if(IS_ARITHMETIC(result_type)) deres = TFHEpp::tlweSymInt32Decrypt<Lvl1>(cres, pow(2., 31), sk.key.lvl1);
        else deres = TFHEpp::tlweSymDecrypt<Lvl1>(cres, sk.key.lvl1);
        if (eres != deres) error_time[4] += 1;
    }
    std::cout << "Test greater than, greater than or equal, less than, less than or equal , equal to , and , add" << std::endl;
    for (size_t i = 0; i < 7; i++)
    {
        std::cout << "Error time" << i << " : " << error_time[i] << std::endl;
        std::cout << "Time" << i << " : " << comparison_time[i] / num_test << "ms" << std::endl;
    }
    
    
}

void msb_related_test(uint32_t plain_bits, int num_test)
{
    uint32_t error_time = 0;
     
    for (int test_time = 0; test_time < num_test; test_time++) 
    {
        std::cout << "test_time : " << test_time << std::endl;
        std::cout << "------ Test of msb_related tttttt------" << std::endl;
        std::cout << "Test Time : " << num_test << std::endl;
        std::cout << "Plain bits : " << plain_bits << std::endl;
        std::random_device seed_gen;
        std::default_random_engine engine(seed_gen());
        using P = Lvl1;
        TFHESecretKey sk;
        TFHEEvalKey ek;
        
        ek.emplacebkfft<Lvl01>(sk);
        ek.emplacebkfft<Lvl02>(sk);
        ek.emplaceiksk<Lvl20>(sk);
        ek.emplaceiksk<Lvl10>(sk);
        ek.emplaceiksk<Lvl21>(sk);

        uint32_t scale_bits;
        
        std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
        std::uniform_int_distribution<typename P::T> type(0, 1);
        scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
        
        std::cout << "scale_bits : " << scale_bits << std::endl;
        
        //Lvl1::T a = message(engine);
        //Lvl1::T b = message(engine);
        Lvl1::T a = 0;
        Lvl2::T b = 0;

        Lvl1::T da, db, d_subs, d_res;
        Lvl2::T db2 ;
        std::cout << "a : " << a << std::endl;
        
        TLWELvl1 ca = TFHEpp::tlweSymInt32Encrypt<P>(a, P::α, pow(2., scale_bits), sk.key.get<P>());
        TLWELvl2 cb = TFHEpp::tlweSymInt32Encrypt<Lvl2>(b, Lvl2::α, pow(2., scale_bits), sk.key.get<Lvl2>());
        TLWELvl1 c_sub, res;
        TLWELvl2 res2;
        
        da = TFHEpp::tlweSymInt32Decrypt_print<P>(ca, pow(2., scale_bits ), sk.key.lvl1);
        //db2 = TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(cb, pow(2., scale_bits ), sk.key.lvl2);
        
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     c_subs[i] = ca - cb;
        // }
        
        // d_subs = TFHEpp::tlweSymInt32Decrypt_print<P>(c_subs, pow(2., scale_bits ), sk.key.lvl1);
        
        my_MSBGateBootstrapping_2(res, ca, ek, ARITHMETIC,3);
        
        //std::cout << "-----------MSB(ca) Decrypt------------ : " << std::endl;
        //d_res = TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res2, pow(2., 25), sk.key.lvl2);
        d_res = TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.lvl1);
        
        if( d_res != 0 )
        {
            error_time ++ ;
            std::cout << "test_time : " << test_time << std::endl;
            // da = TFHEpp::tlweSymInt32Decrypt_print<P>(ca, pow(2., scale_bits ), sk.key.lvl1);
            // std::cout << "-----------MSB(ca) Decrypt------------ : " << std::endl;
            // d_res = TFHEpp::tlweSymInt32Decrypt_print<P>(res, pow(2., 31), sk.key.lvl1);
        }
    }
    
    std::cout << "error_time : " << error_time << std::endl;
}

void euqual_to_test(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of equal to------" << std::endl;
    std::cout << "Test Time : " << num_test << std::endl;
    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    uint32_t scale_bits;
    
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
        
    for (int test = 0; test < num_test; test++) 
        {   
            // p0 = message(engine);
            // p1 = message(engine);
            uint32_t p0 = 0;
            uint32_t p1 = 1;
            
            
            TLWELvl1 c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
            TLWELvl1 c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
            TLWELvl1 cres;
            //equal to
            
            equal<P>(c0, c1, cres, plain_bits, ek, LOGIC);
            
            uint32_t cres_de = TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(cres, pow(2., 29), sk.key.get<Lvl1>());
            

        }
}

int main()
{
    int num_test = 10000;
    euqual_to_test(6, 5);
    //my_comparison_test(5,1);
    //n_and_test(32,100);
    //add_replace_and_test(4, num_test);
    //msb_related_test(5, num_test);
    //cmux_test(5, 20);
    //comparison_modify_test(4,10);
    
    //test(4,1);
    //tlwelvl1_comparison_test(6, num_test);
    //tlwelvl2_comparison_test(16, num_test);
    //tlwelvl2_comparison_test(32, num_test);
}