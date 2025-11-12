#include <iostream>
#include <chrono>
#include <random>
#include<bitset>
#include "HEDB/comparison/comparison.h"
#include "HEDB/utils/utils.h"

#include "HEDB/comparison/tfhepp_utils.h"
#include "HEDB/comparison/extract_msb.h"
#include <gatebootstrapping.hpp>
#include "detwfa.hpp"

#include <cassert>
#include <tfhe++.hpp>


#include <bits/stdint-uintn.h>
#include <tfhe++.hpp>



using namespace std;
using namespace HEDB;
using namespace TFHEpp;


void encrypt_add_test(uint32_t plain_bits, int num_test)
{
    //测试三个1，加密之后相加解密是否正确
    std::cout << "------ Test of add Function tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = 1;
        p1 = 1;
        p2 = 1;
        
        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        p0 = TFHEpp::my_tlweSymInt32Decrypt_print<P>(c0, pow(2., scale_bits), sk.key.lvl1);


        // Add
        ////////////////////////////////////
        TLWELvl1 c_add;
       
        typename std::make_signed<typename P::T>::type d_add;
        
        start = std::chrono::system_clock::now();
       
        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_add[i] = c0[i] + c1[i] + c2[i];
        }
        end = std::chrono::system_clock::now();
        comparison_time[6] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        std::cout<<"d_add-----Decrypt----------------------------:"<< std::endl;
        
        // if(IS_ARITHMETIC(result_type)) d_add = TFHEpp::my_tlweSymInt32Decrypt<P>(c_add, pow(2., 31), sk.key.lvl1, 4);
        // else d_add = TFHEpp::my_tlweSymDecrypt<P>(c_add, sk.key.lvl1);
        d_add = TFHEpp::my_tlweSymInt32Decrypt_print<P>(c_add, pow(2., scale_bits), sk.key.lvl1);

        std::cout<<"d_add:"<< d_add <<std::endl;

    }
}

void my_comparation_test(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of " << plain_bits << "bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ;
    double comparison_time_modified =0; 

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        // p1 = 0;
        // p2 = 0;
        
        int ans = 0;
        if(p2 >= p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        //p0 = TFHEpp::tlweSymInt32Decrypt_print<P>(c0, pow(2., scale_bits), sk.key.lvl1);
        //std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;

        // Add
        ////////////////////////////////////
        TLWELvl1 c_boot;
        TLWELvl1 c_sub;
        TLWELvl1 c_sub2;
        TLWELvl1 c_ans;

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }

        // std::cout<<"-------sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub, pow(2., 27), sk.key.get<Lvl1>());

        //std::cout<<"-------my_comparation-------:"<< std::endl;
        

        ////////////////////my_ExtractMSB9

        // TLWELvl1 shift_tlwe, sign_tlwe5, res;
        // uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits-1;
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     shift_tlwe[i] = c_sub[i] << (plain_bits+1 - 5);
        // }
        // std::cout<<"-------shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe, pow(2., 31), sk.key.get<Lvl1>());

        // std::cout<<"-------1-----boot---------:"<< std::endl;
        // my_MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC, plain_bits+1 - 5);
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe5, pow(2., 31), sk.key.get<Lvl1>());
        
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     res[i] = c_sub[i] - sign_tlwe5[i];
        // }
        // std::cout<<"-------1-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());
        // my_MSBGateBootstrapping(res, res, ek, LOGIC, plain_bits+1 - 5);
        // std::cout<<"-------boot -2---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());

        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than_equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
        
        TLWELvl1 greater_tlwe, less_tlwe, c_add6, c_sub6;
        start = std::chrono::system_clock::now();
        // my_greater_than_equal<P>(c1, c2, greater_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // my_less_than_equal<P>(c1, c2, less_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // //my_HomAND(res, greater_tlwe, less_tlwe, c_2,  ek, ARITHMETIC, k);
        
        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_add6[i] = greater_tlwe[i] + less_tlwe[i];

        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_sub6[i] = c_1[i] - c_add6[i];

        // my_MSBGateBootstrapping(c_ans,c_sub6,ek,ARITHMETIC,3);


        //my_not_equal<Lvl1>(c2,c1,c_ans, c_1, plain_bits, ek, ARITHMETIC, 28);
        my_greater_than_equal<Lvl1>(c2,c1,c_ans, plain_bits, ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        //std::cout<<"d_ans :  "<< d_ans << std::endl;

        if(d_ans != ans) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans: "<< ans << " d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << " p2: "<< p2 <<std::endl;
            std::cout<<"----------------de_greater--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(greater_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------de_less--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(less_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------add--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_add6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------sub--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------cans--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}   

void my_comparation_14bit_test(uint32_t plain_bits, int num_test)
{
    std::cout << "------ Test of add Function tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, ((1 << (plain_bits - 1)) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original , comparison_time_modified; 

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        //p1 = (1 << 11)-1;
        // p1 = 310;
        // p2 = 346;
        
        int ans = 0;
        if(p2 > p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());


        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;

        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }

        // std::cout<<"-------sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

        //std::cout<<"-------my_comparation-------:"<< std::endl;
        

        TFHEpp::TLWE<P> shift_tlwe, sign_tlwe5, res;
        TFHEpp::TLWE<Lvl1> shift_tlwe2, sign_tlwe1, res2 ,c_sub2, sign_tlwe3;
        uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits-1;
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            shift_tlwe[i] = c_sub[i] << (plain_bits+1 - 6);
        }
        // std::cout<<"-----1--shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe, pow(2., 61), sk.key.get<P>());

        //std::cout<<"-------1-----boot---------:"<< std::endl;
        my_MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC, plain_bits+1 - 6);
        //TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe5, pow(2., 61), sk.key.get<P>());
        
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res[i] = c_sub[i] - sign_tlwe5[i];
        }
        // std::cout<<"-------1-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res, pow(2., 61), sk.key.get<Lvl2>());

        //std::cout<<"-------key-switch---------:"<< std::endl;
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res2, res, *ek.iksklvl21);
        //TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            shift_tlwe2[i] = res2[i] << 2;
        }
        // std::cout<<"-------2---shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe2, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(sign_tlwe1, shift_tlwe2, ek, ARITHMETIC, 2);
        // std::cout<<"-------boot -2---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe1, pow(2., 31), sk.key.get<Lvl1>());

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub2[i] = res2[i] - sign_tlwe1[i];
        }
        // std::cout<<"------2-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub2, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(sign_tlwe3, c_sub2, ek, ARITHMETIC, 0);
        // std::cout<<"-------boot -3---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe3, pow(2., 31), sk.key.get<Lvl1>());

        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        start = std::chrono::system_clock::now();
        my_greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC, 0);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 29), sk.key.get<Lvl1>());
        //std::cout<<"d_ans :  "<< d_ans << std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
            std::cout<<"-------sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());
            std::cout<<"-----1--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe, pow(2., 61), sk.key.get<P>());
            std::cout<<"-------1-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe5, pow(2., 61), sk.key.get<P>());
            std::cout<<"-------1-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res, pow(2., 61), sk.key.get<Lvl2>());
            std::cout<<"-------key-switch---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"-------2---shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe2, pow(2., 31), sk.key.get<Lvl1>());
            
            
            std::cout<<"-------boot -2--inner---------:"<< std::endl;
            shift_tlwe2[lvl1param::k * lvl1param::n] += 1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
            std::cout<<"-------inner--add--------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe2, pow(2., 31), sk.key.get<Lvl1>());
            TLWE<lvl0param> tlwelvl0;
            IdentityKeySwitch<lvl10param>(tlwelvl0, shift_tlwe2, *ek.iksklvl10);
            lvl1param::T μ = lvl1param::μ;
            μ = (μ << 1) >> 2;
            GateBootstrappingTLWE2TLWEFFT<lvl01param>(res2, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
            std::cout<<"-------inner--boot--------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());
            res2[lvl1param::k * lvl1param::n] += (μ);
            std::cout<<"-------inner--out--------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-------boot -2---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe1, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"------2-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub2, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"-------boot -3---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe3, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"-------d_ans---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans, pow(2., 29), sk.key.get<Lvl1>());
            
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/1000 << "s"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/1000  << "s"<< std::endl;
}

void my_comparation_19bit_test(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of" << plain_bits<< " bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ; 
    double comparison_time_modified =0 ;

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        // p1 = 0;
        // p2 = 0;
        
        // int ans_great_than = 0;
        // int ans_great_than_equal = 0;
        // int ans_less_than = 0;
        // int ans_less_than_equal =0;
        // if(p2 > p1) {
        //     ans_great_than = 1;
        // }
        // if(p2 >= p1) {
        //     ans_great_than_equal = 1;
        // }
        // if(p2 < p1) {
        //     ans_less_than = 1;
        // }
        // if(p2 <= p1) {
        //     ans_less_than_equal = 1;
        // }

        int ans = 0;
        if(p2 > p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;
        TFHEpp::TLWE<P> shift_tlwe19, sign_tlwe19, res19,  shift_tlwe14, sign_tlwe14, res14;
        TFHEpp::TLWE<Lvl1>  res14_switch ,shift_tlwe9, sign_tlwe9 , c_sub9, res9;
        
    /*
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }

        std::cout<<"-------sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

        std::cout<<"-------my_comparation-------:"<< std::endl;
    
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            shift_tlwe19[i] = c_sub[i] << (plain_bits+1 - 6);
        }
        std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
        std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

        
        my_MSBGateBootstrapping(sign_tlwe19, shift_tlwe19, ek, ARITHMETIC, plain_bits+1 - 6);
        std::cout<<"-------19-----boot---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
        
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res19[i] = c_sub[i] - sign_tlwe19[i];
        }
        std::cout<<"-------19-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

        std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
             shift_tlwe14[i] = res19[i] << (plain_bits+1-5 - 6);
        }
        std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

        my_MSBGateBootstrapping(sign_tlwe14, shift_tlwe14, ek, ARITHMETIC, plain_bits+1-5 - 6);
        std::cout<<"-------14-----boot---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());
        
    
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res14[i] = res19[i] - sign_tlwe14[i];
        }

        std::cout<<"------14-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

        
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res14_switch, res14, *ek.iksklvl21);
        std::cout<<"-------key-switch---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

        std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;
        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            shift_tlwe9[i] = res14_switch[i] << (plain_bits+1-5-5 -5);
        }
        std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(sign_tlwe9, shift_tlwe9, ek, ARITHMETIC, (plain_bits+1-5-5 -5));
        std::cout<<"-------9---boot----------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub9[i] = res14_switch[i] - sign_tlwe9[i];
        }
        std::cout<<"------9-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(res9, c_sub9, ek, ARITHMETIC, 0);
        std::cout<<"------9-----boot2---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
    */
        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, P::α, pow(2., 28), sk.key.get<Lvl1>());
        start = std::chrono::system_clock::now();
        my_greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC, 28);
        //my_not_equal<Lvl2>(c2,c1,c_ans,c_1,plain_bits, ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        
        // std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // std::cout<<"//////////////////////////////////////////////"<< std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
            std::cout<<"-------sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------19-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
            
            std::cout<<"-------19-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

            std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

            std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------14-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"------14-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

            std::cout<<"-------key-switch---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;

            std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-------9---boot----------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----boot2---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}

void my_comparation_24bit_test(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of" << plain_bits<< " bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ; 
    double comparison_time_modified =0 ;

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        // p1 = 0;
        // p2 = 0;
        
        // int ans_great_than = 0;
        // int ans_great_than_equal = 0;
        // int ans_less_than = 0;
        // int ans_less_than_equal =0;
        // if(p2 > p1) {
        //     ans_great_than = 1;
        // }
        // if(p2 >= p1) {
        //     ans_great_than_equal = 1;
        // }
        // if(p2 < p1) {
        //     ans_less_than = 1;
        // }
        // if(p2 <= p1) {
        //     ans_less_than_equal = 1;
        // }

        int ans = 0;
        if(p2 > p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
        
        
        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;
        TFHEpp::TLWE<P> shift_tlwe19, sign_tlwe19, res19,  shift_tlwe14, sign_tlwe14, res14;
        TFHEpp::TLWE<Lvl1>  res14_switch ,shift_tlwe9, sign_tlwe9 , c_sub9, res9;
        
    /*
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }

        std::cout<<"-------sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

        std::cout<<"-------my_comparation-------:"<< std::endl;
    
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            shift_tlwe19[i] = c_sub[i] << (plain_bits+1 - 6);
        }
        std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
        std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

        
        my_MSBGateBootstrapping(sign_tlwe19, shift_tlwe19, ek, ARITHMETIC, plain_bits+1 - 6);
        std::cout<<"-------19-----boot---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
        
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res19[i] = c_sub[i] - sign_tlwe19[i];
        }
        std::cout<<"-------19-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

        std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
             shift_tlwe14[i] = res19[i] << (plain_bits+1-5 - 6);
        }
        std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

        my_MSBGateBootstrapping(sign_tlwe14, shift_tlwe14, ek, ARITHMETIC, plain_bits+1-5 - 6);
        std::cout<<"-------14-----boot---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());
        
    
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res14[i] = res19[i] - sign_tlwe14[i];
        }

        std::cout<<"------14-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

        
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res14_switch, res14, *ek.iksklvl21);
        std::cout<<"-------key-switch---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

        std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;
        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            shift_tlwe9[i] = res14_switch[i] << (plain_bits+1-5-5 -5);
        }
        std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(sign_tlwe9, shift_tlwe9, ek, ARITHMETIC, (plain_bits+1-5-5 -5));
        std::cout<<"-------9---boot----------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub9[i] = res14_switch[i] - sign_tlwe9[i];
        }
        std::cout<<"------9-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(res9, c_sub9, ek, ARITHMETIC, 0);
        std::cout<<"------9-----boot2---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
    */
        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, P::α, pow(2., 28), sk.key.get<Lvl1>());
        start = std::chrono::system_clock::now();
        my_greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC, 28);
        //my_not_equal<Lvl2>(c2,c1,c_ans,c_1,plain_bits, ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        
        // std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // std::cout<<"//////////////////////////////////////////////"<< std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
            std::cout<<"-------sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------19-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
            
            std::cout<<"-------19-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

            std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

            std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------14-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"------14-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

            std::cout<<"-------key-switch---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;

            std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-------9---boot----------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----boot2---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}

void n_and_test( uint32_t k, int num_test){
    std::cout << "n_and_test: " <<  std::endl;
    std::cout << "        n: " << k <<  std::endl;
    std::cout << " num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 2;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0;
    double my_and_time = 0;
    std::chrono::system_clock::time_point start, end;

    TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<P>( 7-1 , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<P>( k%7 -1 , P::α, pow(2., 28), sk.key.get<P>());
    //
    //TLWE<lvl1param> k_times = TFHEpp::tlweSymInt32Encrypt<P>( k/8  , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<P>( 0  , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 1;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);
            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            if(a[i]>b[i]){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }

            ans = ans * arr[i];
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
        }

        

        start = std::chrono::system_clock::now();
        N_HomAND(arr_ciphers ,res , k_rest, k_7, k_0, ek ,29);
        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 29), sk.key.get<P>());
            // for(int i =0 ; i< k ; i++){
            //     std::cout << "-------arr_ciphers[" << i << "]-------"<< std::endl;
            //     TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(arr_ciphers[i], pow(2., 28), sk.key.get<P>());
            // }
            // std::cout << "-------k_rest-------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(k_rest, pow(2., 28), sk.key.get<P>());

            // std::size_t num = arr_ciphers.size();
            // TLWE<lvl1param> sum  = arr_ciphers[0];
            // TLWE<lvl1param> sub ;

            // for(int32_t i = 1 ; i < num ; i++)
            // {
            //     for (size_t j = 0; j <= lvl1param :: n; j++)
            //     {
            //         sum[j] = sum[j] + arr_ciphers[i][j];
            //     }
            // }
            
            // sub = sum;

            // for (size_t i = 0; i <= lvl1param::k * lvl1param::n; i++)
            // {
            //     sub[i] = k_rest[i] - sum[i];
            // }
            
            // my_MSBGateBootstrapping(res, sub, ek, ARITHMETIC, 30 - 29 +1);

            // std::cout << "----------------d_sum-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sum, pow(2., 28), sk.key.get<P>());

            // std::cout << "----------------k_rest-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(k_rest, pow(2., 28), sk.key.get<P>());

            // std::cout << "----------------d_sub-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sub, pow(2., 28), sk.key.get<P>());

            // std::cout << "----------------d_res-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 29), sk.key.get<P>());
        }

        start = std::chrono::system_clock::now();
        TLWELvl1 temp;
        HomAND(arr_ciphers[0], arr_ciphers[1], temp, ek, ARITHMETIC);
        for(int i = 2 ; i< k ; i++){
            HomAND(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "original_and_time: " << original_and_time / num_test << "ms" <<  std::endl;
    std::cout << "my_and_time: " << my_and_time / num_test << "ms" <<  std::endl;
}

void compare_n_and_test( uint32_t k, int num_test){
    std::cout << "n_and_test: " <<  std::endl;
    std::cout << "        n: " << k <<  std::endl;
    std::cout << " num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 4;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0.;
    double my_and_time = 0.;
    double compare_time = 0.;
    std::chrono::system_clock::time_point start, end;

    TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<P>( 7-1 , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<P>( k%7 -1 , P::α, pow(2., 28), sk.key.get<P>());
    //
    //TLWE<lvl1param> k_times = TFHEpp::tlweSymInt32Encrypt<P>( k/8  , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<P>( 0  , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 1;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);
            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            if(a[i]>b[i]){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }

            ans = ans * arr[i];

            start = std::chrono::system_clock::now();
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
            end = std::chrono::system_clock::now();
            compare_time = compare_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        }

        

        start = std::chrono::system_clock::now();
        N_HomAND(arr_ciphers ,res , k_rest, k_7, k_0, ek ,29);
        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            
        }

        start = std::chrono::system_clock::now();
        TLWELvl1 temp;
        HomAND(arr_ciphers[0], arr_ciphers[1], temp, ek, ARITHMETIC);
        for(int i = 2 ; i< k ; i++){
            HomAND(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "compare_time + original_and_time: " << (original_and_time + compare_time) / num_test << "ms" <<  std::endl;
    std::cout << "compare_time + my_and_time: " << (my_and_time + compare_time) / num_test << "ms" <<  std::endl;
}


void q7_and_or_test( uint32_t k, int num_test){
    std::cout << "n_and_test: " <<  std::endl;
    std::cout << "        n: " << k <<  std::endl;
    std::cout << " num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 4;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0.;
    double my_and_time = 0.;
    double compare_time = 0.;
    std::chrono::system_clock::time_point start, end;

    TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<P>( 7-1 , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<P>( k%7 -1 , P::α, pow(2., 28), sk.key.get<P>());
    //
    //TLWE<lvl1param> k_times = TFHEpp::tlweSymInt32Encrypt<P>( k/8  , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<P>( 0  , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 1;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);
            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            if(a[i]>b[i]){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }

            if(i == 3){
                ans = ans*((arr[0] + arr[1] > 1) || (arr[2] + arr[3] > 1) );
            }
            else if(i > 3){
                ans = ans * arr[i];
            }

            

            start = std::chrono::system_clock::now();
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
            end = std::chrono::system_clock::now();
            compare_time = compare_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        }

        std::vector<TLWELvl1> temp_and(2);
        std::vector<TLWELvl1> arr_ciphers2(7); 

        for(int i = 1 ;i < 7; i++){
            arr_ciphers2[i] = arr_ciphers[i+3];
        }

        start = std::chrono::system_clock::now();
        my_HomAND(temp_and[0], arr_ciphers[0], arr_ciphers[1], c_1, ek, ARITHMETIC, 28);
        my_HomAND(temp_and[1], arr_ciphers[2], arr_ciphers[3], c_1, ek, ARITHMETIC, 28);
        N_HomOR(temp_and ,arr_ciphers2[0] , k_0,  ek ,28);

        N_HomAND(arr_ciphers2 ,res , k_rest, k_7, k_0, ek ,29);
        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            
        }


        TLWELvl1 temp;
        TLWELvl1 temp_or;
        std::vector<TLWELvl1> temp_and2(2);

        start = std::chrono::system_clock::now();
    
        HomAND(temp_and2[0], arr_ciphers[0], arr_ciphers[1],  ek, ARITHMETIC);
        HomAND(temp_and2[1], arr_ciphers[2], arr_ciphers[3], ek, ARITHMETIC);
        HomOR(temp_or ,temp_and2[0] ,temp_and2[1] , ek, ARITHMETIC);

        HomAND(temp, temp_or, arr_ciphers[4], ek, ARITHMETIC);
        for(int i = 5 ; i< k ; i++){
            HomAND(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "original_and_time: " << original_and_time  / num_test << "ms" <<  std::endl;
    std::cout << "my_and_time: " << my_and_time  / num_test << "ms" <<  std::endl;
}


void q19_and_or_test( uint32_t k, int num_test){
    k = 21;
    std::cout << "n_and_test: " <<  std::endl;
    std::cout << "         n: " << k <<  std::endl;
    std::cout << "  num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 4;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int ans1 = 1;
    int ans2 = 1;
    int ans3 = 1;

    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0.;
    double my_and_time = 0.;
    double compare_time = 0.;
    std::chrono::system_clock::time_point start, end;

    TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<P>( 7-1 , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<P>( k%7 -1 , P::α, pow(2., 28), sk.key.get<P>());
    //
    //TLWE<lvl1param> k_times = TFHEpp::tlweSymInt32Encrypt<P>( k/8  , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<P>( 0  , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 1;
        ans1 = 1;
        ans2 = 1;
        ans3 = 1;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);
            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            if(a[i]>b[i]){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }

            if(i < 7){
                ans1 = ans1 * arr[i];
            }
            else if(i >= 7 && i < 14){
                ans2 = ans2 * arr[i];
            }
            else if(i >= 14 && i < 21){
                ans3 = ans3 * arr[i];
            }

            start = std::chrono::system_clock::now();
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
            end = std::chrono::system_clock::now();
            compare_time = compare_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }

        ans = ((ans1 + ans2 + ans3)>1);

        std::vector<TLWELvl1> res_part(3); 
        std::vector<TLWELvl1> part1(arr_ciphers.begin(), arr_ciphers.begin() + 7);
        std::vector<TLWELvl1> part2(arr_ciphers.begin() + 7, arr_ciphers.begin() + 2 * 7);
        std::vector<TLWELvl1> part3(arr_ciphers.begin() + 2 * 7, arr_ciphers.end());
       

        start = std::chrono::system_clock::now();
        
        N_HomAND(part1 ,res_part[0] , k_rest, k_7, k_0, ek ,28);
        N_HomAND(part2 ,res_part[1] , k_rest, k_7, k_0, ek ,28);
        N_HomAND(part3 ,res_part[2] , k_rest, k_7, k_0, ek ,28);
        
        N_HomOR(res_part ,res , k_0,  ek ,29);

        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            
        }


        TLWELvl1 temp;
        TLWELvl1 temp_or;
        std::vector<TLWELvl1> temp_and2(2);

        start = std::chrono::system_clock::now();
    
        
        for(int i = 0 ; i< 18 ; i++){
            HomAND(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        for(int i = 0 ; i< 2 ; i++){
            HomOR(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "original_and_time: " << original_and_time  / num_test << "ms" <<  std::endl;
    std::cout << "my_and_time: " << my_and_time  / num_test << "ms" <<  std::endl;
}


void qs1_and_or_test( uint32_t k, int num_test){
    k = 14;
    std::cout << "n_and_test: " <<  std::endl;
    std::cout << "         n: " << k <<  std::endl;
    std::cout << "  num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 4;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    std::vector<bool> d(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int ans1 = 1;
    int ans2 = 1;
    int ans3 = 1;

    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0.;
    double my_and_time = 0.;
    double compare_time = 0.;
    std::chrono::system_clock::time_point start, end;

    TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<P>( 7-1 , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<P>( 6%7 -1 , P::α, pow(2., 28), sk.key.get<P>());
    //
    //TLWE<lvl1param> k_times = TFHEpp::tlweSymInt32Encrypt<P>( k/8  , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<P>( 0  , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 1;
        ans1 = 1;
        ans2 = 1;
        ans3 = 1;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);
            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            a[i] = (a >b );
            

            start = std::chrono::system_clock::now();
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
            end = std::chrono::system_clock::now();
            compare_time = compare_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }

        ans = 
                ((d[0] || d[1] || d[2]) && d[3])
                &&
                ((d[4] || d[5]) && d[6])
                &&
                d[7] 
                && 
                d[8]
                &&
                ((d[9] || d[10]) 
                && 
                (d[11] || d[12] || d[13]))
            ;
             

        std::vector<TLWELvl1> res_part(6); 
        std::vector<TLWELvl1> part1(arr_ciphers.begin(), arr_ciphers.begin() + 3);
        std::vector<TLWELvl1> part2(arr_ciphers.begin() + 4, arr_ciphers.begin() + 6);
        std::vector<TLWELvl1> part3(arr_ciphers.begin() + 9, arr_ciphers.begin() + 11);
        std::vector<TLWELvl1> part4(arr_ciphers.begin() + 11, arr_ciphers.end());

        TLWELvl1 temp_or1;
        TLWELvl1 temp_or2;
        TLWELvl1 temp_or3;
        TLWELvl1 temp_or4;
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());

        start = std::chrono::system_clock::now();
        
        N_HomOR(part1 ,temp_or1 , k_0,  ek  ,28);
        my_HomAND(res_part[0], temp_or1, arr_ciphers[3], c_1, ek, ARITHMETIC, 28);

        N_HomOR(part2 ,temp_or2 , k_0,  ek  ,28);
        my_HomAND(res_part[1], temp_or2, arr_ciphers[6], c_1, ek, ARITHMETIC, 28);

        N_HomOR(part3 ,res_part[2] , k_0,  ek  ,28);
        N_HomOR(part4 ,res_part[3] , k_0,  ek  ,28);
        

        res_part[4] = arr_ciphers[7];
        res_part[5] = arr_ciphers[8];

        N_HomAND(res_part ,res , k_rest, k_7, k_0, ek ,29);
        

        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            
        }


        TLWELvl1 temp;
        TLWELvl1 temp_or;
        std::vector<TLWELvl1> temp_and2(2);

        start = std::chrono::system_clock::now();
    
        
        for(int i = 0 ; i< 7 ; i++){
            HomAND(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        for(int i = 0 ; i< 6 ; i++){
            HomOR(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "original_and_time: " << original_and_time  / num_test << "ms" <<  std::endl;
    std::cout << "my_and_time: " << my_and_time  / num_test << "ms" <<  std::endl;
}

void qs2_and_or_test( uint32_t k, int num_test){
    k = 15;
    std::cout << "n_and_test: " <<  std::endl;
    std::cout << "         n: " << k <<  std::endl;
    std::cout << "  num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 4;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int ans1 = 1;
    int ans2 = 1;
    int ans3 = 1;

    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0.;
    double my_and_time = 0.;
    double compare_time = 0.;
    std::chrono::system_clock::time_point start, end;

    TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<P>( 7-1 , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<P>( 5%7 -1 , P::α, pow(2., 28), sk.key.get<P>());
    //
    //TLWE<lvl1param> k_times = TFHEpp::tlweSymInt32Encrypt<P>( k/8  , P::α, pow(2., 28), sk.key.get<P>());
    TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<P>( 0  , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 1;
        ans1 = 1;
        ans2 = 1;
        ans3 = 1;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);
            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            if(a[i]>b[i]){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }

            if(i < 5){
                ans1 = ans1 * arr[i];
            }
            else if(i >= 5 && i < 10){
                ans2 = ans2 * arr[i];
            }
            else if(i >= 10 && i < 15){
                ans3 = ans3 * arr[i];
            }

            start = std::chrono::system_clock::now();
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
            end = std::chrono::system_clock::now();
            compare_time = compare_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }

        ans = ((ans1 + ans2 + ans3)>0);

        std::vector<TLWELvl1> res_part(3); 
        std::vector<TLWELvl1> part1(arr_ciphers.begin(), arr_ciphers.begin() + 5);
        std::vector<TLWELvl1> part2(arr_ciphers.begin() + 5, arr_ciphers.begin() + 2 * 5);
        std::vector<TLWELvl1> part3(arr_ciphers.begin() + 2 * 5, arr_ciphers.end());
       

        start = std::chrono::system_clock::now();
        
        N_HomAND(part1 ,res_part[0] , k_rest, k_7, k_0, ek ,28);
        N_HomAND(part2 ,res_part[1] , k_rest, k_7, k_0, ek ,28);
        N_HomAND(part3 ,res_part[2] , k_rest, k_7, k_0, ek ,28);
        
        N_HomOR(res_part ,res , k_0,  ek ,29);

        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            
        }


        TLWELvl1 temp;
        TLWELvl1 temp_or;
        std::vector<TLWELvl1> temp_and2(2);

        start = std::chrono::system_clock::now();
    
        
        for(int i = 0 ; i< 12 ; i++){
            HomAND(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        for(int i = 0 ; i< 2 ; i++){
            HomOR(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "original_and_time: " << original_and_time  / num_test << "ms" <<  std::endl;
    std::cout << "my_and_time: " << my_and_time  / num_test << "ms" <<  std::endl;
}

void n_or_test( uint32_t k, int num_test){
    std::cout << "n_or_test: " <<  std::endl;
    std::cout << "        n: " << k <<  std::endl;
    std::cout << " num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 2;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0;
    double my_and_time = 0;
    std::chrono::system_clock::time_point start, end;

    TLWE<lvl1param> k_1 = TFHEpp::tlweSymInt32Encrypt<P>( 0 , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 0;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);

            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            if(a[i]>b[i]){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }

            ans = ans + arr[i];
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
        }
        ans = (ans > 0 ? 1 : 0);

        

        start = std::chrono::system_clock::now();
        N_HomOR(arr_ciphers ,res , k_1,  ek ,29);
        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 29), sk.key.get<P>());
            // for(int i =0 ; i< k ; i++){
            //     std::cout << "-------arr_ciphers[" << i << "]-------"<< std::endl;
            //     TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(arr_ciphers[i], pow(2., 28), sk.key.get<P>());
            // }

            // std::size_t num = arr_ciphers.size();
            // TLWE<lvl1param> sum  = arr_ciphers[0];
            // TLWE<lvl1param> sub ;

            // for(int32_t i = 1 ; i < num ; i++)
            // {
            //     for (size_t j = 0; j <= lvl1param :: n; j++)
            //     {
            //         sum[j] = sum[j] + arr_ciphers[i][j];
            //     }
            // }
            
            // sub = sum;

            // for (size_t i = 0; i <= lvl1param::k * lvl1param::n; i++)
            // {
            //     sub[i] = k_1[i] - sum[i];
            // }
            
            // my_MSBGateBootstrapping(res, sub, ek, ARITHMETIC, 30 - 29 +1);

            // std::cout << "----------------d_sum-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sum, pow(2., 28), sk.key.get<P>());

            // std::cout << "----------------k_1-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(k_1, pow(2., 28), sk.key.get<P>());

            // std::cout << "----------------d_sub-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sub, pow(2., 28), sk.key.get<P>());

            // std::cout << "----------------d_res-----------------"<< std::endl;
            // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 29), sk.key.get<P>());
        }

        start = std::chrono::system_clock::now();
        TLWELvl1 temp;
        HomOR(arr_ciphers[0], arr_ciphers[1], temp, ek, ARITHMETIC);
        for(int i = 2 ; i< k ; i++){
            HomOR(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "original_or_time: " << original_and_time / num_test << "ms" <<  std::endl;
    std::cout << "my_or_time: " << my_and_time / num_test << "ms" <<  std::endl;
}

void compare_n_or_test( uint32_t k, int num_test){
    std::cout << "n_or_test: " <<  std::endl;
    std::cout << "        n: " << k <<  std::endl;
    std::cout << " num_test: " << num_test <<  std::endl;
    using P = Lvl1;
    uint32_t plain_bits = 4;
    uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits - 1;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    TLWELvl1 res ;
    std::vector<TLWELvl1> ca(k);
    std::vector<TLWELvl1> cb(k);

    int ans = 1;
    int d_res = 0;
    std::vector<int> a(k);
    std::vector<int> b(k);
    double original_and_time = 0.;
    double my_and_time = 0.;
    double compare_time = 0.;
    std::chrono::system_clock::time_point start, end;

    TLWE<lvl1param> k_1 = TFHEpp::tlweSymInt32Encrypt<P>( 0 , P::α, pow(2., 28), sk.key.get<P>());

    std::cout << "k/8 : "<< k/8  << std::endl; 
    for (int test = 0; test < num_test; test++) {
        ans = 0;
        for(int i =0 ; i < k; i++){
            
            a[i] = message(engine);
            b[i] = message(engine);

            ca[i] = TFHEpp::tlweSymInt32Encrypt<P>(a[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            cb[i] = TFHEpp::tlweSymInt32Encrypt<P>(b[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            if(a[i]>b[i]){
                arr[i] = 1;
            }
            else{
                arr[i] = 0;
            }

            ans = ans + arr[i];

            start = std::chrono::system_clock::now();
            my_greater_than<P>(ca[i], cb[i], arr_ciphers[i], plain_bits, ek, ARITHMETIC, 28);
            end = std::chrono::system_clock::now();
            compare_time = compare_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        }
        ans = (ans > 0 ? 1 : 0);

        

        start = std::chrono::system_clock::now();
        N_HomOR(arr_ciphers ,res , k_1,  ek ,29);
        end = std::chrono::system_clock::now();

        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "a:[" << i << "]:"<< a[i] << "  a:[" << i << "]:"<< b[i] << std::endl; 
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;
            
        }

        start = std::chrono::system_clock::now();
        TLWELvl1 temp;
        HomOR(arr_ciphers[0], arr_ciphers[1], temp, ek, ARITHMETIC);
        for(int i = 2 ; i< k ; i++){
            HomOR(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "compare_time + original_or_time: " << (original_and_time + compare_time) / num_test << "ms" <<  std::endl;
    std::cout << "compare_time + my_or_time: " << (my_and_time + compare_time) / num_test << "ms" <<  std::endl;
}

void test(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of add Function tttttt------" << std::endl;
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
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original , comparison_time_modified; 

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = 26720;
        p1 = message(engine) ;
        p2 = message(engine);
    
        int ans = 0;
        if(p2 >= p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        std::cout<<"p1: "<< p1  << std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c0, pow(2., scale_bits), sk.key.get<Lvl2>());

        TLWELvl2 c_ans;
        for(int i = 0 ; i < 16 ; i++){
            my_MSBGateBootstrapping(c_ans, c0, ek, LOGIC, i);
            std::cout<<"c_ans---boot"<< i << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_ans, pow(2., scale_bits), sk.key.get<Lvl2>());
        }
        
        ////////////////////////////////////

        // if(d_ans != ans ) {
        //     error_time[0]++;
            
        // }
       
        // typename std::make_signed<typename P::T>::type d_add;
        
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    
}

void mul_test( uint32_t plain_bits, int num_test){
    
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = TFHEpp::lvl1param;
    using privksP = TFHEpp::lvl11param;

    TFHEpp::SecretKey *sk = new TFHEpp::SecretKey();
    
    TFHEpp::relinKeyFFT<P> relinkeyfft =TFHEpp::relinKeyFFTgen<P>(sk->key.get<P>());
    
    TFHEpp::PrivateKeySwitchingKey<privksP> privksk;
    TFHEpp::privkskgen<privksP>(privksk, {1}, *sk);

    uint32_t scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;

    for (int test = 0; test < num_test; test++) {
        std::uniform_int_distribution<typename P::T> message(0, (1 << plain_bits) - 1);

        typename P::T p0, p1, pres, ptrue;
        // p0 = message(engine);
        // p1 = message(engine);

        p0 = 1;
        p1 = 0;
        ptrue = (p0 * p1) % P::plain_modulus;

        // TFHEpp::TLWE<P> c0 = TFHEpp::tlweSymIntEncrypt<P>(p0, P::α, sk->key.get<P>());
        // TFHEpp::TLWE<P> c1 = TFHEpp::tlweSymIntEncrypt<P>(p1, P::α, sk->key.get<P>());

        TFHEpp::TLWE<P> c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk->key.get<P>());
        TFHEpp::TLWE<P> c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk->key.get<P>());

        TFHEpp::TLWE<P> cres;
        TFHEpp::TLWEMult<privksP>(cres, c0, c1, relinkeyfft, privksk);
        //pres = TFHEpp::tlweSymIntDecrypt<P>(cres, sk->key.get<P>());
        pres = TFHEpp::tlweSymInt32Decrypt<Lvl1>(cres, pow(2., scale_bits), sk->key.get<P>());

        std::cout << "p0:"  << p0 << " p1:" << p1 << " pres:" << pres << " ptrue:" << ptrue<< std::endl;
        assert(pres == ptrue);
        
    }
    std::cout << "Passed" << std::endl;
}

void mul_and_test( uint32_t k, int num_test){

    using P = Lvl1;
    int32_t scale_bits = 31 - static_cast<int>(std::log2(k)) -1 ;
    
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplacebkfft<Lvl02>(sk);
    ek.emplaceiksk<Lvl20>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    ek.emplaceiksk<Lvl21>(sk);

    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<typename P::T> message(0, 1);

    std::vector<int> arr(k);
    int error_times = 0;
    std::vector<TLWELvl1> arr_ciphers(k);
    TLWELvl1 res;
    int ans = 1;
    int d_res = 0;
    double original_and_time = 0;
    double my_and_time = 0;
    std::chrono::system_clock::time_point start, end;

    TLWE<lvl1param> k_1 = TFHEpp::tlweSymInt32Encrypt<P>( (k-1) , P::α, pow(2., scale_bits), sk.key.get<P>());
    for (int test = 0; test < num_test; test++) {
        ans = 1;
        for(int i =0 ; i < k; i++){
            
            arr[i] = message(engine);
            ans = ans * arr[i];
            arr_ciphers[i] = TFHEpp::tlweSymInt32Encrypt<P>(arr[i], P::α, pow(2., scale_bits), sk.key.get<P>());
        }

        //Multi_HomAND(arr_ciphers ,res , k_1 ,ek ,29);

        start = std::chrono::system_clock::now();
        std::size_t num = arr_ciphers.size();
        TLWE<lvl1param> sum  = arr_ciphers[0];
        TLWE<lvl1param> sub;
        for(int32_t i = 1 ; i < num ; i++)
        {
            for (size_t j = 0; j <= lvl1param :: n; j++)
            {
                sum[j] = sum[j] + arr_ciphers[i][j];
            }
        }
        
        sub = sum;
        
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub[i] = k_1[i] - sum[i];
        }
        
        my_MSBGateBootstrapping(res, sub, ek, ARITHMETIC, 30 - 29 +1);

        end = std::chrono::system_clock::now();
        my_and_time = my_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        d_res = TFHEpp::tlweSymInt32Decrypt<Lvl1>(res, pow(2., 29), sk.key.get<P>());

        if(ans != d_res){
            error_times ++;
            for(int i =0 ; i< k ; i++){
                std::cout << "arr[" << i << "]:" << arr[i] << std::endl; 
            }
            std::cout << "ans:" << ans << std::endl;
            std::cout << "d_res:" << d_res << std::endl;

            std::cout << "----------------d_sum-----------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sum, pow(2., scale_bits), sk.key.get<P>());

            std::cout << "----------------k_1-----------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(k_1, pow(2., scale_bits), sk.key.get<P>());

            std::cout << "----------------d_sub-----------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sub, pow(2., scale_bits), sk.key.get<P>());
        }

        start = std::chrono::system_clock::now();
        TLWELvl1 temp;
        HomAND(arr_ciphers[0], arr_ciphers[1], temp, ek, ARITHMETIC);
        for(int i = 2 ; i< k ; i++){
            HomAND(temp, arr_ciphers[i], temp, ek, ARITHMETIC);
        }
        end = std::chrono::system_clock::now();
        original_and_time = original_and_time + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }
    std::cout << "error_times: " << error_times <<  std::endl;
    std::cout << "original_and_time: " << original_and_time / 1000.0 << "s" <<  std::endl;
    std::cout << "my_and_time: " << my_and_time / 1000.0 << "s" <<  std::endl;
}


void my_comparation_8bit_all_test(uint32_t plain_bits)
{
    
    std::cout << "------ Test of " << plain_bits << "bits comparation tttttt------" << std::endl;

    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    using P = Lvl1;
    TFHESecretKey sk;
    TFHEEvalKey ek;
    ek.emplacebkfft<Lvl01>(sk);
    ek.emplaceiksk<Lvl10>(sk);
    uint32_t scale_bits;
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ;
    double comparison_time_modified =0; 

    uint32_t num_test = 1ULL << plain_bits;
    std::cout<<"num_test: "<< num_test << std::endl;

    for (int test = 0; test < num_test; test++) 
    {   
        // p0 = message(engine);
        // p1 = message(engine) ;
        // p2 = message(engine);

        p1 = 0;
        p2 = test;
        
        int ans = 0;
        if(p2 > p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        //p0 = TFHEpp::tlweSymInt32Decrypt_print<P>(c0, pow(2., scale_bits), sk.key.lvl1);
        //std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;


        // Add
        ////////////////////////////////////
        TLWELvl1 c_boot;
        TLWELvl1 c_sub;
        TLWELvl1 c_sub2;
        TLWELvl1 c_ans;

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }

        // std::cout<<"-------sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub, pow(2., 27), sk.key.get<Lvl1>());

        //std::cout<<"-------my_comparation-------:"<< std::endl;
        

        ////////////////////my_ExtractMSB9

        // TLWELvl1 shift_tlwe, sign_tlwe5, res;
        // uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits-1;
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     shift_tlwe[i] = c_sub[i] << (plain_bits+1 - 5);
        // }
        // std::cout<<"-------shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe, pow(2., 31), sk.key.get<Lvl1>());

        // std::cout<<"-------1-----boot---------:"<< std::endl;
        // my_MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC, plain_bits+1 - 5);
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe5, pow(2., 31), sk.key.get<Lvl1>());
        
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     res[i] = c_sub[i] - sign_tlwe5[i];
        // }
        // std::cout<<"-------1-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());
        // my_MSBGateBootstrapping(res, res, ek, LOGIC, plain_bits+1 - 5);
        // std::cout<<"-------boot -2---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());

        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
        start = std::chrono::system_clock::now();
        TLWELvl1 greater_tlwe, less_tlwe, c_add6, c_sub6;
        // my_greater_than_equal<P>(c1, c2, greater_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // my_less_than_equal<P>(c1, c2, less_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // //my_HomAND(res, greater_tlwe, less_tlwe, c_2,  ek, ARITHMETIC, k);
        
        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_add6[i] = greater_tlwe[i] + less_tlwe[i];

        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_sub6[i] = c_1[i] - c_add6[i];

        // my_MSBGateBootstrapping(c_ans,c_sub6,ek,ARITHMETIC,3);


        //my_not_equal<Lvl1>(c2,c1,c_ans, c_1, plain_bits, ek, ARITHMETIC, 28);
        my_greater_than<Lvl1>(c2,c1,c_ans, plain_bits, ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        //std::cout<<"d_ans :  "<< d_ans << std::endl;

        if(d_ans != ans) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans: "<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "p2: "<< p2 <<std::endl;
            std::cout<<"----------------de_greater--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(greater_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------de_less--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(less_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------add--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_add6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------sub--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------cans--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}

void my_comparation_16bit_all_test(uint32_t plain_bits)
{
    
    std::cout << "------ Test of" << plain_bits<< " bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ; 
    double comparison_time_modified =0 ;

    uint32_t num_test = 1ULL << plain_bits ;
    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        p1 = test;
        p2 = 0;
        
        // int ans_great_than = 0;
        // int ans_great_than_equal = 0;
        // int ans_less_than = 0;
        // int ans_less_than_equal =0;
        // if(p2 > p1) {
        //     ans_great_than = 1;
        // }
        // if(p2 >= p1) {
        //     ans_great_than_equal = 1;
        // }
        // if(p2 < p1) {
        //     ans_less_than = 1;
        // }
        // if(p2 <= p1) {
        //     ans_less_than_equal = 1;
        // }

        int ans = 0;
        if(p2 > p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
        
        
        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;
        TFHEpp::TLWE<P> shift_tlwe19, sign_tlwe19, res19,  shift_tlwe14, sign_tlwe14, res14;
        TFHEpp::TLWE<Lvl1>  res14_switch ,shift_tlwe9, sign_tlwe9 , c_sub9, res9;
        
    /*
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }

        std::cout<<"-------sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

        std::cout<<"-------my_comparation-------:"<< std::endl;
    
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            shift_tlwe19[i] = c_sub[i] << (plain_bits+1 - 6);
        }
        std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
        std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

        
        my_MSBGateBootstrapping(sign_tlwe19, shift_tlwe19, ek, ARITHMETIC, plain_bits+1 - 6);
        std::cout<<"-------19-----boot---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
        
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res19[i] = c_sub[i] - sign_tlwe19[i];
        }
        std::cout<<"-------19-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

        std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
             shift_tlwe14[i] = res19[i] << (plain_bits+1-5 - 6);
        }
        std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

        my_MSBGateBootstrapping(sign_tlwe14, shift_tlwe14, ek, ARITHMETIC, plain_bits+1-5 - 6);
        std::cout<<"-------14-----boot---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());
        
    
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res14[i] = res19[i] - sign_tlwe14[i];
        }

        std::cout<<"------14-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

        
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res14_switch, res14, *ek.iksklvl21);
        std::cout<<"-------key-switch---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

        std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;
        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            shift_tlwe9[i] = res14_switch[i] << (plain_bits+1-5-5 -5);
        }
        std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(sign_tlwe9, shift_tlwe9, ek, ARITHMETIC, (plain_bits+1-5-5 -5));
        std::cout<<"-------9---boot----------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub9[i] = res14_switch[i] - sign_tlwe9[i];
        }
        std::cout<<"------9-----sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(res9, c_sub9, ek, ARITHMETIC, 0);
        std::cout<<"------9-----boot2---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
    */

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, P::α, pow(2., 28), sk.key.get<Lvl1>());
        start = std::chrono::system_clock::now();
        my_greater_than<Lvl2>(c2,c1,c_ans, plain_bits, ek, ARITHMETIC, 28);
        //my_not_equal<Lvl2>(c2,c1,c_ans,c_1,plain_bits, ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        
        // std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // std::cout<<"//////////////////////////////////////////////"<< std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
            std::cout<<"-------sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------19-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
            
            std::cout<<"-------19-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

            std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

            std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------14-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"------14-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

            std::cout<<"-------key-switch---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;

            std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-------9---boot----------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----boot2---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}

void bwtween_and_test_8bit(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of " << plain_bits << "bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ;
    double comparison_time_modified =0; 

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        // p1 = 0;
        // p2 = 0;
        
        int ans = 0;
        if(p1 >= p0 && p1 <= p2) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        //p0 = TFHEpp::tlweSymInt32Decrypt_print<P>(c0, pow(2., scale_bits), sk.key.lvl1);
        //std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;

        // Add
        ////////////////////////////////////
        TLWELvl1 c_boot;
        TLWELvl1 c_sub;
        TLWELvl1 c_sub2;

        TLWELvl1 c_ans;
        TLWELvl1 c_ans1;
        TLWELvl1 c_ans2;

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }
  
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than_equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
        
        TLWELvl1 greater_tlwe, less_tlwe, c_add6, c_sub6;
        start = std::chrono::system_clock::now();
        


        //my_not_equal<Lvl1>(c2,c1,c_ans, c_1, plain_bits, ek, ARITHMETIC, 28);
        my_greater_than_equal<Lvl1>(c1,c0,c_ans1, plain_bits, ek, ARITHMETIC, 28);
        my_less_than_equal<Lvl1>(c1,c2,c_ans2, plain_bits, ek, ARITHMETIC, 28);
        my_HomAND(c_ans, c_ans1, c_ans2, c_1,  ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        //std::cout<<"d_ans :  "<< d_ans << std::endl;

        if(d_ans != ans) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans: "<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "p2: "<< p2 <<std::endl;
            std::cout<<"----------------de_greater--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(greater_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------de_less--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(less_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------add--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_add6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------sub--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------cans--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}

void bwtween_and_test_16bit(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of" << plain_bits<< " bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ; 
    double comparison_time_modified =0 ;

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        int ans = 0;
        if(p0 <= p1 && p1 <= p2) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;
        TLWELvl1 c_ans1;
        TLWELvl1 c_ans2;

        TFHEpp::TLWE<P> shift_tlwe19, sign_tlwe19, res19,  shift_tlwe14, sign_tlwe14, res14;
        TFHEpp::TLWE<Lvl1>  res14_switch ,shift_tlwe9, sign_tlwe9 , c_sub9, res9;
        
        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        start = std::chrono::system_clock::now();
        my_greater_than_equal<Lvl2>(c1,c0,c_ans1,plain_bits, ek, ARITHMETIC, 28);
        //my_not_equal<Lvl2>(c2,c1,c_ans,c_1,plain_bits, ek, ARITHMETIC, 28);

        my_less_than_equal<Lvl2>(c1,c2,c_ans2,plain_bits, ek, ARITHMETIC, 28);
        my_HomAND(c_ans, c_ans1, c_ans2, c_1,  ek, ARITHMETIC, 28);

        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        
        // std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // std::cout<<"//////////////////////////////////////////////"<< std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
            std::cout<<"-------sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------19-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
            
            std::cout<<"-------19-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

            std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

            std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------14-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"------14-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

            std::cout<<"-------key-switch---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;

            std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-------9---boot----------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----boot2---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}

void in_test_16bit(uint32_t plain_bits, int num_test, int k)
{
    
    std::cout << "------ Test of" << plain_bits<< " bits comparation tttttt------" << std::endl;
    std::cout << "Test Time : " << num_test << std::endl;
    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::cout << "k : " << k << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);
    std::vector<uint32_t> arr(k, 0);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    std::vector<TLWELvl2> arr_c(k);
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ; 
    double comparison_time_modified =0 ;

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        int ans = 0;

        for(int i = 0; i< k; i++){
            arr[i] = message(engine);
            if(p0 == arr[i]){
                ans = 1;
            }
        }
        

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        for(int i = 0; i< k; i++){
            arr_c[i] = TFHEpp::tlweSymInt32Encrypt<P>(arr[i], P::α, pow(2., scale_bits), sk.key.get<P>());
        }

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;
        TLWELvl1 c_ans1;
        TLWELvl1 c_ans2;
        std::vector<TLWELvl1> c_ans_part(k);

        TFHEpp::TLWE<P> shift_tlwe19, sign_tlwe19, res19,  shift_tlwe14, sign_tlwe14, res14;
        TFHEpp::TLWE<Lvl1>  res14_switch ,shift_tlwe9, sign_tlwe9 , c_sub9, res9;
        
        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<Lvl1>( 7-1 , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<Lvl1>( k%7 -1 , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<Lvl1>( 0  , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());

        start = std::chrono::system_clock::now();
        //my_not_equal<Lvl2>(c2,c1,c_ans,c_1,plain_bits, ek, ARITHMETIC, 28);
        for(int i = 0; i< k; i++){
            my_equal<Lvl2>(c0,arr_c[i],c_ans_part[i],c_1,plain_bits, ek, ARITHMETIC, 28);;
        }
        N_HomOR(c_ans_part ,c_ans , k_0, ek ,28);

        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        
        // std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // std::cout<<"//////////////////////////////////////////////"<< std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----my_ImExtractMSB19---------:"<< std::endl;
            std::cout<<"-------sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

            std::cout<<"-----19--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe19, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------19-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe19, pow(2., 61), sk.key.get<P>());
            
            std::cout<<"-------19-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res19, pow(2., 61), sk.key.get<Lvl2>());

            std::cout<<"-----my_ImExtractMSB14---------:"<< std::endl;

            std::cout<<"-----14--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"-------14-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe14, pow(2., 61), sk.key.get<P>());

            std::cout<<"------14-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res14, pow(2., 62), sk.key.get<Lvl2>());

            std::cout<<"-------key-switch---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res14_switch, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-----my_ImExtractMSB9---------:"<< std::endl;

            std::cout<<"-------9---shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-------9---boot----------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub9, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"------9-----boot2---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res9, pow(2., 31), sk.key.get<Lvl1>());
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}

void in_test_8bit(uint32_t plain_bits, int num_test, int k)
{
    
    std::cout << "------ Test of" << plain_bits<< " bits comparation tttttt------" << std::endl;
    std::cout << "Test Time : " << num_test << std::endl;
    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::cout << "k : " << k << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);
    std::vector<uint32_t> arr(k, 0);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    std::vector<TLWELvl1> arr_c(k);
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ; 
    double comparison_time_modified =0 ;

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        int ans = 0;

        for(int i = 0; i< k; i++){
            arr[i] = message(engine);
            if(p0 == arr[i]){
                ans = 1;
            }
        }
        

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        for(int i = 0; i< k; i++){
            arr_c[i] = TFHEpp::tlweSymInt32Encrypt<P>(arr[i], P::α, pow(2., scale_bits), sk.key.get<P>());
        }

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;
        TLWELvl1 c_ans1;
        TLWELvl1 c_ans2;
        std::vector<TLWELvl1> c_ans_part(k);

        TFHEpp::TLWE<P> shift_tlwe19, sign_tlwe19, res19,  shift_tlwe14, sign_tlwe14, res14;
        TFHEpp::TLWE<Lvl1>  res14_switch ,shift_tlwe9, sign_tlwe9 , c_sub9, res9;
        
        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<Lvl1>( 7-1 , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<Lvl1>( k%7 -1 , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<Lvl1>( 0  , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());

        start = std::chrono::system_clock::now();
        //my_not_equal<Lvl2>(c2,c1,c_ans,c_1,plain_bits, ek, ARITHMETIC, 28);
        for(int i = 0; i< k; i++){
            my_equal<P>(c0,arr_c[i],c_ans_part[i],c_1,plain_bits, ek, ARITHMETIC, 28);
        }
        N_HomOR(c_ans_part ,c_ans , k_0, ek ,28);

        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        
        // std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // std::cout<<"//////////////////////////////////////////////"<< std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "   d_ans: "<< d_ans <<std::endl;
            std::cout<<"p0: "<< p0  << std::endl;

            std::cout<<"arr: " << std::endl;
            for(int i = 0; i< k ; i++){
                std::cout<< arr[i] << ", " ;
            }
            std::cout<<" "<< std::endl;

            std::cout<<"arr[i]== p0 : " << std::endl;
            for(int i = 0; i< k ; i++){
                d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans_part[i], pow(2., 28), sk.key.get<Lvl1>());
                std::cout<< d_ans << ", " ;
            }
            std::cout<<" "<< std::endl;
            
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}


void not_in_test_8bit(uint32_t plain_bits, int num_test, int k)
{
    
    std::cout << "------ Test of" << plain_bits<< " bits comparation tttttt------" << std::endl;
    std::cout << "Test Time : " << num_test << std::endl;
    std::cout << "Plain bits : " << plain_bits << std::endl;
    std::cout << "k : " << k << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);
    std::vector<uint32_t> arr(k, 0);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits ) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    std::vector<TLWELvl1> arr_c(k);
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ; 
    double comparison_time_modified =0 ;

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        p1 = message(engine) ;
        p2 = message(engine);

        int ans = 1;

        for(int i = 0; i< k; i++){
            arr[i] = message(engine);
            if(p0 == arr[i]){
                ans = 0;
            }
        }
        

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        for(int i = 0; i< k; i++){
            arr_c[i] = TFHEpp::tlweSymInt32Encrypt<P>(arr[i], P::α, pow(2., scale_bits), sk.key.get<P>());
        }

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;
        TLWELvl1 c_ans1;
        TLWELvl1 c_ans2;
        std::vector<TLWELvl1> c_ans_part(k);

        TFHEpp::TLWE<P> shift_tlwe19, sign_tlwe19, res19,  shift_tlwe14, sign_tlwe14, res14;
        TFHEpp::TLWE<Lvl1>  res14_switch ,shift_tlwe9, sign_tlwe9 , c_sub9, res9;
        
        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        greater_than<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<Lvl1>(1, Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_7 = TFHEpp::tlweSymInt32Encrypt<Lvl1>( 7-1 , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_rest = TFHEpp::tlweSymInt32Encrypt<Lvl1>( k%7 -1 , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_0 = TFHEpp::tlweSymInt32Encrypt<Lvl1>( 0  , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());
        TLWE<lvl1param> k_times = TFHEpp::tlweSymInt32Encrypt<Lvl1>( k/8  , Lvl1::α, pow(2., 28), sk.key.get<Lvl1>());

        start = std::chrono::system_clock::now();
        //my_not_equal<Lvl2>(c2,c1,c_ans,c_1,plain_bits, ek, ARITHMETIC, 28);
        for(int i = 0; i< k; i++){
            my_not_equal<P>(c0,arr_c[i],c_ans_part[i],c_1,plain_bits, ek, ARITHMETIC, 28);
        }
        N_HomAND(c_ans_part ,c_ans ,k_rest, k_7, k_times, ek ,28);

        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        
        // std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // std::cout<<"//////////////////////////////////////////////"<< std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "   d_ans: "<< d_ans <<std::endl;
            std::cout<<"p0: "<< p0  << std::endl;

            std::cout<<"arr: " << std::endl;
            for(int i = 0; i< k ; i++){
                std::cout<< arr[i] << ", " ;
            }
            std::cout<<" "<< std::endl;

            std::cout<<"arr[i]== p0 : " << std::endl;
            for(int i = 0; i< k ; i++){
                d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans_part[i], pow(2., 28), sk.key.get<Lvl1>());
                std::cout<< d_ans << ", " ;
            }
            std::cout<<" "<< std::endl;
            
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}


void my_comparation_14bit_negative256_test(uint32_t plain_bits, int num_test)
{
    std::cout << "------ Test of add Function tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    //std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original , comparison_time_modified; 

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        // p1 = message(engine) ;
        // p2 = message(engine);

        //p1 = (1 << 11)-1;
        p1 = 0;
        p2 = 256;
        
        int ans = 0;
        if(p2 > p1) {
            ans = 1;
        }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        // std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());


        
        ////////////////////////////////////
        TFHEpp::TLWE<P> c_boot;
        TFHEpp::TLWE<P> c_sub;
        
        TLWELvl1 c_ans;

        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            c_sub[i] = c1[i] - c2[i];
        }

        std::cout<<"-------sub---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());

        //std::cout<<"-------my_comparation-------:"<< std::endl;
        

        TFHEpp::TLWE<P> shift_tlwe, sign_tlwe5, res;
        TFHEpp::TLWE<Lvl1> shift_tlwe2, sign_tlwe1, res2 ,c_sub2, sign_tlwe3;
        uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits-1;
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            shift_tlwe[i] = c_sub[i] << (plain_bits+1 - 6);
        }
        // std::cout<<"-----1--shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe, pow(2., 61), sk.key.get<P>());

        //std::cout<<"-------1-----boot---------:"<< std::endl;
        my_MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC, plain_bits+1 - 6);
        //TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe5, pow(2., 61), sk.key.get<P>());
        
        for (size_t i = 0; i <= Lvl2 :: n; i++)
        {
            res[i] = c_sub[i] - sign_tlwe5[i];
        }
        // std::cout<<"-------1-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res, pow(2., 61), sk.key.get<Lvl2>());

        //std::cout<<"-------key-switch---------:"<< std::endl;
        TFHEpp::IdentityKeySwitch<TFHEpp::lvl21param>(res2, res, *ek.iksklvl21);
        //TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            shift_tlwe2[i] = res2[i] << 2;
        }
        // std::cout<<"-------2---shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe2, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(sign_tlwe1, shift_tlwe2, ek, ARITHMETIC, 2);
        // std::cout<<"-------boot -2---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe1, pow(2., 31), sk.key.get<Lvl1>());

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub2[i] = res2[i] - sign_tlwe1[i];
        }
        // std::cout<<"------2-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub2, pow(2., 31), sk.key.get<Lvl1>());

        my_MSBGateBootstrapping(sign_tlwe3, c_sub2, ek, ARITHMETIC, 0);
        // std::cout<<"-------boot -3---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe3, pow(2., 31), sk.key.get<Lvl1>());

        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        //greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        //改进后的比较 计时记录
        start = std::chrono::system_clock::now();
        my_greater_than<Lvl2>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC, 1);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); 

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 29), sk.key.get<Lvl1>());
        //std::cout<<"d_ans :  "<< d_ans << std::endl;

        if(d_ans != ans ) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans:"<< ans << "d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c1, pow(2., 61), sk.key.get<P>());
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c2, pow(2., 61), sk.key.get<P>());
            std::cout<<"-------sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(c_sub, pow(2., 61), sk.key.get<P>());
            std::cout<<"-----1--shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(shift_tlwe, pow(2., 61), sk.key.get<P>());
            std::cout<<"-------1-----boot---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(sign_tlwe5, pow(2., 61), sk.key.get<P>());
            std::cout<<"-------1-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl2>(res, pow(2., 61), sk.key.get<Lvl2>());
            std::cout<<"-------key-switch---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"-------2---shift_tlwe---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe2, pow(2., 31), sk.key.get<Lvl1>());
            
            
            std::cout<<"-------boot -2--inner---------:"<< std::endl;
            shift_tlwe2[lvl1param::k * lvl1param::n] += 1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
            std::cout<<"-------inner--add--------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe2, pow(2., 31), sk.key.get<Lvl1>());
            TLWE<lvl0param> tlwelvl0;
            IdentityKeySwitch<lvl10param>(tlwelvl0, shift_tlwe2, *ek.iksklvl10);
            lvl1param::T μ = lvl1param::μ;
            μ = (μ << 1) >> 2;
            GateBootstrappingTLWE2TLWEFFT<lvl01param>(res2, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
            std::cout<<"-------inner--boot--------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());
            res2[lvl1param::k * lvl1param::n] += (μ);
            std::cout<<"-------inner--out--------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res2, pow(2., 31), sk.key.get<Lvl1>());

            std::cout<<"-------boot -2---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe1, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"------2-----sub---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub2, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"-------boot -3---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe3, pow(2., 31), sk.key.get<Lvl1>());
            std::cout<<"-------d_ans---------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans, pow(2., 29), sk.key.get<Lvl1>());
            
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time :"<< error_time[0]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/1000 << "s"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/1000  << "s"<< std::endl;
}

void my_comparation_8bit_negative256_test(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of " << plain_bits << "bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ;
    double comparison_time_modified =0; 

    for (int test = 0; test < num_test; test++) 
    {   
        p0 = message(engine);
        // p1 = message(engine) ;
        // p2 = message(engine);
        p0 = 256;
        p1 = 255;
        p2 = 0;
        
        int ans = 1;
        // if(p2 >= p1) {
        //     ans = 1;
        // }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits+1), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        //p0 = TFHEpp::tlweSymInt32Decrypt_print<P>(c0, pow(2., scale_bits), sk.key.lvl1);
        //std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;

        // Add
        ////////////////////////////////////
        TLWELvl1 c_boot;
        TLWELvl1 c_sub;
        TLWELvl1 c_sub2;
        TLWELvl1 c_ans;
        TLWELvl1 c_ans2;

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub[i] = c2[i] - c1[i] ;
        }
        c_sub[Lvl1 :: n] = c_sub[Lvl1 :: n] - (1 << 23);
        // std::cout<<"-------c0---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c0, pow(2., 27), sk.key.get<Lvl1>());
        
        my_ExtractMSB9(c_ans,c_sub, plain_bits , ek, ARITHMETIC, 28);
        ExtractMSB9(c_ans2, c_sub, plain_bits , ek, ARITHMETIC);

        // my_ExtractMSB9(c_ans, c_sub, plain_bits + 1 , ek, ARITHMETIC, 28);
        // ExtractMSB9(c_ans2, c_sub, plain_bits + 1, ek, ARITHMETIC);
        //std::cout<<"-------my_comparation-------:"<< std::endl;
        

        ////////////////////my_ExtractMSB9

        // TLWELvl1 shift_tlwe, sign_tlwe5, res;
        // uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits-1;
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     shift_tlwe[i] = c_sub[i] << (plain_bits+1 - 5);
        // }
        // std::cout<<"-------shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe, pow(2., 31), sk.key.get<Lvl1>());

        // std::cout<<"-------1-----boot---------:"<< std::endl;
        // my_MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC, plain_bits+1 - 5);
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe5, pow(2., 31), sk.key.get<Lvl1>());
        
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     res[i] = c_sub[i] - sign_tlwe5[i];
        // }
        // std::cout<<"-------1-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());
        // my_MSBGateBootstrapping(res, res, ek, LOGIC, plain_bits+1 - 5);
        // std::cout<<"-------boot -2---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());

        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        //greater_than_equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
        
        TLWELvl1 greater_tlwe, less_tlwe, c_add6, c_sub6;
        start = std::chrono::system_clock::now();
        // my_greater_than_equal<P>(c1, c2, greater_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // my_less_than_equal<P>(c1, c2, less_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // //my_HomAND(res, greater_tlwe, less_tlwe, c_2,  ek, ARITHMETIC, k);
        
        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_add6[i] = greater_tlwe[i] + less_tlwe[i];

        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_sub6[i] = c_1[i] - c_add6[i];

        // my_MSBGateBootstrapping(c_ans,c_sub6,ek,ARITHMETIC,3);


        //my_not_equal<Lvl1>(c2,c1,c_ans, c_1, plain_bits, ek, ARITHMETIC, 28);

        //my_greater_than_equal<Lvl1>(c2,c1,c_ans, plain_bits, ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        int d_ans2 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans2, pow(2., 31), sk.key.get<Lvl1>());
        //std::cout<<"d_ans :  "<< d_ans << std::endl;

        if(d_ans != ans) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans: "<< ans << " d_ans: "<< d_ans <<std::endl;
            std::cout<<"----------------de_ans--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
            
        }

        if(d_ans2 != ans) {
            error_time[1]++;
            std::cout<<"----------------de_ans2--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans2, pow(2., 31), sk.key.get<Lvl1>());
            
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time_0 :"<< error_time[0]<< std::endl;
    std::cout<<"error time_1 :"<< error_time[1]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
} 


void my_comparation_5bit_negative256_test(uint32_t plain_bits, int num_test)
{
    
    std::cout << "------ Test of " << plain_bits << "bits comparation tttttt------" << std::endl;
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
    
    
    std::vector<uint32_t> error_time(3, 0);
    std::vector<double> comparison_time(7, 0.);

    // For simplicity，the input range is [0, 2^(p-1) -1]
    std::uniform_int_distribution<typename P::T> message(0, (1 << (plain_bits - 1) - 1));
    std::uniform_int_distribution<typename P::T> type(0, 1);
    scale_bits = std::numeric_limits<P::T>::digits - plain_bits - 1;
    
    std::cout << "scale_bits : " << scale_bits << std::endl;
    
    typename P::T p0, p1, p2, p3, gres, geres, lres, leres, eres, dgres, dgeres, dlres, dleres, deres, ands, d_ands;
    TFHEpp::TLWE<P> c0, c1, c2, c3, c;

    TLWELvl1 cres,cres1, cres2, c_ands;
    
    std::chrono::system_clock::time_point start, end;
    double comparison_time_original =0 ;
    double comparison_time_modified =0; 

    for (int test = 0; test < num_test; test++) 
    {   
        // p0 = message(engine);
        // p1 = message(engine) ;
        // p2 = message(engine);

        p0 = 16;
        p1 = 255;
        p2 = 0;
        
        int ans = 1;
        // if(p2 >= p1) {
        //     ans = 1;
        // }

        bool result_type = ARITHMETIC;
       
        c0 = TFHEpp::tlweSymInt32Encrypt<P>(p0, P::α, pow(2., scale_bits+1), sk.key.get<P>());
        c1 = TFHEpp::tlweSymInt32Encrypt<P>(p1, P::α, pow(2., scale_bits), sk.key.get<P>());
        c2 = TFHEpp::tlweSymInt32Encrypt<P>(p2, P::α, pow(2., scale_bits), sk.key.get<P>());

        //p0 = TFHEpp::tlweSymInt32Decrypt_print<P>(c0, pow(2., scale_bits), sk.key.lvl1);
        //std::cout<<"p1: "<< p1 << "  p2:"<< p2 << std::endl;

        // Add
        ////////////////////////////////////
        TLWELvl1 c_boot;
        TLWELvl1 c_sub;
        TLWELvl1 c_sub2;
        TLWELvl1 c_ans;
        TLWELvl1 c_ans2;

        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_sub[i] = c2[i] - c1[i] ;
        }
        c_sub[Lvl1 :: n] = c_sub[Lvl1 :: n] - (1 << 23);
        std::cout<<"-------c0---------:"<< std::endl;
        TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c0, pow(2., 27), sk.key.get<Lvl1>());
        
        my_MSBGateBootstrapping(c_ans, c0, ek,  ARITHMETIC, 3);
        MSBGateBootstrapping(c_ans2, c0, ek, ARITHMETIC);

        //my_ExtractMSB9(c0, c_ans, plain_bits + 1 , ek, ARITHMETIC, 28);
        // ExtractMSB9(c0, c_ans2, plain_bits + 1, ek, ARITHMETIC);
        //std::cout<<"-------my_comparation-------:"<< std::endl;
        

        ////////////////////my_ExtractMSB9

        // TLWELvl1 shift_tlwe, sign_tlwe5, res;
        // uint32_t scale_bits = std::numeric_limits<Lvl1::T>::digits - plain_bits-1;
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     shift_tlwe[i] = c_sub[i] << (plain_bits+1 - 5);
        // }
        // std::cout<<"-------shift_tlwe---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(shift_tlwe, pow(2., 31), sk.key.get<Lvl1>());

        // std::cout<<"-------1-----boot---------:"<< std::endl;
        // my_MSBGateBootstrapping(sign_tlwe5, shift_tlwe, ek, ARITHMETIC, plain_bits+1 - 5);
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(sign_tlwe5, pow(2., 31), sk.key.get<Lvl1>());
        
        // for (size_t i = 0; i <= Lvl1 :: n; i++)
        // {
        //     res[i] = c_sub[i] - sign_tlwe5[i];
        // }
        // std::cout<<"-------1-----sub---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());
        // my_MSBGateBootstrapping(res, res, ek, LOGIC, plain_bits+1 - 5);
        // std::cout<<"-------boot -2---------:"<< std::endl;
        // TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(res, pow(2., 31), sk.key.get<Lvl1>());

        
        //原始比较 计时记录
        start = std::chrono::system_clock::now();
        //greater_than_equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //equal<Lvl1>(c2,c1,c_ans,plain_bits, ek, ARITHMETIC);
        //HomNOT(c_ans,c_ans);
        end = std::chrono::system_clock::now();
        comparison_time_original = comparison_time_original + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //改进后的比较 计时记录
        TLWELvl1 c_1 = TFHEpp::tlweSymInt32Encrypt<P>(1, P::α, pow(2., 28), sk.key.get<P>());
        
        TLWELvl1 greater_tlwe, less_tlwe, c_add6, c_sub6;
        start = std::chrono::system_clock::now();
        // my_greater_than_equal<P>(c1, c2, greater_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // my_less_than_equal<P>(c1, c2, less_tlwe, plain_bits, ek, ARITHMETIC, 28);
        // //my_HomAND(res, greater_tlwe, less_tlwe, c_2,  ek, ARITHMETIC, k);
        
        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_add6[i] = greater_tlwe[i] + less_tlwe[i];

        // for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
        //     c_sub6[i] = c_1[i] - c_add6[i];

        // my_MSBGateBootstrapping(c_ans,c_sub6,ek,ARITHMETIC,3);


        //my_not_equal<Lvl1>(c2,c1,c_ans, c_1, plain_bits, ek, ARITHMETIC, 28);

        //my_greater_than_equal<Lvl1>(c2,c1,c_ans, plain_bits, ek, ARITHMETIC, 28);
        end = std::chrono::system_clock::now();
        comparison_time_modified = comparison_time_modified + std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();  ;

        //std::cout<<"-------d_ans---------:"<< std::endl;
        int d_ans = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        int d_ans2 = TFHEpp::tlweSymInt32Decrypt<Lvl1>(c_ans2, pow(2., 31), sk.key.get<Lvl1>());
        //std::cout<<"d_ans :  "<< d_ans << std::endl;

        if(d_ans != ans) {
            error_time[0]++;
            std::cout<<"//////////////////////////////////////////////"<< std::endl;
            std::cout<<"ans: "<< ans << " d_ans: "<< d_ans <<std::endl;
            std::cout<<"p1: "<< p1 << " p2: "<< p2 <<std::endl;
            std::cout<<"----------------de_greater--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(greater_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------de_less--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(less_tlwe, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------add--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_add6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------sub--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_sub6, pow(2., 28), sk.key.get<Lvl1>());
            std::cout<<"----------------cans--------------------"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<Lvl1>(c_ans, pow(2., 28), sk.key.get<Lvl1>());
        }

        if(d_ans2 != ans) {
            error_time[1]++;
            
        }
       
        typename std::make_signed<typename P::T>::type d_add;
        
        
    
    }
    std::cout<<"error time_0 :"<< error_time[0]<< std::endl;
    std::cout<<"error time_1 :"<< error_time[1]<< std::endl;
    std::cout<<" comparison_time_original: "<< comparison_time_original/num_test << "ms"<< std::endl;
    std::cout<<" comparison_time_modified: "<< comparison_time_modified/num_test  << "ms"<< std::endl;
}




int main()
{
    std::chrono::system_clock::time_point start, end;
    int num_test = 100;
    start = std::chrono::system_clock::now();

    // my_comparation_test(8, 20000);
    my_comparation_19bit_test(16,100);
    //my_comparation_14bit_negative256_test(9,1);
    //my_comparation_8bit_negative256_test(9,10000);
    //my_comparation_5bit_negative256_test(5,10);

    end = std::chrono::system_clock::now();
    int total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout<<" total_time: "<< total_time/1000  << "s"<< std::endl; 

    
}


