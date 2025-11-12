#include <iostream>
#include <chrono>
#include <random>
#include<bitset>
#include "src/HEDB/comparison/comparison.h"
#include "src/HEDB/utils/utils.h"

#include "src/HEDB/comparison/tfhepp_utils.h"
#include <gatebootstrapping.hpp>
#include "detwfa.hpp"

#include <cassert>
#include <tfhe++.hpp>

using namespace std;
using namespace HEDB;
using namespace TFHEpp;


void add_replace_7and_test(uint32_t plain_bits, int num_test)
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
    std::vector<uint32_t> error_time_list(num +3, 0);
    std::vector<double> times(2, 0.);
    
    std::vector<TLWELvl1> c0_list(num);
    std::vector<TLWELvl1> c1_list(num);
    std::vector<TLWELvl1> cres_list(num);
    
    uint32_t u_right_move_bit = 2;
    bool result_type = LOGIC;
    std::chrono::system_clock::time_point start, end;
    
    for (int test_time = 0; test_time < num_test; test_time++) 
    {
        int32_t add = 0;
        for(int i = 0; i < num; i++)
        {
            p0_list[i] = message(engine);
            p1_list[i] = message(engine);
            
            if(p0_list[i] >= p1_list[i]) gres_list[i] = 1;
            else gres_list[i] = 0;
            
            c0_list[i] = TFHEpp::tlweSymInt32Encrypt<P>(p0_list[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            c1_list[i] = TFHEpp::tlweSymInt32Encrypt<P>(p1_list[i], P::α, pow(2., scale_bits), sk.key.get<P>());
            
            //std::cout<< "gres["<< i << "]: "<< gres_list[i]<<"   p0["<< i << "]:"<< p0_list[i] << "   p1["<< i << "]:"<< p1_list[i]  << std::endl;
            
            my_greater_than_equal<P>(c0_list[i], c1_list[i], cres_list[i], plain_bits, ek, ARITHMETIC, 28);
            
            //std::cout<<"Greater than-----Decrypt--------------------:"<< std::endl;
            dgres_list[i] = TFHEpp:: tlweSymInt32Decrypt<P>(cres_list[i], pow(2., 28), sk.key.lvl1);
        
            
            //std::cout<<"dgres["<< i << "]: "<< dgres_list[i] << std::endl;            
            if (gres_list[i] != dgres_list[i]){
                error_time_list[i] += 1;
            }
            
            //cres_list[i][Lvl1 :: n] +=  P::μ >> u_right_move_bit;  
            

        }
        
        
        TLWELvl1 c_add;
        //typename std::make_signed<typename P::T>::type d_add;
        P::T d_add;
        
        for (size_t i = 0; i < num; i++)
        {
            if (gres_list[i] < 1) add -= 0;
            else add +=1;
        }
        
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i <= Lvl1 :: n; i++)
        {
            c_add[i] = 0;
            for (size_t j = 0; j < num; j++)
            {
                c_add[i] += cres_list[j][i];
            }
            
        }
        end = std::chrono::system_clock::now();
        times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
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
        //uint32_t add_plain_bits = std::numeric_limits<typename P::T>::digits - (29 - u_right_move_bit);
        P::T add_dres;
        
        //uint32_t  threshold = (1ULL << (add_plain_bits-1)) + filter -1;
        uint32_t  threshold = filter-1;
    
        TLWELvl1 c_threshold = TFHEpp::tlweSymInt32Encrypt<P>(threshold, P::α, pow(2., 29 - u_right_move_bit +1), sk.key.get<P>());
        
        //c_add[Lvl1 :: n] += 1ULL << (std::numeric_limits<typename P::T>::digits -1);
        
        // std::cout<<"--------c_add += 2^31:---------- "<<std::endl;
        // TFHEpp::my_tlweSymDecrypt<P>(c_add, sk.key.lvl1);
        
        start = std::chrono::system_clock::now();
        TLWELvl1 c_subs;
        for (size_t i = 0; i <= Lvl1 :: n; i++)
            {
                c_subs[i] = c_threshold[i] - c_add[i];
            }
        // std::cout<<"c_threshold[i] - c_add[i]+++++++++ Decrypt +++++++++++++++++:"<<std::endl;
        // TFHEpp::tlweSymInt32Decrypt<P>(c_subs, pow(2., scale_bits), sk.key.lvl1);
        
        my_greater_than<P>(c_add, c_threshold, add_cres, plain_bits, ek, ARITHMETIC,31);
        TLWELvl1 add_cres1 = add_cres;
        
        //GateBootstrapping(add_cres, add_cres, ek);

        end = std::chrono::system_clock::now();
        times[0] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        //std::cout<<"Threshold Greater than-----Decrypt--------------------:"<< std::endl;
        add_dres = TFHEpp::tlweSymInt32Decrypt<P>(add_cres, pow(2., 31), sk.key.lvl1);
        
        if(!(add == 7 && add_dres == 1 || add < 7 && add_dres == 0))
        {
            std::cout<<"--------c_add += 2^31:---------- "<<std::endl;
            TFHEpp::tlweSymDecrypt_print<P>(c_add, sk.key.lvl1);
            std::cout<<"c_threshold[i] - c_add[i]+++++++++ Decrypt +++++++++++++++++:"<<std::endl;
            TFHEpp::tlweSymDecrypt_print<P>(c_subs, sk.key.lvl1);
            //TFHEpp::tlweSymInt32Decrypt<P>(c_subs, pow(2., scale_bits), sk.key.lvl1);
            
            std::cout<<"Threshold Greater than-----Decrypt--------------------:"<< std::endl;
            TFHEpp::tlweSymInt32Decrypt_print<P>(add_cres1, pow(2., 31), sk.key.lvl1);
            std::cout<<"d_add: "<< d_add <<std::endl;
            std::cout<<"add_dres: "<< add_dres <<std::endl;
            error_time_list[num +1] += 1;
        }

        //直接用普通的and运算
        TLWELvl1 and_cres;
        P::T and_dres;

        and_dres = TFHEpp::tlweSymDecrypt<P>(and_cres, sk.key.lvl1);
        start = std::chrono::system_clock::now();

        HomAND(and_cres, cres_list[0], cres_list[1], ek, LOGIC);
        and_dres = TFHEpp::tlweSymDecrypt<P>(and_cres, sk.key.lvl1);
        //std::cout<< "and_dres[1]:" << and_dres <<std::endl;

        for(int32_t i =2 ;i< num ;i++){
            HomAND(and_cres, cres_list[i], and_cres, ek, LOGIC);
            and_dres = TFHEpp::tlweSymDecrypt<P>(and_cres, sk.key.lvl1);
            //std::cout<< "and_dres[" << i << "]:"<< and_dres <<std::endl;

        }
        end = std::chrono::system_clock::now();
        times[1] += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
    }

    for (size_t i = 0; i < num +3; i++)
    {
        std::cout << "Error time" << i << " : " << error_time_list[i] << std::endl;
    }
    std::cout << "++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << "7 times and time" << " : " << times[1]/num_test << std::endl;
    std::cout << "add_replace_and time" << " : " << times[0]/num_test << std::endl;

}

    int main()
{
    int num_test = 100;
    add_replace_7and_test(4, num_test);
}