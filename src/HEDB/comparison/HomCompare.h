#pragma once
#include "extract_msb.h"
#include "operators.h"

//tfhepp_utils.h 里面还要加my_MSBGateBootstrapping的声明
#include "tfhepp_utils.h"

/***
 * A Greater than B
 * A Greater than or equal to B
 * A less than B
 * A less than or equal to B
 * A equal to B
***/
namespace HEDB
{
    // cipher1 > cipher2 <=> msb(cipher2 - cipher1) == 1
    template <typename P>
    void greater_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        }
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }

    //*******************改的greater_than*******************///
    template <typename P>
    void my_greater_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type, uint32_t k)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        }
        //HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
        if(plain_bits <= 4){
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_MSBGateBootstrapping(res, sub_tlwe, ek, result_type, (31-k));  //my_MSBGateBootstrapping  in  src/comparison/tfhepp_utils.cpp
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if (plain_bits <= 8){
            
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_ExtractMSB9(res, sub_tlwe, plain_bits + 1 , ek, result_type, k);
            }
            else if constexpr (std::is_same_v<P, TFHEpp::lvl2param>){
                my_ImExtractMSB9(res, sub_tlwe, plain_bits + 1 , ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        else if(plain_bits <= 13){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB14
                my_ImExtractMSB14(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if(plain_bits <= 18){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB19(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }else if(plain_bits <= 23){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB24(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 28){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB29(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 33){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB34(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        

    }
    //////////////////////////////////////////////////////////////////
    
    // cipher1 >= cipher2 <=> NOT(msb(cipher1 - cipher2)) == 1
    template <typename P>
    void greater_than_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        }
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) TFHEpp::LOG_to_ARI(res, res, ek);
    }
    
    //测试使用 —— Rainyz ——
    template <typename P>
    void my_greater_than_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type, uint32_t k,TFHESecretKey &sk)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        }

        //
            // int decrypted = TFHEpp::tlweSymInt32Decrypt_print<P>(sub_tlwe, pow(2., 24), sk.key.get<Lvl1>() );
            // std::cout << "Decrypted sub_tlwe (cipher1 - cipher2) = " << decrypted << std::endl;
        //

        //HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        if(plain_bits <= 4){
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_MSBGateBootstrapping(res, sub_tlwe, ek, result_type, (31-k));  //my_MSBGateBootstrapping  in  src/comparison/tfhepp_utils.cpp
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if (plain_bits <= 8){    
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_ExtractMSB9(res, sub_tlwe, plain_bits + 1 , ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        else if(plain_bits <= 13){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                
                my_ImExtractMSB14(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if(plain_bits <= 18){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB19(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        } else if(plain_bits <= 23){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB24(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 28){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB29(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 33){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB34(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        
        HomNOT<Lvl1>(res, res);
        res[TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] += (1ULL << k); // @TODO: 这里的作用是什么？  Rainyz——2025.8.9
        //if (IS_ARITHMETIC(result_type)) TFHEpp::LOG_to_ARI(res, res, ek);
    }



    //*******************改动的greater_than_equal*******************//
    template <typename P>
    void my_greater_than_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type, uint32_t k)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        }


        //HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        if(plain_bits <= 4){
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_MSBGateBootstrapping(res, sub_tlwe, ek, result_type, (31-k));  //my_MSBGateBootstrapping  in  src/comparison/tfhepp_utils.cpp
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if (plain_bits <= 8){    
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_ExtractMSB9(res, sub_tlwe, plain_bits + 1 , ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        else if(plain_bits <= 13){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                
                my_ImExtractMSB14(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if(plain_bits <= 18){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB19(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }else if(plain_bits <= 23){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB24(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 28){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB29(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 33){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB34(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        
        HomNOT<Lvl1>(res, res);
        res[TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] += (1ULL << k);  //@TODO 这一行不知道作用是什么。 --RAINYZ--
        //if (IS_ARITHMETIC(result_type)) TFHEpp::LOG_to_ARI(res, res, ek); //这一行本来就是注释的，不变
    }
    ///////////////////////////////////////////////////////////////
    

    ///*******************改的less_than*******************///
    template <typename P>
    void my_less_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type, uint32_t k)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        }
        //HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);

        if(plain_bits <= 4){
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_MSBGateBootstrapping(res, sub_tlwe, ek, result_type, (31-k));  //my_MSBGateBootstrapping  in  src/comparison/tfhepp_utils.cpp
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if (plain_bits <= 8){
            
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_ExtractMSB9(res, sub_tlwe, plain_bits + 1 , ek, result_type, k);
            }
            else if constexpr (std::is_same_v<P, TFHEpp::lvl2param>){
                my_ImExtractMSB9(res, sub_tlwe, plain_bits + 1 , ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        else if(plain_bits <= 13){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB14
                my_ImExtractMSB14(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if(plain_bits <= 18){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB19(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }else if(plain_bits <= 23){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB24(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 28){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB29(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }else if(plain_bits <= 33){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB34(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////
    
    

    // cipher1 < cipher2 <=> msb(cipher1 - cipher2) == 1
    template <typename P>
    void less_than(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher1[i] - cipher2[i];
        }
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
    }
    
    //********************改的less_than_equal*****************//
    template <typename P>
    void my_less_than_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type, uint32_t k)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        }
        //HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        
        if(plain_bits <= 4){
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_MSBGateBootstrapping(res, sub_tlwe, ek, result_type, (31-k));  //my_MSBGateBootstrapping  in  src/comparison/tfhepp_utils.cpp
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if (plain_bits <= 8){
            
            if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
                
                my_ExtractMSB9(res, sub_tlwe, plain_bits + 1 , ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
        }
        else if(plain_bits <= 13){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB14
                my_ImExtractMSB14(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        }
        else if(plain_bits <= 18){
            if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
                // 如果 P 是 TFHEpp::lvl2param，可以直接调用 my_ImExtractMSB19
                my_ImExtractMSB19(res, sub_tlwe, plain_bits +1, ek, result_type, k);
            }
            else{
                HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
            }
            
        } 
        else{
            
            HomMSB(res, sub_tlwe, plain_bits + 1, ek, result_type);
        }
        
        HomNOT<Lvl1>(res, res);
        res[TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] += (1ULL << k);
        //if (IS_ARITHMETIC(result_type)) TFHEpp::LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void my_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, TLWELvl1 &c_2, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type , uint32_t k)
    {
        TLWELvl1 greater_tlwe, less_tlwe;
        my_greater_than_equal<P>(cipher1, cipher2, greater_tlwe, plain_bits, ek, result_type, k);
        my_less_than_equal<P>(cipher1, cipher2, less_tlwe, plain_bits, ek, result_type, k);
        my_HomAND(res, greater_tlwe, less_tlwe, c_2,  ek, result_type, k);
    }

    template <typename P>
    void my_not_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res,TLWELvl1 &c_2, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type , uint32_t k)
    {
        TLWELvl1 greater_tlwe, less_tlwe;
        my_equal<P>(cipher1, cipher2, res, c_2, plain_bits, ek, result_type, k);
        HomNOT<Lvl1>(res, res);
        res[TFHEpp::lvl1param::k * TFHEpp::lvl1param::n] += (1ULL << k);
    }
    ////////////////////////////////////////////////////////////////////////
    
    // cipher1 <= cipher2 <=> NOT(msb(cipher2 - cipher1)) == 1
    template <typename P>
    void less_than_equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type)
    {
        TFHEpp::TLWE<P> sub_tlwe; 
        for (size_t i = 0; i <= P::k * P::n; i++)
        {
            sub_tlwe[i] = cipher2[i] - cipher1[i];
        }
        HomMSB(res, sub_tlwe, plain_bits + 1, ek, LOGIC);
        HomNOT<Lvl1>(res, res);
        if (IS_ARITHMETIC(result_type)) TFHEpp::LOG_to_ARI(res, res, ek);
    }

    template <typename P>
    void equal(TFHEpp::TLWE<P> &cipher1, TFHEpp::TLWE<P> &cipher2, TLWELvl1 &res, uint32_t plain_bits, TFHEEvalKey &ek, bool result_type)
    {
        TLWELvl1 greater_tlwe, less_tlwe;
        greater_than_equal<P>(cipher1, cipher2, greater_tlwe, plain_bits, ek, LOGIC);
        less_than_equal<P>(cipher1, cipher2, less_tlwe, plain_bits, ek, LOGIC);
        HomAND(res, greater_tlwe, less_tlwe, ek, result_type);
    }
} // namespace HEDB