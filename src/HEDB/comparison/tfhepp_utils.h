#pragma once
#include <cmath>
#include <limits>
#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"
#include "gatebootstrapping.hpp"

#include "circuitbootstrapping.hpp"

#include<bitset>
using namespace std;

namespace TFHEpp
{
    // c = (a, a * s + m + e)
    template <class P>
    TLWE<P> tlweSymInt32Encrypt(const typename P::T p, const double α, const double scale, const Key<P> &key)
    {
        std::uniform_int_distribution<typename P::T> Torusdist(0, std::numeric_limits<typename P::T>::max());
        TLWE<P> res = {};
        res[P::k * P::n] =
            ModularGaussian<P>(static_cast<typename P::T>(p * scale), α);
        // res[P::k * P::n] = static_cast<typename P::T>(p * scale);
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++) {
                res[k * P::n + i] = Torusdist(generator);
                // res[k * P::n + i] = 1;
                res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
            }
        return res;
    }

////////////////////////////////////////////////////////////////////////////自己加的

    template <class P>
    pair<TLWE<P> , uint32_t> new_tlweSymInt32Encrypt(const typename P::T p, const double α, const double scale, const Key<P> &key)
    {
        std::uniform_int_distribution<typename P::T> Torusdist(0, std::numeric_limits<typename P::T>::max());
        TLWE<P> res = {};
        uint64_t temp = 0;
        uint32_t extra_first_bit ;
        
        res[P::k * P::n] = ModularGaussian<P>(static_cast<typename P::T>(p * scale), α);
        temp = res[P::k * P::n];
        
        std::cout<<"res[n]:"<< bitset<32>(res[P::k * P::n ]) << std::endl;
        std::cout<<"temp:"<< bitset<64>(temp) << std::endl;
        // res[P::k * P::n] = static_cast<typename P::T>(p * scale);
        for (int k = 0; k < P::k; k++){
            for (int i = 0; i < P::n; i++) {
                res[k * P::n + i] = Torusdist(generator);
                
                // res[k * P::n + i] = 1;
                res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
                temp += res[k * P::n + i] * key[k * P::n + i];
            }
            
            std::cout<<"res[n]:0000000000000000000000000000000@"<< bitset<32>(res[P::k * P::n ]) << std::endl;
            std::cout<<"temp  :"<< bitset<64>(temp) << std::endl;
            
            extra_first_bit = ((temp - res[P::k * P::n]) >> 32)%2;
            std::cout<<"extra_first_bit:"<< bitset<32>(extra_first_bit) << std::endl;
        }
        return make_pair(res, extra_first_bit);
    }
    
    template <class P>
    pair<typename P::T , uint64_t> new_tlweSymInt32Decrypt(const TLWE<P> &c, const double scale, const Key<P> &key)
    {
        uint64_t temp = 0 ;
        
        typename P::T phase = c[P::k * P::n];
        std::cout<<"c[n]:"<< phase << std::endl;
        
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++){
            for (int i = 0; i < P::n; i++){
                phase -= c[k * P::n + i] * key[k * P::n + i];
                
                temp += c[k * P::n + i] * key[k * P::n + i];
            }
        }
        ////////////////////////////////////////////////////////////////////////
        
        std::cout<<"plain_modulus:"<< plain_modulus << std::endl;
        std::cout<<"phase:"<< phase << std::endl;
        std::cout<<"phase_ComplementBits:"<< bitset<32>(phase) << std::endl;
        //std::cout<<"static_cast...:"<< static_cast<typename P::T>(std::round(phase / scale)) << std::endl;
        
        ////////////////////////////////////////////////////////////////////////
        
        typename P::T res = 
        static_cast<typename P::T>(std::round(phase / scale)) % plain_modulus;
        
        temp += phase;
        
        std::cout<<"解密得 :"<< res << std::endl;
        std::cout<<"对照看    :"<< bitset<32>(0UL) << std::endl;
        std::cout<<"无溢出c[n]:"<< bitset<64>(temp) << std::endl;
        
        return make_pair(res, temp );;
    }
    
    

    template <class P>
    typename std::make_signed<typename P::T>::type my_tlweSymInt32Decrypt_print(const TLWE<P> &c, const double scale, const Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
                
        
        typename std::make_signed<typename P::T>::type phase1 = 
        static_cast<typename std::make_signed<typename P::T>::type>(phase);
        
        typename std::make_signed<typename P::T>::type plain_modulus1 = 
        static_cast<typename std::make_signed<typename P::T>::type>(plain_modulus);
        
        ///////////////////////////////////////////////////////////////
        std::cout<<"plain_modulus:"<< plain_modulus << std::endl;
        std::cout<<"phase1:"<< phase1 << std::endl;
        std::cout<<"phase1:"<< bitset<32>(phase1) << std::endl;
        std::cout<<"phase1 / scale:"<< phase1 / scale << std::endl;
        std::cout<<"std::round...:"<< std::round(phase1 / scale) << std::endl;
        std::cout<<"-2 mod plain_modulus:"<< -2 % plain_modulus1 << std::endl;
        std::cout<<"std::round mod plain_modulus:"<< static_cast<typename std::make_signed<typename P::T>::type>(std::round(phase1 / scale)) % plain_modulus1 << std::endl;
        //////////////////////////////////////////////////////////////////
        
        // typename P::T res = 
        // static_cast<typename P::T>(std::round(phase / scale)) % plain_modulus;
        // return res;
        
        typename std::make_signed<typename P::T>::type res = static_cast<typename std::make_signed<typename P::T>::type>(std::round(phase1 / scale)) % plain_modulus1;
        return res;
    }
    
    //  my_tlweSymInt32Decrypt 不输出打印信息版
    template <class P>
    typename std::make_signed<typename P::T>::type my_tlweSymInt32Decrypt(const TLWE<P> &c, const double scale, const Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
                
        
        typename std::make_signed<typename P::T>::type phase1 = 
        static_cast<typename std::make_signed<typename P::T>::type>(phase);
        
        typename std::make_signed<typename P::T>::type plain_modulus1 = 
        static_cast<typename std::make_signed<typename P::T>::type>(plain_modulus);
        
        // typename P::T res = 
        // static_cast<typename P::T>(std::round(phase / scale)) % plain_modulus;
        // return res;
        
        typename std::make_signed<typename P::T>::type res = static_cast<typename std::make_signed<typename P::T>::type>(std::round(phase1 / scale)) % plain_modulus1;
        return res;
    }
    
    
    template <class P>
    bool tlweSymDecrypt_print(const TLWE<P> &c, const Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        /////////////////////////////////////////
        std::cout<<"phase:"<< phase << std::endl;
        std::cout<<"phase_ComplementBits:"<< bitset<32>(phase) << std::endl;
        std::cout<<"make_signed phase:"<< static_cast<typename std::make_signed<typename P::T>::type>(phase) << std::endl;
        std::cout<<"make_signed phase/2^28 :"<< static_cast<typename std::make_signed<typename P::T>::type>(phase) / pow(2., 28)<< std::endl;
        /////////////////////////////////////////
        bool res =
            static_cast<typename std::make_signed<typename P::T>::type>(phase) > 0;
        return res;
    }
    
    template <class P>
    typename P::T tlweSymInt32Decrypt_print(const TLWE<P> &c, const double scale, const Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        
        ////////////////////////////////////////////////////////////////////////
        
        std::cout<<"plain_modulus:"<< plain_modulus << std::endl;
        std::cout<<"phase:"<< phase << std::endl;
        if constexpr (std::is_same_v<P, TFHEpp::lvl1param>) {
            std::cout<<"phase_ComplementBits:"<< bitset<32>(phase) << std::endl;
        }
        if constexpr (std::is_same_v<P, TFHEpp::lvl2param>) {
            std::cout<<"phase_ComplementBits:"<< bitset<64>(phase) << std::endl;
        }
        
        std::cout<<"phase / scale:"<< phase / scale <<"  ";
        std::cout<<"std::round(phase / scale):"<< std::round(phase / scale) <<"  ";
        std::cout<<"static_cast(std::round(phase / scale)):"<< static_cast<typename P::T>(std::round(phase / scale)) << std::endl;
        
        ////////////////////////////////////////////////////////////////////////
        
        typename P::T res = 
        static_cast<typename P::T>(std::round(phase / scale)) % plain_modulus;
        return res;
    }
    
////////////////////////////////////////////////////////////////////////////////

    template <class P>
    typename P::T tlweSymInt32Decrypt(const TLWE<P> &c, const double scale, const Key<P> &key)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        
        typename P::T res = 
        static_cast<typename P::T>(std::round(phase / scale)) % plain_modulus;
        return res;
    }
    
    

    // For debug
    template <class P>
    typename P::T tlweSymInt32Decrypt(const TLWE<P> &c, const double scale, const Key<P> &key, double &noise)
    {
        typename P::T phase = c[P::k * P::n];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        typename P::T res = 
        static_cast<typename P::T>(std::round(phase / scale)) % plain_modulus;
        noise = (std::round(phase / scale) * scale > phase) ? std::round(phase / scale) * scale - phase : phase - std::round(phase / scale) * scale;
        return res;
    }

    template <class P>
    typename P::T tlweSymDecryptNoise(const TLWE<P> &c, typename P::T plain, const Key<P> &key, double &noise)
    {
        typename P::T phase = c[P::k * P::n];
        for (int k = 0; k < P::k; k++)
            for (int i = 0; i < P::n; i++)
                phase -= c[k * P::n + i] * key[k * P::n + i];
        noise = (phase > plain)? phase - plain : plain - phase;
        bool res =
        static_cast<typename std::make_signed<typename P::T>::type>(phase) > 0;
        return res;
    }

    // c = (a, as + m +e)
    template <class P>
    TRLWE<P> trlweSymInt32Encrypt(const array<typename P::T, P::n> &p, const double α, const double scale, const Key<P> &key)
    {
        TRLWE<P> c = trlweSymEncryptZero<P>(α, key);
        for (int i = 0; i < P::n; i++)
            c[P::k][i] += static_cast<typename P::T>(scale* p[i]);
        return c;
    }

    template <class P>
    Polynomial<P> trlweSymInt32Decrypt(const TRLWE<P> &c, double scale, const Key<P> &key)
    {
        Polynomial<P> phase = c[P::k];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++) {
            Polynomial<P> mulres;
            std::array<typename P::T, P::n> partkey;
            for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
            PolyMul<P>(mulres, c[k], partkey);
            for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
        }

        Polynomial<P> p;
        for (int i = 0; i < P::n; i++)
            p[i] = static_cast<typename P::T>(std::round(phase[i] / scale)) % plain_modulus;
        return p;
    }

    template <class P>
    Polynomial<P> trlweSymInt32Decrypt_print(const TRLWE<P> &c, double scale, const Key<P> &key)
    {
        Polynomial<P> phase = c[P::k];
        typename P::T plain_modulus = (1ULL << (std::numeric_limits<typename P::T>::digits -1)) / scale;
        plain_modulus *= 2;
        for (int k = 0; k < P::k; k++) {
            Polynomial<P> mulres;
            std::array<typename P::T, P::n> partkey;
            for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
            PolyMul<P>(mulres, c[k], partkey);
            for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
            
        }

        std::cout << "phase: " << std::endl;
        for(int i = 0 ; i<P::n ; i++){
            std::cout << phase[i] << ' ';
        }
        std::cout << std::endl;
    
        Polynomial<P> p;
        for (int i = 0; i < P::n; i++)
            p[i] = static_cast<typename P::T>(std::round(phase[i] / scale)) % plain_modulus;
        return p;
    }
    
    template <class P>
    constexpr Polynomial<P> μ_polygen(typename P::T μ)
    {
        Polynomial<P> poly;
        for (typename P::T &p : poly) p = -μ;
        return poly;
    }

    template <class P>
    Polynomial<P> gpolygen(uint32_t plain_bits, uint32_t scale_bits)
    {
        Polynomial<P> poly;
        uint32_t padding_bits = P :: nbit - plain_bits;
        for (int i = 0; i< P::n; i++) poly[i] = (1ULL << scale_bits) * (i >> padding_bits);
        return poly;
    }


    template <class P>
    void TLWEAdd(TLWE<P> &ca, TLWE<P> &cb, TLWE<P> &res)
    {
        for (int i = 0; i <= P::k * P::n; i++) res[i] = ca[i] + cb[i];
        
    }

    template <class P>
    void TLWESub(TLWE<P> &ca, TLWE<P> &cb, TLWE<P> &res)
    {
        for (int i = 0; i <= P::k * P::n; i++) res[i] = ca[i] - cb[i];
        
    }

    template <class P>
    void TRLWEAdd(TRLWE<P> &ca, TRLWE<P> &cb, TRLWE<P> &res)
    {
        for (int k = 0; k <= P :: k; k++)
        {
            for (int i = 0; i < P::n; i++) res[k][i] = ca[k][i] + cb[k][i];
        }     
    }

    template <class P>
    void TRLWESub(TRLWE<P> &ca, TRLWE<P> &cb, TRLWE<P> &res)
    {
        for (int k = 0; k <= P :: k; k++)
        {
            for (int i = 0; i < P::n; i++) res[k][i] = ca[k][i] - cb[k][i];
        }
        
    }

    template <class P>
    void TLWEToTRLWE(TLWE<P> &c, TRLWE<P> &res)
    {
        res[0][0] = c[0];
        for (int i = 1; i < P::n; i++)
        {
            res[0][i] = -c[P::n - i];
        }
        res[P::k][0] = c[P::n];
    }

    void ARI_to_LOG(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek);

    void LOG_to_ARI(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek);

    void log_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, uint32_t scale_bits, const EvalKey &ek);

    void ari_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, uint32_t scale_bits, const EvalKey &ek);

    void MSBGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek, bool result_type);
    
    ////////////加的///////////
    void my_MSBGateBootstrapping(TLWE<lvl1param> &res,
                                 const TLWE<lvl1param> &tlwe,
                                 const EvalKey &ek, bool result_type, uint32_t k,
                                 uint32_t plain_bits_eff = 5,
                                 bool use_gap_offset = false);
    void my_MSBGateBootstrapping(TLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, const EvalKey &ek, bool result_type, uint32_t k);
    void my_MSBGateBootstrapping_2(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek, bool result_type, uint32_t k);
    void Multi_HomAND(std::vector<TLWE<lvl1param>> &arr_ciphers, TLWE<lvl1param> &res, TLWE<lvl1param> &k_1, const EvalKey &ek, uint32_t out_bit );
    void N_HomAND(std::vector<TLWE<lvl1param>> &arr_ciphers, TLWE<lvl1param> &res, TLWE<lvl1param> &k_rest, TLWE<lvl1param> &k_7, TLWE<lvl1param> &k_times, const EvalKey &ek, uint32_t out_bit );
    void Multi_HomOR(std::vector<TLWE<lvl1param>> &arr_ciphers, TLWE<lvl1param> &res, TLWE<lvl1param> &k_1, const EvalKey &ek, uint32_t out_bit );
    void N_HomOR(std::vector<TLWE<lvl1param>> &arr_ciphers, TLWE<lvl1param> &res, TLWE<lvl1param> &k_1, const EvalKey &ek, uint32_t out_bit );

    void MSBGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl2param> &tlwe, const EvalKey &ek, bool result_type);
    
    void MSBGateBootstrapping(TLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, const EvalKey &ek, bool result_type);

    void IdeGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, uint32_t scale_bits, const EvalKey &ek);

    void IdeGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl2param> &tlwe, uint32_t scale_bits, const EvalKey &ek);

    void IdeGateBootstrapping(TLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, uint32_t scale_bits, const EvalKey &ek);
} // namespace TFHEpp
