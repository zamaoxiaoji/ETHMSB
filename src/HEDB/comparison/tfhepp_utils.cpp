#include "tfhepp_utils.h"
#include "HEDB/utils/types.h"

namespace TFHEpp
{

    ////////////////////
    void Multi_HomAND(std::vector<TLWE<lvl1param>> &arr_ciphers, TLWE<lvl1param> &res, TLWE<lvl1param> &k_1, const EvalKey &ek, uint32_t out_bit )
    {
        std::size_t num = arr_ciphers.size();
        TLWE<lvl1param> sum  = arr_ciphers[0];
        TLWE<lvl1param> sub ;

        for(int32_t i = 1 ; i < num ; i++)
        {
            for (size_t j = 0; j <= lvl1param :: n; j++)
            {
                sum[j] = sum[j] + arr_ciphers[i][j];
            }
        }
        
        sub = sum;

        for (size_t i = 0; i <= lvl1param::k * lvl1param::n; i++)
        {
            sub[i] = k_1[i] - sum[i];
        }
        
        my_MSBGateBootstrapping(res, sub, ek, ARITHMETIC, 30 - out_bit +1);
    }

    void Multi_HomOR(std::vector<TLWE<lvl1param>> &arr_ciphers, TLWE<lvl1param> &res, TLWE<lvl1param> &c_0, const EvalKey &ek, uint32_t out_bit )
    {
        std::size_t num = arr_ciphers.size();
        TLWE<lvl1param> sum  = arr_ciphers[0];
        TLWE<lvl1param> sub ;

        for(int32_t i = 1 ; i < num ; i++)
        {
            for (size_t j = 0; j <= lvl1param :: n; j++)
            {
                sum[j] = sum[j] + arr_ciphers[i][j];
            }
        }
        
        sub = sum;
        
        for (size_t i = 0; i <= lvl1param::k * lvl1param::n; i++)
        {
            sub[i] = c_0[i] - sum[i];
        }
        
        my_MSBGateBootstrapping(res, sub, ek, ARITHMETIC, 30 - out_bit +1);
    }

    // In arithmetic mode, 1 -> 1 *scale, 0 -> 0 * scale
    // In logic mode, 1 -> q/8, 0 -> -q/8
    void ARI_to_LOG(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek)
    {
        lvl1param::T μ = lvl1param::μ;
        TLWE<lvl1param> tlweoffset = tlwe;
        tlweoffset[lvl1param::k * lvl1param::n] += μ;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
    }

    void LOG_to_ARI(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek)
    {
        lvl1param::T μ = lvl1param::μ;
        μ = μ << 1;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(-μ));
        res[lvl1param::k * lvl1param::n] += lvl1param::μ;
    }

    void log_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, uint32_t scale_bits, const EvalKey &ek)
    {
        lvl1param::T μ = 1ULL << (scale_bits - 1);
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(-μ));
        res[lvl1param::k * lvl1param::n] += μ;
    }

    void ari_rescale(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, uint32_t scale_bits, const EvalKey &ek)
    {
        lvl1param::T μ = 1ULL << (scale_bits - 1);
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
        TLWE<lvl1param> tlweoffset = tlwe;
        tlweoffset[lvl1param::k * lvl1param::n] += offset;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
        res[lvl1param::k * lvl1param::n] += μ;
    }
    
    ////////////////////////////////////////////////////////////////////////
    void my_MSBGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek, bool result_type, uint32_t k)
    {
        
        lvl1param::T μ = 1U << 29;
        if(IS_LOGIC(result_type)) μ = (μ << 2) >> k;
        if (IS_ARITHMETIC(result_type)) μ = (μ << 1) >> k;
        //μ = (μ << 1) >> k;
        uint64_t offset = 1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
        offset = (3*offset)/4 ;
        
        TLWE<lvl1param> tlweoffset = tlwe;
        tlweoffset[lvl1param::k * lvl1param::n] += offset;
        
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
        
        if (IS_ARITHMETIC(result_type)) res[lvl1param::k * lvl1param::n] += (μ);
    }

    /// 测试验证用的
    void my_MSBGateBootstrapping_2(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek, bool result_type, uint32_t k)
    {
        
        uint32_t μ = 1U << 29;
        //if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        μ = (μ << 1) >> k;
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
        TLWE<lvl1param> tlweoffset = tlwe;
        tlweoffset[lvl1param::k * lvl1param::n] += offset;
        
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(-μ));
        
        if (IS_ARITHMETIC(result_type)) res[lvl1param::k * lvl1param::n] += (μ);
    }

    void my_MSBGateBootstrapping(TLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, const EvalKey &ek, bool result_type, uint32_t k)
    {
        uint64_t μ = 1ULL << 61;
        if(IS_LOGIC(result_type)) μ = (μ << 2) >> k;
        if (IS_ARITHMETIC(result_type)) μ = (μ << 1) >> k;
        //μ = (μ << 1) >> k;
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl2param::T>::digits - 7);
        TLWE<lvl2param> tlweoffset = tlwe;
        tlweoffset[lvl2param::k * lvl2param::n] += offset;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        GateBootstrappingTLWE2TLWEFFT<lvl02param>(res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<lvl2param>(μ));
        if (IS_ARITHMETIC(result_type)) res[lvl2param::k * lvl2param::n] += μ;
    }
    ////////////////////////////////////////////////////////////////////////

    // ARI 0 / 1/2, LOG 1/8 , -1/8
    void MSBGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, const EvalKey &ek, bool result_type)
    {
        uint32_t μ = 1U << 29;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;  
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
        TLWE<lvl1param> tlweoffset = tlwe;
        tlweoffset[lvl1param::k * lvl1param::n] += offset;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
        if (IS_ARITHMETIC(result_type)) res[lvl1param::k * lvl1param::n] += μ;
    }

    void MSBGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl2param> &tlwe, const EvalKey &ek, bool result_type)
    {
        uint32_t μ = 1U << 29;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl2param::T>::digits - 6);
        TLWE<lvl2param> tlweoffset = tlwe;
        tlweoffset[lvl2param::k * lvl2param::n] += offset;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, μ_polygen<lvl1param>(μ));
        if (IS_ARITHMETIC(result_type)) res[lvl1param::k * lvl1param::n] += μ;
    }
    
    void MSBGateBootstrapping(TLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, const EvalKey &ek, bool result_type)
    {
        uint64_t μ = 1ULL << 61;
        if (IS_ARITHMETIC(result_type)) μ = μ << 1;
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl2param::T>::digits - 7);
        TLWE<lvl2param> tlweoffset = tlwe;
        tlweoffset[lvl2param::k * lvl2param::n] += offset;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        GateBootstrappingTLWE2TLWEFFT<lvl02param>(res, tlwelvl0, *ek.bkfftlvl02, μ_polygen<lvl2param>(μ));
        if (IS_ARITHMETIC(result_type)) res[lvl2param::k * lvl2param::n] += μ;
    }

    void IdeGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe, uint32_t scale_bits, const EvalKey &ek)
    {
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl1param::T>::digits - 6);
        TLWE<lvl1param> tlweoffset = tlwe;
        tlweoffset[lvl1param::k * lvl1param::n] += offset;
        constexpr uint32_t plain_bits = 4;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl10param>(tlwelvl0, tlweoffset, *ek.iksklvl10);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, gpolygen<lvl1param>(plain_bits, scale_bits));
    }

    void IdeGateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl2param> &tlwe, uint32_t scale_bits, const EvalKey &ek)
    {
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl2param::T>::digits - 6);
        TLWE<lvl2param> tlweoffset = tlwe;
        tlweoffset[lvl2param::k * lvl2param::n] += offset;
        constexpr uint32_t plain_bits = 4;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01, gpolygen<lvl1param>(plain_bits, scale_bits - 32));
    }

    void IdeGateBootstrapping(TLWE<lvl2param> &res, const TLWE<lvl2param> &tlwe, uint32_t scale_bits, const EvalKey &ek)
    {
        constexpr uint64_t offset = 1ULL << (std::numeric_limits<lvl2param::T>::digits - 7);
        TLWE<lvl2param> tlweoffset = tlwe;
        tlweoffset[lvl2param::k * lvl2param::n] += offset;
        constexpr uint32_t plain_bits = 5;
        TLWE<lvl0param> tlwelvl0;
        IdentityKeySwitch<lvl20param>(tlwelvl0, tlweoffset, *ek.iksklvl20);
        GateBootstrappingTLWE2TLWEFFT<lvl02param>(res, tlwelvl0, *ek.bkfftlvl02, gpolygen<lvl2param>(plain_bits, scale_bits));
    }

    void N_HomAND(std::vector<TLWE<lvl1param>> &arr_ciphers, TLWE<lvl1param> &res, TLWE<lvl1param> &k_rest, TLWE<lvl1param> &k_7, TLWE<lvl1param> &c_0, const EvalKey &ek, uint32_t out_bit )
    {

        std::size_t total_num = arr_ciphers.size();

        if(total_num < 7){
            Multi_HomAND(arr_ciphers, res , k_rest, ek , out_bit);
            return ;
        }
        else if(total_num == 7){
            Multi_HomAND(arr_ciphers, res , k_7, ek , out_bit);
            return ;
        }

        std::size_t loop_times = total_num / 7;
        if( total_num %7 != 0){
            loop_times += 1;
        }

        std::vector<TLWE<lvl1param>> boot_res(loop_times );

        if(loop_times > 1){
            for(int k = 0 ;k < loop_times; k++){
                if( (k+1)*7 <= total_num){
                    std::vector<TLWE<lvl1param>> new_ciphers(arr_ciphers.begin()+ k * 7, arr_ciphers.begin() + (k+1)*7);
                    Multi_HomAND(new_ciphers, boot_res[k] , k_7, ek , 28);
                }
                else{
                    if( k * 7 == total_num -1){
                        boot_res[k] = arr_ciphers[total_num -1];
                    }
                    else{
                        std::vector<TLWE<lvl1param>> new_ciphers(arr_ciphers.begin()+ k * 7, arr_ciphers.end());
                        Multi_HomAND(new_ciphers, boot_res[k] , k_rest, ek , 28);
                    }
                    
                }
            }
        }
        

        // if(total_num < 7){
        //     Multi_HomAND(arr_ciphers, res , k_rest, ek , out_bit);
        // }
        // else if(total_num == 7){
        //     Multi_HomAND(arr_ciphers, res , k_7, ek , out_bit);
        // }
        // else{
        //     Multi_HomAND(boot_res, res , k_times, ek , out_bit);
        // }
        
        TLWE<lvl1param> c_loop_rest = c_0;
        uint32_t loop_rest_scale = (1ULL << 28) * ((loop_times %7) -1) ;
        c_loop_rest[lvl1param::k * lvl1param::n] += loop_rest_scale;

        // if(loop_times < 7){
        //     Multi_HomAND(boot_res, res , c_loop_rest, ek , out_bit);
        // }
        // else if(loop_times == 7){
        //     Multi_HomAND(boot_res, res , k_7, ek , out_bit);
        // }
        // else{
        //     N_HomAND(boot_res, res , c_loop_rest, k_7, c_0, ek , out_bit);
        // }
        N_HomAND(boot_res, res , c_loop_rest, k_7, c_0, ek , out_bit);
    
    }

    void N_HomOR(std::vector<TLWE<lvl1param>> &arr_ciphers,TLWE<lvl1param> &res,  TLWE<lvl1param> &c_0, const EvalKey &ek, uint32_t out_bit )
    {

        std::size_t total_num = arr_ciphers.size();
        if(total_num <= 7){
            Multi_HomOR(arr_ciphers, res , c_0, ek , out_bit);
            return ;
        }
        
        std::size_t loop_times = total_num / 7;
        if( total_num %7 != 0){
            loop_times += 1;
        }
        
        std::vector<TLWE<lvl1param>> boot_res(loop_times );

        if(loop_times > 1){
            for(int k = 0 ;k < loop_times; k++){
                if( (k+1)*7 <= total_num){
                    std::vector<TLWE<lvl1param>> new_ciphers(arr_ciphers.begin()+ k * 7, arr_ciphers.begin() + (k+1)*7);
                    Multi_HomOR(new_ciphers, boot_res[k] , c_0, ek , 28);
                }
                else{
                    if( k * 7 == total_num -1){
                        boot_res[k] = arr_ciphers[total_num -1];
                    }
                    else{
                        std::vector<TLWE<lvl1param>> new_ciphers(arr_ciphers.begin()+ k * 7, arr_ciphers.end());
                        Multi_HomOR(new_ciphers, boot_res[k] , c_0, ek , 28);
                    }
                    
                }
            }
        }
        
        // if(total_num <= 7){
        //     Multi_HomOR(arr_ciphers, res , c_0, ek , out_bit);
        // }
        // else{
        //     Multi_HomOR(boot_res, res , c_0, ek , out_bit);
        // }

        // if(loop_times <= 7){
        //     Multi_HomOR(boot_res, res , c_0, ek , out_bit);
        // }
        // else{
        //     N_HomOR(boot_res, res ,c_0, ek , out_bit);
        // }
        N_HomOR(boot_res, res ,c_0, ek , out_bit);
        
    }


} // namespace TFHEpp