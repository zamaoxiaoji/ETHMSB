#include "operators.h"
#include "tfhepp_utils.h"

namespace HEDB
{
    // c0 + c1 - 1/8

    void HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb, const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T offset = Lvl1::μ;
        if(IS_ARITHMETIC(result_type)) offset = (offset << 1);
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            res[i] = ca[i] + cb[i];
        res[Lvl1::k * Lvl1::n] -= Lvl1::μ >> 1;   // - 1/8
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, res, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<Lvl1>(-offset));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += offset;
        
    }

    void my_HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb, TLWELvl1 &c_1, const TFHEEvalKey &ek, bool result_type, uint32_t k)
    {
        //Lvl1::T offset = Lvl1::μ;
        uint32_t offset = 1ULL << k;
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            res[i] = ca[i] + cb[i];

        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            res[i] = c_1[i] - res[i];

        my_MSBGateBootstrapping(res,res,ek,result_type, 31-k);
    }

    // void my_HomAND(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb,  const TFHEEvalKey &ek, bool result_type, uint32_t k)
    // {
    //     //Lvl1::T offset = Lvl1::μ;
    //     uint32_t offset = 1ULL << k;
    //     for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
    //         res[i] = ca[i] + cb[i];
    //     res[Lvl1::k * Lvl1::n] = res[Lvl1::k * Lvl1::n] - offset;   // - 1/8
    //     my_MSBGateBootstrapping(res,res,ek,result_type,k);
    // }

    // c0 + c1 + 1/8
    void HomOR(TLWELvl1 &res, const TLWELvl1 &ca, const TLWELvl1 &cb, const TFHEEvalKey &ek, bool result_type)
    {
        Lvl1::T offset = Lvl1::μ;
        if(IS_ARITHMETIC(result_type)) offset = (offset << 1);
        for (int i = 0; i <= Lvl1::k * Lvl1::n; i++)
            res[i] = ca[i] + cb[i];
        res[Lvl1::k * Lvl1::n] += (Lvl1::μ >> 1);   // + 1/8
        TLWELvl0 tlwelvl0;
        TFHEpp::IdentityKeySwitch<Lvl10>(tlwelvl0, res, *ek.iksklvl10);
        TFHEpp::GateBootstrappingTLWE2TLWEFFT<Lvl01>(res, tlwelvl0, *ek.bkfftlvl01, TFHEpp::μ_polygen<Lvl1>(-offset));
        if (IS_ARITHMETIC(result_type)) res[Lvl1::k * Lvl1::n] += offset;
    }
    
} // namespace HEDB