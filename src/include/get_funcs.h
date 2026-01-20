
#ifndef GET_FUNCTIONS_H
#define GET_FUNCTIONS_H
#include <vector>
#include "real_type.h"
# include "data_structs.h"



void get_ueinv(const Isol& isol, Real* ueinv);

void get_cp(Post& post, const Oper& oper, const Param& param, const Real edgeVelocity[Ncoords]);

Real get_uk(const Real& incompSpeed,const Param&param,Real&dCompSpeed_dIncompSpeed);

Real get_Mach2(const Real&edgeVel,const Param&param,Real*dMsqrd_dState);

Real get_H(const Real th, const Real ds, Real* H_U);

Real get_Hw(const Real th, const Real wgap, Real* Hw_U);

Real get_Hk(const Real th, const Real ds, const Real ue, const Param&param, Real*Hk_U);

Real get_Hss(const Real th, const Real ds, const Real ue, const Param&param, Real*Hss_U);

Real get_de(const Real th, const Real ds, const Real ue, const Param&param,Real* de_U);

Real get_Ret(const Real th, const Real ds, const Real ue, const Param&param, Real* Ret_U);

Real get_cf(const Real th, const Real ds, const Real sa, const Real ue,const bool turb,const bool wake, const Param& param, Real* cf_U);

Real get_cfxt(const Real th, const Real ds, const Real sa, const Real ue, const Real dist, const bool turb,const bool wake, const Param& param, Real*cfxt_U,Real& cfxt_x);

Real get_cteq(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Param&param,Real (&cteq_U)[4]);

Real get_damp(const Real th,const Real ds,const Real sa,const Real ue, const Param& param, Real (&damp_U)[4]);

Real get_damp_forced(const Real th,const Real ds,const Real sa,const Real ue, const Param& param, const Real&ncrit, Real (&damp_U)[4]);

Real get_cttr(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const Param&param,Real (&cttr_U)[4]);

Real get_Hs(
    const Real th, const Real ds, const Real sa, const Real ue,
    const Param& param, const bool turb, const bool wake,
    Real* Hs_U  // output (length 4)
);

Real get_Us(const Real th, const Real ds, const Real sa, const Real ue, const Param& param,
    const bool turb, const bool wake, Real* Us_U);


Real get_uq(const Real ds,const Real (&ds_U)[8],
const Real cf,const Real (&cf_U)[8],
Real Hk, Real (&Hk_U)[8],
const Real Ret,const Real (&Ret_U)[8],
const bool wake,
const Param&param,
Real (&uq_U)[8]);

Real get_cDi_turbwall(
    const Real th, const Real ds, const Real sa, const Real ue,const bool turb, const bool wake,
    const Param& param,
    Real* cDi_U  // output: length-4
);

Real get_cDi_lam(const Real th,const Real ds,const Real sa,const Real ue,const Param& param,Real (&cDi_U)[4]);

Real get_cDi_lamwake(const Real th,const Real ds,const Real sa,const Real ue, const bool turb, const bool wake, const Param& param,Real (&cDi_U)[4]);

Real get_cDi_outer(const Real th, const Real ds,const Real sa,const Real ue, const bool turb,const bool wake,const Param& param, Real (&cDi_U)[4]);

Real get_cDi_lamstress(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Param& param, Real (&cDi_U)[4]) ;

Real get_cDi(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Param& param,
    Real (&cDi_U)[4]);

Real get_cDixt(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Real dist, const Param& param,
    Real (&cDixt_U)[4],Real& cDixt_x);

Real get_upw(const Real th1,const Real ds1,const Real sa1,const Real ue1,
    const Real th2,const Real ds2,const Real sa2,const Real ue2, const bool wake, const Param&param,
    Real (&upw_U)[8]);
#endif