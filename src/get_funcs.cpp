#include <iostream>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "main_func.h"




/*

For optimisation, a lot of functions involve re-calculating parametrs
that are already found in earlier sections of the code. Can reduce 
computation 

*/

void get_ueinv(const Isol& isol, Real* ueinv){

    for (int i = 0; i < Ncoords; ++i) {
        ueinv[i] = isol.edgeVelSign[i]*isol.gammas[i];
    }

    for (int i = 0; i < Nwake; ++i) {
        ueinv[Ncoords+i] = isol.uewi[i] ;
    }

    ueinv[Ncoords] = ueinv[Ncoords-1]; // continuity of upper surface and wake velocity
};


void get_cp(Post& post, const Oper& oper, const Param& param, const Real edgeVelocity[Ncoords+Nwake]) {
    constexpr int Nsys = Ncoords+Nwake;
    if (oper.Ma > 0.0) {
        for (int i = 0; i < Nsys; ++i) {
            post.cp[i] = (1 - (edgeVelocity[i] / oper.Vinf) * (edgeVelocity[i] / oper.Vinf));

            Real den = param.KTb + 0.5 * param.KTl * (1 + param.KTb) * post.cp[i];
            post.cp[i] /= den;
        }
    } else {
        for (int i = 0; i < Nsys; ++i) {
            post.cp[i] = 1 - (edgeVelocity[i] / oper.Vinf) * (edgeVelocity[i] / oper.Vinf);
        }
    }
}


Real get_uk(const Real& incompSpeed,const Param&param,Real&dCompSpeed_dIncompSpeed){

    // Finds compressible speed using corrections

    Real compSpeed = incompSpeed;
    dCompSpeed_dIncompSpeed = 1.0;
    
    if (param.Minf > 0.0){
        Real l = param.KTl;
        Real Vinf = param.Vinf;

        Real den = 1.0-l*((incompSpeed/Vinf)*(incompSpeed/Vinf));
        Real den_u = -2.0*l*incompSpeed/(Vinf*Vinf);

        compSpeed = incompSpeed*(1-l)/den ;
        dCompSpeed_dIncompSpeed = (1-l)/den - (compSpeed/den)*den_u ;
    }

    return compSpeed ;
}

Real get_Mach2(const Real&edgeVel,const Param&param,Real*dMsqrd_dState){

    // Finding squared Mach number, doing compressible correction if needed
    // Returns sqaured mach by value, and change wrt state U as reference.
    // Initialise dMsqrd_sState as 4 zeros !!
    
    Real M2 = 0.0;
    dMsqrd_dState[0] = 0;
    dMsqrd_dState[1] = 0;
    dMsqrd_dState[2] = 0;
    dMsqrd_dState[3] = 0;

    if (param.Minf > 0.0){
        Real H0 = param.H0;
        Real g = param.gam;

        Real compSpeed = 0.0;
        Real dCompSpeed_dIncompSpeed = 0.0;
        
        compSpeed = get_uk(edgeVel,param,dCompSpeed_dIncompSpeed);
        


        Real c2 = (g-1)*(H0-0.5*compSpeed*compSpeed);
        Real c2_compSpeed = (g-1)*(-compSpeed); 

        M2 = (compSpeed*compSpeed)/c2;
        Real M2_uk = 2*compSpeed/c2 - M2/c2*c2_compSpeed ;

        dMsqrd_dState[3] = M2_uk*dCompSpeed_dIncompSpeed;

        return M2;


    }
    else
    
    { return M2;}
}

Real get_H(const Real th, const Real ds, Real* H_U){
    
    Real H = ds / th;

    H_U[0] = -H / th;
    H_U[1] = 1.0 / th;
    H_U[2] = 0.0;
    H_U[3] = 0.0;

    return H;
}

Real get_Hw(const Real th, const Real wgap, Real* Hw_U){
    
    const Real Hw = wgap/th ;

    Hw_U[0] = -Hw / th;
    Hw_U[1] = 0.0;
    Hw_U[2] = 0.0;
    Hw_U[3] = 0.0;

    return Hw;
}

Real get_Hk(const Real th, const Real ds, const Real ue, const Param&param, Real*Hk_U){

    // gets the kinematic shape parameter and also change wrt state U.
    
    
    Real H_U[4] = {0.0};
    Real H = get_H(th,ds,H_U);

    Real Hk = H;

    if (param.Minf >0.0){

        Real M2_U[4]= {0.0};
        Real M2 = get_Mach2(ue,param,M2_U);
        Real den = (1+0.113*M2); 
        Real den_M2 = 0.113;
        Hk = (H-0.29*M2)/den;
        
        for (int i=0;i<4;++i){
            Hk_U[i] = (H_U[i]-0.29*M2_U[i])/den - Hk/den*den_M2*M2_U[i];
        }
    }
    else {

        for (int i=0;i<4;++i){
            Hk_U[i] = H_U[i];
        }
    }

    return Hk;
}


Real get_Hss(const Real th, const Real ds, const Real ue, const Param&param, Real*Hss_U){

    // gets density shape parameter and change wrt state U

    Real M2_U[4] = {0.0};
    Real Hk_U[4] = {0.0};

    Real M2 = get_Mach2(ue,param,M2_U);
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);

    Real num = 0.064/(Hk-0.8) + 0.251;
    Real num_U[4]={0} ;

    for (int i=0;i<4;++i){
        num_U[i] = -0.064/((Hk-0.8)*(Hk-0.8))*Hk_U[i];
    }
    
    Real Hss = M2*num;
    for (int i=0;i<4;++i){
        Hss_U[i] = M2_U[i]*num + M2*num_U[i];
    }

    return Hss ;
}


Real get_de(const Real th, const Real ds, const Real ue, const Param&param,Real* de_U){

    Real Hk_U[4] = {0.0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U) ;

    Real aa = 3.15 + 1.72/(Hk-1);
    
    Real aa_U[4]={0};

    for (int i=0;i<4;++i){
        aa_U[i] = -1.72/((Hk-1)*(Hk-1))*Hk_U[i] ;
    }

    Real de = th*aa + ds;
    
    Real temp[4] = {aa,1.0,0.0,0.0};

    for (int i=0;i<4;++i){
        de_U[i] = temp[i] + th*aa_U[i] ;
    }
    
    Real dmx = 12;
    if (de > dmx*th){
         
        de = dmx*th;
        for (int i=0;i<4;++i){ de_U[i] = 0.0 ;}
        de_U[0] = dmx ;
    }

    return de;
}


Real get_Ret(const Real th, const Real ds, const Real ue, const Param&param, Real* Ret_U){


    Real Ret = param.rho0*th*ue/param.mu0 ;
    if (param.Minf > 0.0){

    
        Real M2_U[4] = {0.0};
        Real uk_u = 0.0;

        Real M2 = get_Mach2(ue,param,M2_U);
        Real uk = get_uk(ue,param,uk_u);

        Real H0 = param.H0;
        Real gmi = param.gam-1; 
        Real Ts =  param.Tsrat ;


        Real Tr = 1-0.5*(uk*uk)/H0; 
        Real Tr_uk = -uk/H0;  // edge/stagnation temperature ratio

        Real f = (std::pow(Tr,1.5))*(1+Ts)/(Tr+Ts); 
        Real f_Tr = 1.5*f/Tr-f/(Tr+Ts); // Sutherland's ratio
        Real mu = param.mu0*f;
        Real mu_uk = param.mu0*f_Tr*Tr_uk; // local dynamic viscosity
        Real den = 1+0.5*gmi*M2;
        Real den_M2 = 0.5*gmi;
        Real rho = param.rho0/(std::pow(den,1/gmi));
        Real rho_U[4]={0}; // density
        
        for (int i=0;i<4;++i){
            rho_U[i] = (-1/gmi)*rho/den*den_M2*M2_U[i];
        }

        Ret =rho*uk*th/mu ;

        Real const1 = rho*th/mu-Ret/mu*mu_uk ;
        Real const2 = rho*uk/mu ;

        Real array1[4] = {0,0,0,uk_u};
        Real array2[4] = {1,0,0,0};

        for (int i=0;i<4;++i){
            Ret_U[i] = rho_U[i]*uk*th/mu + const1*array1[i] + const2*array2[i] ;
        }
    }
    else {
        Ret_U[0] = ue/param.mu0;
        Ret_U[1] = 0;
        Ret_U[2] = 0;
        Ret_U[3] = th/param.mu0 ;
    }

    return Ret ;
}

Real get_cf(const Real th, const Real ds, const Real sa, const Real ue,const bool turb,const bool wake, const Param& param, Real* cf_U)  // output: cf linearisation w.r.t. th, ds, sa, ue
{
    for (int i = 0; i < 4; ++i) cf_U[i] = 0.0;

    if (wake) return 0.0;

    Real Hk_U[4] = {0.0}, Ret_U[4] = {0.0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real Ret = get_Ret(th,ds,ue,param,Ret_U);

    if (turb) {
        Real M2_U[4]={0};
        Real M2 = get_Mach2(ue, param, M2_U);

        const Real g1 = 0.5 * (param.gam - 1.0);
        Real Fc = std::sqrt(1.0 + g1 * M2);
        Real Fc_U[4]={0};
        Real invFc = 1.0 / Fc;
        for (int i = 0; i < 4; ++i)
            Fc_U[i] = g1 * 0.5 * M2_U[i] * invFc;

        Real aa = -1.33 * Hk;
        Real aa_U[4]={0};
        for (int i = 0; i < 4; ++i)
            aa_U[i] = -1.33 * Hk_U[i];

        if (aa < -17.0) {
            Real expterm = std::exp((aa + 17.0) / 3.0);
            aa = -20.0 + 3.0 * expterm;
            Real scale = (aa + 20.0) / 3.0;
            for (int i = 0; i < 4; ++i)
                aa_U[i] *= scale;
        }

        Real bb = std::log(Ret / Fc);
        Real bb_U[4]={0};
        for (int i = 0; i < 4; ++i)
            bb_U[i] = Ret_U[i] / Ret - Fc_U[i] / Fc;

        if (bb < 3.0) {
            bb = 3.0;
            for (int i = 0; i < 4; ++i) bb_U[i] = 0.0;
        }

        const Real log10inv = 1.0 / std::log(10.0);
        bb *= log10inv;
        for (int i = 0; i < 4; ++i)
            bb_U[i] *= log10inv;

        Real cc = -1.74 - 0.31 * Hk;
        Real cc_U[4]={0};
        for (int i = 0; i < 4; ++i)
            cc_U[i] = -0.31 * Hk_U[i];

        Real dd = std::tanh(4.0 - Hk / 0.875);
        Real dd_U[4]={0};
        Real sech2 = 1.0 - dd * dd;
        for (int i = 0; i < 4; ++i)
            dd_U[i] = sech2 * (-Hk_U[i] / 0.875);

        Real log_bb = std::log(bb);
        Real bb_pow = std::pow(bb, cc);
        Real exp_aa = std::exp(aa);
        Real cf0 = 0.3 * exp_aa * bb_pow;
        Real cf0_U[4]={0};
        for (int i = 0; i < 4; ++i)
            cf0_U[i] = cf0 * aa_U[i]
                    + 0.3 * exp_aa * cc * std::pow(bb, cc-1) * bb_U[i]
                    + cf0 * log_bb * cc_U[i];

        Real cf = (cf0 + 1.1e-4 * (dd - 1.0)) * invFc;
        for (int i = 0; i < 4; ++i)
            cf_U[i] = (cf0_U[i] + 1.1e-4 * dd_U[i]) * invFc - cf * invFc * Fc_U[i];

        return cf;
    } else {
        Real num, num_Hk;
        if (Hk < 5.5) {
            Real dh = 5.5 - Hk;
            num = 0.0727 * dh * dh * dh / (Hk + 1.0) - 0.07;
            num_Hk = 0.0727 * (3.0 * dh * dh / (Hk + 1.0) * (-1.0)
                      - dh * dh * dh / ((Hk + 1.0) * (Hk + 1.0)));
        } else {
            Real denom = Hk - 4.5;
            Real factor = 1.0 - 1.0 / denom;
            num = 0.015 * factor * factor - 0.07;
            num_Hk = 0.015 * 2.0 * factor / (denom * denom);
        }

        Real cf = num / Ret;
        for (int i = 0; i < 4; ++i)
            cf_U[i] = (num_Hk * Hk_U[i] - num * Ret_U[i] / Ret) / Ret;

        return cf;
    }
}


Real get_cfxt(const Real th, const Real ds, const Real sa, const Real ue, const Real dist, const bool turb,const bool wake, const Param& param, Real*cfxt_U,Real& cfxt_x){

    Real cf_U[4] = {0.0};

    Real cf = get_cf(th,ds,sa,ue,turb,wake,param,cf_U); 
    Real cfxt = cf*dist/th ;
    
    
    for (int i=0; i<4;++i){
        cfxt_U[i] = cf_U[i]*dist/th;
    }
    cfxt_U[0] -= cfxt/th ;
    
    cfxt_x = cf/th;

    return cfxt ;
}



Real get_Hs(
    const Real th, const Real ds, const Real sa, const Real ue,
    const Param& param, const bool turb, const bool wake,
    Real* Hs_U  // output (length 4)
) {
    
    Real Hs;
    
    for (int i = 0; i < 4; ++i) Hs_U[i] = 0.0;

    Real Hk_U[4]={0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);

    // Limit Hk
    if (wake && Hk < 1.00005) {
        Hk = 1.00005;
        for (int i = 0; i < 4; ++i) Hk_U[i] = 0.0;
    }
    if (!wake && Hk < 1.05) {
        Hk = 1.05;
        for (int i = 0; i < 4; ++i) Hk_U[i] = 0.0;
    }

    if (turb) {
        const Real Hsmin = 1.5;
        const Real dHsinf = 0.015;

        Real Ret_U[4]={0};
        Real Ret = get_Ret(th,ds,ue,param,Ret_U);
        Real Ho = 4.0;
        Real Ho_U[4] = {0.0};
        if (Ret > 400.0) {
            Ho = 3.0 + 400.0 / Ret;
            for (int i = 0; i < 4; ++i)
                Ho_U[i] = -400.0 / (Ret * Ret) * Ret_U[i];
        }

        Real Reb = Ret;
        Real Reb_U[4]={0};
        for (int i = 0; i < 4; ++i) Reb_U[i] = Ret_U[i];
        if (Ret < 200.0) {
            Reb = 200.0;
            for (int i = 0; i < 4; ++i) Reb_U[i] = 0.0;
        }

        if (Hk < Ho) {
            Real Hr = (Ho - Hk) / (Ho - 1.0);
            Real Hr_U[4]={0};
            for (int i = 0; i < 4; ++i)
                Hr_U[i] = (Ho_U[i] - Hk_U[i]) / (Ho - 1.0)
                        - (Ho - Hk) / ((Ho - 1.0) * (Ho - 1.0)) * Ho_U[i];

            Real invReb = 1.0 / Reb;
            Real aa = (2.0 - Hsmin - 4.0 * invReb) * Hr * Hr;
            Real aa_U[4]={0};
            for (int i = 0; i < 4; ++i)
                aa_U[i] = (4.0 / (Reb * Reb)) * Reb_U[i] * Hr * Hr +
                          (2.0 - Hsmin - 4.0 * invReb) * 2.0 * Hr * Hr_U[i];

            Real denom = Hk + 0.5;
            Hs = Hsmin + 4.0 * invReb + aa * 1.5 / denom;
            for (int i = 0; i < 4; ++i)
                Hs_U[i] = -4.0 / (Reb * Reb) * Reb_U[i]
                        + aa_U[i] * 1.5 / denom
                        - aa * 1.5 / (denom * denom) * Hk_U[i];
        } 
        else {
            Real lrb = std::log(Reb);
            Real lrb_U[4]={0};
            for (int i = 0; i < 4; ++i)
                lrb_U[i] = (1.0 / Reb) * Reb_U[i];

            Real invLrb2 = 1.0 / (lrb * lrb);
            Real aa = Hk - Ho + 4.0 / lrb;
            Real aa_U[4]={0};
            for (int i = 0; i < 4; ++i)
                aa_U[i] = Hk_U[i] - Ho_U[i] - 4.0 * invLrb2 * lrb_U[i];

            Real bb = 0.007 * lrb / (aa * aa) + dHsinf / Hk;
            Real bb_U[4]={0};
            for (int i = 0; i < 4; ++i)
                bb_U[i] = 0.007 * (lrb_U[i] / (aa * aa) - 2.0 * lrb * aa_U[i] / (aa * aa * aa))
                        - dHsinf / (Hk * Hk) * Hk_U[i];

            Real deltaH = Hk - Ho;
            Hs = Hsmin + 4.0 / Reb + deltaH * deltaH * bb;
            for (int i = 0; i < 4; ++i)
                Hs_U[i] = -4.0 / (Reb * Reb) * Reb_U[i]
                        + 2.0 * deltaH * (Hk_U[i] - Ho_U[i]) * bb
                        + deltaH * deltaH * bb_U[i];
        }

        // Mach correction
        Real M2_U[4]={0};
        Real M2 = get_Mach2(ue,param,M2_U);
        Real den = 1.0 + 0.014 * M2;
        Real den_U[4]={0};
        for (int i = 0; i < 4; ++i) den_U[i] = 0.014 * M2_U[i];

        Real num = Hs + 0.028 * M2;
        for (int i = 0; i < 4; ++i)
            Hs_U[i] = (Hs_U[i] + 0.028 * M2_U[i]) / den - num / (den * den) * den_U[i];
        Hs = num / den;
    } 
    else {
        Real a = Hk - 4.35;
        Real Hs_Hk, num;
        if (Hk < 4.35) {
            num = 0.0111 * a * a - 0.0278 * a * a * a;
            Real denom = Hk + 1.0;
            Hs = num / denom + 1.528 - 0.0002 * (a * Hk) * (a * Hk);
            Hs_Hk = (0.0111 * 2 * a - 0.0278 * 3 * a * a) / denom
                    - num / (denom * denom)
                    - 0.0002 * 2.0 * a * Hk * (Hk + a);
        } else {
            Hs = 0.015 * a * a / Hk + 1.528;
            Hs_Hk = 0.015 * 2.0 * a / Hk - 0.015 * a * a / (Hk * Hk);
        }
        for (int i = 0; i < 4; ++i)
            Hs_U[i] = Hs_Hk * Hk_U[i];
    }

    return Hs;
}



Real get_Us(const Real th, const Real ds, const Real sa, const Real ue, const Param& param,
    const bool turb, const bool wake, Real* Us_U) {
    
    for (int i = 0; i < 4; ++i) Us_U[i] = 0.0;

    Real Hs_U[4]={0}, Hk_U[4]={0}, H_U[4]={0};
    Real Hs = get_Hs(th,ds,sa,ue,param,turb,wake,Hs_U);
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real H  = get_H(th, ds, H_U);

    // limit Hk
    if (wake && Hk < 1.00005) {
        Hk = 1.00005;
        for (int i = 0; i < 4; ++i) Hk_U[i] = 0.0;
    }
    if (!wake && Hk < 1.05) {
        Hk = 1.05;
        for (int i = 0; i < 4; ++i) Hk_U[i] = 0.0;
    }

    Real beta = param.GB;
    Real bi = 1.0 / beta;
    Real Hk_minus_1 = Hk - 1.0;
    Real invH = 1.0 / H;
    Real invH2 = invH * invH;
    Real factor = 1.0 - bi * Hk_minus_1 * invH;

    Real Us = 0.5 * Hs * factor;

    for (int i = 0; i < 4; ++i) {
        Real dHk = Hk_U[i];
        Real dHs = Hs_U[i];
        Real dH  = H_U[i];

        Us_U[i] = 0.5 * dHs * factor
                + 0.5 * Hs * (-bi * dHk * invH + bi * Hk_minus_1 * dH * invH2);
    }

    // limit Us
    if (!wake && Us > 0.95) {
        Us = 0.98;
        for (int i = 0; i < 4; ++i) Us_U[i] = 0.0;
    }
    if (!wake && Us > 0.99995) {
        Us = 0.99995;
        for (int i = 0; i < 4; ++i) Us_U[i] = 0.0;
    }

    return Us;
}


Real get_uq(const Real ds,const Real (&ds_U)[8],
const Real cf,const Real (&cf_U)[8],
Real Hk, Real (&Hk_U)[8],
const Real Ret,const Real (&Ret_U)[8],
const bool wake,
const Param&param,
Real (&uq_U)[8])
{
    Real beta = param.GB,A=param.GA,C=param.GC;
    if (wake){A *= param.Dlr, C = 0;}

    if (wake && Hk<1.00005){
        Hk = 1.00005;
        for (int i =0;i<8;++i){
            Hk_U[i] = 0;
        }
    }

    if (!wake && Hk<1.05){
        Hk = 1.05;
        for (int i =0;i<8;++i){
            Hk_U[i] = 0;
        }
    }

    Real Hkc = Hk-1 - C/Ret ;
    Real Hkc_U[8]={0};
    for (int i =0;i<8;++i){Hkc_U[i] = Hk_U[i] + C/(Ret*Ret)*Ret_U[i];}

    if (Hkc < 0.01){
        Hkc = 0;
        for (int i =0;i<8;++i){
            Hkc_U[i] = 0;
        }
    }

    Real AxHk = (A*Hk) ;
    Real ut = 0.5*cf - (Hkc/(AxHk))*(Hkc/(AxHk));
    Real ut_U[8]={0};
    for (int i =0;i<8;++i){
        ut_U[i] = 0.5*cf_U[i] - 2*(Hkc/(AxHk))*(Hkc_U[i]/(AxHk) - Hkc/(AxHk*Hk)*Hk_U[i]);
    }

    Real  uq = ut/(beta*ds);
    for (int i =0;i<8;++i){
        uq_U[i] = ut_U[i]/(beta*ds) - uq/ds * ds_U[i] ;
    }
    
    return uq ;





}



Real get_cDi_turbwall(
    const Real th, const Real ds, const Real sa, const Real ue,const bool turb, const bool wake,
    const Param& param,
    Real* cDi_U  // output: length-4
) {
    for (int i = 0; i < 4; ++i) cDi_U[i] = 0.0;

    if (wake)
        return 0.0;

    // Get required quantities and derivatives
    Real cf_U[4]={0}, Hk_U[4]={0}, Hs_U[4]={0}, Us_U[4]={0}, Ret_U[4]={0};
    Real cf = get_cf(th, ds, sa, ue,turb,wake, param, cf_U);
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real Hs = get_Hs(th,ds,sa,ue,param,turb,wake,Hs_U);
    Real Us = get_Us(th, ds, sa, ue, param, turb, wake, Us_U);
    Real Ret = get_Ret(th,ds,ue,param,Ret_U);

    // Compute log(Ret) and its derivative
    Real lr = std::log(Ret);
    Real lr_U[4]={0};
    for (int i = 0; i < 4; ++i)
        lr_U[i] = Ret_U[i] / Ret;

    // Hmin = 1 + 2.1 / lr
    Real Hmin = 1.0 + 2.1 / lr;
    Real Hmin_U[4]={0};
    for (int i = 0; i < 4; ++i)
        Hmin_U[i] = -2.1 / (lr * lr) * lr_U[i];

    // fac = 0.5 + 0.5 * tanh((Hk - 1)/(Hmin - 1))
    Real num = Hk - 1.0;
    Real denom = Hmin - 1.0;
    Real x = num / denom;
    Real tanhx = std::tanh(x);
    Real fac = 0.5 + 0.5 * tanhx;

    Real fac_U[4]={0};
    Real sech2 = 1.0 - tanhx * tanhx;
    for (int i = 0; i < 4; ++i)
        fac_U[i] = 0.5 * sech2 *
                   (Hk_U[i] / denom - num / (denom * denom) * Hmin_U[i]);

    // cDi = 0.5 * cf * Us * (2 / Hs) * fac
    Real inv_Hs = 1.0 / Hs;
    Real two_over_Hs = 2.0 * inv_Hs;
    Real cDi = 0.5 * cf * Us * two_over_Hs * fac;

    for (int i = 0; i < 4; ++i) {
        cDi_U[i] =
            cf_U[i] * Us * inv_Hs * fac +
            cf * Us_U[i] * inv_Hs * fac -
            cDi * Hs_U[i] / Hs +
            cf * Us * inv_Hs * fac_U[i];
    }

    return cDi;
}

Real get_cDi_lam(const Real th,const Real ds,const Real sa,const Real ue,const Param& param,Real (&cDi_U)[4]) {

    Real Hk_U[4]={0}, Ret_U[4]={0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real Ret = get_Ret(th,ds,ue,param,Ret_U);

    Real num, num_Hk;
    if (Hk < 4.0) {
    Real dHk = 4.0 - Hk;
    num = 0.00205 * std::pow(dHk, 5.5) + 0.207;
    num_Hk = -0.00205 * 5.5 * std::pow(dHk, 4.5);
    } else {
    Real Hk1 = Hk - 4.0;
    Real Hk1_2 = Hk1 * Hk1;
    Real denom = 1.0 + 0.02 * Hk1_2;
    num = -0.0016 * Hk1_2 / denom + 0.207;
    num_Hk = -0.0016 * (2.0 * Hk1 / denom - Hk1_2 * 0.04 * Hk1 / (denom * denom));
    }

    Real invRet = 1.0 / Ret;
    Real invRet2 = invRet * invRet;
    for (int i = 0; i < 4; ++i)
    cDi_U[i] = num_Hk * Hk_U[i] * invRet - num * invRet2 * Ret_U[i];

    return num * invRet;
}

Real get_cDi_lamwake(const Real th,const Real ds,const Real sa,const Real ue, const bool turb, const bool wake, const Param& param,Real (&cDi_U)[4]) {


    //turb = false;
    Real Hk_U[4]={0}, Hs_U[4]={0}, Ret_U[4]={0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real Hs = get_Hs(th,ds,sa,ue,param,false,wake,Hs_U);
    Real Ret = get_Ret(th, ds, ue, param, Ret_U);

    Real invHk = 1.0 / Hk;
    Real oneMinusInvHk = 1.0 - invHk;
    Real num = 2.2 * oneMinusInvHk * oneMinusInvHk * invHk;

    Real invHk2 = invHk * invHk;
    Real num_Hk = 2.2 * (2.0 * oneMinusInvHk * invHk2 * invHk - oneMinusInvHk * oneMinusInvHk * invHk2);

    Real HsRet = Hs * Ret;
    Real invHsRet = 1.0 / HsRet;
    Real invHsRet2 = invHsRet * invHsRet;

    for (int i = 0; i < 4; ++i){
        cDi_U[i] = num_Hk * Hk_U[i] * invHsRet - num * invHsRet2 * (Hs_U[i] * Ret + Hs * Ret_U[i]);
    }

    return num * invHsRet;
}


Real get_cDi_outer(const Real th, const Real ds,const Real sa,const Real ue, const bool turb,const bool wake,const Param& param, Real (&cDi_U)[4]) {
    
    if (!turb) {
        for (int i = 0; i < 4; ++i) cDi_U[i] = 0.0;
        return 0.0;
    }

    Real Hs_U[4]={0}, Us_U[4]={0};
    Real Hs = get_Hs(th,ds,sa,ue,param,turb,wake,Hs_U);
    Real Us = get_Us(th,ds,sa,ue,param,turb,wake,Us_U);

    Real ct = sa * sa;
    Real ct_U[4] = {0.0, 0.0, 2.0 * sa, 0.0};
    Real factor = (0.995 - Us)*2.0;

    for (int i = 0; i < 4; ++i){
        cDi_U[i] = ct_U[i]*factor/Hs + ct*(-Us_U[i])*2/Hs - ct*factor/(Hs*Hs)*Hs_U[i];
    }

    return ct*factor/Hs;
}


Real get_cDi_lamstress(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Param& param, Real (&cDi_U)[4]) {
    
    Real Hs_U[4]={0}, Us_U[4]={0}, Ret_U[4]={0};
    Real Hs = get_Hs(th,ds,sa,ue,param,turb,wake,Hs_U);
    Real Us = get_Us(th,ds,sa,ue,param,turb,wake,Us_U);
    Real Ret = get_Ret(th, ds, ue, param, Ret_U);

    Real oneMinusUs = 0.995 - Us;
    Real num = 0.15 * oneMinusUs * oneMinusUs * 2.0;
    Real num_Us = -0.6 * oneMinusUs;

    Real HsRet = Hs * Ret;
    Real invHsRet = 1.0 / HsRet;
    Real invHsRet2 = invHsRet * invHsRet;

    for (int i = 0; i < 4; ++i){
        cDi_U[i] = num_Us * Us_U[i] * invHsRet - num * invHsRet2 * (Hs_U[i] * Ret + Hs * Ret_U[i]);
    }

    return num * invHsRet;
}



Real get_cDi(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Param& param,
    Real (&cDi_U)[4])
{

    if (turb){

        Real cDi = 0.0;
        Real cDil_U[4]={0};
        Real cDil = 0.0;


        if(!wake){

            Real cDi0_U[4] = {0};
            for (int i=0;i<4;++i){cDi_U[i] = 0.0;}
            Real cDi0 = get_cDi_turbwall(th,ds,sa,ue,turb,wake,param,cDi0_U);
            cDi += cDi0;

            for (int i=0;i<4;++i){cDi_U[i] += cDi0_U[i];}
            cDil = get_cDi_lam(th,ds,sa,ue,param,cDil_U);
        }
        else{
            cDil = get_cDi_lamwake(th,ds,sa,ue,turb,wake,param,cDil_U);
        }


        Real cDi0_Uouter[4] = {0};
        Real cDi0outer = get_cDi_outer(th,ds,sa,ue,turb,wake,param,cDi0_Uouter);

        cDi += cDi0outer;
        for (int i=0;i<4;++i){cDi_U[i] += cDi0_Uouter[i];}

        Real cDi0_Ulam[4] = {0};
        Real cDi0lam = get_cDi_lamstress(th,ds,sa,ue,turb,wake,param,cDi0_Ulam);

        cDi += cDi0lam;
        for (int i=0;i<4;++i){cDi_U[i] += cDi0_Ulam[i];}

        if (cDil >cDi){

            cDi = cDil;
            for (int i=0;i<4;++i){cDi_U[i] = cDil_U[i];}
        }

        if (wake){

            cDi *= 2;
            for (int i=0;i<4;++i){cDi_U[i] *= 2;}
        }
    
        return cDi;
    }
    else{

        Real cDi = get_cDi_lam(th,ds,sa,ue,param,cDi_U);
        return cDi;
    }
}


Real get_cDixt(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Real dist, const Param& param,
    Real (&cDixt_U)[4],Real& cDixt_x)
{

    Real cDi_U[4]={0};
    Real cDi = get_cDi(th,ds,sa,ue,turb,wake,param,cDi_U);

    Real temp = dist/th ;
    Real cDixt = cDi*temp ;

    for (int i=0;i<4;++i){cDixt_U[i] = cDi_U[i]*temp;}
    cDixt_U[0] -= cDixt/th; 

    cDixt_x = cDi/th;

    return cDixt;
}


Real get_upw(const Real th1,const Real ds1,const Real sa1,const Real ue1,
    const Real th2,const Real ds2,const Real sa2,const Real ue2, const bool wake, const Param&param,
    Real (&upw_U)[8]){

    Real Hk_U1[4] = {0};
    Real Hk1 = get_Hk(th1,ds1,ue1,param,Hk_U1);
    Real Hk_U2[4] = {0};
    Real Hk2 = get_Hk(th2,ds2,ue2,param,Hk_U2);

    Real Hut = 1.0;
    Real Z[4] = {0};
    Real C = 5.0;
    if (wake){C = 1.0;}

    Real Huc = C*Hut/(Hk2*Hk2); // only depends on U2
    
    Real Huc_U[8] = {0};
    Real temp = -2*Huc/Hk2 ;
    for (int i=4;i<8;++i){Huc_U[i] = temp*Hk_U2[i-4];}

    Real aa = (Hk2-1.)/(Hk1-1.);

    Real sign = 1.0;

    if (aa==0.0){sign = 0.0;}
    else if (aa<0.0){sign = -1.0 ;}

    Real la = std::log(sign*aa);

    Real la_U[8]={0} ;
    for (int i=0;i<4;++i){la_U[i] = -1./(Hk1-1.)*Hk_U1[i];}
    for (int i=4;i<8;++i){la_U[i] = 1./(Hk2-1.)*Hk_U2[i-4];}

    Real Hls = la*la;
    Real Hls_U[8]={0};
    temp = 2*la ;
    for (int i=0;i<8;++i){Hls_U[i] = temp*la_U[i];}

    if (Hls > 15){
        Hls =  15;
        for (int i=0;i<8;++i){Hls_U[i] = 0;}
    }

    temp = -0.5*std::exp(-Hls*Huc);
    Real upw = 1.0+temp;

    for (int i=0;i<8;++i){upw_U[i] = temp*(-Hls_U[i]*Huc-Hls*Huc_U[i]);}

    return upw ;
}


Real get_cteq(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const bool wake, const Param&param,Real (&cteq_U)[4])
{
    Real CC = 0.5/(param.GA*param.GA*param.GB), C = param.GC ;
    Real Hk_U[4] = {0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real Hs_U[4] = {0};
    Real Hs = get_Hs(th,ds,sa,ue,param,turb,wake,Hs_U);
    Real H_U[4]  = {0};
    Real H = get_H(th,ds,H_U);
    Real Ret_U[4] = {0};
    Real Ret = get_Ret(th,ds,ue,param,Ret_U);
    Real Us_U[4] = {0};
    Real Us = get_Us(th,ds,sa,ue,param,turb,wake,Us_U);

    Real Hkc_U[4]={0};
    Real Hkc ;
    if (wake){

        if (Hk<1.00005){
            
            Hk = 1.00005;
            for (int i=0; i<4;++i){Hk_U[i]=0;}
        }
        Hkc = Hk-1;
        for (int i=0; i<4;++i){Hkc_U[i]=Hk_U[i];}
    }
    else{

        if (Hk<1.05){
            Hk = 1.05;
            for (int i=0; i<4;++i){Hk_U[i]=0;}
        }
        Hkc = Hk-1 - C/Ret;
        for (int i=0; i<4;++i){Hkc_U[i]=Hk_U[i]+ C/(Ret*Ret)*Ret_U[i];}
        
        if (Hkc < 0.01){
            Hkc = 0.01;
            for (int i=0; i<4;++i){Hkc_U[i]=0;}
        }
    }

    Real num = CC*Hs*(Hk-1)*Hkc*Hkc;
    Real num_U[4]={0};
    for (int i=0; i<4;++i){
        num_U[i] = CC*(Hs_U[i]*(Hk-1)*Hkc*Hkc + Hs*Hk_U[i]*Hkc*Hkc + Hs*(Hk-1)*2*Hkc*Hkc_U[i]);
    }
    
    Real den = (1-Us)*H*Hk*Hk;
    Real den_U[4]={0};
    for (int i=0; i<4;++i){
        den_U[i] = (-Us_U[i])*H*Hk*Hk + (1-Us)*H_U[i]*Hk*Hk + (1-Us)*H*2*Hk*Hk_U[i];
    }
    Real cteq = std::sqrt(num/den);
    for (int i=0; i<4;++i){
        cteq_U[i] = 0.5/cteq*(num_U[i]/den - num/(den*den)*den_U[i]);
    }

    return cteq ;
}


Real get_damp(const Real th,const Real ds,const Real sa,const Real ue, const Param& param, Real (&damp_U)[4]) {
    
    
    
    Real Hk_U[4]={0}, Ret_U[4]={0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real Ret = get_Ret(th,ds,ue,param,Ret_U);

    // Default no amplification
    for (int i = 0; i < 4; ++i) {damp_U[i] = 0.0;}
    Real damp = 0.0;

    // Limit Hk
    if (Hk < 1.05) {
        Hk = 1.05;
        for (int i = 0; i < 4; ++i) Hk_U[i] = 0.0;
    }

    // Compute amplification
    Real Hmi = 1.0 / (Hk - 1.0);
    Real Hmi_U[4]={0};
    for (int i = 0; i < 4; ++i){
        Hmi_U[i] = -Hmi * Hmi * Hk_U[i];    
    }

    Real aa = 2.492 * std::pow(Hmi, 0.43);
    Real aa_U[4]={0};
    for (int i = 0; i < 4; ++i){
        aa_U[i] = 0.43 * aa / Hmi * Hmi_U[i];
    }
    Real bb = std::tanh(14.0 * Hmi - 9.24);
    Real bb_U[4]={0};
    for (int i = 0; i < 4; ++i){
        bb_U[i] = (1.0 - bb * bb) * 14.0 * Hmi_U[i];
    }
    Real lrc = aa + 0.7 * (bb + 1.0);
    Real lrc_U[4]={0};
    for (int i = 0; i < 4; ++i){
        lrc_U[i] = aa_U[i] + 0.7 * bb_U[i];
    }
    Real lr = std::log(Ret) / std::log(10.0);
    Real lr_U[4]={0};
    for (int i = 0; i < 4; ++i){
        lr_U[i] = Ret_U[i] / (Ret * std::log(10.0));
    }
    Real dl = 0.1;

    if (lr >= lrc - dl) {
        Real rn = (lr - (lrc - dl)) / (2.0 * dl);
        Real rn_U[4]={0};
        for (int i = 0; i < 4; ++i)
            rn_U[i] = (lr_U[i] - lrc_U[i]) / (2.0 * dl);

        Real rf, rf_U[4]={0};
        if (rn >= 1.0) {
            rf = 1.0;
            for (int i = 0; i < 4; ++i) rf_U[i] = 0.0;
        } else {
            rf = 3.0 * rn * rn - 2.0 * rn * rn * rn;
            for (int i = 0; i < 4; ++i)
                rf_U[i] = (6.0 * rn - 6.0 * rn * rn) * rn_U[i];
        }

        Real ar = 3.87 * Hmi - 2.52;
        Real ar_U[4]={0};
        for (int i = 0; i < 4; ++i)
            ar_U[i] = 3.87 * Hmi_U[i];

        Real ex = std::exp(-ar * ar);
        Real ex_U[4]={0};
        for (int i = 0; i < 4; ++i)
            ex_U[i] = ex * (-2.0 * ar * ar_U[i]);

        Real da = 0.028 * (Hk - 1.0) - 0.0345 * ex;
        Real da_U[4]={0};
        for (int i = 0; i < 4; ++i)
            da_U[i] = 0.028 * Hk_U[i] - 0.0345 * ex_U[i];

        Real af = -0.05 + 2.7 * Hmi - 5.5 * Hmi * Hmi + 3.0 * Hmi * Hmi * Hmi + 0.1 * std::exp(-20.0 * Hmi);
        Real af_U[4]={0};
        for (int i = 0; i < 4; ++i)
            af_U[i] = (2.7 - 11.0 * Hmi + 9.0 * Hmi * Hmi - std::exp(-20.0 * Hmi)) * Hmi_U[i];

        damp = rf * af * da / th;
        for (int i = 0; i < 4; ++i)
            damp_U[i] = (rf_U[i] * af * da + rf * af_U[i] * da + rf * af * da_U[i]) / th - (damp / th) * (i == 0 ? 1.0 : 0.0);
    }

    // Additional amplification
    Real ncrit = param.ncrit;
    Real nx = 5.0 * (sa - ncrit);
    Real eex = 1.0 + std::tanh(nx);
    Real ed = eex * 0.001 / th;

    Real nx_U[4] = {0.0, 0.0, 5.0, 0.0};
    Real eex_U[4]={0}, ed_U[4]={0};
    Real tanh_val = std::tanh(nx);
    for (int i = 0; i < 4; ++i) {
        eex_U[i] = (1.0 - tanh_val * tanh_val) * nx_U[i];
        ed_U[i] = eex_U[i] * 0.001 / th - ed / th * (i == 0 ? 1.0 : 0.0);
        damp_U[i] += ed_U[i];
    }

    damp += ed;
    return damp;
}


Real get_damp_forced(const Real th,const Real ds,const Real sa,const Real ue, const Param& param, const Real&ncrit, Real (&damp_U)[4]) {
    
    
    
    Real Hk_U[4]={0}, Ret_U[4]={0};
    Real Hk = get_Hk(th,ds,ue,param,Hk_U);
    Real Ret = get_Ret(th,ds,ue,param,Ret_U);

    // Default no amplification
    for (int i = 0; i < 4; ++i) {damp_U[i] = 0.0;}
    Real damp = 0.0;

    // Limit Hk
    if (Hk < 1.05) {
        Hk = 1.05;
        for (int i = 0; i < 4; ++i) Hk_U[i] = 0.0;
    }

    // Compute amplification
    Real Hmi = 1.0 / (Hk - 1.0);
    Real Hmi_U[4]={0};
    for (int i = 0; i < 4; ++i){
        Hmi_U[i] = -Hmi * Hmi * Hk_U[i];    
    }

    Real aa = 2.492 * std::pow(Hmi, 0.43);
    Real aa_U[4]={0};
    for (int i = 0; i < 4; ++i){
        aa_U[i] = 0.43 * aa / Hmi * Hmi_U[i];
    }
    Real bb = std::tanh(14.0 * Hmi - 9.24);
    Real bb_U[4]={0};
    for (int i = 0; i < 4; ++i){
        bb_U[i] = (1.0 - bb * bb) * 14.0 * Hmi_U[i];
    }
    Real lrc = aa + 0.7 * (bb + 1.0);
    Real lrc_U[4]={0};
    for (int i = 0; i < 4; ++i){
        lrc_U[i] = aa_U[i] + 0.7 * bb_U[i];
    }
    Real lr = std::log(Ret) / std::log(10.0);
    Real lr_U[4]={0};
    for (int i = 0; i < 4; ++i){
        lr_U[i] = Ret_U[i] / (Ret * std::log(10.0));
    }
    Real dl = 0.1;

    if (lr >= lrc - dl) {
        Real rn = (lr - (lrc - dl)) / (2.0 * dl);
        Real rn_U[4]={0};
        for (int i = 0; i < 4; ++i)
            rn_U[i] = (lr_U[i] - lrc_U[i]) / (2.0 * dl);

        Real rf, rf_U[4]={0};
        if (rn >= 1.0) {
            rf = 1.0;
            for (int i = 0; i < 4; ++i) rf_U[i] = 0.0;
        } else {
            rf = 3.0 * rn * rn - 2.0 * rn * rn * rn;
            for (int i = 0; i < 4; ++i)
                rf_U[i] = (6.0 * rn - 6.0 * rn * rn) * rn_U[i];
        }

        Real ar = 3.87 * Hmi - 2.52;
        Real ar_U[4]={0};
        for (int i = 0; i < 4; ++i)
            ar_U[i] = 3.87 * Hmi_U[i];

        Real ex = std::exp(-ar * ar);
        Real ex_U[4]={0};
        for (int i = 0; i < 4; ++i)
            ex_U[i] = ex * (-2.0 * ar * ar_U[i]);

        Real da = 0.028 * (Hk - 1.0) - 0.0345 * ex;
        Real da_U[4]={0};
        for (int i = 0; i < 4; ++i)
            da_U[i] = 0.028 * Hk_U[i] - 0.0345 * ex_U[i];

        Real af = -0.05 + 2.7 * Hmi - 5.5 * Hmi * Hmi + 3.0 * Hmi * Hmi * Hmi + 0.1 * std::exp(-20.0 * Hmi);
        Real af_U[4]={0};
        for (int i = 0; i < 4; ++i)
            af_U[i] = (2.7 - 11.0 * Hmi + 9.0 * Hmi * Hmi - std::exp(-20.0 * Hmi)) * Hmi_U[i];

        damp = rf * af * da / th;
        for (int i = 0; i < 4; ++i)
            damp_U[i] = (rf_U[i] * af * da + rf * af_U[i] * da + rf * af * da_U[i]) / th - (damp / th) * (i == 0 ? 1.0 : 0.0);
    }

    // Additional amplification
    //Real ncrit = param.ncrit;
    Real nx = 5.0 * (sa - ncrit);
    Real eex = 1.0 + std::tanh(nx);
    Real ed = eex * 0.001 / th;

    Real nx_U[4] = {0.0, 0.0, 5.0, 0.0};
    Real eex_U[4]={0}, ed_U[4]={0};
    Real tanh_val = std::tanh(nx);
    for (int i = 0; i < 4; ++i) {
        eex_U[i] = (1.0 - tanh_val * tanh_val) * nx_U[i];
        ed_U[i] = eex_U[i] * 0.001 / th - ed / th * (i == 0 ? 1.0 : 0.0);
        damp_U[i] += ed_U[i];
    }

    damp += ed;
    return damp;
}

Real get_cttr(const Real th,const Real ds,const Real sa,const Real ue,const bool turb,const Param&param,Real (&cttr_U)[4]){

    Real cteq,cteq_U[4]={0};
    Real Hk,Hk_U[4]={0};
    cteq = get_cteq(th,ds,sa,ue,turb,false,param,cteq_U);
    Hk = get_Hk(th,ds,ue,param,Hk_U);

    if (Hk<1.05){
        Hk=1.05;
        for (int i=0;i<4;++i){Hk_U[i]=0;};
    }
    
    Real C=param.CtauC, E=param.CtauE;

    Real c = C*std::exp(-E/(Hk-1.));
    Real c_U[4]={0};
    for (int i=0;i<4;++i){c_U[i] = c*E/((Hk-1)*(Hk-1))*Hk_U[i];}

    Real cttr = c*cteq;

    for (int i=0;i<4;++i){cttr_U[i] = c_U[i]*cteq + c*cteq_U[i];}

    return cttr;
}
