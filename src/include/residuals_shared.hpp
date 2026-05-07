#pragma once

#include <cmath>
// Requires get_funcs.h (fwd) or get_funcs.hpp (AD) and vector_ops.hpp
// to be included before this header.

template<typename Real>
Real upwind(
    const Real upw,
    const Real (&upw_U)[8],
    const Real f1,
    const Real (&f1_U1)[4],
    const Real f2,
    const Real (&f2_U2)[4],
    Real (&f_U)[8])
{
    Real f = (1-upw)*f1 + upw*f2 ;

    Real temp[8]={0};
    for (int i=0;i<4;++i){temp[i] = (1-upw)*f1_U1[i];}
    for (int i=4;i<8;++i){temp[i] = upw*f2_U2[i-4];}

    for (int i=0;i<8;++i){f_U[i] = (-upw_U[i])*f1 + upw_U[i]*f2 + temp[i];}
    
    return f ;
}

template<typename Real>
Real upwind_half(
    
    // Sets upw to 0.5 and upw_U to 0 
    const Real f1,
    const Real (&f1_U1)[4],
    const Real f2,
    const Real (&f2_U2)[4],
    Real (&f_U)[8])
{
    Real f = 0.5*f1 + 0.5*f2 ;

    for (int i=0;i<4;++i){f_U[i] = 0.5*f1_U1[i];}
    for (int i=4;i<8;++i){f_U[i] = 0.5*f2_U2[i-4];}

    return f ;
}

template<bool ComputeJacobian=true, typename Real, typename ParamT>
void residual_station(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Real aux1,
    const Real aux2,
    const bool wake,
    const bool turb,
    const bool simi,
    const ParamT& param,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6])
{   

    // Extract elements of BL States
    Real th1 = U1[0],th2 = U2[0];
    Real ds1 = U1[1]-aux1,ds2 = U2[1]-aux2;
    Real sa1 = U1[2],sa2 = U2[2] ;
    Real ue1 = U1[3],ue2 = U2[3];

    // Compressibility correction on edge velocity
    Real uk1_u={0},uk2_u={0};
    Real uk1 = get_uk(ue1,param,uk1_u);
    Real uk2 = get_uk(ue2,param,uk2_u);

    // Log changes
    Real thlog = std::log(th2 / th1);
    Real thlog_U[8] = {0};
    thlog_U[0] = -1.0 / th1;
    thlog_U[4] =  1.0 / th2;

    Real uelog = std::log(uk2 / uk1);
    Real uelog_U[8] = {0};
    uelog_U[3] = -uk1_u/uk1;
    uelog_U[7] =  uk2_u/uk2;

    Real xlog = std::log(x2 / x1);
    Real xlog_x[2] = {-1/x1,1/x2} ;
    
    Real dx = x2-x1;
    Real dx_x[2] = {-1,1};

    // upwinding factor
    Real upw_U[8] = {0};
    Real upw = get_upw(th1,ds1,sa1,ue1,th2,ds2,sa2,ue2,wake,param,upw_U);

    // shape parameter
    Real H1_U1[4] = {0},H2_U2[4] = {0} ;
    Real H1 = get_H(th1,ds1,H1_U1), H2= get_H(th2,ds2,H2_U2);
    Real H = 0.5*(H1+H2);
    Real H_U[8] = {0.5*H1_U1[0], 0.5*H1_U1[1], 0, 0, 0.5*H2_U2[0], 0.5*H2_U2[1], 0, 0};

    // KE shape parameter, averaged across nodes
    Real Hs1_U1[4]={0},Hs2_U2[4]={0};
    Real Hs1 = get_Hs(th1,ds1,sa1,ue1,param,turb,wake,Hs1_U1), Hs2 = get_Hs(th2,ds2,sa2,ue2,param,turb,wake,Hs2_U2);
    Real Hs_U[8]= {0};
    Real Hs = upwind_half(Hs1,Hs1_U1,Hs2,Hs2_U2,Hs_U);

    // log changes in KE shape parameter
    Real Hslog = std::log(Hs2/Hs1);
    Real Hslog_U[8]={0};
    for (int i=0;i<4;++i){Hslog_U[i]= -1/Hs1*Hs1_U1[i];}
    for (int i=4;i<8;++i){Hslog_U[i]=  1/Hs2*Hs2_U2[i-4];}

    if (simi){

        thlog=0;
        Hslog=0;
        uelog=1;
        xlog =1;
        dx = 0.5*(x1+x2);

        for (int i=0;i<8;++i){
            thlog_U[i] = 0;
            Hslog_U[i] = 0;
            uelog_U[i] = 0;
        }
        xlog_x[0]=0; xlog_x[1]=0;
        dx_x[0]=0.5; dx_x[1]=0.5;
    }

    // wake shape parameter
    Real Hw1_U1[4]={0}, Hw2_U2[4]={0};
    Real Hw1 = get_Hw(th1,aux1,Hw1_U1), Hw2 = get_Hw(th2,aux2,Hw2_U2);
    Real Hw_U[8]={0};
    Real Hw = upwind_half(Hw1,Hw1_U1,Hw2,Hw2_U2,Hw_U);


    Real Rlag,Rlag_U[8]={0},Rlag_x[2]={0};
    if (turb){

        Real de1, de2, de;
        Real Us1, Us2, Us;
        Real Hk1, Hk2, Hk;
        Real Ret1, Ret2, Ret;
        Real cf1, cf2, cf;
        Real uq;

        Real de1_U1[4] = {0}, de2_U2[4] = {0}, de_U[8] = {0};
        Real Us1_U1[4] = {0}, Us2_U2[4] = {0}, Us_U[8] = {0};
        Real Hk1_U1[4] = {0}, Hk2_U2[4] = {0}, Hk_U[8] = {0};
        Real Ret1_U1[4] = {0}, Ret2_U2[4] = {0}, Ret_U[8] = {0};
        Real cf1_U1[4] = {0}, cf2_U2[4] = {0}, cf_U[8] = {0};
        Real uq_U[8]     = {0};


        Real salog = std::log(sa2/sa1);
        Real salog_U[8] ={0,0,-1./sa1,0, 0,0,1./sa2,0};

        // TODO: remove un-needed repeated calcs of Hk and other params etc
        // Be careful with modifying values if they are needed later in different form

        // --- BL thickness measure
        de1 = get_de(th1,ds1,ue1,param,de1_U1);
        de2 = get_de(th2,ds2,ue2,param,de2_U2);
        de  = upwind_half(de1, de1_U1, de2, de2_U2, de_U);

        // --- Normalized slip velocity
        Us1 = get_Us(th1,ds1,sa1,ue1,param,true,wake,Us1_U1);
        Us2 = get_Us(th2,ds2,sa2,ue2,param,true,wake,Us2_U2);
        Us  = upwind_half(Us1, Us1_U1, Us2, Us2_U2, Us_U);

        // --- Hk, upwinded
        Hk1 = get_Hk(th1,ds1,ue1,param,Hk1_U1);
        Hk2 = get_Hk(th2,ds2,ue2,param,Hk2_U2);
        Hk  = upwind(upw, upw_U, Hk1, Hk1_U1, Hk2, Hk2_U2, Hk_U);

        // --- Re_theta, averaged
        Ret1 = get_Ret(th1,ds1,ue1,param,Ret1_U1);
        Ret2 = get_Ret(th2,ds2,ue2,param,Ret2_U2);
        Ret  = upwind_half(Ret1, Ret1_U1, Ret2, Ret2_U2, Ret_U);

        // --- Skin friction, upwinded
        cf1 = get_cf(th1,ds1,sa1,ue1,true,wake,param,cf1_U1);
        cf2 = get_cf(th2,ds2,sa2,ue2,true,wake,param,cf2_U2);
        cf  = upwind(upw, upw_U, cf1, cf1_U1, cf2, cf2_U2, cf_U);


        // displacement thickness, averaged
        Real dsa = 0.5*(ds1 + ds2);
        Real dsa_U[8] ={0,0.5,0,0, 0,0.5,0,0};
        uq = get_uq(dsa, dsa_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, wake, param, uq_U);

        // cteq = root equilibrium wake layer shear coeficient: (ctau eq)^.5
        Real cteq1_U1[4]={0},cteq2_U2[4]={0};
        Real cteq1 = get_cteq(th1,ds1,sa1,ue1,turb,wake,param,cteq1_U1);
        Real cteq2 = get_cteq(th2,ds2,sa2,ue2,turb,wake,param,cteq2_U2);
        Real cteq_U[8] ;
        Real cteq = upwind(upw, upw_U, cteq1, cteq1_U1, cteq2, cteq2_U2,cteq_U);

        // root of shear coefficient (a state), upwinded
        Real f1_U[4] = {0,0,1,0};
        Real saa_U[8]={0};
        Real saa = upwind(upw,upw_U,sa1,f1_U,sa2,f1_U,saa_U);

        // Lag coeff
        Real Clag = param.SlagK/param.GB*1/(1+Us);
        Real Clag_U[8]={0};
        for (int i=0;i<8;++i){Clag_U[i] = -Clag/(1.+Us)*Us_U[i];}

        // extra dissipation in wake
        Real ald = 1.0;
        if (wake){ ald = param.Dlr;}


        // shear lag equation
        Rlag = Clag*(cteq-ald*saa)*dx - 2*de*salog + 2*de*(uq*dx-uelog)*param.Cuq;
        for (int i=0;i<8;++i){

            Rlag_U[i] = Clag_U[i]*(cteq-ald*saa)*dx + Clag*(cteq_U[i]-ald*saa_U[i])*dx \
            - 2*de_U[i]*salog - 2*de*salog_U[i] \
            + 2*de_U[i]*(uq*dx-uelog)*param.Cuq + 2*de*(uq_U[i]*dx-uelog_U[i])*param.Cuq ;
        }
        for (int i=0;i<2;++i){Rlag_x[i] = Clag*(cteq-ald*saa)*dx_x[i] + 2*de*uq*dx_x[i];}
    }
    else{
        // laminar, amplification factor equation

        if (simi){

            Rlag = sa2-sa1;
            Rlag_U[2] = 1, Rlag_U[6] = 1;
            //Rlag_x is already zeros
        }
        else{

            Real damp1_U1[4]={0}, damp2_U2[4]={0},damp_U[8]={0};

            Real damp1 = get_damp(th1,ds1,sa1,ue1,param,damp1_U1);
            Real damp2 = get_damp(th2,ds2,sa2,ue2,param,damp2_U2);
            Real damp = upwind_half(damp1,damp1_U1,damp2,damp2_U2,damp_U);

            Rlag = sa2-sa1 - damp*dx;
            Rlag_U[2] = -1,Rlag_U[6]=1 ;
            for (int i=0;i<8;++i){Rlag_U[i] -= damp_U[i]*dx;}
            for (int i=0;i<2;++i){Rlag_x[i] = -damp*dx_x[i];}
        }
    }

    Real Ms1, Ms2, Ms, Ms1_U[4]={0}, Ms2_U[4]={0}, Ms_U[8]={0};
    Ms1 = get_Mach2(ue1,param,Ms1_U);
    Ms2 = get_Mach2(ue2,param,Ms2_U);
    Ms = upwind_half(Ms1, Ms1_U, Ms2, Ms2_U,Ms_U);

    Real cfxt1, cfxt2, cfxtm;
    Real cfxt1_U[4]={0}, cfxt2_U[4]={0}, cfxtm_U[4]={0};
    Real cfxt1_x, cfxt2_x, cfxtm_x;
    cfxt1 = get_cfxt(th1,ds1,sa1,ue1,x1,turb,wake,param, cfxt1_U, cfxt1_x);
    cfxt2 = get_cfxt(th2,ds2,sa2,ue2,x2,turb,wake,param, cfxt2_U, cfxt2_x);
    Real th_avg = 0.5*(th1+th2), ds_avg = 0.5*(ds1+ds2);
    Real sa_avg = 0.5*(sa1+sa2), ue_avg = 0.5*(ue1+ue2);
    Real x_avg  = 0.5*(x1+x2);
    cfxtm = get_cfxt(th_avg, ds_avg, sa_avg, ue_avg, x_avg,
                     turb,wake,param,cfxtm_U,cfxtm_x);
    
    Real cfxt = 0.25 * cfxt1 + 0.5 * cfxtm + 0.25 * cfxt2;
    Real cfxt_U[8]={0}, cfxt_x[2]={0};
    for (int i = 0; i < 4; ++i) {
        cfxt_U[i]     = 0.25 * (cfxt1_U[i] + cfxtm_U[i]);
        cfxt_U[i + 4] = 0.25 * (cfxtm_U[i] + cfxt2_U[i]);
    }
    cfxt_x[0] = 0.25 * (cfxt1_x + cfxtm_x);
    cfxt_x[1] = 0.25 * (cfxtm_x + cfxt2_x);

    // Momentum residual
    Real Rmom = thlog + (2. + H + Hw - Ms) * uelog - 0.5 * xlog * cfxt;
    Real Rmom_U[8]={0},Rmom_x[2]={0};
    for (int i = 0; i < 8; ++i) {
        Rmom_U[i] = thlog_U[i] + (H_U[i]+Hw_U[i]-Ms_U[i])*uelog + 
                    (2+H+Hw-Ms)*uelog_U[i] - 0.5*xlog*cfxt_U[i];
    }
    Rmom_x[0] = -0.5 * xlog_x[0] * cfxt - 0.5 * xlog * cfxt_x[0];
    Rmom_x[1] = -0.5 * xlog_x[1] * cfxt - 0.5 * xlog * cfxt_x[1];

    // Dissipation cDi
    Real cDixt1, cDixt2, cDixt, cDixt1_U[4]={0}, cDixt2_U[4]={0}, cDixt_U[8]={0};
    Real cDixt1_x, cDixt2_x, cDixt_x[2]={0};
    cDixt1 = get_cDixt(th1,ds1,sa1,ue1,turb,wake,x1,param,cDixt1_U,cDixt1_x);
    cDixt2 = get_cDixt(th2,ds2,sa2,ue2,turb,wake,x2,param,cDixt2_U,cDixt2_x);
    cDixt = upwind(upw, upw_U, cDixt1, cDixt1_U, cDixt2, cDixt2_U, cDixt_U);
    cDixt_x[0] = (1. - upw) * cDixt1_x;
    cDixt_x[1] = upw * cDixt2_x;

    // cfxt upwinded
    Real cfxtu, cfxtu_U[8]={0}, cfxtu_x[2]={0};
    cfxtu = upwind(upw, upw_U, cfxt1, cfxt1_U, cfxt2, cfxt2_U, cfxtu_U);
    cfxtu_x[0] = (1. - upw) * cfxt1_x;
    cfxtu_x[1] = upw * cfxt2_x;

    // Hss
    Real Hss1, Hss2, Hss, Hss1_U[4]={0}, Hss2_U[4]={0}, Hss_U[8]={0};
    Hss1 = get_Hss(th1,ds1,ue1,param,Hss1_U);
    Hss2 = get_Hss(th2,ds2,ue2,param,Hss2_U);
    Hss = upwind_half(Hss1, Hss1_U, Hss2, Hss2_U, Hss_U);

    // Shape residual
    Real Rshape = Hslog + (2. * Hss / Hs+1.0-H-Hw)*uelog + xlog*(0.5*cfxtu - cDixt);
    Real Rshape_U[8]={0},Rshape_x[2]={0};

    for (int i = 0; i < 8; ++i) {
        Rshape_U[i] = Hslog_U[i] +
                    (2. * Hss_U[i] / Hs - 2. * Hss / (Hs * Hs) * Hs_U[i] -
                     H_U[i] - Hw_U[i]) *
                        uelog +
                    (2. * Hss / Hs + 1. - H - Hw) * uelog_U[i] +
                    xlog * (0.5 * cfxtu_U[i] - cDixt_U[i]);
    }
    Rshape_x[0] = xlog_x[0] * (0.5 * cfxtu - cDixt) + xlog * (0.5 * cfxtu_x[0] - cDixt_x[0]);
    Rshape_x[1] = xlog_x[1] * (0.5 * cfxtu - cDixt) + xlog * (0.5 * cfxtu_x[1] - cDixt_x[1]);

    // Final residual assembly
    R[0] = Rmom;
    R[1] = Rshape;
    R[2] = Rlag;
    if constexpr (ComputeJacobian) {
        
        int step = 0;
        for (int col = 0; col < 8; ++col){
            R_U[step] = Rmom_U[col];
            step += 3;
        }
        step = 0;
        for (int col = 0; col < 8; ++col){
            R_U[step+1] = Rshape_U[col];
            step += 3;
        }
        step = 0;
        for (int col = 0; col < 8; ++col){
            R_U[step+2] = Rlag_U[col];
            step += 3;
        }
        step = 0;
        
        R_x[0]=Rmom_x[0],  R_x[3]=Rmom_x[1];
        R_x[1]=Rshape_x[0],R_x[4]=Rshape_x[1];
        R_x[2]=Rlag_x[0],  R_x[5]=Rlag_x[1];
    }
}

// Convenience overload: no Jacobian output (used by the AD solver path).
template<typename Real, typename ParamT>
void residual_station(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Real aux1, const Real aux2,
    const bool wake, const bool turb, const bool simi,
    const ParamT& param,
    Real (&R)[3])
{
    Real _dummy_U[24]={0}, _dummy_x[6]={0};
    residual_station<false>(U1, U2, x1, x2, aux1, aux2, wake, turb, simi,
                            param, R, _dummy_U, _dummy_x);
}

template<bool ComputeJacobian=true, typename Real, typename ParamT>
void residual_station_forced(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Real aux1,
    const Real aux2,
    const bool wake,
    const bool turb,
    const bool simi,
    const ParamT& param,
    const Real& ncrit,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6])
{   

    // Extract elements of BL States
    Real th1 = U1[0],th2 = U2[0];
    Real ds1 = U1[1]-aux1,ds2 = U2[1]-aux2;
    Real sa1 = U1[2],sa2 = U2[2] ;
    Real ue1 = U1[3],ue2 = U2[3];

    // Compressibility correction on edge velocity
    Real uk1_u={0},uk2_u={0};
    Real uk1 = get_uk(ue1,param,uk1_u);
    Real uk2 = get_uk(ue2,param,uk2_u);

    // Log changes
    Real thlog = std::log(th2 / th1);
    Real thlog_U[8] = {0};
    thlog_U[0] = -1.0 / th1;
    thlog_U[4] =  1.0 / th2;

    Real uelog = std::log(uk2 / uk1);
    Real uelog_U[8] = {0};
    uelog_U[3] = -uk1_u/uk1;
    uelog_U[7] =  uk2_u/uk2;

    Real xlog = std::log(x2 / x1);
    Real xlog_x[2] = {-1/x1,1/x2} ;
    
    Real dx = x2-x1;
    Real dx_x[2] = {-1,1};

    // upwinding factor
    Real upw_U[8] = {0};
    Real upw = get_upw(th1,ds1,sa1,ue1,th2,ds2,sa2,ue2,wake,param,upw_U);

    // shape parameter
    Real H1_U1[4] = {0},H2_U2[4] = {0} ;
    Real H1 = get_H(th1,ds1,H1_U1), H2= get_H(th2,ds2,H2_U2);
    Real H = 0.5*(H1+H2);
    Real H_U[8] = {0.5*H1_U1[0], 0.5*H1_U1[1], 0, 0, 0.5*H2_U2[0], 0.5*H2_U2[1], 0, 0};

    // KE shape parameter, averaged across nodes
    Real Hs1_U1[4]={0},Hs2_U2[4]={0};
    Real Hs1 = get_Hs(th1,ds1,sa1,ue1,param,turb,wake,Hs1_U1), Hs2 = get_Hs(th2,ds2,sa2,ue2,param,turb,wake,Hs2_U2);
    Real Hs_U[8]= {0};
    Real Hs = upwind_half(Hs1,Hs1_U1,Hs2,Hs2_U2,Hs_U);

    // log changes in KE shape parameter
    Real Hslog = std::log(Hs2/Hs1);
    Real Hslog_U[8]={0};
    for (int i=0;i<4;++i){Hslog_U[i]= -1/Hs1*Hs1_U1[i];}
    for (int i=4;i<8;++i){Hslog_U[i]=  1/Hs2*Hs2_U2[i-4];}

    if (simi){

        thlog=0;
        Hslog=0;
        uelog=1;
        xlog =1;
        dx = 0.5*(x1+x2);

        for (int i=0;i<8;++i){
            thlog_U[i] = 0;
            Hslog_U[i] = 0;
            uelog_U[i] = 0;
        }
        xlog_x[0]=0; xlog_x[1]=0;
        dx_x[0]=0.5; dx_x[1]=0.5;
    }

    // wake shape parameter
    Real Hw1_U1[4]={0}, Hw2_U2[4]={0};
    Real Hw1 = get_Hw(th1,aux1,Hw1_U1), Hw2 = get_Hw(th2,aux2,Hw2_U2);
    Real Hw_U[8]={0};
    Real Hw = upwind_half(Hw1,Hw1_U1,Hw2,Hw2_U2,Hw_U);


    Real Rlag,Rlag_U[8]={0},Rlag_x[2]={0};
    if (turb){

        Real de1, de2, de;
        Real Us1, Us2, Us;
        Real Hk1, Hk2, Hk;
        Real Ret1, Ret2, Ret;
        Real cf1, cf2, cf;
        Real uq;

        Real de1_U1[4] = {0}, de2_U2[4] = {0}, de_U[8] = {0};
        Real Us1_U1[4] = {0}, Us2_U2[4] = {0}, Us_U[8] = {0};
        Real Hk1_U1[4] = {0}, Hk2_U2[4] = {0}, Hk_U[8] = {0};
        Real Ret1_U1[4] = {0}, Ret2_U2[4] = {0}, Ret_U[8] = {0};
        Real cf1_U1[4] = {0}, cf2_U2[4] = {0}, cf_U[8] = {0};
        Real uq_U[8]     = {0};


        Real salog = std::log(sa2/sa1);
        Real salog_U[8] ={0,0,-1./sa1,0, 0,0,1./sa2,0};

        // TODO: remove un-needed repeated calcs of Hk and other params etc
        // Be careful with modifying values if they are needed later in different form

        // --- BL thickness measure
        de1 = get_de(th1,ds1,ue1,param,de1_U1);
        de2 = get_de(th2,ds2,ue2,param,de2_U2);
        de  = upwind_half(de1, de1_U1, de2, de2_U2, de_U);

        // --- Normalized slip velocity
        Us1 = get_Us(th1,ds1,sa1,ue1,param,true,wake,Us1_U1);
        Us2 = get_Us(th2,ds2,sa2,ue2,param,true,wake,Us2_U2);
        Us  = upwind_half(Us1, Us1_U1, Us2, Us2_U2, Us_U);

        // --- Hk, upwinded
        Hk1 = get_Hk(th1,ds1,ue1,param,Hk1_U1);
        Hk2 = get_Hk(th2,ds2,ue2,param,Hk2_U2);
        Hk  = upwind(upw, upw_U, Hk1, Hk1_U1, Hk2, Hk2_U2, Hk_U);

        // --- Re_theta, averaged
        Ret1 = get_Ret(th1,ds1,ue1,param,Ret1_U1);
        Ret2 = get_Ret(th2,ds2,ue2,param,Ret2_U2);
        Ret  = upwind_half(Ret1, Ret1_U1, Ret2, Ret2_U2, Ret_U);

        // --- Skin friction, upwinded
        cf1 = get_cf(th1,ds1,sa1,ue1,true,wake,param,cf1_U1);
        cf2 = get_cf(th2,ds2,sa2,ue2,true,wake,param,cf2_U2);
        cf  = upwind(upw, upw_U, cf1, cf1_U1, cf2, cf2_U2, cf_U);


        // displacement thickness, averaged
        Real dsa = 0.5*(ds1 + ds2);
        Real dsa_U[8] ={0,0.5,0,0, 0,0.5,0,0};
        uq = get_uq(dsa, dsa_U, cf, cf_U, Hk, Hk_U, Ret, Ret_U, wake, param, uq_U);

        // cteq = root equilibrium wake layer shear coeficient: (ctau eq)^.5
        Real cteq1_U1[4]={0},cteq2_U2[4]={0};
        Real cteq1 = get_cteq(th1,ds1,sa1,ue1,turb,wake,param,cteq1_U1);
        Real cteq2 = get_cteq(th2,ds2,sa2,ue2,turb,wake,param,cteq2_U2);
        Real cteq_U[8] ;
        Real cteq = upwind(upw, upw_U, cteq1, cteq1_U1, cteq2, cteq2_U2,cteq_U);

        // root of shear coefficient (a state), upwinded
        Real f1_U[4] = {0,0,1,0};
        Real saa_U[8]={0};
        Real saa = upwind(upw,upw_U,sa1,f1_U,sa2,f1_U,saa_U);

        // Lag coeff
        Real Clag = param.SlagK/param.GB*1/(1+Us);
        Real Clag_U[8]={0};
        for (int i=0;i<8;++i){Clag_U[i] = -Clag/(1.+Us)*Us_U[i];}

        // extra dissipation in wake
        Real ald = 1.0;
        if (wake){ ald = param.Dlr;}


        // shear lag equation
        Rlag = Clag*(cteq-ald*saa)*dx - 2*de*salog + 2*de*(uq*dx-uelog)*param.Cuq;
        for (int i=0;i<8;++i){

            Rlag_U[i] = Clag_U[i]*(cteq-ald*saa)*dx + Clag*(cteq_U[i]-ald*saa_U[i])*dx \
            - 2*de_U[i]*salog - 2*de*salog_U[i] \
            + 2*de_U[i]*(uq*dx-uelog)*param.Cuq + 2*de*(uq_U[i]*dx-uelog_U[i])*param.Cuq ;
        }
        for (int i=0;i<2;++i){Rlag_x[i] = Clag*(cteq-ald*saa)*dx_x[i] + 2*de*uq*dx_x[i];}
    }
    else{
        // laminar, amplification factor equation

        if (simi){

            Rlag = sa2-sa1;
            Rlag_U[2] = 1, Rlag_U[6] = 1;
            //Rlag_x is already zeros
        }
        else{

            Real damp1_U1[4]={0}, damp2_U2[4]={0},damp_U[8]={0};

            Real damp1 = get_damp_forced(th1,ds1,sa1,ue1,param,ncrit,damp1_U1);
            Real damp2 = get_damp_forced(th2,ds2,sa2,ue2,param,ncrit,damp2_U2);
            Real damp = upwind_half(damp1,damp1_U1,damp2,damp2_U2,damp_U);

            Rlag = sa2-sa1 - damp*dx;
            Rlag_U[2] = -1,Rlag_U[6]=1 ;
            for (int i=0;i<8;++i){Rlag_U[i] -= damp_U[i]*dx;}
            for (int i=0;i<2;++i){Rlag_x[i] = -damp*dx_x[i];}
        }
    }

    Real Ms1, Ms2, Ms, Ms1_U[4]={0}, Ms2_U[4]={0}, Ms_U[8]={0};
    Ms1 = get_Mach2(ue1,param,Ms1_U);
    Ms2 = get_Mach2(ue2,param,Ms2_U);
    Ms = upwind_half(Ms1, Ms1_U, Ms2, Ms2_U,Ms_U);

    Real cfxt1, cfxt2, cfxtm;
    Real cfxt1_U[4]={0}, cfxt2_U[4]={0}, cfxtm_U[4]={0};
    Real cfxt1_x, cfxt2_x, cfxtm_x;
    cfxt1 = get_cfxt(th1,ds1,sa1,ue1,x1,turb,wake,param, cfxt1_U, cfxt1_x);
    cfxt2 = get_cfxt(th2,ds2,sa2,ue2,x2,turb,wake,param, cfxt2_U, cfxt2_x);
    Real th_avg = 0.5*(th1+th2), ds_avg = 0.5*(ds1+ds2);
    Real sa_avg = 0.5*(sa1+sa2), ue_avg = 0.5*(ue1+ue2);
    Real x_avg  = 0.5*(x1+x2);
    cfxtm = get_cfxt(th_avg, ds_avg, sa_avg, ue_avg, x_avg,
                     turb,wake,param,cfxtm_U,cfxtm_x);
    
    Real cfxt = 0.25 * cfxt1 + 0.5 * cfxtm + 0.25 * cfxt2;
    Real cfxt_U[8]={0}, cfxt_x[2]={0};
    for (int i = 0; i < 4; ++i) {
        cfxt_U[i]     = 0.25 * (cfxt1_U[i] + cfxtm_U[i]);
        cfxt_U[i + 4] = 0.25 * (cfxtm_U[i] + cfxt2_U[i]);
    }
    cfxt_x[0] = 0.25 * (cfxt1_x + cfxtm_x);
    cfxt_x[1] = 0.25 * (cfxtm_x + cfxt2_x);

    // Momentum residual
    Real Rmom = thlog + (2. + H + Hw - Ms) * uelog - 0.5 * xlog * cfxt;
    Real Rmom_U[8]={0},Rmom_x[2]={0};
    for (int i = 0; i < 8; ++i) {
        Rmom_U[i] = thlog_U[i] + (H_U[i]+Hw_U[i]-Ms_U[i])*uelog + 
                    (2+H+Hw-Ms)*uelog_U[i] - 0.5*xlog*cfxt_U[i];
    }
    Rmom_x[0] = -0.5 * xlog_x[0] * cfxt - 0.5 * xlog * cfxt_x[0];
    Rmom_x[1] = -0.5 * xlog_x[1] * cfxt - 0.5 * xlog * cfxt_x[1];

    // Dissipation cDi
    Real cDixt1, cDixt2, cDixt, cDixt1_U[4]={0}, cDixt2_U[4]={0}, cDixt_U[8]={0};
    Real cDixt1_x, cDixt2_x, cDixt_x[2]={0};
    cDixt1 = get_cDixt(th1,ds1,sa1,ue1,turb,wake,x1,param,cDixt1_U,cDixt1_x);
    cDixt2 = get_cDixt(th2,ds2,sa2,ue2,turb,wake,x2,param,cDixt2_U,cDixt2_x);
    cDixt = upwind(upw, upw_U, cDixt1, cDixt1_U, cDixt2, cDixt2_U, cDixt_U);
    cDixt_x[0] = (1. - upw) * cDixt1_x;
    cDixt_x[1] = upw * cDixt2_x;

    // cfxt upwinded
    Real cfxtu, cfxtu_U[8]={0}, cfxtu_x[2]={0};
    cfxtu = upwind(upw, upw_U, cfxt1, cfxt1_U, cfxt2, cfxt2_U, cfxtu_U);
    cfxtu_x[0] = (1. - upw) * cfxt1_x;
    cfxtu_x[1] = upw * cfxt2_x;

    // Hss
    Real Hss1, Hss2, Hss, Hss1_U[4]={0}, Hss2_U[4]={0}, Hss_U[8]={0};
    Hss1 = get_Hss(th1,ds1,ue1,param,Hss1_U);
    Hss2 = get_Hss(th2,ds2,ue2,param,Hss2_U);
    Hss = upwind_half(Hss1, Hss1_U, Hss2, Hss2_U, Hss_U);

    // Shape residual
    Real Rshape = Hslog + (2. * Hss / Hs+1.0-H-Hw)*uelog + xlog*(0.5*cfxtu - cDixt);
    Real Rshape_U[8]={0},Rshape_x[2]={0};

    for (int i = 0; i < 8; ++i) {
        Rshape_U[i] = Hslog_U[i] +
                    (2. * Hss_U[i] / Hs - 2. * Hss / (Hs * Hs) * Hs_U[i] -
                     H_U[i] - Hw_U[i]) *
                        uelog +
                    (2. * Hss / Hs + 1. - H - Hw) * uelog_U[i] +
                    xlog * (0.5 * cfxtu_U[i] - cDixt_U[i]);
    }
    Rshape_x[0] = xlog_x[0] * (0.5 * cfxtu - cDixt) + xlog * (0.5 * cfxtu_x[0] - cDixt_x[0]);
    Rshape_x[1] = xlog_x[1] * (0.5 * cfxtu - cDixt) + xlog * (0.5 * cfxtu_x[1] - cDixt_x[1]);

    // Final residual assembly
    R[0] = Rmom;
    R[1] = Rshape;
    R[2] = Rlag;
    
    if constexpr (ComputeJacobian) {
        int step = 0;
        for (int col = 0; col < 8; ++col){
            R_U[step] = Rmom_U[col];
            step += 3;
        }
        step = 0;
        for (int col = 0; col < 8; ++col){
            R_U[step+1] = Rshape_U[col];
            step += 3;
        }
        step = 0;
        for (int col = 0; col < 8; ++col){
            R_U[step+2] = Rlag_U[col];
            step += 3;
        }
        step = 0;
        
        R_x[0]=Rmom_x[0],  R_x[3]=Rmom_x[1];
        R_x[1]=Rshape_x[0],R_x[4]=Rshape_x[1];
        R_x[2]=Rlag_x[0],  R_x[5]=Rlag_x[1];
    }
}

// Convenience overload: no Jacobian output (used by the AD solver path).
template<typename Real, typename ParamT>
void residual_station_forced(
    const Real* U1, const Real* U2,
    const Real x1, const Real x2,
    const Real aux1, const Real aux2,
    const bool wake, const bool turb, const bool simi,
    const ParamT& param,
    const Real& ncrit,
    Real (&R)[3])
{
    Real _dummy_U[24]={0}, _dummy_x[6]={0};
    residual_station_forced<false>(U1, U2, x1, x2, aux1, aux2, wake, turb, simi,
                                   param, ncrit, R, _dummy_U, _dummy_x);
}
