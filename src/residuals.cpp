#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "vector_ops.hpp"




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



void wake_sys(const Vsol& vsol, const Foil& foil, const Glob& glob, const Param& param,
    Real (&R)[3], Real (&R_U)[36], int (&J)[3]) {

    const Real* U = glob.U;

    // --- Get trailing edge states ---
    int il = 0;
    int iu = Ncoords-1;
    int iw = Ncoords;

    Real Ul[4]={0}, Uu[4]={0}, Uw[4]={0};
    for (int i = 0; i < 4; ++i) {
        Ul[i] = U[colMajorIndex(i,il,4)];
        Uu[i] = U[colMajorIndex(i,iu,4)];
        Uw[i] = U[colMajorIndex(i,iw,4)];
    }

    // --- Compute wake shear stress (ctw) ---
    Real ctl, ctu, ctl_Ul[4]={0}, ctu_Uu[4]={0};

    if (vsol.turb[il]) {
        ctl = Ul[2];
        ctl_Ul[0] = 0; ctl_Ul[1] = 0; ctl_Ul[2] = 1; ctl_Ul[3] = 0;
    } else {
        ctl = get_cttr(Ul[0],Ul[1],Ul[2],Ul[3],true,param,ctl_Ul);
    }

    if (vsol.turb[iu]) {
        ctu = Uu[2];
        ctu_Uu[2] = 1;
    } else {
        ctu = get_cttr(Uu[0],Uu[1],Uu[2],Uu[3],true,param,ctu_Uu);
    }

    Real thsum = Ul[0] + Uu[0];
    Real ctw = (ctl * Ul[0] + ctu * Uu[0]) / thsum;

    

    Real ctw_Ul[4]={0}, ctw_Uu[4]={0};

    Real zero = 0;
    for (int i = 0; i < 4; ++i) {
        ctw_Ul[i] = (ctl_Ul[i] * Ul[0] + (i == 0 ? ctl - ctw : zero)) / thsum;
        ctw_Uu[i] = (ctu_Uu[i] * Uu[0] + (i == 0 ? ctu - ctw : zero)) / thsum;
    }

    // --- Residual vector ---
    R[0] = Uw[0] - (Ul[0] + Uu[0]);
    R[1] = Uw[1] - (Ul[1] + Uu[1] + foil.te.hTE);
    R[2] = Uw[2] - ctw;

    
    // --- Jacobian block indices ---
    J[0] = il;
    J[1] = iu;
    J[2] = iw;

    // --- Residual Jacobian matrix: 3x12 ---
    // R_Ul (cols 0-3): -I (2 rows) and -ctw_Ul
    R_U[0] = -1.0;
    R_U[colMajorIndex(1,1,3)] = -1.0;
    for (int i = 0; i < 4; ++i){
        R_U[colMajorIndex(2,i,3)] = -ctw_Ul[i];
    }

    // R_Uu (cols 4-7): -I (2 rows) and -ctw_Uu
    R_U[colMajorIndex(0,4,3)] = -1.0;
    R_U[colMajorIndex(1,5,3)] = -1.0;
    for (int i = 0; i < 4; ++i){
        R_U[colMajorIndex(2,i+4,3)] = -ctw_Uu[i];
    }
    // R_Uw (cols 8-11): identity
    R_U[colMajorIndex(0,8,3)]  = 1.0;
    R_U[colMajorIndex(1,9,3)]  = 1.0;
    R_U[colMajorIndex(2,10,3)] = 1.0;
}

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
    const Param&param,
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
    cfxtm = get_cfxt(0.5*(th1+th2), 0.5*(ds1+ds2), 0.5*(sa1+sa2), 0.5*(ue1+ue2),
                     0.5*(x1+x2),turb,wake,param,cfxtm_U,cfxtm_x);
    
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
    const Param&param,
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
    cfxtm = get_cfxt(0.5*(th1+th2), 0.5*(ds1+ds2), 0.5*(sa1+sa2), 0.5*(ue1+ue2),
                     0.5*(x1+x2),turb,wake,param,cfxtm_U,cfxtm_x);
    
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


void residual_transition(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Real aux1,
    const Real aux2,
    const Param& param,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6]
) {
    const int nNewton = 20;
    const Real dx = x2 - x1;
    
    Real ncrit = param.ncrit;


    // Extract elements of BL States
    Real th1 = U1[0],th2 = U2[0];
    Real ds1 = U1[1],ds2 = U2[1];
    Real sa1 = U1[2],sa2 = U2[2];
    Real ue1 = U1[3],ue2 = U2[3];


    Real damp1, damp1_U1[4];
    
    damp1 = get_damp(th1,ds1,sa1,ue1,param,damp1_U1); // is constant in newton loop
    
    
    Real dampt, dampt_Ut[4]={0};
    Real Rxt, Rxt_xt, dxt;
    Real Ut[4]={0},Ut_xt[4]={0};
    Real w1,w2;
    Real xt = x1 + 0.5 * dx; // initial guess
    
    

    // Not forcing transition, iteratively finding transition location within the panel
    for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
        
        // weights
        w2 = (xt - x1)/dx, w1 =1-w2;

        // State at transition point xt 
        for (int i = 0; i < 4; ++i) {
            Ut[i] = w1*U1[i] + w2*U2[i];
            Ut_xt[i] = (U2[i]-U1[i])/dx;
        }

        Ut[2] = ncrit;   // amplification at transition
        Ut_xt[2] = 0.0;

        dampt = get_damp(Ut[0],Ut[1],Ut[2],Ut[3],param,dampt_Ut);
        dampt_Ut[2] = 0.0;

        Rxt = ncrit - sa1 - 0.5 * (damp1+dampt)*(xt-x1);
        Rxt_xt = -0.5 * (damp1 + dampt);
        for (int i = 0; i < 4; ++i){Rxt_xt -= 0.5*(xt-x1) * dampt_Ut[i]*Ut_xt[i];}

        dxt = -Rxt / Rxt_xt;
        Real dmax = 0.2 * dx * (1.1 - Real(iNewton) / nNewton);
        if (std::abs(dxt) > dmax){ dxt *= dmax / std::abs(dxt);}
        if (std::abs(Rxt) < 1e-10) break;

        xt += dxt;
    }
    
    
    Real Rxt_U[8] = {0};
    for (int i = 0; i < 4; ++i) {
        Rxt_U[i]   = -0.5 * (xt-x1) * (damp1_U1[i] + dampt_Ut[i]*w1);
        Rxt_U[i+4] = -0.5 * (xt-x1) * (dampt_Ut[i]*w2);
    }
    Rxt_U[2] -= 1.0;

    Real Ut_x1[4]={0}, Ut_x2[4]={0};
    for (int i = 0; i < 4; ++i) {
        Ut_x1[i] = (U2[i]-U1[i]) * (w2 - 1) / dx;
        Ut_x2[i] = (U2[i]-U1[i]) * (-w2) / dx;
    }
    
    
    // for natural transition the amp is fixed at Ut
    Ut_x1[2] = 0;
    Ut_x2[2] = 0;
    
    
    Real Rxt_x1 = 0.5 * (damp1 + dampt);    
    Real Rxt_x2 = 0.0;                    
    for (int i = 0; i < 4; ++i) {
        Rxt_x1 -= 0.5 * (xt - x1) * dampt_Ut[i] * Ut_x1[i];
        Rxt_x2 -= 0.5 * (xt - x1) * dampt_Ut[i] * Ut_x2[i];
    }

    // sensitivity of xt w.r.t. U,x from Rxt(xt,U,x) = 0 constraint
    Real xt_U[8]={0}, xt_U1[4]={0}, xt_U2[4]={0};
    Real xt_x1 = -Rxt_x1 / Rxt_xt;
    Real xt_x2 = -Rxt_x2 / Rxt_xt;
    
    for (int i = 0; i < 8; ++i){xt_U[i] = -Rxt_U[i] / Rxt_xt;}

    for (int i = 0; i < 4; ++i) {
        xt_U1[i] = xt_U[i];
        xt_U2[i] = xt_U[i+4];
    }

    for (int i = 0; i < 4; ++i) {
        Ut_x1[i] += Ut_xt[i] * xt_x1;
        Ut_x2[i] += Ut_xt[i] * xt_x2;
    }

    // include derivatives w.r.t. xt in Ut_x1 and Ut_x2
    Real Ut_U1[16]={0}, Ut_U2[16]={0};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            Ut_U1[i + 4*j] = (i == j ? w1 : 0.0) + (U2[i]-U1[i]) * xt_U1[j] / dx;  // the xt_U1 term goes to zero for fixed
            Ut_U2[i + 4*j] = (i == j ? w2 : 0.0) + (U2[i]-U1[i]) * xt_U2[j] / dx;
        }

    // Copy to Utl, Utt (laminar and turbulent transition states)
    Real Utl[4] = {0}, Utt[4] = {0}, Utl_x1[4] = {0}, Utl_x2[4] = {0}, Utt_x1[4] = {0}, Utt_x2[4] = {0};
    Real Utl_U1[16] = {0}, Utl_U2[16] = {0}, Utt_U1[16] = {0}, Utt_U2[16] = {0};

    
    for (int i = 0; i < 4; ++i) {
        Utl[i] = Ut[i]; Utt[i] = Ut[i];
        Utl_x1[i] = Ut_x1[i]; Utl_x2[i] = Ut_x2[i];
        Utt_x1[i] = Ut_x1[i]; Utt_x2[i] = Ut_x2[i];
        for (int j = 0; j < 4; ++j) {
            Utl_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utl_U2[i + 4*j] = Ut_U2[i + 4*j];
            Utt_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utt_U2[i + 4*j] = Ut_U2[i + 4*j];
        }
    }
    Utl[2] = ncrit;
    Utl_x1[2] = 0; Utl_x2[2] = 0;
    for (int j = 0; j < 4; ++j) {
        Utl_U1[2 + 4*j] = 0.0;
        Utl_U2[2 + 4*j] = 0.0;
    }

    Real cttr, cttr_Ut[4]={0};
    cttr = get_cttr(Ut[0],Ut[1],Ut[2],Ut[3],true,param,cttr_Ut);
    Utt[2] = cttr;
    for (int j = 0; j < 4; ++j) {
        Real sumU1 = 0.0, sumU2 = 0.0;
        for (int i = 0; i < 4; ++i) {
            sumU1 += cttr_Ut[i] * Ut_U1[i + 4*j];
            sumU2 += cttr_Ut[i] * Ut_U2[i + 4*j];
        }
        Utt_U1[2 + 4*j] = sumU1;
        Utt_U2[2 + 4*j] = sumU2;
    }
    Utt_x1[2] = 0.0;
    Utt_x2[2] = 0.0;
    for (int i = 0; i < 4; ++i) {
        Utt_x1[2] += cttr_Ut[i] * Ut_x1[i];
        Utt_x2[2] += cttr_Ut[i] * Ut_x2[i];
    }

    //param.turb = false;
    Real Rl[3]={0}, Rl_U[24]={0}, Rl_x[6]={0};
    residual_station(U1, Utl, x1, xt, aux1, aux2, false, false, false, param, Rl, Rl_U, Rl_x);
    
    //param.turb = true;
    Real Rt[3]={0}, Rt_U[24]={0}, Rt_x[6]={0};
    residual_station(Utt, U2, xt, x2, aux1, aux2, false, true, false, param, Rt, Rt_U, Rt_x);

    for (int i = 0; i < 3; ++i){
        R[i] = Rl[i] + Rt[i];
    }
    // R_U = dRl/dU1 + dRl/dU2 + dRt/dU1 + dRt/dU2
    // Compute contributions using chain rule with intermediate derivatives

    Real Rl_U1[12]={0}, Rl_Utl[12]={0};
    Real Rt_Utt[12]={0}, Rt_U2[12]={0};
    
    // do with pointers ??
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 4; ++c){
            Rl_U1[colMajorIndex(r,c,3)]    = Rl_U[colMajorIndex(r, c, 3)];
            Rl_Utl[colMajorIndex(r,c,3)]   = Rl_U[colMajorIndex(r, c+4, 3)];
            Rt_Utt[colMajorIndex(r,c,3)]   = Rt_U[colMajorIndex(r, c, 3)];
            Rt_U2[colMajorIndex(r,c,3)]    = Rt_U[colMajorIndex(r,c+4,3)];
        }
    }

    // R_x = dRl/dx1 + dRl/dx2 + dRt/dx1 + dRt/dx2
    for (int i = 0; i < 3; ++i) {
        R_x[i] = Rl_x[i] + Rt_x[i] * xt_x1;        // dR/dx1
        R_x[i+3] = Rl_x[i+3] + Rt_x[i] * xt_x2;    // dR/dx2
    }


    Real R_U1[12] = {0},R_U2[12] = {0},tmp1[12]={0};

    // Calculate R_U1
    cnp::add_inplace<12>(R_U1,Rl_U1);
    cnp::matmat_mul<3,4,4>(Rl_Utl, Utl_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    const Real* Rl_xCol1 = Rl_x + 3;  // Rl_x[:,1]
    cnp::outer_product<3,4>(Rl_xCol1,xt_U1,tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    cnp::outer_product<3,4>(Rt_x,xt_U1,tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    
    // Calculate R_U2
    cnp::add_inplace<12>(R_U2,Rt_U2);
    cnp::matmat_mul<3,4,4>(Rl_Utl,Utl_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    const Real* Rt_xCol1 = Rt_x + 3;  // Rl_x[:,1]
    cnp::outer_product<3,4>(Rl_xCol1,xt_U2,tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    cnp::outer_product<3,4>(Rt_x,xt_U2,tmp1); cnp::add_inplace<12>(R_U2,tmp1);

    cnp::hstack<12,12>(R_U1,R_U2,R_U); // R_U

    // do R_x: 

    Real R_x1[3]={0},R_x2[3]={0};
    Real tmp2[3]={0};

    cnp::add_inplace<3>(R_x1,Rl_x);
    cnp::scalar_mul<3>(Rl_xCol1,xt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::scalar_mul<3>(Rt_x,xt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);

    cnp::add_inplace<3>(R_x2,Rt_xCol1);
    cnp::scalar_mul<3>(Rl_xCol1,xt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::scalar_mul<3>(Rt_x,xt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);

    cnp::hstack<3,3>(R_x1,R_x2,R_x) ;
}


void residual_transition_forced(
    const Real* U1,
    const Real* U2,
    const Real x1,
    const Real x2,
    const Param& param,
    const Real transPos,
    Real (&R)[3],
    Real (&R_U)[24],
    Real (&R_x)[6]
) {

    const Real dx = x2 - x1;
    Real ncrit = param.ncrit;

    // Extract elements of BL States
    Real th1 = U1[0],th2 = U2[0];
    Real ds1 = U1[1],ds2 = U2[1];
    Real sa1 = U1[2],sa2 = U2[2];
    Real ue1 = U1[3],ue2 = U2[3];
    Real Ut[4]={0}; // state at transiton (initially laminar)
    

    

    // position is fixed, giving fixed weightings 
    Real xt = transPos;
   
    // Fixed weights and interpolated state
    Real w2 = (xt - x1) / dx;
    Real w1 = 1.0 - w2;


    for (int i = 0; i < 4; ++i){
        Ut[i] = w1*U1[i] + w2*U2[i];
    }
    
    // gradient of how transition states change wrt starting and ending locations (linear variation across panel)
    Real Ut_x1[4]={0}, Ut_x2[4]={0};
    for (int i = 0; i < 4; ++i) {
        Ut_x1[i] = (U2[i]-U1[i]) * (w2 - 1) / dx;
        Ut_x2[i] = (U2[i]-U1[i]) * (-w2) / dx;
    }
    

    // gradient of transition state wrt start and end states (again linear interp across)
    Real Ut_U1[16]={0}, Ut_U2[16]={0};
    for (int i = 0; i < 4; ++i) {
        int idx = colMajorIndex(i,i,4);  // Diagonal element in column-major
        Ut_U1[idx] = w1;
        Ut_U2[idx] = w2;
    }

    // Copy to Utl, Utt (laminar and turbulent transition states)
    Real Utl[4] = {0}, Utt[4] = {0}, Utl_x1[4] = {0}, Utl_x2[4] = {0}, Utt_x1[4] = {0}, Utt_x2[4] = {0};
    Real Utl_U1[16] = {0}, Utl_U2[16] = {0}, Utt_U1[16] = {0}, Utt_U2[16] = {0};

    // coptying lam/turb states at the transition point
    for (int i = 0; i < 4; ++i) {
        Utl[i] = Ut[i];
        Utt[i] = Ut[i];
        Utl_x1[i] = Ut_x1[i];
        Utl_x2[i] = Ut_x2[i];
        Utt_x1[i] = Ut_x1[i];
        Utt_x2[i] = Ut_x2[i];
        for (int j = 0; j < 4; ++j) {
            Utl_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utl_U2[i + 4*j] = Ut_U2[i + 4*j];
            Utt_U1[i + 4*j] = Ut_U1[i + 4*j];
            Utt_U2[i + 4*j] = Ut_U2[i + 4*j];
        }
    }
    

    // We dont care about the laminar amplicifcation growth at transition as we force it, so remove 
    // residuals for that in Utl
    Utl_x1[2] = 0.0;
    Utl_x2[2] = 0.0;
    for (int j = 0; j < 4; ++j){
        Utl_U1[colMajorIndex(2,j,4)] = 0.0;
        Utl_U2[colMajorIndex(2,j,4)] = 0.0;
    }

    // initialise the turb state with value of shear stress coeff (note Ut[2] has no effect)
    Real cttr, cttr_Ut[4]={0};
    cttr = get_cttr(Ut[0],Ut[1],Ut[2],Ut[3],true,param,cttr_Ut);
    Utt[2] = cttr;
    

    // doing that dUtt/dU1 [2,:] = dot(dcttr/dUt, dUt/dU1)
    // doing that dUtt/dU2 [2,:] = dot(dcttr/dUt, dUt/dU2)
    // but dUt/dU1 is just diaginal with value of w1 so code is simpler: 
    for (int j = 0; j < 4; ++j) {
        Utt_U1[colMajorIndex(2,j,4)] = cttr_Ut[j] * w1;
        Utt_U2[colMajorIndex(2,j,4)] = cttr_Ut[j] * w2;
    }
    
    

    // using dUtt/dx1 [2,:] = dot(dcttr/dUt, dUt/dx1)
    // using dUtt/dx2 [2,:] = dot(dcttr/dUt, dUt/dx2)
    Utt_x1[2] = 0.0;
    Utt_x2[2] = 0.0;
    for (int i = 0; i < 4; ++i) {
        Utt_x1[2] += cttr_Ut[i] * Ut_x1[i];
        Utt_x2[2] += cttr_Ut[i] * Ut_x2[i];
    }

    // residual for laminar section (remember that dont care about amplifcation residual here)
    Real Rl[3]={0}, Rl_U[24]={0}, Rl_x[6]={0};
    residual_station(U1, Utl, x1, xt, 0.0, 0.0, false, false, false, param, Rl, Rl_U, Rl_x);
    
    Rl[2] = 0.0;
    for (int i = 0; i < 8; ++i) {
        Rl_U[colMajorIndex(2,i,3)] = 0.0;
    }
     for (int i = 0; i < 2; ++i) {
        Rl_x[colMajorIndex(2,i,3)] = 0.0;
    }

    // residual for the turbulent section, here residual for shear stress matters
    Real Rt[3]={0}, Rt_U[24]={0}, Rt_x[6]={0};
    residual_station(Utt, U2, xt, x2, 0.0, 0.0, false, true, false, param, Rt, Rt_U, Rt_x);

    // sum residuals for total panel residuals
    for (int i = 0; i < 3; ++i){
        R[i] = Rl[i] + Rt[i];
    }

    Real Rl_U1[12]={0}, Rl_Utl[12]={0};
    Real Rt_Utt[12]={0}, Rt_U2[12]={0};
    
    // seperating Rl_U correctly and Rt_U
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 4; ++c){
            Rl_U1[colMajorIndex(r,c,3)]    = Rl_U[colMajorIndex(r, c, 3)];
            Rl_Utl[colMajorIndex(r,c,3)]   = Rl_U[colMajorIndex(r, c+4, 3)];
            Rt_Utt[colMajorIndex(r,c,3)]   = Rt_U[colMajorIndex(r, c, 3)];
            Rt_U2[colMajorIndex(r,c,3)]    = Rt_U[colMajorIndex(r,c+4,3)];
        }
    }

    // R_x = dRl/dx1 + dRl/dx2 + dRt/dx1 + dRt/dx2
    for (int i = 0; i < 3; ++i) {
        R_x[i] = Rl_x[i];       // dR/dx1
        R_x[i+3] = Rl_x[i+3];  // dR/dx2
    }


    Real R_U1[12] = {0},R_U2[12] = {0},tmp1[12]={0};

    // Calculate R_U1
    cnp::add_inplace<12>(R_U1,Rl_U1);
    cnp::matmat_mul<3,4,4>(Rl_Utl, Utl_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U1, tmp1); cnp::add_inplace<12>(R_U1,tmp1);
    
    // Calculate R_U2
    cnp::add_inplace<12>(R_U2,Rt_U2);
    cnp::matmat_mul<3,4,4>(Rl_Utl,Utl_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    cnp::matmat_mul<3,4,4>(Rt_Utt,Utt_U2, tmp1); cnp::add_inplace<12>(R_U2,tmp1);
    
    cnp::hstack<12,12>(R_U1,R_U2,R_U); // combine into single Residual jacobian

    // do R_x: 
    Real R_x1[3]={0},R_x2[3]={0};
    Real tmp2[3]={0};

    const Real* Rl_xCol1 = Rl_x + 3;  // Rl_x[:,1]
    const Real* Rt_xCol1 = Rt_x + 3;  // Rt_x[:,1]

    cnp::add_inplace<3>(R_x1,Rl_x); // add Rl_x[:,0]
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x1,tmp2); cnp::add_inplace<3>(R_x1,tmp2);

    cnp::add_inplace<3>(R_x2,Rt_xCol1);
    cnp::matmat_mul<3,4,1>(Rl_Utl,Utl_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);
    cnp::matmat_mul<3,4,1>(Rt_Utt,Utt_x2,tmp2); cnp::add_inplace<3>(R_x2,tmp2);

    cnp::hstack<3,3>(R_x1,R_x2,R_x) ;
}