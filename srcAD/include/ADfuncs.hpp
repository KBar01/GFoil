#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.hpp"
#include "data_structs.hpp"
#include "spline.hpp"
#include "solve_inv.hpp"
#include "panel_funcs.hpp"
#include "calc_ue_m.hpp"
#include "main_func.hpp"

#include "extract_BL_TE.hpp"
#include "sound.hpp"
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cassert>
#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;


template<typename Real>
double partialOutputspartialInputs(

    // geometry parameters
    const Real nCrit, const Real Ufac, const Real TEfac, const Real chordScaling, 
    const double (&inXcoords)[Nin], const Real Re, const Real Ma, const Real rhoInf, 
    const Real kinViscInf,
    const std::string model, const Real sampleTE, const Real sampleTE_hi,
    const Real* obsX, const Real* obsY, const Real* obsZ, int nObs,
    const Real S,
    
    // inputs x ---------------------------------------------------
    Real (&inYcoords)[Nin],
    Real alphad,
    // state vector omega (or U in this case) ---------------------
    Real (&states)[RVdimension],
    const int (&turb)[Ncoords+Nwake],

    // Outputs ----------------------------------------------------
    double (&jacobianCL_y)[Nin], 
    double (&jacobianOASPL_y)[Nin],
    double& jacobianCL_alf,
    double& jacobianOASPL_alf,
    double (&jacobianCL_states)[RVdimension],
    double (&jacobianOASPL_states)[RVdimension],
    int aWeighting = 0,
    double f_min = 200.0,
    double f_max = 20000.0
    ){
    
    using Tape = typename Real::Tape;
    // set up AD inputs (ycoords and alphad)
    Tape& tape = Real::getTape();
    tape.setActive();

    for (int i = 0; i < Nin; ++i) {
        tape.registerInput(inYcoords[i]);
    }
    tape.registerInput(alphad);

    for (int i = 0; i < RVdimension; ++i) {
        tape.registerInput(states[i]);
    }

    Real alpha = (alphad/180)*M_PI;
    Oper oper(alpha,Re,Ma);
    oper.rho = rhoInf;
    Geom<Real> geom;
    Real flattenedCoords[2*Ncoords]={0};
    Real inCoords[2*Nin]={0};
    for (int i=0;i<Nin;++i){
        inCoords[colMajorIndex(0,i,2)] = inXcoords[i];
        inCoords[colMajorIndex(1,i,2)] = inYcoords[i];
    }
    make_panels(inCoords,flattenedCoords,Ufac,TEfac); // does spline to redist nodes over aerofoil for fixed number of 200 nodes

    Foil foil(flattenedCoords);

    Param<Real> param;
    param.ncrit = nCrit;

    init_thermo(oper,param,geom);

    Glob<Real> glob;
    for (int i=0;i<RVdimension;++i){
        glob.U[i] = states[i] ;
    }


    Post<Real> post;
    calc_force(oper,geom,param,foil,glob,post);

    const Real Uinf = (Re*kinViscInf)/(chordScaling) ;

    Real topsurf[7],botsurf[7];
    Real xcoords[Ncoords]={0};
    for (int i=0;i<Ncoords;++i){
        xcoords[i] = flattenedCoords[colMajorIndex(0,i,2)];
    }

    interpolate_at_95_both_surfaces(xcoords,glob.U,post.cp,oper,turb,param,topsurf,botsurf,Uinf,sampleTE,sampleTE_hi,chordScaling);
    Real OASPL = calc_OASPL<Real>(botsurf,topsurf,chordScaling,Uinf,obsX,obsY,obsZ,nObs,S,kinViscInf,rhoInf,model,f_min,f_max,aWeighting,alpha);

    Real outputs[2] = {post.cl,OASPL} ;

    constexpr int jacobianHeight = 2;
    for (int i=0;i<jacobianHeight;++i){
        tape.registerOutput(outputs[i]);
    }
    tape.setPassive();
    for (int i=0;i<jacobianHeight;++i){outputs[i].gradient()[i] = 1.0 ;}
    tape.evaluate();

    for (int i = 0; i < Nin; ++i) {
        jacobianCL_y[i] = (inYcoords[i].getGradient()[0]);
        jacobianOASPL_y[i] = (inYcoords[i].getGradient()[1]);
    }

    jacobianCL_alf = (alphad.getGradient()[0]);
    jacobianOASPL_alf = (alphad.getGradient()[1]);

    for (int i = 0; i < RVdimension; ++i) {
        jacobianCL_states[i] = (states[i].getGradient()[0]);
        jacobianOASPL_states[i] = (states[i].getGradient()[1]);
    }

    tape.resetHard();
    return post.cl.getValue();

};

template<typename Real>
void partialRpartialx(


    // geometry parameters
    const Real nCrit, const Real Ufac, const Real TEfac,
    const double (&inXcoords)[Nin], const Real Re, const Real Ma, const Real rhoInf,

    int (&currStag)[2],
    // forced transition x/c per surface (0=lower, 1=upper; 1.0 = no forcing)
    const double xft_xc[2],
    
    // ADJOINT VECTOR size  4*(Ncoords+Nwake)
    const double (&adlambdaCL)[RVdimension],
    const double (&adlambdaCD)[RVdimension],
    const double (&adlambdaOASPL)[RVdimension],
    // inputs x 
    Real (&inYcoords)[Nin],
    Real alphad,
    // state vector omega (or U in this case)
    const double (&states)[RVdimension],
    const int (&turb)[Ncoords+Nwake],

    // OUTPUT — gradient of h wrt geometry coords + alpha
    double (&dgCLdy)[Nin], double (&dgCDdy)[Nin], double (&dgOASPLdy)[Nin],
    double& dgCLdalpha, double& dgCDdalpha, double& dgOASPLdalpha
    ){
    
    using Tape = typename Real::Tape;
    Tape& tape = Real::getTape();
    tape.setActive();

    for (int i = 0; i < Nin; ++i){
        tape.registerInput(inYcoords[i]);
    }
    tape.registerInput(alphad);

    Real alpha = (alphad/180)*M_PI;
    Oper<Real> oper(alpha,Re,Ma);
    oper.rho = rhoInf;
    Geom<Real> geom;
    Real flattenedCoords[2*Ncoords]={0};
    Real inCoords[2*Nin]={0};
    for (int i=0;i<Nin;++i){
        inCoords[colMajorIndex(0,i,2)] = inXcoords[i];
        inCoords[colMajorIndex(1,i,2)] = inYcoords[i];
    }

    make_panels(inCoords,flattenedCoords,Ufac,TEfac); // does spline to redist nodes over aerofoil for fixed number of 200 nodes


    Foil<Real> foil(flattenedCoords);
    Isolc<Real> isolc;
    Isolv<Real> isol_pre;
    Vsol<Real> vsol ;
    Param<Real> param;
    param.ncrit = nCrit;
    Wake<Real> wake;

    Glob<Real> glob;
    for (int i=0;i<RVdimension;++i){glob.U[i] = states[i];}
    for (int i=0;i<(Ncoords+Nwake);++i){vsol.turb[i] = turb[i];}

    build_gamma_codi<Real>(isolc,foil,oper);
    init_thermo(oper,param,geom);
    build_wake(foil,geom,oper,isolc,wake);
    stagpoint_find(isolc,isol_pre,foil,wake);
    identify_surfaces(isol_pre,vsol);
    set_wake_gap(foil,isol_pre,vsol);
    calc_ue_m<Real>(foil,wake,isolc,vsol);

    Isolv<Real> isol_final;
    stagpoint_move_AD(isol_final,glob,foil,wake,vsol,currStag);

    // Compute vsol.forcet[si] and vsol.xift[si] from xft_xc and foil geometry.
    // Must run after identify_surfaces (populates vsol.Is) AND after
    // stagpoint_move_AD, so that xift is measured against isol_FINAL.distFromStag
    // — the SAME arc-length array build_glob_RV_AD uses for the residual station
    // coordinates. This matches the forward solver order (coupled.cpp:
    // stagpoint_move -> update_transition -> build_glob_RV) and is essential for
    // correctness: the residual depends on the transition position RELATIVE to
    // its stations (xt - x1), so a stag-point shift (e.g. from changing alpha)
    // must cancel between xift and the stations. Using isol_pre here would leave
    // a spurious dxift/dalpha that breaks that cancellation.
    //
    // The PANEL BRACKET containing the x/c forcing station is an integer/branch
    // selection — kept passive via .getValue() comparisons. The interpolation
    // that locates xift WITHIN that bracket is continuous in the (taped) node
    // coordinates foil.x and arc-lengths distFromStag, so it is computed in Real
    // arithmetic. This restores d(xift)/d(geometry) on the tape, which the
    // adjoint needs as dR/dxift * dxift/dy (see bench/results/GRAD_VERIFY.md).
    for (int si = 0; si < 2; ++si) {
        vsol.forcet[si] = false;
        vsol.xift[si]   = Real(0.0);
        if (xft_xc[si] < 1.0 - 1e-9 && si < (int)vsol.Is.size()) {
            const std::vector<int>& Is_si = vsol.Is[si];
            int nSP = (int)Is_si.size();
            // x_max is the (taped) max node abscissa — the chord scale moves with
            // the geometry, so keep the value taped while selecting the index passively.
            Real x_max = foil.x[0];
            for (int k = 1; k < Ncoords; ++k)
                if (foil.x[2*k].getValue() > x_max.getValue()) x_max = foil.x[2*k];
            Real xft_abs = xft_xc[si] * x_max;
            for (int k = 1; k < nSP; ++k) {
                Real x_prev = foil.x[2 * Is_si[k-1]];
                Real x_curr = foil.x[2 * Is_si[k  ]];
                if ((x_prev.getValue() - xft_abs.getValue()) *
                    (x_curr.getValue() - xft_abs.getValue()) <= 0.0) {
                    Real xi_prev = isol_final.distFromStag[Is_si[k-1]];
                    Real xi_curr = isol_final.distFromStag[Is_si[k  ]];
                    Real frac = (x_curr.getValue() == x_prev.getValue()) ? Real(0.0) :
                                  (xft_abs - x_prev) / (x_curr - x_prev);
                    vsol.xift[si]   = xi_prev + (xi_curr - xi_prev) * frac;
                    vsol.forcet[si] = true;
                    break;
                }
            }
        }
    }

    build_glob_RV_AD(foil,vsol,isol_final,glob,param);
    finishdRdU_AD(foil,isolc,isol_final,glob,vsol,oper);

    Real hCL = 0.0;
    Real hCD = 0.0;
    Real hOASPL = 0.0;

    for (int i = 0; i < RVdimension; ++i){
        hCL -= adlambdaCL[i] * glob.R[i];
        hCD -= adlambdaCD[i] * glob.R[i];
        hOASPL -= adlambdaOASPL[i] * glob.R[i];
    }

    Real out[3] = {hCL,hCD,hOASPL} ;
    constexpr int jacobianHeight = 3;
    for (int i=0;i<jacobianHeight;++i){
        tape.registerOutput(out[i]);
    }
    tape.setPassive();
    for (int i=0;i<jacobianHeight;++i){out[i].gradient()[i] = 1.0 ;}
    tape.evaluate();

    for (int i = 0; i < Nin; ++i){
        dgCLdy[i] = inYcoords[i].getGradient()[0];
        dgCDdy[i] = inYcoords[i].getGradient()[1];
        dgOASPLdy[i] = inYcoords[i].getGradient()[2];
    }
    dgCLdalpha = alphad.getGradient()[0];
    dgCDdalpha = alphad.getGradient()[1];
    dgOASPLdalpha = alphad.getGradient()[2];

    tape.resetHard();
}


