#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.hpp"
#include "data_structs.hpp"
#include "spline.hpp"
#include "solve_inv.hpp"
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


//# ifdef USE_CODIPACK
// helper for column-major indexing

/*
// 1. Type aliases for sparse Eigen
template<typename T>
using MatrixSparse = Eigen::SparseMatrix<T>;
template<typename T>
using Vector      = Eigen::Matrix<T, Eigen::Dynamic, 1>;

// 2. Your own solver for the numeric step
template<typename T>
void sparseSolveFunc(MatrixSparse<T> const& A, Vector<T> const& rhs, Vector<T>& sol) {
    // choose your factorization; SparseLU works for general unsymmetric
    Eigen::SparseLU<MatrixSparse<T>> lu;
    lu.compute(A);
    sol = lu.solve(rhs);
}

// 3. Wrap in CoDiPack's sparse linear system
template<typename Number>
struct SparseEigenSolver
  : public codi::SparseEigenLinearSystem<Number, MatrixSparse, Vector> {

    using Base       = codi::SparseEigenLinearSystem<Number, MatrixSparse, Vector>;
    using MatrixReal = typename Base::MatrixReal;  // numeric (Real) matrix
    using VectorReal = typename Base::VectorReal;  // numeric (Real) vector

    void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
        sparseSolveFunc(*A, *b, *x);  // just delegate to your numeric routine
    }
};

// 4. Your driver
template<typename Real>
void solve_sys_lambda(Real (&vals)[119700], int (&rows)[119700], int (&cols)[119700], int latestVal, Real (&dOutdStates)[RVdimension], Real (&adlambda)[RVdimension]) {
    constexpr int Nsize = RVdimension;

    // build sparse matrix from your glob arrays
    std::vector<Eigen::Triplet<Real>> triplets;
    triplets.reserve(latestVal);
    for (int k = 0; k < latestVal; ++k) {
        triplets.emplace_back(rows[k],
                              cols[k],
                              vals[k]);
    }

    MatrixSparse<Real> A(Nsize, Nsize);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Vector<Real> rhs(Nsize);
    for (int i = 0; i < Nsize; ++i) rhs(i) = dOutdStates[i];

    Vector<Real> sol(Nsize);

    // 5. Let CoDiPack handle AD by calling its wrapper:
    codi::solveLinearSystem(SparseEigenSolver<Real>(), A, rhs, sol);

    // 6. Write solution back into your state
    for (int i = 0; i < Nsize; ++i) {
        adlambda[i] = sol(i);
    }
}
*/

template<typename Real>
double partialOutputspartialInputs(

    // geometry parameters
    const Real nCrit, const Real Ufac, const Real TEfac, const Real chordScaling, 
    const double (&inXcoords)[Nin], const Real Re, const Real Ma, const Real rhoInf, 
    const Real kinViscInf,
    const std::string model, const Real sampleTE, const Real X, const Real Y, const Real Z,
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
    double (&jacobianOASPL_states)[RVdimension]
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

    interpolate_at_95_both_surfaces(xcoords,glob.U,post.cp,oper,turb,param,topsurf,botsurf,Uinf,sampleTE,chordScaling);
    //Real OASPL = calc_OASPL_AD(botsurf,topsurf,chordScaling,Uinf,X,Y,Z,S,kinViscInf,rhoInf,model);

    const bool doRotation = true ;
    Real OASPL = 0.0;
    if (doRotation){

        Real obsX = -Z*std::sin(alpha);
        Real obsZ = Z*std::cos(alpha);
        OASPL = calc_OASPL(botsurf,topsurf,chordScaling,Uinf,obsX,Y,obsZ,S,kinViscInf,rhoInf,doCps,model);
    }
    else{
        OASPL = calc_OASPL(botsurf,topsurf,chordScaling,Uinf,X,Y,Z,S,kinViscInf,rhoInf,doCps,model);
    }

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

    return post.cl.getValue();
    tape.reset();

};

template<typename Real>
void partialRpartialx(

    
    // geometry parameters
    const Real nCrit, const Real Ufac, const Real TEfac, 
    const double (&inXcoords)[Nin], const Real Re, const Real Ma, const Real rhoInf, 
    const Real &topTransPos, const Real &botTransPos,
    
    int (&currStag)[2],
    
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
   
    // finding node positions to force transition ------------------------
    Real xTransTop = topTransPos * geom.chord ;
    Real xTransBot = botTransPos * geom.chord ;
    
    int idx_closest_bot = 0;
    int idx_closest_top = Ncoords - 1;  // top TE node
    
    Trans<Real> tdata;

    tdata.transNode[0] = idx_closest_bot;
    tdata.transNode[1] = idx_closest_top;
    tdata.transPos[0] = xTransBot ;
    tdata.transPos[1] = xTransTop ;

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

    build_gamma_codi(isolc,foil,oper);
    init_thermo(oper,param,geom);
    build_wake(foil,geom,oper,isolc,wake);
    stagpoint_find(isolc,isol_pre,foil,wake);
    identify_surfaces(isol_pre,vsol);
    set_wake_gap(foil,isol_pre,vsol);
    calc_ue_m(foil,wake,isolc,vsol);
    Isolv<Real> isol_final;
    stagpoint_move_AD(isol_final,glob,foil,wake,vsol,currStag);
    build_glob_RV_AD(foil,vsol,isol_final,glob,param,tdata);
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

    tape.reset();
}


