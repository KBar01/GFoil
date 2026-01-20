#include <iostream>
#include <cmath>
#include <codi.hpp>
#include <Eigen/Dense>
#include "real_type.h"
#include "panel_funcs.h"
#include "data_structs.h"
#include "get_funcs.h"
#include "residuals.h"

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;


void thwaites_init(const Real&stagConstant, const Param&param,Real& momThickness,Real&dispThickness){

    Real nu = param.mu0/param.rho0 ;
    momThickness = std::sqrt(0.45*nu/(6.0*stagConstant));
    dispThickness = 2.2*momThickness;
}

/*
Thought on optimising: calcuolating residual currently involves the calculation of parameters
on both nodes, however you only iteratively change the end node in the newton steps, the start
node remains unchanged, so its parameters are fixed for that iteration, dont need to recompute.
*/ 


#ifndef USE_CODIPACK

void solve_linear_system(const Real* A, const Real* RHS, Real*xOut, const int matDim,const int outDim) {


    // Wrap the flattened column-major AIC vector into an Eigen matrix
    Eigen::Map<const Eigen::Matrix <Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> 
        A_eigen(A, matDim, matDim);

    // Wrap the flattened column-major RHS vector into an Eigen matrix
    // Eigen-map is class that maps raw data (like c++ array) to eigen matrix/vector
    // WITHOUT copying data.
    // Inputs to this class are datatype (Real in this case), dynamic means non-fixed at compile time
    // TODO: eigen sizes should be static.
    // ColMajor specifies how data is stored, matching arrays

    Eigen::Map<const Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> 
        rhs_eigen(RHS, matDim, 1);
    
    // Solve Ax=b using Eigen's solver

    // Note: x here is a Eigen matrix 
    Eigen::MatrixXd x = A_eigen.colPivHouseholderQr().solve(rhs_eigen);

    // Ensure gamref is correctly mapped

    // This is taking the array of isol struct and temporarily mapping an
    // Eigen matrix to it, which allows easy update using solution as done in 
    // next step
    Eigen::Map<Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> 
        xOut_eigen(xOut, outDim, 1);


    if (matDim == outDim){
        // Store solution in xOUT, ignoring the
        xOut_eigen  = x;
    }
    else{
        xOut_eigen.block(0,0,matDim,1) = x.block(0,0,matDim,1);
    }


}
#else



template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
 
template<typename Type>
void func(Matrix<Type> const& A, Vector<Type> const& rhs, Vector<Type>& sol) {
  sol = A.colPivHouseholderQr().solve(rhs);
}
 
template<typename Number>
struct EigenSolver : public codi::EigenLinearSystem<Number, Matrix, Vector> {
  public:
 
    using Base = codi::EigenLinearSystem<Number, Matrix, Vector>;  
    using MatrixReal = typename Base::MatrixReal;                  
    using VectorReal = typename Base::VectorReal;                  
 
    void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
        func(*A, *b, *x);
    }
};


void solve_linear_system(const Real* A, const Real* RHS, Real*xOut, const int matDim,const int outDim) {


    // Wrap the flattened column-major AIC vector into an Eigen matrix
   
    // Wrap the flattened column-major RHS vector into an Eigen matrix
    // Eigen-map is class that maps raw data (like c++ array) to eigen matrix/vector
    // WITHOUT copying data.
    // Inputs to this class are datatype (Real in this case), dynamic means non-fixed at compile time
    // TODO: eigen sizes should be static.
    // ColMajor specifies how data is stored, matching arrays
    Matrix<Real> AEigen(matDim, matDim);
    Vector<Real> rhsEigen(matDim);
    Vector<Real> sol(matDim);
    
    for (int i = 0; i < matDim; ++i) {
        for (int j = 0; j < matDim; ++j) {
            AEigen(i, j) = A[colMajorIndex(i, j, matDim)];
        }
        rhsEigen(i) = RHS[i];
    }

    // Solve Ax=b using Eigen's solver

    
    // Note: x here is a Eigen matrix 
    codi::solveLinearSystem(EigenSolver<Real>(), AEigen, rhsEigen, sol);

    // Ensure gamref is correctly mapped

    // This is taking the array of isol struct and temporarily mapping an
    // Eigen matrix to it, which allows easy update using solution as done in 
    // next step


    if (matDim == outDim){
        // Store solution in xOUT, ignoring the
        for (int i=0; i<matDim; ++i){
            xOut[i] = sol(i);
        }
    }
    else{
        for (int i=0; i<matDim; ++i){
            xOut[i] = sol(i);
        }
    }

    


}
#endif


void wake_init(const Vsol& vsol, const Foil& foil, const Glob& glob, const Param& param, Real ue, Real (&Uw)[4]) {

    // Get first wake index
    int iw = vsol.Is[2].front();

    // Extract Uw from global U
    for (int i = 0; i < 4; ++i) {Uw[i] = glob.U[colMajorIndex(i,iw,4)];}

    // Compute residual system
    Real R[3]={0}, R_U[36]={0};
    int J[3]={0};
    wake_sys(vsol, foil, glob, param, R, R_U, J);


    // Solve the residual system by subtracting R
    for (int i = 0; i < 3; ++i) {Uw[i] -= R[i];}
    
    Uw[3] = ue;
}






void init_boundary_layer(const Oper&oper, const Foil&foil, Param&param, Isol&isol, Vsol&vsol, Glob&glob, Trans&tdata, const bool force) {
    
    constexpr int Nsys = Ncoords + Nwake;
    const Real Hmaxl = 3.8;
    const Real Hmaxt = 2.5;

    Real ue[Nsys]={0};
    get_ueinv(isol,ue);  // get inviscid edge velocity

    for (int surf = 0; surf < 3; ++surf) {

        const std::vector<int> indexList = vsol.Is[surf] ;
        int N = indexList.size(); // How many nodes are on this surface

        // Ensure edge velocities are not tiny fofr each surface
        Real uemax = 0.0;
        for (int i = 0; i < N; ++i){
            uemax = std::max(uemax, std::abs(ue[indexList[i]]));
        }
        
        for (int i = 0; i < N; ++i){
            ue[indexList[i]] = std::max(ue[indexList[i]], 1e-8*uemax);
        }

        
        bool turb = false, wake = false,simi=false;
        
        int i0 = 0;
        if (surf < 2) {
            
            // Solving at stagnation node, initialising viscous properties using Thwaites
            bool hitstag = false;
            Real K = (isol.distFromStag[indexList[0]] < 1e-8*isol.distFromStag[indexList[N-1]]) ? (ue[indexList[1]] / isol.distFromStag[indexList[1]]) : (ue[indexList[0]] / isol.distFromStag[indexList[0]]);
            hitstag = (isol.distFromStag[indexList[0]] < 1e-8 * isol.distFromStag[indexList[N-1]]);

            Real th, ds;
            thwaites_init(K, param, th, ds);

            Real xst = 1e-6;
            Real Ust[4] = {th, ds, 0, K*xst};   // The stagnation state initial guess

            // Do convergence using BL residuals on the initial stagnation state guess
            constexpr int nNewton = 20;
            Real R[3]={0}, R_U[24]={0}, R_x[6]={0};  // Residuals init

            for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
                
                //param.turb = false;
                //param.simi = true;

                // Find residuals, given its the stagnation point (controlled by bool flags)
                residual_station(Ust,Ust,xst,xst,0,0,false,false,true,param,R,R_U,R_x);

                //param.simi = false;

                if (norm2_3D(R) < 1e-10) break;   // Convergence check

                Real A[9], b[3];
                for (int i = 0; i < 3; ++i) {
                    b[i] = -R[i];
                    for (int j = 0; j < 3; ++j){
                        A[colMajorIndex(i,j,3)] = R_U[colMajorIndex(i,j+4,3)] + R_U[colMajorIndex(i,j,3)];
                    }
                }

                Real dU[4]={0};
                solve_linear_system(A, b, dU, 3, 4);        // solves for updates to state (0 update on edge vel)
                
                // under-relaxation 
                Real dm = std::max(std::abs(dU[0] / Ust[0]), std::abs(dU[1] / Ust[1]));
                Real omega = 1.0;
                if (dm >= 0.2) omega = 0.2/dm;

                for (int i = 0; i < 2; ++i) {Ust[i] += omega*dU[i];}
            }

            if (hitstag) {
                for (int j=0; j<4; ++j){ glob.U[colMajorIndex(j,indexList[0],4)] = Ust[j];}
                glob.U[colMajorIndex(3,indexList[0],4)] = ue[indexList[0]]; // Reset edge vel
                i0 = 1;
            }

            for (int j=0; j<3; ++j){ glob.U[colMajorIndex(j,indexList[i0],4)] = Ust[j];}
            glob.U[colMajorIndex(3,indexList[i0],4)] = ue[indexList[i0]];
        } 
        
        else { // initialising wake viscous states
            Real Uw[4]={0};
            wake_init(vsol,foil,glob,param,ue[indexList[0]],Uw);

            //check_for_nans(ue,230,"ue vals");

            for (int j=0; j<3; ++j){ glob.U[colMajorIndex(j,indexList[0],4)] = Uw[j];}
            glob.U[colMajorIndex(3,indexList[0],4)] = ue[indexList[0]];
            turb = true;
            wake = true;
            vsol.turb[indexList[0]] = true;

            //check_for_nans(glob.U,920,"wake init section");
        }

        // Looping over rest of point on given surface, using prevbious node as intial guess
        bool tran = false;
        int i = i0 + 1;
        
        
        Real prevState[4];
        for (int r=0;r<4;++r){
            prevState[r] = glob.U[colMajorIndex(r,indexList[i0],4)]; // prev state carry-over array
        }
        
        while (i < N) {
            
            
            // index of curr and prev node on given surface
            int prevNode = indexList[i-1];
            int currNode = indexList[i];

            // create initial guess using previsous state + current node edge velocity
            Real currState[4] = {prevState[0],prevState[1],prevState[2],ue[currNode]} ;

            // If forcing transition, and on node before tran location, switch to solving transition over the panel
            if ((!turb && !tran && force) && prevNode == tdata.transNode[surf]){
                tran = true;
                tdata.isForced[surf] = 1;
            }

            if (tran) {
                Real ct, ct_U[4]={0};
                ct = get_cttr(currState[0],currState[1],currState[2],currState[3],turb,param,ct_U);
                currState[2] = ct;
            }

            vsol.turb[currNode] = (tran || turb);
            bool direct = true;
            constexpr int nNewton = 30;
            const int iNswitch = 12;
            Real Hktgt={0} ;
            // Main newton loop for convergence of BL properties
            for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
                Real R[3]={0}, R_U[24]={0}, R_x[6]={0};
                
                if (tran) {
                    
                    // do linear interp to find y value
                    Real x1 = foil.x[colMajorIndex(0,prevNode,2)],y1 = foil.x[colMajorIndex(1,prevNode,2)];
                    Real x2 = foil.x[colMajorIndex(0,currNode,2)],y2 = foil.x[colMajorIndex(1,currNode,2)];

                    Real yt = y1 + ((tdata.transPos[surf] - x1) / (x2 - x1)) * (y2 - y1);
                    
                    Real distFromStagTrans =  isol.distFromStag[prevNode] + std::sqrt((tdata.transPos[surf]-x1)*(tdata.transPos[surf]-x1)  + (yt-y1)*(yt-y1));
                    if (!tdata.isForced[surf]){
                        residual_transition(prevState,currState,isol.distFromStag[prevNode],isol.distFromStag[currNode],0.0,0.0,param,R,R_U,R_x);
                    }
                    else{
                        residual_transition_forced(prevState,currState,isol.distFromStag[prevNode],isol.distFromStag[currNode],param,distFromStagTrans,R,R_U,R_x);
                    }
                } 
                else {
    
                    Real aux1=0,aux2=0;
                    if (wake){aux1=vsol.wgap[prevNode-Ncoords];aux2 = vsol.wgap[currNode-Ncoords] ;}
                    residual_station(prevState,currState,isol.distFromStag[prevNode],isol.distFromStag[currNode],aux1,aux2,wake,turb,false,param,R,R_U,R_x);
                }

                if (norm2_3D(R) < 1e-10) break;

                Real dU[4]={0};
                if (direct) {
                    Real A[9]={0}, b[3]={0};
                    for (int r = 0; r < 3; ++r) {
                        b[r] = -R[r];
                        for (int c = 0; c < 3; ++c){
                            A[colMajorIndex(r,c,3)] = R_U[colMajorIndex(r,c+4,3)];}
                    }
                    solve_linear_system(A, b, dU, 3, 4);
                } 
                else {
                    // Handle inverse mode
                    Real Hk, Hk_U[4]={0};
                    Hk = get_Hk(currState[0],currState[1],currState[3],param,Hk_U);

                    Real A[16]={0}, b[4]={0};
                    for (int r = 0; r < 3; ++r) {
                        b[r] = -R[r];
                        for (int c = 0; c < 4; ++c)
                            A[colMajorIndex(r,c,4)] = R_U[colMajorIndex(r,c+4,3)];
                    }

                    // edge velocity gets swapped out in state for shape param, and this is used instead
                    for (int c = 0; c < 4; ++c){A[colMajorIndex(3,c,4)] = Hk_U[c];}
                    b[3] = Hktgt - Hk;
                    
                    solve_linear_system(A, b, dU, 4, 4);
                }
                
                // under relaxation :
                Real dm = std::max(std::abs(dU[0]/prevState[0]), std::abs(dU[1]/prevState[1]));
                if (!direct){ dm = std::max(dm,std::abs(dU[3]/prevState[3]));}
                if (turb)   { dm = std::max(dm,std::abs(dU[2]/prevState[2]));}
                else if (direct){ dm = std::max(dm,std::abs(dU[2]/10));}
                

                
                Real omega = 1;
                if (dm>0.3){omega = 0.3/dm;}

                if (tran){

                    // Tentative update
                    Real Ui_test2 = currState[2] + omega * dU[2];

                    // If U[2] would become too small or negative, scale omega
                    if (Ui_test2 < 0.00001) {
                        omega = (0.00001 - currState[2]) / dU[2];
                        omega = std::max(omega, 1e-6); // prevent divide-by-zero or too small
                    }
                }

                for (int r=0;r<4;++r){dU[r] *= omega;}


                // intermediate updated state
                Real Ui[4] = {
                    currState[0]+dU[0],
                    currState[1]+dU[1],
                    currState[2]+dU[2],
                    currState[3]+dU[3]
                };
                

                if (turb){ Ui[2] = std::max(std::min(Ui[2], 0.3), 1e-7);} //clip values

                // Check for near seperation
                Real Hmax = turb ? Hmaxt : Hmaxl;
                Real Hk_U[4]={0};
                Real Hk = get_Hk(Ui[0],Ui[1],Ui[3],param,Hk_U);

                
                if (direct && ((Hk>Hmax) || (iNewton>iNswitch))){
                    // Not taking update, switching to inverse mode insteasd and retry
                    
                    
                    direct = false;
                    for (int k=0; k<4;++k){Hk_U[k]=0;}
                    Hk = get_Hk(prevState[0],prevState[1],prevState[3],param,Hk_U);
                    Real Hkr = (isol.distFromStag[currNode]-isol.distFromStag[prevNode])/prevState[0] ;
                    if (wake){
                        
                        Real H2=Hk;
                        for (int k=0; k<6;++k){ H2 -= (H2+.03*Hkr*std::pow(H2-1,3)-Hk)/(1+.09*Hkr*(H2-1)*(H2-1));}
                        Hktgt = std::max(H2,1.01);
                    }
                    else if (turb) {Hktgt = Hk - 0.15*Hkr;}
                    else {Hktgt = Hk + 0.03*Hkr;}

                    if (! wake){ Hktgt = std::max(Hktgt, Hmax);}
                    
                    // Re-intialise the guess state 
                    if (iNewton > iNswitch){
                        
                        for (int r=0;r<3;++r){currState[r] = prevState[r];}
                        currState[3] = ue[currNode];
                    }
                }
                else {
                    // Taking the update
                    for (int r=0;r<4;++r){currState[r] = Ui[r];}
                }
            }

            // Have converged BL state for given node, now need to check for transition based on amplification value
            
            // HERE!!! add an additional statement inside the if, saying if not tran & ncrit exceeded OR forced trans node
            // set tran = true 
            
            if (!turb && ((!tran && currState[2]>param.ncrit))){
                tran = true;
                tdata.isForced[surf] = 0; // in this case transition is natural so dont need to force 
                
                continue ; // amplification exceeds ncrit, redo node position with trans= true
            }

            //if ((!turb && !tran && force) && (prevNode==tdata.transNode[surf])){
            //    tran = true;
            //    tdata.isForced[surf] = 1;  // using the forced approach, changes the residual calc buy fixing x trans
            //    continue ;
            //}

            if (tran){turb = true; tran=false;} // after tranistion, all nodes are turbulent
            
            for (int r=0;r<4;++r){   // store in global struct
                glob.U[colMajorIndex(r,currNode,4)] = currState[r];
                prevState[r] = currState[r];
            }
            i++; // move on to next node    
        }
    }
}



void init_boundary_layer_from_xfoil(const Oper&oper, const Foil&foil, const Param&param, Isol&isol, Vsol&vsol, Glob&glob) {
    
    constexpr int Nsys = Ncoords + Nwake;
    const Real Hmaxl = 3.8;
    const Real Hmaxt = 2.5;

    Real ue[Nsys]={0};
    get_ueinv(isol,ue);  // get inviscid edge velocity

    // Read in the input json with xfoil data in it
    // Open JSON file
    
    #ifndef USE_CODIPACK
    std::ifstream infile("xfoilstates.json");
    if (!infile) {
        std::cerr << "Failed to open xfoilstates.json\n";
    }

    // Parse the JSON
    json j;
    infile >> j;

    for (int i = 0; i < (Ncoords*4); ++i) {
        glob.U[i] = j["states"][i];       // Xfoil states on foil (col major )
    }
    for (int i = 0; i < (Ncoords); ++i){
        vsol.turb[i] = j["turb"][i].get<int>();
    }

    #endif
    int surf = 2;
    const std::vector<int> indexList = vsol.Is[surf] ;
    int N = indexList.size(); // How many nodes are on this surface

    // Ensure edge velocities are not tiny fofr each surface
    Real uemax = 0.0;
    for (int i = 0; i < N; ++i){
        uemax = std::max(uemax, std::abs(ue[indexList[i]]));
    }
    
    for (int i = 0; i < N; ++i){
        ue[indexList[i]] = std::max(ue[indexList[i]], 1e-8*uemax);
    }

    
    bool turb = false, wake = false,simi=false;
    
    int i0 = 0;
    
    
    // initialising wake viscous states
    Real Uw[4]={0};
    wake_init(vsol,foil,glob,param,ue[indexList[0]],Uw);

    //check_for_nans(ue,230,"ue vals");

    for (int j=0; j<3; ++j){ glob.U[colMajorIndex(j,indexList[0],4)] = Uw[j];}
    glob.U[colMajorIndex(3,indexList[0],4)] = ue[indexList[0]];
    turb = true;
    wake = true;
    vsol.turb[indexList[0]] = true;

    // Looping over rest of point on given surface, using prevbious node as intial guess
    bool tran = false;
    int i = i0 + 1;
    
    
    Real prevState[4];
    for (int r=0;r<4;++r){
        prevState[r] = glob.U[colMajorIndex(r,indexList[i0],4)]; // prev state carry-over array
    }
    
    while (i < N) {
        
        
        // index of curr and prev node on given surface
        int prevNode = indexList[i-1];
        int currNode = indexList[i];

        // create initial guess using previsous state + current node edge velocity
        Real currState[4] = {prevState[0],prevState[1],prevState[2],ue[currNode]} ;


        if (tran) {
            Real ct, ct_U[4]={0};
            ct = get_cttr(currState[0],currState[1],currState[2],currState[3],turb,param,ct_U);
            currState[2] = ct;
        }

        vsol.turb[currNode] = (tran || turb);
        bool direct = true;
        constexpr int nNewton = 30;
        const int iNswitch = 12;
        
        Real Hktgt={0} ;
        // Main newton loop for convergence of BL properties
        for (int iNewton = 0; iNewton < nNewton; ++iNewton) {
            Real R[3]={0}, R_U[24]={0}, R_x[6]={0};
            
            if (tran) {
                
                residual_transition(prevState,currState,isol.distFromStag[prevNode],isol.distFromStag[currNode],0.0,0.0,param,R,R_U,R_x);
            } 
            else {

                Real aux1=0,aux2=0;
                if (wake){aux1=vsol.wgap[prevNode-Ncoords];aux2 = vsol.wgap[currNode-Ncoords] ;}
                residual_station(prevState,currState,isol.distFromStag[prevNode],isol.distFromStag[currNode],aux1,aux2,wake,turb,false,param,R,R_U,R_x);
            }

            if (norm2_3D(R) < 1e-10) break;

            Real dU[4]={0};
            if (direct) {
                Real A[9]={0}, b[3]={0};
                for (int r = 0; r < 3; ++r) {
                    b[r] = -R[r];
                    for (int c = 0; c < 3; ++c){
                        A[colMajorIndex(r,c,3)] = R_U[colMajorIndex(r,c+4,3)];}
                }
                solve_linear_system(A, b, dU, 3, 4);
            } 
            else {
                // Handle inverse mode
                Real Hk, Hk_U[4]={0};
                Hk = get_Hk(currState[0],currState[1],currState[3],param,Hk_U);

                Real A[16]={0}, b[4]={0};
                for (int r = 0; r < 3; ++r) {
                    b[r] = -R[r];
                    for (int c = 0; c < 4; ++c)
                        A[colMajorIndex(r,c,4)] = R_U[colMajorIndex(r,c+4,3)];
                }

                // edge velocity gets swapped out in state for shape param, and this is used instead
                for (int c = 0; c < 4; ++c){A[colMajorIndex(3,c,4)] = Hk_U[c];}
                b[3] = Hktgt - Hk;
                
                solve_linear_system(A, b, dU, 4, 4);
            }
            
            // under relaxation :
            Real dm = std::max(std::abs(dU[0]/prevState[0]), std::abs(dU[1]/prevState[1]));
            if (!direct){ dm = std::max(dm,std::abs(dU[3]/prevState[3]));}
            if (turb)   { dm = std::max(dm,std::abs(dU[2]/prevState[2]));}
            else if (direct){ dm = std::max(dm,std::abs(dU[2]/10));}
            
            Real omega = 1;
            if (dm>0.3){omega = 0.3/dm;}
            for (int r=0;r<4;++r){dU[r] *= omega;}

            // intermediate updated state
            Real Ui[4] = {
                currState[0]+dU[0],
                currState[1]+dU[1],
                currState[2]+dU[2],
                currState[3]+dU[3]
            };
            

            if (turb){ Ui[2] = std::max(std::min(Ui[2], 0.3), 1e-7);} //clip values

            // Check for near seperation
            Real Hmax = turb ? Hmaxt : Hmaxl;
            Real Hk_U[4]={0};
            Real Hk = get_Hk(Ui[0],Ui[1],Ui[3],param,Hk_U);

            
            if (direct && ((Hk>Hmax) || (iNewton>iNswitch))){
                // Not taking update, switching to inverse mode insteasd and retry
                
                
                direct = false;
                for (int k=0; k<4;++k){Hk_U[k]=0;}
                Hk = get_Hk(prevState[0],prevState[1],prevState[3],param,Hk_U);
                Real Hkr = (isol.distFromStag[currNode]-isol.distFromStag[prevNode])/prevState[0] ;
                if (wake){
                    
                    Real H2=Hk;
                    for (int k=0; k<6;++k){ H2 -= (H2+.03*Hkr*std::pow(H2-1,3-Hk))/(1+.09*Hkr*(H2-1)*(H2-1));}
                    Hktgt = std::max(H2,1.01);
                }
                else if (turb) {Hktgt = Hk - 0.15*Hkr;}
                else {Hktgt = Hk + 0.03*Hkr;}

                if (! wake){ Hktgt = std::max(Hktgt, Hmax);}
                
                // Re-intialise the guess state 
                if (iNewton > iNswitch){
                    
                    for (int r=0;r<3;++r){currState[r] = prevState[r];}
                    currState[3] = ue[currNode];
                }
            }
            else {
                // Taking the update
                for (int r=0;r<4;++r){currState[r] = Ui[r];}
            }
        }

        // Have converged BL state for given node, now need to check for transition based on amplification value
        if (!turb && (!tran && currState[2]>param.ncrit)){
            tran = true;
            continue ; // amplification exceeds ncrit, redo node position with trans= true
        }

        if (tran){turb = true; tran=false;} // after tranistion, all nodes are turbulent
        
        for (int r=0;r<4;++r){   // store in global struct
            glob.U[colMajorIndex(r,currNode,4)] = currState[r];
            prevState[r] = currState[r];
        }
        i++; // move on to next node    
    }

}

