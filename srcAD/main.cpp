#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.hpp"
#include "ADfuncs.hpp"
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

using MatrixSparse = Eigen::SparseMatrix<double>;
using VectorD      = Eigen::Matrix<double, Eigen::Dynamic, 1>;

void solve_sys(int* R_V_cols, int* R_V_rows, double* RV_vals, const int RVlatest,
    double (&dCLdStates)[RVdimension],
    double (&dCDdStates)[RVdimension],
    double (&dOASPLdStates)[RVdimension],
    double (&adlambdaCL)[RVdimension],
    double (&adlambdaCD)[RVdimension],
    double (&adlambdaOASPL)[RVdimension])
{
    
    // === 1. Build triplets ===
    std::vector<Eigen::Triplet<double>> triplets;

    triplets.reserve(RVlatest);

    for (int k = 0; k < RVlatest; ++k) {
        int row = R_V_rows[k];
        int col = R_V_cols[k];
        double val = RV_vals[k];

        //if (val != double(0))
        triplets.emplace_back(row, col, val);
    }

    // === 2. Create sparse matrix ===
    MatrixSparse A_sparse(RVdimension, RVdimension);
    A_sparse.setFromTriplets(triplets.begin(), triplets.end());

    // === 3. Map RHS vectors ===
    VectorD rhs1(RVdimension);
    VectorD rhs2(RVdimension);
    VectorD rhs3(RVdimension);
    for (int i = 0; i < RVdimension; ++i){
        rhs1(i) = dCLdStates[i];
        rhs2(i) = dCDdStates[i];
        rhs3(i) = dOASPLdStates[i];
    }

    //Eigen::Map<const Vec> rhs1(dCLdStates);
    //Eigen::Map<const Vec> rhs2(dCDdStates);
    //Eigen::Map<const Vec> rhs3(dOASPLdStates);

    VectorD sol1(RVdimension);
    VectorD sol2(RVdimension);
    VectorD sol3(RVdimension);
    // === 4. Factorization ===
    Eigen::SparseLU<MatrixSparse> sparse_solver;
    sparse_solver.compute(A_sparse);

    if (sparse_solver.info() != Eigen::Success) {
        std::cerr << "Sparse LU factorization failed!\n";
        return;
    }

    // === 5. Solve for each RHS ===
    sol1 = sparse_solver.solve(rhs1);
    sol2 = sparse_solver.solve(rhs2);
    sol3 = sparse_solver.solve(rhs3);

    // === 6. Copy results back to output arrays ===

    for (int i=0;i<RVdimension;++i){
        adlambdaCL[i] = sol1(i);
        adlambdaCD[i] = sol2(i);
        adlambdaOASPL[i] = sol3(i);

    }
}

void solve_sys_single(int* R_V_cols, int* R_V_rows, double* RV_vals, const int RVlatest,
    double (&dCLdStates)[RVdimension],
    double (&adlambdaCL)[RVdimension]) 
    {
    
    constexpr int Nsize = 4 * (Ncoords + Nwake);

    // === 1. Build triplets from glob arrays ===
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(RVlatest);
    for (int k = 0; k < RVlatest; ++k) {
        int row = R_V_rows[k];
        int col = R_V_cols[k];
        double val = RV_vals[k];
        if (val != 0.0) {
            triplets.emplace_back(row, col, val);
        }
    }

    // === 2. Fill sparse matrix from triplets ===
    Eigen::SparseMatrix<double> A_sparse(Nsize, Nsize);
    A_sparse.setFromTriplets(triplets.begin(), triplets.end());
    
    Eigen::Map<const Eigen::Matrix<double, RVdimension, 1, Eigen::ColMajor>>
    rhs_eigen(dCLdStates, Nsize, 1);

    // Use SparseLU solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>> sparse_solver;
    sparse_solver.compute(A_sparse);

    if(sparse_solver.info() != Eigen::Success) {
        std::cerr << "Sparse solver failed during factorization!\n";
        return; // or handle the error appropriately
    }
    
    Eigen::Matrix<double, RVdimension, 1> x = sparse_solver.solve(rhs_eigen);

    // Map the solution to output vector
    Eigen::Map<Eigen::Matrix<double, RVdimension, 1, Eigen::ColMajor>>
        x_eigen(adlambdaCL, Nsize, 1);
    x_eigen = x;
}

int main(){


    // Open JSON file
    std::ifstream infile("input.json");
    if (!infile) {
        std::cerr << "Failed to open input.json\n";
        return 1;
    }

    // Parse the JSON
    json j;
    infile >> j;

    //using Realfwd = codi::RealForward ;
    using RealVec2 = codi::RealReverseVec<2> ;
    using RealRev = codi::RealReverseVec<3> ;
    using Realfwd = codi::RealForward;
   
    RealVec2 inYcoords_2[Nin]={0};
    RealRev inYcoords_Rev[Nin]={0};
    double inXcoords_d[Nin] = {0};
    for (int i = 0; i < Nin; ++i) {   
        inXcoords_d[i] = j["xcoords"][i];  
        inYcoords_2[i] = j["ycoords"][i];
        inYcoords_Rev[i] = j["ycoords"][i];
    }

    // Read input variables
    RealVec2 targetAlphaDeg = j["alpha_degrees"].get<double>();
    RealVec2 Re = j["Re"].get<double>();
    RealVec2 Ma = j["Ma"].get<double>();
    RealVec2 rhoInf = j["rho"].get<double>();
    RealVec2 nuInf = j["nu"].get<double>();
    RealVec2 custChord = j["chord"].get<double>();
    RealVec2 sampleTE = j["sampleTE"].get<double>();
    const RealVec2 X = j["X"].get<double>();
    const RealVec2 Y = j["Y"].get<double>();
    const RealVec2 Z = j["Z"].get<double>();
    const RealVec2 S = j["S"].get<double>();
    const RealVec2 Ncrit = j["ncrit"].get<double>();
    const RealVec2 Ufac = j["Ufac"].get<double>();
    const RealVec2 TEfac = j["TEfac"].get<double>();
    // forcing transition variables
    const bool force = j["forcetrans"].get<int>();
    const RealVec2 topTransPos = j["toptrans"].get<double>();
    const RealVec2 botTransPos = j["bottrans"].get<double>();
    const std::string model = j["model"].get<std::string>();

    // Read input variables
    
    RealRev targetAlphaDeg_r = j["alpha_degrees"].get<double>();
    RealRev Re_r = j["Re"].get<double>();
    RealRev Ma_r = j["Ma"].get<double>();
    RealRev rhoInf_r = j["rho"].get<double>();
    const RealRev Ncrit_r = j["ncrit"].get<double>();
    const RealRev Ufac_r = j["Ufac"].get<double>();
    const RealRev TEfac_r = j["TEfac"].get<double>();

    // forcing transition variables
    const RealRev topTransPos_r = j["toptrans"].get<double>();
    const RealRev botTransPos_r = j["bottrans"].get<double>();

    // Open JSON file
    std::ifstream restartfile("restart.json");
    if (!restartfile) {
        std::cerr << "Failed to open restart.json\n";
        return 1;
    }
    // Parse the JSON
    json jr;
    restartfile >> jr;
    static RealVec2 states[RVdimension] = {0};
    static double states_d[RVdimension] = {0};
    int turb[Ncoords+Nwake] = {0};
    for (int i = 0; i < RVdimension; ++i) {
        states[i] = jr["states"][i]; 
        states_d[i] = jr["states"][i];
    }
    for (int i = 0; i < (Ncoords+Nwake); ++i) {
        turb[i] = jr["turb"][i].get<int>(); 
    }  

    const double theta = jr["states"][RVdimension-4];
    const double disp = jr["states"][RVdimension-3];
    const double ue = jr["states"][RVdimension-1];
    
    int currStag[2] = {jr["stag"][0], jr["stag"][1]} ;

    static double dRdU_vals[119700] = {0};
 
    static int dRdU_rows[119700] = {0};
    static int dRdU_cols[119700] = {0};  
    int latetsEntry = jr["RVnz"] ;

    for (int i=0;i<latetsEntry;++i){

        dRdU_rows[i] = jr["RVcols"][i].get<int>() ;
        dRdU_cols[i] = jr["RVrows"][i].get<int>() ;
        dRdU_vals[i] = jr["RVvals"][i] ;
    }
    
    double d_CL_d_y[Nin];
    double d_OASPL_d_y[Nin];
    double d_CL_dalpha ;
    double d_OASPL_dalpha;
    static double d_CL_d_States[RVdimension] = {0};
    static double d_CD_d_States[RVdimension] = {0};
    static double d_OASPL_d_States[RVdimension] = {0};

    // Exponent: (5/2 + disp/(2*theta))
    Realfwd exponent = 2.5 + disp / (2.0 * theta);
    // Common term: ue^(exponent)
    Realfwd ue_pow = std::pow(ue, exponent);
    // dC/d cmom
    Realfwd dCdcmom = ue_pow * (2.0 * theta - disp * std::log(ue)) / theta;
    // Cdddisp = log(ue) * ue^(exponent)
    Realfwd Cdddisp = std::log(ue) * ue_pow;
    // Cddue = (2*theta) * exponent * ue^(exponent - 1)
    Realfwd Cddue = (2.0 * theta) * exponent * std::pow(ue, exponent - 1.0);

    double OASPL = partialOutputspartialInputs<RealVec2>(Ncrit,Ufac,TEfac,custChord,inXcoords_d,Re,Ma,rhoInf,nuInf,
        model,sampleTE,X,Y,Z,S,inYcoords_2,targetAlphaDeg,states,turb,
        d_CL_d_y,d_OASPL_d_y,d_CL_dalpha,d_OASPL_dalpha,d_CL_d_States,d_OASPL_d_States
    );

    static double adlambda_CL[RVdimension] = {0.0}; 
    static double adlambda_CD[RVdimension] = {0.0};
    static double adlambda_OASPL[RVdimension] = {0.0};

    d_CD_d_States[RVdimension-4] = dCdcmom.getValue();
    d_CD_d_States[RVdimension-3] = Cdddisp.getValue();
    d_CD_d_States[RVdimension-1] = Cddue.getValue();
   
    solve_sys(dRdU_cols,dRdU_rows,dRdU_vals,latetsEntry,
        d_CL_d_States,d_CD_d_States,d_OASPL_d_States,
        adlambda_CL,adlambda_CD,adlambda_OASPL
    );
    
    static double dgdy_CL[Nin] = {0};
    static double dgdy_CD[Nin] = {0};
    static double dgdy_OASPL[Nin] = {0};
    double dgdalpha_CL = {0};
    double dgdalpha_CD = {0};
    double dgdalpha_OASPL = {0};
    
    partialRpartialx<RealRev>(Ncrit_r,Ufac_r,TEfac_r,inXcoords_d,Re_r,Ma_r,rhoInf_r,topTransPos_r,botTransPos_r,currStag,
        adlambda_CL,adlambda_CD,adlambda_OASPL,
        inYcoords_Rev,targetAlphaDeg_r,states_d,turb,
        dgdy_CL,dgdy_CD,dgdy_OASPL,dgdalpha_CL,dgdalpha_CD,dgdalpha_OASPL
    );

    // total gradient
    double totalDerivative_CL[Nin] = {0.0};
    double totalDerivative_CD[Nin] = {0.0};
    double totalDerivative_OASPL[Nin] = {0.0};

    for(int i=0;i<Nin;i++){
        totalDerivative_CL[i] = d_CL_d_y[i] + dgdy_CL[i];
        totalDerivative_CD[i] = dgdy_CD[i];
        totalDerivative_OASPL[i] = d_OASPL_d_y[i] + dgdy_OASPL[i];
    }

    double totalDerivativeCLAlpha = d_CL_dalpha + dgdalpha_CL ;
    double totalDerivativeCDAlpha = dgdalpha_CD ;
    double totalDerivativeOASPLAlpha = d_OASPL_dalpha + dgdalpha_OASPL ;
    
    json impGrad;
    impGrad["d cl / d ycoords"] = totalDerivative_CL ;
    impGrad["d cl / d alpha"] = totalDerivativeCLAlpha;
    impGrad["d cd / d ycoords"] = totalDerivative_CD ;
    impGrad["d cd / d alpha"] =  totalDerivativeCDAlpha ;
    impGrad["d OASPL / d ycoords"] = totalDerivative_OASPL ;
    impGrad["d OASPL / d alpha"] = totalDerivativeOASPLAlpha ;
    std::ofstream fout("ad_gradients.json");
    fout << impGrad.dump(4); 

    return 1;
};