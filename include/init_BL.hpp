#pragma once
#include <codi.hpp>
#include <Eigen/Dense>
#include "real_type.hpp"
#include "data_structs.hpp"
#include "get_funcs.hpp"
#include "newton_residuals.hpp"

// Helper types at file scope (required — cannot define templates inside functions in C++)
template<typename T>
using BLMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using BLVec = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename Number>
struct BLLinSolver : public codi::EigenLinearSystem<Number, BLMat, BLVec> {
    using Base = codi::EigenLinearSystem<Number, BLMat, BLVec>;
    using MatrixReal = typename Base::MatrixReal;
    using VectorReal = typename Base::VectorReal;
    void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
        *x = A->colPivHouseholderQr().solve(*b);
    }
};

template<typename Real>
void solve_linear_system_local(const Real* A, const Real* RHS, Real* xOut, const int matDim, const int outDim)
{
    BLMat<Real> AEig(matDim, matDim);
    BLVec<Real> rhsEig(matDim);
    BLVec<Real> sol(matDim);
    for (int i=0;i<matDim;++i) {
        for (int j=0;j<matDim;++j) AEig(i,j) = A[colMajorIndex(i,j,matDim)];
        rhsEig(i) = RHS[i];
    }
    codi::solveLinearSystem(BLLinSolver<Real>(), AEig, rhsEig, sol);
    for (int i=0;i<matDim && i<outDim;++i) xOut[i] = sol(i);
}

template<typename Real>
void thwaites_init(const Real& stagConstant, const Param<Real>& param, Real& momThickness, Real& dispThickness)
{
    Real nu = param.mu0 / param.rho0;
    momThickness  = std::sqrt(0.45*nu / (6.0*stagConstant));
    dispThickness = 2.2*momThickness;
}

template<typename Real>
void wake_init(const Vsol<Real>& vsol, const Foil<Real>& foil, const Glob<Real>& glob,
    const Param<Real>& param, Real ue, Real (&Uw)[4])
{
    int iw = vsol.Is[2].front();
    for (int i=0;i<4;++i) Uw[i] = glob.U[colMajorIndex(i,iw,4)];

    Real R[3]={0}, R_U[36]={0};
    int J[3]={0};
    wake_sys_full<Real>(vsol,foil,glob,param,R,R_U,J);
    for (int i=0;i<3;++i) Uw[i] -= R[i];
    Uw[3] = ue;
}

template<typename Real>
void init_boundary_layer(const Oper<Real>& oper, const Foil<Real>& foil, Param<Real>& param,
    const Isolc<Real>& isolc, Isolv<Real>& isolv, Vsol<Real>& vsol, Glob<Real>& glob,
    Trans<Real>& tdata, const bool force)
{
    const int Nsys = oper.Vinf > 0 ? glob.nc + glob.nw : Ncoords + Nwake; // runtime
    (void)Nsys; // Nsys not directly used here, Ncoords/Nwake from structs
    const Real Hmaxl = 3.8;
    const Real Hmaxt = 2.5;

    Real ue[Nsys]={0};
    get_ueinv<Real>(isolc, isolv, ue);

    for (int surf=0; surf<3; ++surf) {
        const std::vector<int>& indexList = vsol.Is[surf];
        int N = indexList.size();

        Real uemax = 0.0;
        for (int i=0;i<N;++i) uemax = std::max(uemax, std::abs(ue[indexList[i]]));
        for (int i=0;i<N;++i) ue[indexList[i]] = std::max(ue[indexList[i]], 1e-8*uemax);

        bool turb=false, wake=false;
        int i0=0;

        if (surf < 2) {
            bool hitstag = (isolv.distFromStag[indexList[0]] < 1e-8 * isolv.distFromStag[indexList[N-1]]);
            Real K = hitstag ? (ue[indexList[1]] / isolv.distFromStag[indexList[1]])
                             : (ue[indexList[0]] / isolv.distFromStag[indexList[0]]);

            Real th, ds;
            thwaites_init<Real>(K, param, th, ds);

            Real xst = 1e-6;
            Real Ust[4] = {th, ds, 0, K*xst};

            constexpr int nNewton=20;
            for (int iNewton=0; iNewton<nNewton; ++iNewton) {
                Real R[3]={0}, R_U[24]={0}, R_x[6]={0};
                residual_station_full<Real>(Ust,Ust,xst,xst,0,0,false,false,true,param,R,R_U,R_x);

                Real norm = std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
                if (norm < 1e-10) break;

                Real A[9]={0}, b[3]={0};
                for (int r=0;r<3;++r) {
                    b[r]=-R[r];
                    for (int c=0;c<3;++c) A[colMajorIndex(r,c,3)] = R_U[colMajorIndex(r,c+4,3)] + R_U[colMajorIndex(r,c,3)];
                }
                Real dU[4]={0};
                solve_linear_system_local<Real>(A,b,dU,3,4);
                Real dm = std::max(std::abs(dU[0]/Ust[0]), std::abs(dU[1]/Ust[1]));
                Real omega = (dm >= 0.2) ? 0.2/dm : Real(1.0);
                for (int r=0;r<2;++r) Ust[r] += omega*dU[r];
            }

            if (hitstag) {
                for (int j=0;j<4;++j) glob.U[colMajorIndex(j,indexList[0],4)] = Ust[j];
                glob.U[colMajorIndex(3,indexList[0],4)] = ue[indexList[0]];
                i0=1;
            }
            for (int j=0;j<3;++j) glob.U[colMajorIndex(j,indexList[i0],4)] = Ust[j];
            glob.U[colMajorIndex(3,indexList[i0],4)] = ue[indexList[i0]];
        } else {
            Real Uw[4]={0};
            wake_init<Real>(vsol,foil,glob,param,ue[indexList[0]],Uw);
            for (int j=0;j<3;++j) glob.U[colMajorIndex(j,indexList[0],4)] = Uw[j];
            glob.U[colMajorIndex(3,indexList[0],4)] = ue[indexList[0]];
            turb=true; wake=true; vsol.turb[indexList[0]]=1;
        }

        bool tran=false;
        int i = i0+1;
        Real prevState[4];
        for (int r=0;r<4;++r) prevState[r] = glob.U[colMajorIndex(r,indexList[i0],4)];

        while (i < N) {
            int prevNode = indexList[i-1];
            int currNode = indexList[i];
            Real currState[4] = {prevState[0], prevState[1], prevState[2], ue[currNode]};

            if (!turb && !tran && force && prevNode == tdata.transNode[surf]) {
                tran = true; tdata.isForced[surf] = 1;
            }
            if (tran) {
                Real ct_U[4]={0};
                currState[2] = get_cttr(currState[0],currState[1],currState[2],currState[3],turb,param,ct_U);
            }
            vsol.turb[currNode] = (tran || turb) ? 1 : 0;

            bool direct=true;
            constexpr int nNewton=30;
            const int iNswitch=12;
            Real Hktgt=0;

            for (int iNewton=0; iNewton<nNewton; ++iNewton) {
                Real R[3]={0}, R_U[24]={0}, R_x[6]={0};

                if (tran) {
                    Real x1=foil.x[colMajorIndex(0,prevNode,2)], y1=foil.x[colMajorIndex(1,prevNode,2)];
                    Real x2=foil.x[colMajorIndex(0,currNode,2)], y2=foil.x[colMajorIndex(1,currNode,2)];
                    Real yt = y1 + ((tdata.transPos[surf]-x1)/(x2-x1))*(y2-y1);
                    Real distTrans = isolv.distFromStag[prevNode] +
                        std::sqrt((tdata.transPos[surf]-x1)*(tdata.transPos[surf]-x1)+(yt-y1)*(yt-y1));
                    if (!tdata.isForced[surf])
                        residual_transition_full<Real>(prevState,currState,
                            isolv.distFromStag[prevNode],isolv.distFromStag[currNode],0,0,param,R,R_U,R_x);
                    else
                        residual_transition_forced_full<Real>(prevState,currState,
                            isolv.distFromStag[prevNode],isolv.distFromStag[currNode],param,distTrans,R,R_U,R_x);
                } else {
                    Real aux1=0, aux2=0;
                    if (wake) { aux1=vsol.wgap[prevNode-glob.nc]; aux2=vsol.wgap[currNode-glob.nc]; }
                    residual_station_full<Real>(prevState,currState,
                        isolv.distFromStag[prevNode],isolv.distFromStag[currNode],
                        aux1,aux2,wake,turb,false,param,R,R_U,R_x);
                }

                Real norm = std::sqrt(R[0]*R[0]+R[1]*R[1]+R[2]*R[2]);
                if (norm < 1e-10) break;

                Real dU[4]={0};
                if (direct) {
                    Real A[9]={0}, b[3]={0};
                    for (int r=0;r<3;++r) {
                        b[r]=-R[r];
                        for (int c=0;c<3;++c) A[colMajorIndex(r,c,3)] = R_U[colMajorIndex(r,c+4,3)];
                    }
                    solve_linear_system_local<Real>(A,b,dU,3,4);
                } else {
                    Real Hk_U[4]={0};
                    Real Hk = get_Hk(currState[0],currState[1],currState[3],param,Hk_U);
                    Real A[16]={0}, b[4]={0};
                    for (int r=0;r<3;++r) {
                        b[r]=-R[r];
                        for (int c=0;c<4;++c) A[colMajorIndex(r,c,4)] = R_U[colMajorIndex(r,c+4,3)];
                    }
                    for (int c=0;c<4;++c) A[colMajorIndex(3,c,4)] = Hk_U[c];
                    b[3] = Hktgt-Hk;
                    solve_linear_system_local<Real>(A,b,dU,4,4);
                }

                Real dm = std::max(std::abs(dU[0]/prevState[0]), std::abs(dU[1]/prevState[1]));
                if (!direct) dm = std::max(dm, std::abs(dU[3]/prevState[3]));
                if (turb)    dm = std::max(dm, std::abs(dU[2]/prevState[2]));
                else if (direct) dm = std::max(dm, std::abs(dU[2]/10));

                Real omega = (dm > 0.3) ? 0.3/dm : Real(1.0);

                if (tran) {
                    Real Ui2_test = currState[2] + omega*dU[2];
                    if (Ui2_test < 0.00001)
                        omega = std::max((0.00001-currState[2])/dU[2], Real(1e-6));
                }

                for (int r=0;r<4;++r) dU[r] *= omega;

                Real Ui[4] = {currState[0]+dU[0], currState[1]+dU[1], currState[2]+dU[2], currState[3]+dU[3]};
                if (turb) Ui[2] = std::max(std::min(Ui[2], Real(0.3)), Real(1e-7));

                Real Hmax = turb ? Hmaxt : Hmaxl;
                Real Hk2_U[4]={0};
                Real Hk2 = get_Hk(Ui[0],Ui[1],Ui[3],param,Hk2_U);

                if (direct && (Hk2 > Hmax || iNewton > iNswitch)) {
                    direct = false;
                    Real Hk_p[4]={0};
                    Real Hkp = get_Hk(prevState[0],prevState[1],prevState[3],param,Hk_p);
                    Real Hkr = (isolv.distFromStag[currNode]-isolv.distFromStag[prevNode])/prevState[0];
                    if (wake) {
                        Real H2=Hkp;
                        for (int k=0;k<6;++k) H2 -= (H2+0.03*Hkr*std::pow(H2-1,3)-Hkp)/(1+0.09*Hkr*(H2-1)*(H2-1));
                        Hktgt = std::max(H2, Real(1.01));
                    } else if (turb) {
                        Hktgt = Hkp - 0.15*Hkr;
                    } else {
                        Hktgt = Hkp + 0.03*Hkr;
                    }
                    if (!wake) Hktgt = std::max(Hktgt, Hmax);
                    if (iNewton > iNswitch) {
                        for (int r=0;r<3;++r) currState[r]=prevState[r];
                        currState[3]=ue[currNode];
                    }
                } else {
                    for (int r=0;r<4;++r) currState[r]=Ui[r];
                }
            }

            if (!turb && !tran && currState[2] > param.ncrit) {
                tran=true; tdata.isForced[surf]=0; continue;
            }
            if (tran) { turb=true; tran=false; }

            for (int r=0;r<4;++r) {
                glob.U[colMajorIndex(r,currNode,4)] = currState[r];
                prevState[r] = currState[r];
            }
            ++i;
        }
    }
}
