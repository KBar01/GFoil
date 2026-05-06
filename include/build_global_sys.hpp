#pragma once
// Builds the global residual vector R and the sparse COO Jacobian dR/dU
// for the Newton iteration. Uses full-Jacobian residual functions.

#include "real_type.hpp"
#include "data_structs.hpp"
#include "vector_ops.hpp"
#include "newton_residuals.hpp"

template<typename Real>
void equate_block_inplace_sparse(
    Glob<Real>& glob,
    int startEntryRow, int startEntryCol,
    const Real* blockToInsert,
    int blockLd,
    int blockStartRow, int blockStartCol,
    int nRowsInsert, int nColsInsert)
{
    for (int j = 0; j < nColsInsert; ++j) {
        const Real* src = blockToInsert + (blockStartCol + j)*blockLd + blockStartRow;
        int globalCol = startEntryCol + j;
        for (int i = 0; i < nRowsInsert; ++i) {
            glob.R_V_rows.push_back(startEntryRow + i);
            glob.R_V_cols.push_back(globalCol);
            glob.R_V_vals.push_back(src[i]);
        }
    }
}

template<typename Real>
void findColumnIndices(const Glob<Real>& glob, int colIndex, std::vector<int>& indicesOut)
{
    indicesOut.clear();
    indicesOut.reserve(6);
    for (int k = 0; k < (int)glob.R_V_vals.size(); ++k) {
        if (glob.R_V_cols[k] == colIndex) indicesOut.push_back(k);
    }
}

template<typename Real>
void addColumnValues(
    Glob<Real>& glob, int colIndex,
    const std::vector<int>& existingIndices,
    const Real* valsToAdd, int nRowsAdd)
{
    for (int i = 0; i < nRowsAdd; ++i) {
        int row = i;
        Real val = valsToAdd[i];
        bool found = false;
        for (int idx : existingIndices) {
            if (glob.R_V_rows[idx] == row) { glob.R_V_vals[idx] += val; found = true; break; }
        }
        if (!found) {
            glob.R_V_rows.push_back(row);
            glob.R_V_cols.push_back(colIndex);
            glob.R_V_vals.push_back(val);
        }
    }
}

template<typename Real>
void stagnation_state_full(const Real* U1, const Real* U2, const Real x1, const Real x2,
    Real (&Ust)[4], Real (&Ust_U)[32], Real (&Ust_x)[8], Real& xst)
{
    Real dx = x2-x1;
    Real dx_x[2] = {-1, 1};
    Real rx = x2/x1;
    Real rx_x[2] = {-rx/x1, 1/x1};

    Real w1 = x2/dx,  w1_x[2] = {-w1/dx*dx_x[0], -w1/dx*dx_x[1] + 1/dx};
    Real w2 = -x1/dx, w2_x[2] = {-w2/dx*dx_x[0] - 1/dx, -w2/dx*dx_x[1]};

    for (int i=0;i<4;++i) Ust[i] = U1[i]*w1 + U2[i]*w2;

    Real wk1 = rx/dx,      wk1_x[2] = {rx_x[0]/dx - wk1/dx*dx_x[0], rx_x[1]/dx - wk1/dx*dx_x[1]};
    Real wk2 = -1/(rx*dx), wk2_x[2] = {-wk2*(rx_x[0]/rx + dx_x[0]/dx), -wk2*(rx_x[1]/rx + dx_x[1]/dx)};
    Real K   = wk1*U1[3] + wk2*U2[3];
    Real K_U[8] = {0,0,0,wk1, 0,0,0,wk2};
    Real K_x[2] = {U1[3]*wk1_x[0]+U2[3]*wk2_x[0], U1[3]*wk1_x[1]+U2[3]*wk2_x[1]};

    xst = 1e-6;
    Ust[3] = K*xst;

    for (int c=0;c<4;++c) for (int r=0;r<4;++r)
        Ust_U[r + c*4] = (r==c && r<3) ? w1 : 0.0;
    for (int c=0;c<4;++c) for (int r=0;r<4;++r)
        Ust_U[r + (c+4)*4] = (r==c && r<3) ? w2 : 0.0;
    for (int c=0;c<8;++c) Ust_U[3 + c*4] = K_U[c]*xst;

    Real temp1[6]={0}, temp2[6]={0};
    cnp::outer_product<3,2>(U1,w1_x,temp1);
    cnp::outer_product<3,2>(U2,w2_x,temp2);
    cnp::add_inplace<6>(temp1,temp2);
    Real botRow[2]; cnp::scalar_mul<2>(K_x,xst,botRow);
    cnp::vstack<3,1,2>(temp1,botRow,Ust_x);
}

template<typename Real>
void build_glob_RV(const Foil<Real>& foil, const Vsol<Real>& vsol,
    const Isolv<Real>& isol, Glob<Real>& glob, Param<Real>& param, Trans<Real>& tdata)
{
    const int nc = foil.nc, nw = vsol.nw;
    const int RXsize = 3*(nc+nw);

    std::vector<Real> R_st(RXsize, 0);
    const Real* xi = isol.distFromStag.data();

    for (int si = 0; si < 3; ++si) {
        const std::vector<int>& Is = vsol.Is[si];
        const int nSurfPoints = Is.size();
        int i0 = ((si < 2) && (xi[Is[0]] < 1e-8 * xi[Is[nSurfPoints-1]])) ? 1 : 0;

        bool turb = false, wake = false;

        if (si < 2) {
            Real R1[3], R1_U[24], R1_x[6];
            Real* U1 = &glob.U[colMajorIndex(0,Is[i0],4)];
            Real* U2 = &glob.U[colMajorIndex(0,Is[i0+1],4)];
            Real Ust[4], Ust_U[32], Ust_x[8], xst;
            stagnation_state_full<Real>(U1,U2,xi[Is[i0]],xi[Is[i0+1]],Ust,Ust_U,Ust_x,xst);

            Real R1_Ut[24]={0};
            residual_station_full<Real>(Ust,Ust,xst,xst,0,0,false,false,true,param,R1,R1_Ut,R1_x);
            Real R1_Ust[12]={0};
            cnp::add<12>(R1_Ut,R1_Ut+12,R1_Ust);
            cnp::matmat_mul<3,4,8>(R1_Ust,Ust_U,R1_U);
            cnp::matmat_mul<3,4,2>(R1_Ust,Ust_x,R1_x);

            int J[2] = {Is[i0], Is[i0+1]};
            int Ig = 3*Is[i0];
            for (int row=0;row<3;++row) glob.R[Ig+row] = R1[row];

            for (int j=0;j<2;++j)
                equate_block_inplace_sparse<Real>(glob,Ig,4*J[j],R1_U,3,0,4*j,3,4);

            for (int row=0;row<3;++row) {
                Real fv = R1_x[colMajorIndex(row,0,3)];
                if (isol.edgeVelSign[J[0]]==-1) R_st[Ig+row]+=fv; else R_st[Ig+row]-=fv;
                Real sv = R1_x[colMajorIndex(row,1,3)];
                if (isol.edgeVelSign[J[1]]==-1) R_st[Ig+row]+=sv; else R_st[Ig+row]-=sv;
            }
        } else {
            Real R1[3], R1_U[36]={0};
            int J[3];
            wake_sys_full<Real>(vsol,foil,glob,param,R1,R1_U,J);
            int Ig = 3*Is[i0];
            for (int row=0;row<3;++row) glob.R[Ig+row] = R1[row];
            for (int j=0;j<3;++j)
                equate_block_inplace_sparse<Real>(glob,Ig,4*J[j],R1_U,3,0,4*j,3,4);
            wake = true; turb = true;
        }

        for (int i = i0+1; i < nSurfPoints; ++i) {
            int prevI = i-1, currI = i;
            bool tran = vsol.turb[Is[prevI]] ^ vsol.turb[Is[currI]];

            Real Ri[3], Ri_U[24]={0}, Ri_x[6]={0};
            Real* Uprev = &glob.U[colMajorIndex(0,Is[prevI],4)];
            Real* Ucurr = &glob.U[colMajorIndex(0,Is[currI],4)];

            if (tran) {
                int isForced = tdata.isForced[si];
                if (!isForced) {
                    residual_transition_full<Real>(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],0,0,param,Ri,Ri_U,Ri_x);
                } else {
                    Real x1 = foil.x[colMajorIndex(0,Is[prevI],2)], y1 = foil.x[colMajorIndex(1,Is[prevI],2)];
                    Real x2 = foil.x[colMajorIndex(0,Is[currI],2)], y2 = foil.x[colMajorIndex(1,Is[currI],2)];
                    Real yt = y1 + ((tdata.transPos[si]-x1)/(x2-x1))*(y2-y1);
                    Real distTrans = isol.distFromStag[Is[prevI]] +
                        std::sqrt((tdata.transPos[si]-x1)*(tdata.transPos[si]-x1) + (yt-y1)*(yt-y1));
                    residual_transition_forced_full<Real>(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],param,distTrans,Ri,Ri_U,Ri_x);
                }
            } else {
                Real aux1=0, aux2=0;
                if (wake) { aux1=vsol.wgap[Is[prevI]-nc]; aux2=vsol.wgap[Is[currI]-nc]; }
                residual_station_full<Real>(Uprev,Ucurr,xi[Is[prevI]],xi[Is[currI]],aux1,aux2,wake,turb,false,param,Ri,Ri_U,Ri_x);
            }

            int Ig = 3*Is[i];
            for (int j=0;j<3;++j) glob.R[Ig+j] += Ri[j];

            equate_block_inplace_sparse<Real>(glob,Ig,4*Is[i-1],Ri_U,3,0,0,3,4);
            equate_block_inplace_sparse<Real>(glob,Ig,4*Is[i  ],Ri_U,3,0,4,3,4);

            for (int row=0;row<3;++row) {
                Real fv = Ri_x[colMajorIndex(row,0,3)];
                if (wake || isol.edgeVelSign[Is[i-1]]==1) R_st[Ig+row]-=fv; else R_st[Ig+row]+=fv;
                Real sv = Ri_x[colMajorIndex(row,1,3)];
                if (wake || isol.edgeVelSign[Is[i  ]]==1) R_st[Ig+row]-=sv; else R_st[Ig+row]+=sv;
            }

            if (tran) turb = true;
        }
    }

    for (int k=0;k<RXsize;++k) R_st[k] *= isol.sstag_ue[0];
    std::vector<int> colIdx;
    findColumnIndices<Real>(glob, 4*isol.stagIndex[0]+3, colIdx);
    addColumnValues<Real>(glob, 4*isol.stagIndex[0]+3, colIdx, R_st.data(), RXsize);

    if (std::abs(isol.sstag_ue[0]) > 1e-14) {
        Real scale = isol.sstag_ue[1] / isol.sstag_ue[0];
        for (int k=0;k<RXsize;++k) R_st[k] *= scale;
        std::vector<int> colIdx2;
        findColumnIndices<Real>(glob, 4*isol.stagIndex[1]+3, colIdx2);
        addColumnValues<Real>(glob, 4*isol.stagIndex[1]+3, colIdx2, R_st.data(), RXsize);
    }
}
