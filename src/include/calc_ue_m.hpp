#pragma once

// Shared wake-influence-matrix assembly and solve.
//
// solve_sys_ue : CoDiPack-aware multi-column block solve.
//   IsolT is duck-typed — accepts Isol (fwd) or Isolc<Real> (AD);
//   both expose infMatrix[].  tape.isActive() selects passive vs
//   external-function path internally, so this one template covers both.
//
// compute_Dw : Dw = Cgam * Bp + Csig  (explicit loop, no cnp:: dependency).
//
// calc_ue_m : full wake influence matrix assembly.
//   All struct params duck-typed for Foil/Foil<Real>, Wake/Wake<Real>,
//   Isol/Isolc<Real>, Vsol/Vsol<Real>.
//
// rebuild_ue_m is AD-specific (uses Isolv<Real>::edgeVelSign) and lives
// in srcAD/include/main_func.hpp.
//
// Callers must include real_type.h/hpp, data_structs, and panel_funcs
// before this header.

#include <Eigen/Dense>
#include <codi.hpp>
#include <vector>
#include "panel_funcs.hpp"
// Callers must include the appropriate panel_funcs declaration header:
//   fwd: panel_funcs.h  (declares dvelocity_dgamma for Foil&)
//   AD:  panel_funcs_foil.hpp  (template dvelocity_dgamma for Foil<Real>&)
// Do NOT include panel_funcs.h here — it pulls in data_structs.h with
// non-template Foil/Isol which would conflict with the AD template structs.

// ── CoDiPack external-function data for the multi-column block solve ──────────

template<typename Active>
struct ImplicitBlockSolveData {
  using Real         = typename Active::Real;
  using Identifier   = typename Active::Identifier;
  using Tape         = typename Active::Tape;

  int n; // Ncoords+1
  int m; // Npanels

  Eigen::MatrixXd A;                 // n x n
  Eigen::MatrixXd X;                 // n x m  (solution Bp_full)
  std::vector<Identifier> A_id;      // n*n
  std::vector<Identifier> B_id;      // n*m
  std::vector<Identifier> X_id;      // (n-1)*m

  ImplicitBlockSolveData(int n_, int m_)
    : n(n_), m(m_), A(n_,n_), X(n_,m_),
      A_id(n_*n_), B_id(n_*m_), X_id((n_-1)*m_) {}
};

template<typename Active>
static void implicit_block_solve_b(typename Active::Tape*, void* d,
    codi::VectorAccessInterface<typename Active::Real,
                                typename Active::Identifier>* adj)
{
  using Real       = typename Active::Real;
  using Identifier = typename Active::Identifier;

  auto* data = static_cast<ImplicitBlockSolveData<Active>*>(d);
  const int n = data->n, m = data->m;

  Eigen::PartialPivLU<Eigen::MatrixXd> luAT(data->A.transpose());
  const size_t maxDim = adj->getVectorSize();

  for (size_t dim = 0; dim < maxDim; ++dim) {
    Eigen::MatrixXd X_b = Eigen::MatrixXd::Zero(n, m);
    for (int col = 0; col < m; ++col)
      for (int row = 0; row < n-1; ++row) {
        int k = colMajorIndex(row, col, n-1);
        Identifier id = data->X_id[k];
        X_b(row, col) = adj->getAdjoint(id, dim);
        adj->resetAdjoint(id, dim);
      }

    Eigen::MatrixXd Lambda = luAT.solve(X_b);

    for (int col = 0; col < m; ++col)
      for (int row = 0; row < n; ++row) {
        int k = colMajorIndex(row, col, n);
        adj->updateAdjoint(data->B_id[k], dim, Real(-Lambda(row, col)));
      }

    Eigen::MatrixXd dA = -(Lambda * data->X.transpose());
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i) {
        int k = colMajorIndex(i, j, n);
        adj->updateAdjoint(data->A_id[k], dim, Real(dA(i, j)));
      }
  }
}

template<typename Active>
static void implicit_block_solve_delete(typename Active::Tape*, void* d)
{
  delete static_cast<ImplicitBlockSolveData<Active>*>(d);
}

// ── solve_sys_ue ──────────────────────────────────────────────────────────────
// IsolT duck-typed: works for Isol (fwd, tape always inactive) and
// Isolc<Real> (AD, tape may be active).

template<typename Real, typename IsolT>
void solve_sys_ue(IsolT& isol, const Real* RHS, Real* Bp)
{
  constexpr int n = Ncoords + 1;
  constexpr int m = Ncoords + Nwake - 2;

  using Active = Real;
  using Tape   = typename Active::Tape;
  Tape& tape   = Active::getTape();

  auto* data = new ImplicitBlockSolveData<Active>(n, m);

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i) {
      int k = colMajorIndex(i, j, n);
      const Active& a = isol.infMatrix[k];
      data->A(i,j)   = a.getValue();
      data->A_id[k]  = a.getIdentifier();
    }

  Eigen::MatrixXd Bval(n, m);
  for (int col = 0; col < m; ++col)
    for (int row = 0; row < n; ++row) {
      int k = colMajorIndex(row, col, n);
      const Active& b = RHS[k];
      Bval(row, col) = b.getValue();
      data->B_id[k]  = b.getIdentifier();
    }

  Eigen::PartialPivLU<Eigen::MatrixXd> lu(data->A);
  data->X = lu.solve(-Bval);

  for (int col = 0; col < m; ++col)
    for (int row = 0; row < Ncoords; ++row) {
      int outk = colMajorIndex(row, col, Ncoords);
      Bp[outk] = data->X(row, col);
    }

  if (!tape.isActive()) { delete data; return; }

  for (int col = 0; col < m; ++col)
    for (int row = 0; row < Ncoords; ++row) {
      int outk = colMajorIndex(row, col, Ncoords);
      (void)tape.registerExternalFunctionOutput(Bp[outk]);
      data->X_id[outk] = Bp[outk].getIdentifier();
    }

  tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(
    &implicit_block_solve_b<Active>, data,
    &implicit_block_solve_delete<Active>));
}

// ── compute_Dw ────────────────────────────────────────────────────────────────
// Dw[Nwake x nPanels] = Cgam[Nwake x Ncoords] * Bp[Ncoords x nPanels] + Csig
// Uses explicit loops to avoid cnp:: template-ordering issues between src and srcAD.

template<typename Real>
void compute_Dw(const Real* Cgam, const Real* Bp, const Real* Csig, Real* Dw)
{
  constexpr int nPanels = Nwake + Ncoords - 2;
  for (int i = 0; i < Nwake; ++i)
    for (int j = 0; j < nPanels; ++j) {
      Real s = 0;
      for (int k = 0; k < Ncoords; ++k)
        s += Cgam[colMajorIndex(i,k,Nwake)] * Bp[colMajorIndex(k,j,Ncoords)];
      Dw[colMajorIndex(i,j,Nwake)] = s + Csig[colMajorIndex(i,j,Nwake)];
    }
}

// ── calc_ue_m ─────────────────────────────────────────────────────────────────
// All struct params duck-typed.  Real must be specified explicitly at the call
// site (e.g. calc_ue_m<Real>(foil, wake, isol, vsol)) because it cannot be
// deduced from the template struct arguments alone.

template<typename Real, typename FoilT, typename WakeT, typename IsolT, typename VsolT>
void calc_ue_m(const FoilT& foil, const WakeT& wake, IsolT& isol, VsolT& vsol)
{
    Real Cgam[Nwake * Ncoords] = {0.0};
    for (int i = 0; i < Nwake; ++i) {
        Real V_G[2*Ncoords] = {0.0};
        dvelocity_dgamma(foil,
            wake.x[colMajorIndex(0,i,2)],
            wake.x[colMajorIndex(1,i,2)], V_G);
        for (int j = 0; j < Ncoords; ++j) {
            Real V_Gx = V_G[colMajorIndex(0, j, 2)];
            Real V_Gy = V_G[colMajorIndex(1, j, 2)];
            Cgam[colMajorIndex(i, j, Nwake)] =
                V_Gx*wake.t[colMajorIndex(0,i,2)] +
                V_Gy*wake.t[colMajorIndex(1,i,2)];
        }
    }

    constexpr int nPanels = Ncoords + Nwake - 2;
    Real B[(Ncoords+1) * nPanels] = {0.0};

    PanelInfo<Real> info;
    for (int i = 0; i < Ncoords; ++i) {
        Real foilPoint[2] = {foil.x[colMajorIndex(0,i,2)],
                             foil.x[colMajorIndex(1,i,2)]};

        for (int j = 0; j < Ncoords-1; ++j) {
            Real a;
            panel_constsource_stream(
                foil.x[colMajorIndex(0,j,2)],   foil.x[colMajorIndex(1,j,2)],
                foil.x[colMajorIndex(0,j+1,2)], foil.x[colMajorIndex(1,j+1,2)],
                foilPoint[0], foilPoint[1], info, a);
            B[colMajorIndex(i,j,Ncoords+1)] = a;
        }

        for (int j = 0; j < Nwake-1; ++j) {
            Real wL[2] = {wake.x[colMajorIndex(0,j,2)],   wake.x[colMajorIndex(1,j,2)]};
            Real wR[2] = {wake.x[colMajorIndex(0,j+1,2)], wake.x[colMajorIndex(1,j+1,2)]};
            Real wM[2] = {0.5*(wL[0]+wR[0]), 0.5*(wL[1]+wR[1])};
            if (j == Nwake-2) {
                wR[0] = 2*wR[0]-wM[0];
                wR[1] = 2*wR[1]-wM[1];
            }
            Real a, b;
            panel_linsource_stream(wL[0],wL[1],wM[0],wM[1],foilPoint[0],foilPoint[1],info,a,b);
            if (j > 0) {
                B[colMajorIndex(i,Ncoords-1+j,  Ncoords+1)] += 0.5*a + b;
                B[colMajorIndex(i,Ncoords-1+j-1,Ncoords+1)] += 0.5*a;
            } else {
                B[colMajorIndex(i,Ncoords-1+j,Ncoords+1)] += b;
            }
            panel_linsource_stream(wM[0],wM[1],wR[0],wR[1],foilPoint[0],foilPoint[1],info,a,b);
            B[colMajorIndex(i,Ncoords-1+j,Ncoords+1)] += a + 0.5*b;
            if (j < Nwake-2)
                B[colMajorIndex(i,Ncoords+j,Ncoords+1)] += 0.5*b;
            else
                B[colMajorIndex(i,Ncoords-1+j,Ncoords+1)] += 0.5*b;
        }
    }

    Real Bp[Ncoords * nPanels];
    solve_sys_ue<Real>(isol, B, Bp);

    Real Csig[Nwake * nPanels] = {0.0};
    Real a1,a2,b1,b2;
    for (int i = 0; i < Nwake; ++i) {
        Real xi[2] = {wake.x[colMajorIndex(0,i,2)], wake.x[colMajorIndex(1,i,2)]};
        Real ti[2] = {wake.t[colMajorIndex(0,i,2)], wake.t[colMajorIndex(1,i,2)]};
        int jstart = (i==0) ? 1 : 0;
        int jend   = (i==0) ? Ncoords-2 : Ncoords-1;
        for (int j = jstart; j < jend; ++j) {
            panel_constsource_velocity(
                foil.x[colMajorIndex(0,j,2)],   foil.x[colMajorIndex(1,j,2)],
                foil.x[colMajorIndex(0,j+1,2)], foil.x[colMajorIndex(1,j+1,2)],
                xi[0],xi[1],info,a1,a2);
            Csig[colMajorIndex(i,j,Nwake)] = a1*ti[0] + a2*ti[1];
        }
        for (int j = 0; j < Nwake; ++j) {
            int I0 = std::max(j-1, 0);
            int I1 = j;
            int I2 = std::min(j+1, Nwake-1);
            Real leftX = 0.5*(wake.x[colMajorIndex(0,I0,2)]+wake.x[colMajorIndex(0,I1,2)]);
            Real leftY = 0.5*(wake.x[colMajorIndex(1,I0,2)]+wake.x[colMajorIndex(1,I1,2)]);
            Real midX  =      wake.x[colMajorIndex(0,I1,2)];
            Real midY  =      wake.x[colMajorIndex(1,I1,2)];
            Real rightX = 0.5*(midX + wake.x[colMajorIndex(0,I2,2)]);
            Real rightY = 0.5*(midY + wake.x[colMajorIndex(1,I2,2)]);
            if (j == Nwake-1) { rightX = 2*midX-leftX; rightY = 2*midY-leftY; }
            Real lL[2] = {midX-leftX,   midY-leftY};
            Real rL[2] = {rightX-midX,  rightY-midY};
            Real d1 = norm2(lL);
            Real d2 = norm2(rL);
            if (i == j) {
                if (j == 0) {
                    Real loP[2] = {foil.x[colMajorIndex(0,1,2)]-foil.x[colMajorIndex(0,0,2)],
                                   foil.x[colMajorIndex(1,1,2)]-foil.x[colMajorIndex(1,0,2)]};
                    Real upP[2] = {foil.x[colMajorIndex(0,Ncoords-1,2)]-foil.x[colMajorIndex(0,Ncoords-2,2)],
                                   foil.x[colMajorIndex(1,Ncoords-1,2)]-foil.x[colMajorIndex(1,Ncoords-2,2)]};
                    Real dl = norm2(loP), du = norm2(upP);
                    Csig[colMajorIndex(i,0,Nwake)]         += (0.5/M_PI)*(std::log(dl/d2)+1.0);
                    Csig[colMajorIndex(i,Ncoords-2,Nwake)] += (0.5/M_PI)*(std::log(du/d2)+1.0);
                    Csig[colMajorIndex(i,Ncoords-1,Nwake)] += -0.5/M_PI;
                } else if (j == Nwake-1) {
                    // self effect = 0
                } else {
                    Real aa = (0.25/M_PI)*std::log(d1/d2);
                    Csig[colMajorIndex(i,Ncoords-2+j,Nwake)] += aa + 0.5/M_PI;
                    Csig[colMajorIndex(i,Ncoords-1+j,Nwake)] += aa - 0.5/M_PI;
                }
            } else {
                if (j == 0) {
                    panel_linsource_velocity(midX,midY,rightX,rightY,xi[0],xi[1],info,a1,b1,a2,b2);
                    Real a = a1*ti[0]+a2*ti[1], b = b1*ti[0]+b2*ti[1];
                    Csig[colMajorIndex(i,Ncoords-1,Nwake)] += b;
                    Csig[colMajorIndex(i,0,Nwake)]         += a;
                    Csig[colMajorIndex(i,Ncoords-2,Nwake)] += a;
                } else if (j == Nwake-1) {
                    panel_constsource_velocity(leftX,leftY,rightX,rightY,xi[0],xi[1],info,a1,a2);
                    Csig[colMajorIndex(i,Ncoords+Nwake-3,Nwake)] += (a1*ti[0]+a2*ti[1]);
                } else {
                    panel_linsource_velocity(leftX,leftY,midX,midY,xi[0],xi[1],info,a1,b1,a2,b2);
                    Csig[colMajorIndex(i,Ncoords-2+j,Nwake)] +=
                        (a1*ti[0]+a2*ti[1]) + 0.5*(b1*ti[0]+b2*ti[1]);
                    panel_linsource_velocity(midX,midY,rightX,rightY,xi[0],xi[1],info,a1,b1,a2,b2);
                    Csig[colMajorIndex(i,Ncoords-1+j,Nwake)] +=
                        0.5*(a1*ti[0]+a2*ti[1]) + (b1*ti[0]+b2*ti[1]);
                }
            }
        }
    }

    Real Dw[Nwake * nPanels];
    compute_Dw<Real>(Cgam, Bp, Csig, Dw);
    for (int j = 0; j < nPanels; ++j)
        Dw[colMajorIndex(0,j,Nwake)] = Bp[colMajorIndex(Ncoords-1,j,Ncoords)];

    for (int i = 0; i < Ncoords; ++i)
        for (int j = 0; j < nPanels; ++j)
            vsol.ue_sigma[colMajorIndex(i,j,Ncoords+Nwake)] =
                Bp[colMajorIndex(i,j,Ncoords)];

    for (int i = 0; i < Nwake; ++i)
        for (int j = 0; j < nPanels; ++j)
            vsol.ue_sigma[colMajorIndex(Ncoords+i,j,Ncoords+Nwake)] =
                Dw[colMajorIndex(i,j,Nwake)];
}
