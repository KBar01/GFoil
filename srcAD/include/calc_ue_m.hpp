#pragma once

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <codi.hpp>
#include "real_type.hpp"
#include "data_structs.hpp"
#include "vector_ops.hpp"
#include "panel_funcs.hpp"




template<typename Real> struct Isolc;
template<typename Real> struct Isolv;
template<typename Real> struct Vsol;
template<typename Real> struct Foil;
template<typename Real> struct Param;
template<typename Real> struct Post;
template<typename Real> struct Oper;
template<typename Real> struct Geom;
template<typename Real> struct Wake;
template<typename Real> struct Glob;
template<typename Real> struct Trans;


template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
/*
template<typename Type>
void func_ue(Matrix<Type> const& A, Vector<Type> const& rhs, Vector<Type> & sol) {
  
    auto before = std::chrono::high_resolution_clock::now();
    auto qr = A.colPivHouseholderQr();
    auto afterA = std::chrono::high_resolution_clock::now();
    sol = qr.solve(rhs);
    auto aftersolve = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elap1 = afterA - before;
    std::cout << "qr " << elap1.count() << " seconds\n";
    std::chrono::duration<double> elap2 = aftersolve - afterA;
    std::cout << "solve " << elap2.count() << " seconds\n";
}

template<typename Number>
struct EigenSolver_ue : public codi::EigenLinearSystem<Number, Matrix, Vector> {
  public:
 
    using Base = codi::EigenLinearSystem<Number, Matrix, Vector>;  
    using MatrixReal = typename Base::MatrixReal;                                   
    using VectorReal = typename Base::VectorReal; 

    void doQr(*A)
    
    
    void solveSystem(MatrixReal const* A, VectorReal const* b, VectorReal* x) {
        func_ue(*A, *b, *x);
    }
};

template<typename Real>
void solve_sys_ue(Isolc<Real>&isolc, const Real* RHS, Real*Bp) {
    
    constexpr int Nsize = Ncoords + 1;  // Matrix dimensions
    constexpr int Npanels = Ncoords+Nwake-2;
    
    // Wrap the flattened column-major AIC vector into an Eigen matrix
    Matrix<Real> A(Nsize, Nsize);
    Vector<Real> rhsEigen(Ncoords+1);
    Vector<Real> sol(Ncoords+1);
    
    for (int i = 0; i < Nsize; ++i) {
        for (int j = 0; j < Nsize; ++j) {
            A(i, j) = isolc.infMatrix[colMajorIndex(i, j, Nsize)];
        }
    }

    for (int col = 0; col < Npanels; ++col) {
        
        
        for (int row = 0; row < Nsize; ++row) {
            rhsEigen(row) = RHS[colMajorIndex(row, col, Nsize)];
        }

        codi::solveLinearSystem(EigenSolver_ue<Real>(), A, rhsEigen, sol);

        for (int row = 0; row < Ncoords; ++row) {
            Bp[colMajorIndex(row, col, Ncoords)] =-1*sol(row);
        }
    }
}




template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename Number>
struct CachedEigenQRSystem
  : public codi::EigenLinearSystem<Number, Matrix, Vector>
{
  using Base = codi::EigenLinearSystem<Number, Matrix, Vector>;
  using MatrixReal = typename Base::MatrixReal;
  using VectorReal = typename Base::VectorReal;

  // Separate caches
  mutable std::shared_ptr<Eigen::ColPivHouseholderQR<MatrixReal>> qrA;
  mutable std::shared_ptr<Eigen::ColPivHouseholderQR<MatrixReal>> qrAT;

  mutable MatrixReal const* lastA  = nullptr;
  mutable MatrixReal const* lastAT = nullptr;

 
  void solveSystem(MatrixReal const* A,
                   VectorReal const* b,
                   VectorReal* x)
  {
    // Case 1: this is A
    if (A == lastA || lastA == nullptr) {
      if (!qrA || A != lastA) {
        qrA = std::make_shared<Eigen::ColPivHouseholderQR<MatrixReal>>(*A);
        lastA = A;
      }
      *x = qrA->solve(*b);
      return;
    }

    // Case 2: this is A^T
    if (A == lastAT || lastAT == nullptr) {
      if (!qrAT || A != lastAT) {
        qrAT = std::make_shared<Eigen::ColPivHouseholderQR<MatrixReal>>(*A);
        lastAT = A;
      }
      *x = qrAT->solve(*b);
      return;
    }

    // Fallback (should never happen)
    Eigen::ColPivHouseholderQR<MatrixReal> qrTmp(*A);
    *x = qrTmp.solve(*b);
  }

  MatrixReal* transposeMatrix(MatrixReal* A) {
    return new MatrixReal(A->transpose());
  }

  void solveSystemPrimal(MatrixReal const* A,
                         VectorReal const* b,
                         VectorReal* x)
  {
    solveSystem(A, b, x);
  }
};

template<typename Real>
void solve_sys_ue(Isolc<Real>&isolc, const Real* RHS, Real*Bp) {
    
    constexpr int Nsize = Ncoords + 1;  // Matrix dimensions
    constexpr int Npanels = Ncoords+Nwake-2;
    
    // Wrap the flattened column-major AIC vector into an Eigen matrix
    Matrix<Real> A(Nsize, Nsize);
    Vector<Real> rhsEigen(Ncoords+1);
    Vector<Real> sol(Ncoords+1);
    


    for (int i = 0; i < Nsize; ++i) {
        for (int j = 0; j < Nsize; ++j) {
            A(i, j) = isolc.infMatrix[colMajorIndex(i, j, Nsize)];
        }
    }

    for (int col = 0; col < Npanels; ++col) {
        
        
        for (int row = 0; row < Nsize; ++row) {
            rhsEigen(row) = RHS[colMajorIndex(row, col, Nsize)];
        }
        
        auto befFunc = std::chrono::high_resolution_clock::now();
        
        codi::solveLinearSystem(
            CachedEigenQRSystem<Real>(),
            A,
            rhsEigen,
            sol
        );
        auto aftFunc = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedt7 = aftFunc - befFunc;
        std::cout << "newfunc" << elapsedt7.count() << " seconds\n";

        for (int row = 0; row < Ncoords; ++row) {
            Bp[colMajorIndex(row, col, Ncoords)] =-1*sol(row);
        }
    }
}

*/


template<typename Active>
struct ImplicitBlockSolveData {
  using Real       = typename Active::Real;
  using Identifier = typename Active::Identifier;
  using Tape       = typename Active::Tape;
  using VectorAccess = codi::VectorAccessInterface<Real, Identifier>;

  int n; // Ncoords+1
  int m; // Npanels

  // Stored primal values (passive)
  Eigen::MatrixXd A;   // n x n
  Eigen::MatrixXd X;   // n x m   (solution, i.e. Bp_full)

  // Identifiers for inputs/outputs (needed to read/write adjoints)
  std::vector<Identifier> A_id;   // n*n
  std::vector<Identifier> B_id;   // n*m
  std::vector<Identifier> X_id;   // (Ncoords)*m  (your output storage size)

  ImplicitBlockSolveData(int n_, int m_)
    : n(n_), m(m_), A(n_, n_), X(n_, m_),
      A_id(n_*n_), B_id(n_*m_), X_id((n_-1)*m_) {}
};


// Reverse callback:
// Given adjoints of X (Bp) -> compute adjoints of A and B using implicit rule
template<typename Active>
static void implicit_block_solve_b(typename Active::Tape* /*tape*/, void* d,
                                   codi::VectorAccessInterface<typename Active::Real, typename Active::Identifier>* adj)
{
  using Real       = typename Active::Real;
  using Identifier = typename Active::Identifier;

  auto* data = static_cast<ImplicitBlockSolveData<Active>*>(d);

  const int n = data->n;
  const int m = data->m;

  const size_t maxDim = adj->getVectorSize();

  // Factorization of A^T (passive) - reuse per dimension
  Eigen::PartialPivLU<Eigen::MatrixXd> luAT(data->A.transpose());

  for (size_t dim = 0; dim < maxDim; ++dim) {

    // 1) Gather X_b (n x m) from tape adjoints.
    // Your output only stores first (n-1)=Ncoords rows; last row is treated as 0 adjoint.
    Eigen::MatrixXd X_b = Eigen::MatrixXd::Zero(n, m);

    for (int col = 0; col < m; ++col) {
      for (int row = 0; row < n-1; ++row) {
        const int k = colMajorIndex(row, col, n-1);
        Identifier id = data->X_id[k];
        X_b(row, col) = adj->getAdjoint(id, dim);
        adj->resetAdjoint(id, dim); // x_b = 0 in CoDiPack’s convention
      }
    }

    // 2) Solve adjoint system:  A^T * Lambda = X_b
    Eigen::MatrixXd Lambda = luAT.solve(X_b); // n x m

    // 3) Update adjoints:
    //    B_b += -Lambda   (because A X = -B)
    for (int col = 0; col < m; ++col) {
      for (int row = 0; row < n; ++row) {
        const int k = colMajorIndex(row, col, n);
        Identifier id = data->B_id[k];
        adj->updateAdjoint(id, dim, Real(-Lambda(row, col)));
      }
    }

    //    A_b += -(Lambda * X^T)
    //    (A is n x n, Lambda n x m, X n x m => Lambda*X^T is n x n)
    Eigen::MatrixXd dA = -(Lambda * data->X.transpose());

    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        const int k = colMajorIndex(i, j, n);
        Identifier id = data->A_id[k];
        adj->updateAdjoint(id, dim, Real(dA(i, j)));
      }
    }
  }
}

template<typename Active>
static void implicit_block_solve_delete(typename Active::Tape* /*tape*/, void* d)
{
  delete static_cast<ImplicitBlockSolveData<Active>*>(d);
}


template<typename Real>
void solve_sys_ue(Isolc<Real>& isolc, const Real* RHS, Real* Bp)
{
  constexpr int n = Ncoords + 1;
  constexpr int m = Ncoords + Nwake - 2;  // nPanels

  using Active = Real;
  using Tape   = typename Active::Tape;

  Tape& tape = Active::getTape(); // this exists for CoDiPack ActiveType

  // Allocate external-function data (only used if tape is active)
  auto* data = new ImplicitBlockSolveData<Active>(n, m);

  // 1) Extract primal A (n x n) and store identifiers
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      const int k = colMajorIndex(i, j, n);
      const Active& a = isolc.infMatrix[k];
      data->A(i, j) = a.getValue();
      data->A_id[k] = a.getIdentifier();
    }
  }

  // 2) Extract primal B (n x m) and store identifiers
  Eigen::MatrixXd Bval(n, m);
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < n; ++row) {
      const int k = colMajorIndex(row, col, n);
      const Active& b = RHS[k];
      Bval(row, col) = b.getValue();
      data->B_id[k]  = b.getIdentifier();
    }
  }

  // 3) Passive solve: A * X = -B
  Eigen::PartialPivLU<Eigen::MatrixXd> lu(data->A);
  data->X = lu.solve(-Bval); // n x m

  // 4) Write primal output values into Bp (Ncoords x m)
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < Ncoords; ++row) {
      const int outk = colMajorIndex(row, col, Ncoords);
      Bp[outk] = data->X(row, col); // primal assignment
    }
  }

  // 5) If tape inactive, we're done
  if (!tape.isActive()) {
    delete data;
    return;
  }

  // 6) Register outputs on tape + store their identifiers
  //    IMPORTANT: CoDiPack needs outputs registered as "external function outputs"
  for (int col = 0; col < m; ++col) {
    for (int row = 0; row < Ncoords; ++row) {
      const int outk = colMajorIndex(row, col, Ncoords);
      // registerExternalFunctionOutput mutates the variable's identifier bookkeeping
      (void)tape.registerExternalFunctionOutput(Bp[outk]);
      data->X_id[outk] = Bp[outk].getIdentifier();
    }
  }

  // 7) Push external function with reverse callback
  tape.pushExternalFunction(
    codi::ExternalFunction<Tape>::create(
      /*reverse*/ &implicit_block_solve_b<Active>,
      /*data*/    data,
      /*del*/     &implicit_block_solve_delete<Active>
      // forward/primal callbacks omitted => CoDiPack uses stored primals already
    )
  );
}



template<typename Real>
void compute_Dw(const Real* Cgam, const Real* Bp, const Real* Csig, Real* Dw){

    constexpr int nPanels = Nwake+Ncoords-2 ;
    cnp::matmat_mul<Real,Nwake,Ncoords,nPanels>(Cgam,Bp,Dw);
    cnp::add_inplace<Real,(Nwake*nPanels)>(Dw,Csig);
}

template<typename Real>
void calc_ue_m(const Foil<Real>&foil,const Wake<Real>&wake,Isolc<Real>&isolc,Vsol<Real> &vsol) {


    // Allocate Cgam: Nwake x Ncoords
    Real Cgam[Nwake * Ncoords] = {0.0}; 
    for (int i = 0; i < Nwake; ++i) {
        
        Real V_G[2*Ncoords] = {0.0};
        dvelocity_dgamma(foil,wake.x[colMajorIndex(0,i,2)],wake.x[colMajorIndex(1,i,2)],V_G);

        for (int j = 0; j < Ncoords; ++j) {
            Real V_Gx = V_G[colMajorIndex(0, j, 2)];
            Real V_Gy = V_G[colMajorIndex(1, j, 2)];
            Real Cgamx = V_Gx*wake.t[colMajorIndex(0, i, 2)];
            Real Cgamy = V_Gy*wake.t[colMajorIndex(1, i, 2)];
            Cgam[colMajorIndex(i, j, Nwake)] = Cgamx+Cgamy;
        }
    }

    constexpr int nPanels = Ncoords + Nwake - 2;

    // Allocate B: (Ncoords+1) x npan
    Real B[(Ncoords+1) * nPanels] = {0.0};

    PanelInfo<Real> info;
    for (int i = 0; i < Ncoords; ++i) {
        
        Real foilPoint[2] = {foil.x[colMajorIndex(0, i, 2)], foil.x[colMajorIndex(1, i, 2)]};

        // Airfoil panels
        for (int j = 0; j < Ncoords-1; ++j) {
            
            Real foilPanelx1 = foil.x[colMajorIndex(0, j, 2)];
            Real foilPanely1 = foil.x[colMajorIndex(1, j, 2)];
            Real foilPanelx2 = foil.x[colMajorIndex(0, j+1, 2)];
            Real foilPanely2 = foil.x[colMajorIndex(1, j+1, 2)];
            Real a;
            panel_constsource_stream(foilPanelx1,foilPanely1,foilPanelx2,foilPanely2,foilPoint[0],foilPoint[1],info,a);
            B[colMajorIndex(i, j, Ncoords+1)] = a;
        }

        // Wake panels
        for (int j = 0; j < Nwake-1; ++j) {
            
            Real wakePanelLeft[2] = {wake.x[colMajorIndex(0,j,2)],wake.x[colMajorIndex(1,j,2)]};  // 3 pts: left, mid, right
            Real wakePanelRight[2] = {wake.x[colMajorIndex(0,j+1,2)],wake.x[colMajorIndex(1,j+1,2)]};  // 3 pts: left, mid, right
            Real wakePanelMid[2] = {0.5*(wakePanelLeft[0]+wakePanelRight[0]),0.5*(wakePanelLeft[1]+wakePanelRight[1])};
            
            if (j == Nwake-2) {
    
                wakePanelRight[0] = 2*wakePanelRight[0]-wakePanelMid[0];  // ghost extension
                wakePanelRight[1] = 2*wakePanelRight[1]-wakePanelMid[1];
            }
            
            Real a, b;

            // Left half
            panel_linsource_stream(wakePanelLeft[0],wakePanelLeft[1],wakePanelMid[0],wakePanelMid[1],foilPoint[0],foilPoint[1],info,a,b);
            if (j > 0) {
                B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += 0.5*a + b;
                B[colMajorIndex(i, Ncoords-1+j-1, Ncoords+1)] += 0.5*a;
            } else {
                B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += b;
            }

            // Right half
            panel_linsource_stream(wakePanelMid[0],wakePanelMid[1],wakePanelRight[0],wakePanelRight[1],foilPoint[0],foilPoint[1],info,a,b);
            B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += a + 0.5*b;
            if (j < Nwake-2) {
                B[colMajorIndex(i, Ncoords+j, Ncoords+1)] += 0.5*b;
            } else {
                B[colMajorIndex(i, Ncoords-1+j, Ncoords+1)] += 0.5*b;
            }
        }
    }  

    auto bfsol = std::chrono::high_resolution_clock::now();

    // Solve Bp = -AIC^{-1} * B (note B is (Ncoords+1)xnpan)
    Real Bp[Ncoords * nPanels];
    solve_sys_ue(isolc,B,Bp);
    
    auto afsolv = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedt = afsolv - bfsol;
    std::cout << "solve_sys_ue " << elapsedt.count() << " seconds\n";
    
    //for (int i=0;i<Ncoords;++i){
    //  std::cout << "Bp[" <<i<<"] = " << Bp[i] << std::endl;
    //}

    //std::cout << Bp[colMajorIndex(Ncoords-1,nPanels-1,Ncoords)] << std::endl;

    // Csig: Nwake x npan
    Real Csig[Nwake * nPanels]={0.0};

    Real a1,a2,b1,b2;
    for (int i = 0; i < Nwake; ++i) {
        
        Real xi[2] = {wake.x[colMajorIndex(0, i, 2)],wake.x[colMajorIndex(1, i, 2)]};
        Real ti[2] = {wake.t[colMajorIndex(0, i, 2)],wake.t[colMajorIndex(1, i, 2)]};

        // Constant sources on foil
        int jstart = (i == 0) ? 1 : 0;
        int jend   = (i == 0) ? Ncoords-2 : Ncoords-1;

        for (int j = jstart; j < jend; ++j) {
            
            panel_constsource_velocity(
                foil.x[colMajorIndex(0,j,2)],
                foil.x[colMajorIndex(1,j,2)],
                foil.x[colMajorIndex(0,j+1,2)],
                foil.x[colMajorIndex(1,j+1,2)],
                xi[0],xi[1],info,a1,a2);
            
            Csig[colMajorIndex(i, j, Nwake)] = a1*ti[0] + a2*ti[1];

        }
        // piecewise linear sources across wake panel halves (else singular)
        for (int j = 0; j < Nwake; ++j) {
            
            
            int I0 = std::max(j-1, 0);
            int I1 = j;
            int I2 = std::min(j+1, Nwake-1);
        
            Real leftX = 0.5*(wake.x[colMajorIndex(0,I0,2)] + wake.x[colMajorIndex(0,I1,2)]);
            Real leftY = 0.5*(wake.x[colMajorIndex(1,I0,2)] + wake.x[colMajorIndex(1,I1,2)]);
            Real midX  =      wake.x[colMajorIndex(0,I1,2)];
            Real midY  =      wake.x[colMajorIndex(1,I1,2)];
            Real rightX = 0.5*(midX + wake.x[colMajorIndex(0,I2,2)]);
            Real rightY = 0.5*(midY + wake.x[colMajorIndex(1,I2,2)]);

        
            if (j == Nwake - 1) {
                rightX = 2*midX - leftX;
                rightY = 2*midY - leftY;
            }
            
            Real leftLength[2] = {
                midX - leftX,
                midY - leftY
            };

            Real rightLength[2] = {
                rightX - midX,
                rightY - midY
            };

            Real d1 = norm2(leftLength); // norm between center and mid-left
            Real d2 = norm2(rightLength); // norm between mid-right and center
        
            if (i == j) {
                // first point: special TE system (three panels meet)
                if (j == 0) {
                    
                    Real lowPanelLength[2] = {
                        foil.x[colMajorIndex(0,1,2)]-foil.x[colMajorIndex(0,0,2)],
                        foil.x[colMajorIndex(1,1,2)]-foil.x[colMajorIndex(1,0,2)]
                    };

                    Real upPanelLength[2] = {
                        foil.x[colMajorIndex(0,Ncoords-1,2)]-foil.x[colMajorIndex(0,Ncoords-2,2)],
                        foil.x[colMajorIndex(1,Ncoords-1,2)]-foil.x[colMajorIndex(1,Ncoords-2,2)]
                    };

                    Real dl = norm2(lowPanelLength);
                    Real du = norm2(upPanelLength);
        
                    Csig[colMajorIndex(i, 0, Nwake)]         += (0.5/M_PI) * (std::log(dl/d2) + 1.0);
                    Csig[colMajorIndex(i, Ncoords-2, Nwake)] += (0.5/M_PI) * (std::log(du/d2) + 1.0);
                    Csig[colMajorIndex(i, Ncoords-1, Nwake)] += -0.5/M_PI;
                
                } else if (j == Nwake - 1) {
                    // self effect = 0
                } else {
                    Real aa = (0.25/M_PI) * std::log(d1/d2);
                    Csig[colMajorIndex(i, Ncoords-2+j, Nwake)] += aa + 0.5/M_PI;
                    Csig[colMajorIndex(i, Ncoords-1+j, Nwake)] += aa - 0.5/M_PI;
                }
            } else {
                if (j == 0) {
                    
                    
                    panel_linsource_velocity(midX,midY,rightX,rightY,xi[0],xi[1],info,a1,b1,a2,b2);

                    Real a = a1*ti[0] + a2*ti[1];
                    Real b = b1*ti[0] + b2*ti[1];

                    Csig[colMajorIndex(i, Ncoords-1, Nwake)]   += b;
                    Csig[colMajorIndex(i, 0, Nwake)]           += a;
                    Csig[colMajorIndex(i, Ncoords-2, Nwake)]   += a;

                } else if (j == Nwake-1) {
                    
                    panel_constsource_velocity(leftX,leftY,rightX,rightY,xi[0],xi[1],info,a1,a2); // Xj[:, [0,2]]
                    Csig[colMajorIndex(i, Ncoords+Nwake-3, Nwake)] += (a1*ti[0] + a2*ti[1]);

                } else {

                    panel_linsource_velocity(leftX,leftY,midX,midY,xi[0],xi[1],info,a1,b1,a2,b2);
                    Csig[colMajorIndex(i, Ncoords-2+j, Nwake)]     += (a1*ti[0] + a2*ti[1]) + 0.5 * (b1*ti[0] + b2*ti[1]);
                    
                    panel_linsource_velocity(midX,midY,rightX,rightY,xi[0],xi[1],info,a1,b1,a2,b2);   // Xj[:, [1,2]]
                    Csig[colMajorIndex(i, Ncoords - 1 + j, Nwake)] += 0.5 * (a1*ti[0] + a2*ti[1]) + (b1*ti[0] + b2*ti[1]);
                }
            }
        }
    }

    // Combine Dw = Cgam * Bp + Csig
    Real Dw[Nwake * nPanels];
    //compute_Dw_eigen(Cgam,Bp,Csig,Dw);  // Cgam [Nwake x Ncoords], Bp [Ncoords x npan], Dw [Nwake x npan]

    compute_Dw(Cgam,Bp,Csig,Dw);
    for (int j = 0; j < nPanels; ++j) {
        Dw[colMajorIndex(0, j, Nwake)] = Bp[colMajorIndex(Ncoords-1, j, Ncoords)];
    }


    for (int i=0; i<Ncoords;++i){
        for (int j=0; j<nPanels; ++j){
            vsol.ue_sigma[colMajorIndex(i,j,Ncoords+Nwake)] = Bp[colMajorIndex(i,j,Ncoords)];
        }
    }

    for (int i=0; i<Nwake;++i){
        for (int j=0; j<nPanels; ++j){
            vsol.ue_sigma[colMajorIndex(Ncoords+i,j,Ncoords+Nwake)] = Dw[colMajorIndex(i,j,Nwake)];
        }
    }

}



template<typename Real>
void rebuild_ue_m(const Foil<Real>&foil,const Wake<Real>&wake,const Isolv<Real>&isolv,Vsol<Real>&vsol,bool realloc){


    // Functions builds ue_m from ue_sigma, sigma_m and edgevelocity signs.
    // This looks more complex as taking advantage of known sparsity pattern
    // of sigma_m, therefore removing unnessessary computation.

    // TODO: mayb make iterators in loops more clear etc
    
    Real sigma_m[2*(Ncoords-1)];

    if (realloc){
        for (int i=0;i<(Ncoords+Nwake)*(Ncoords+Nwake);++i){
            vsol.ue_m[i] = 0;
        }
    }

    // Creation of sigma_m, holding non-zero vals only
    for (int i=0; i<Ncoords-1;++i){

        Real ds = foil.s[i+1] - foil.s[i];
        int ind = 2*i;
        sigma_m[ind] = -1.0*isolv.edgeVelSign[i] / ds;
        sigma_m[ind+1] = 1.0*isolv.edgeVelSign[i+1] / ds;
    }

    // first column in ue_m
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[row] = vsol.ue_sigma[row]*sigma_m[0];
    }

    // columns 1 to Ncoords-1
    int ue_mInt = Ncoords+Nwake;
    int sigma_mIndex = 1;
    for (int col=1; col<=Ncoords-2;++col){

        Real sigma_mValue = sigma_m[sigma_mIndex];

        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,col-1,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mValue = sigma_m[sigma_mIndex+1];
        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,col,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mIndex += 2;
        ue_mInt += (Ncoords+Nwake);
    }

    // Ncoords-1 column in ue_m (0 index)
    int start = (Ncoords+Nwake)*(Ncoords-1);
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[colMajorIndex(row,Ncoords-1,Ncoords+Nwake)] = vsol.ue_sigma[colMajorIndex(row,Ncoords-2,Ncoords+Nwake)]*sigma_m[2*(Ncoords-1)-1];
    }


    // Now contributions from wake section

    Real sigma_mWake[2*(Nwake-1)];

    for (int i=0; i<Nwake-1;++i){
        int ind = 2*i;
        Real ds = wake.s[i+1] - wake.s[i];
        sigma_mWake[ind] = -1.0 / ds;
        sigma_mWake[ind+1] = 1.0 / ds;
    }
    
    // Ncoord value

    // Ncoords column
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[colMajorIndex(row,Ncoords,Ncoords+Nwake)] = vsol.ue_sigma[colMajorIndex(row,Ncoords-1,Ncoords+Nwake)]*sigma_mWake[0];
    }


    // columns Ncoords+1 to end-1
    ue_mInt = (Ncoords+Nwake)*(Ncoords) + Ncoords+Nwake;
    sigma_mIndex = 1;
    for (int col=1; col<=Nwake-2;++col){

        Real sigma_mValue = sigma_mWake[sigma_mIndex];

        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,Ncoords+col-2,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mValue = sigma_mWake[sigma_mIndex+1];
        for (int row=0; row<Ncoords+Nwake; ++row){
            vsol.ue_m[ue_mInt+row] += vsol.ue_sigma[colMajorIndex(row,Ncoords+col-1,Ncoords+Nwake)]*sigma_mValue;
        }

        sigma_mIndex += 2;
        ue_mInt += (Ncoords+Nwake);
    }

    // Last coluumn
    for (int row=0; row<Ncoords+Nwake;++row){
        vsol.ue_m[(Ncoords+Nwake)*(Ncoords+Nwake-1)+row] = vsol.ue_sigma[colMajorIndex(row,Ncoords+Nwake-3,Ncoords+Nwake)]*sigma_mWake[2*(Nwake-1)-1];
    }


    // sgnue switching
    for (int row = 0; row < Ncoords; ++row) {
        if (isolv.edgeVelSign[row] == 1.0) {
            continue; // skip row — multiplying by 1 does nothing
        }
    
        for (int col = 0; col < Ncoords+Nwake; ++col) {
            vsol.ue_m[colMajorIndex(row,col,Ncoords+Nwake)] *= -1.0;
        }
    }

}