#ifndef MAIN_FUNCS_H
#define MAIN_FUNCS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.h"

struct Isol;
struct Vsol;
struct Foil;
struct Param;
struct Post;
struct Oper;
struct Geom;
struct Wake;
struct Glob;
struct Trans;

int colMajorIndex(int row, int col, int num_rows);

void make_panels(const Real (&inCoords)[2*Nin], Real (&outCoords)[2*Ncoords],const Real Ufac,const Real TEfac);

void build_gamma_codi(Isol &isol, const Foil& foil, const Oper& oper);

void init_thermo(const Oper& oper,Param& param,const Geom& geom);

void build_wake(const Foil& foil, const Geom& geom, const Oper& op, Isol& isol, Wake& wake);

void stagpoint_find(Isol& isol, const Foil& foil,const Wake&wake);

void identify_surfaces(const Isol& isol, Vsol& vsol);

void set_wake_gap(const Foil&foil,const Isol&isol,Vsol&vsol);

void calc_ue_m(const Foil&foil,const Wake&wake,Isol&isol,Vsol &vsol);

void rebuild_ue_m(const Foil&foil,const Wake&wake,const Isol&isol,Vsol&vsol,bool realloc);

void init_boundary_layer(const Oper&oper, const Foil&foil, Param&param, Isol&isol, Vsol&vsol, Glob&glob, Trans&tdata, const bool force);

void init_boundary_layer_from_xfoil(const Oper&oper, const Foil&foil, const Param&param, Isol&isol, Vsol&vsol, Glob&glob);

void stagpoint_move(Isol& isol,Glob& glob,const Foil& foil,const Wake& wake,Vsol&vsol);

void build_glob_RV(const Foil&foil, const Vsol&vsol,const Isol&isol,Glob&glob, Param&param, Trans&tdata);

void solve_glob(const Foil&foil, const Isol&isol, Glob& glob, Vsol& vsol, const Oper& oper, const int doSolve);

void update_state(const Oper&oper, const Param&param, Glob&glob,Vsol&vsol);

void update_transition(Glob &glob, Vsol &vsol, Isol &isol, Param&param, Trans&tdata,const bool force);

void clear_RV(Glob&glob, const Isol&isol,const Vsol&vsol, const Foil&foil,const Param&param);

bool solve_coupled(const Oper& oper, const Foil& foil, const Wake& wake,
    Param& param, Vsol& vsol, Isol& isol, Glob& glob, Trans&tdata, const bool force);

void calc_force(const Oper&op, const Geom&geom, const Param&par, const Isol&isol,const Foil&foil, const Glob&glob, Post& post);

void turbulent_BL_profile_XFOIL(
    Real theta,        // Momentum thickness [m]
    Real Ue,           // Edge velocity [m/s]
    Real utau,         // Friction velocity [m/s]
    Real nu,           // Kinematic viscosity [m^2/s]
    Real& delta                    // Output: BL thickness (where u/Ue = 0.99)
);


void interpolate_at_95_both_surfaces(const Real* xcoords, const Real* states, const Real*Cps, const Oper&oper, const Vsol&vsol, const Param&param,
    Real (&topBLStates)[7],Real (&botBLStates)[7],const Real Uinf, const Real x_target, const Real chordScale);


Real calc_OASPL(const Real* botStates, const Real* topStates, const Real chordScale, const Real Uinf,
    const Real X,const Real Y,const Real Z, const Real S, const Real nu, const Real rho,
    const int WPSjson,const std::string& model);


#endif