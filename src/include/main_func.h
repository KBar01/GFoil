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

void build_gamma_codi(Isol &isol, const Foil& foil, const Oper& oper);

void init_thermo(const Oper& oper,Param& param,const Geom& geom);

void build_wake(const Foil& foil, const Geom& geom, const Oper& op, Isol& isol, Wake& wake);

void stagpoint_find(Isol& isol, const Foil& foil,const Wake&wake);

void identify_surfaces(const Isol& isol, Vsol& vsol);

void set_wake_gap(const Foil&foil,const Isol&isol,Vsol&vsol);

// rebuild_ue_m → template in solver_funcs.hpp

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

// calc_force → template in solver_funcs.hpp



// interpolate_at_95_both_surfaces is now template<typename Real, typename OperT, typename ParamT>
// in src/include/extract_BL_TE.hpp — include that header directly where needed.
// calc_OASPL is now template<typename Real, bool WriteJSON> in src/include/sound.hpp.
// calc_ue_m is now template<typename Real, ...> in src/include/calc_ue_m.hpp.

#endif