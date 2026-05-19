#ifndef MAIN_FUNCS_H
#define MAIN_FUNCS_H

#include <iostream>
#include <vector>
#include <cmath>
#include "real_type.h"
#include "panel_funcs.hpp"

// Forward declarations — needed because main_func.h is #included inside data_structs.h
// before the struct bodies are defined.
struct Isol; struct Glob;
struct Geom; struct Foil; struct Wake; struct Oper; struct Param;
struct Vsol; struct Post;
// Trans removed: forced transition deleted.

// Functions that remain as separate .cpp translation units:
void init_boundary_layer(const Oper&oper, const Foil&foil, Param&param, Isol&isol, Vsol&vsol, Glob&glob);
void init_boundary_layer_from_xfoil(const Oper&oper, const Foil&foil, const Param&param, Isol&isol, Vsol&vsol, Glob&glob);
void build_glob_RV(const Foil&foil, const Vsol&vsol, const Isol&isol, Glob&glob, Param&param);
void solve_glob(const Foil&foil, const Isol&isol, Glob& glob, Vsol& vsol, const Oper& oper, const int doSolve);
void update_state(const Oper&oper, const Param&param, Glob&glob, Vsol&vsol);
void update_transition(Glob &glob, Vsol &vsol, Isol &isol, Param&param);
bool solve_coupled(const Oper& oper, const Foil& foil, const Wake& wake,
    Param& param, Vsol& vsol, Isol& isol, Glob& glob);

// Everything else is now a template in a shared .hpp:
//   colMajorIndex       → inline in real_type.h
//   build_gamma_codi    → inline wrapper in solve_inv.hpp
//   init_thermo         → solver_funcs.hpp
//   build_wake_impl     → solver_funcs.hpp
//   stagpoint_find_impl → solver_funcs.hpp
//   identify_surfaces   → solver_funcs.hpp
//   set_wake_gap        → solver_funcs.hpp
//   rebuild_ue_m        → solver_funcs.hpp
//   stagpoint_move      → solver_funcs.hpp (AIRFOIL_STRUCTS_H guard)
//   clear_RV            → solver_funcs.hpp (AIRFOIL_STRUCTS_H guard)
//   calc_force          → solver_funcs.hpp

#endif
