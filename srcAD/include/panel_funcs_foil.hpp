// inviscid_velocity and dvelocity_dgamma are now template<typename Real, typename FoilT>
// in src/include/panel_funcs.hpp (accessible to both builds via src/include/ on the
// include path). This file is retained as a thin forwarding header so any existing
// includers don't break.
#pragma once

#include "panel_funcs.hpp"
