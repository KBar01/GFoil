#include <iostream>
#include <cmath>
#include "real_type.h"
#include "data_structs.h"
#include "get_funcs.h"



void calc_force(const Oper&op, const Geom&geom, const Param&par, const Isol&isol,const Foil&foil, const Glob&glob, Post& post) {

    // Compute dynamic pressure
    Real qinf = 0.5*op.rho *(op.Vinf * op.Vinf);
    constexpr int Nsys = Ncoords+Nwake;
    
    // Calculate the pressure coefficient at each node
    Real ue[Nsys]={0};
    for (int i=0;i<Nsys;++i){ue[i] = glob.U[colMajorIndex(3,i,4)];}
    get_cp(post,op,par,ue);     // assigns cp in post struct

    // Pre-allocate variables for memory efficiency
    Real cos_alpha = std::cos(op.alpha);
    Real sin_alpha = std::sin(op.alpha);
    
    // Loop over panels for integration
    for (int i0 = 1; i0 <= Ncoords; ++i0) {
        
        int i = (i0 == Ncoords) ? 0 : i0;
        int ip = (i0 == Ncoords) ? Ncoords-1 : i0-1;

        const Real x1[2] = {foil.x[colMajorIndex(0,ip,2)],foil.x[colMajorIndex(1,ip,2)]}; // points
        const Real x2[2] = {foil.x[colMajorIndex(0,i ,2)],foil.x[colMajorIndex(1,i ,2)]};

        Real dxv[2] = {x2[0] - x1[0], x2[1] - x1[1]};
        Real dx1[2] = {x1[0] - geom.xref[0], x1[1] - geom.xref[1]};
        Real dx2[2] = {x2[0] - geom.xref[0], x2[1] - geom.xref[1]};

        Real dx1nds = dxv[0] * dx1[0] + dxv[1] * dx1[1];
        Real dx2nds = dxv[0] * dx2[0] + dxv[1] * dx2[1];
        //Real dx1nds = (x1[0] - geom.xref[0]) * dxv[1] + (x1[1] - geom.xref[1]) * dxv[0];
        //Real dx2nds = (x2[0] - geom.xref[0]) * dxv[1] - (x2[1] - geom.xref[1]) * dxv[0];

        Real dx = -dxv[0] * cos_alpha - dxv[1] * sin_alpha;
        Real dz =  dxv[1] * cos_alpha - dxv[0] * sin_alpha;

        Real cp1 = post.cp[ip], cp2 = post.cp[i];
        Real cpbar = 0.5 * (cp1 + cp2);

        post.cl += dx * cpbar;
        //post.cl_ue[ip] += dx * 0.5 * cp_ue[ip];
        //post.cl_ue[i]  += dx * 0.5 * cp_ue[i];

        //cl_alpha += cpbar * (sind(alpha) * dxv[0] - cosd(alpha) * dxv[1]) * deg2rad;

        post.cm += cp1 * dx1nds / 3.0 + cp1 * dx2nds / 6.0 + cp2 * dx1nds / 6.0 + cp2 * dx2nds / 3.0;
        //cdpi += dz * cpbar;
    }

    // Normalize by chord
    post.cl /= geom.chord;
    post.cm /= (geom.chord*geom.chord);
    
    int iw = Ncoords+Nwake-1;  // end of wake
    const Real* U = &glob.U[colMajorIndex(0,iw,4)];

    Real H, H_U[4];
    H = get_H(U[0],U[1],H_U);

    Real uk, uk_ue;
    uk = get_uk(U[3],par,uk_ue);

    post.cd = 2.0 * U[0] * pow(uk / op.Vinf, (5.0 + H) / 2.0);
    
}
