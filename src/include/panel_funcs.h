#ifndef PANEL_FUNCTIONS_H
#define PANEL_FUNCTIONS_H
#include "real_type.h"

// Structure to hold panel information
struct PanelInfo {
    Real t[2];
    Real n[2];
    Real x, z, d, r1, r2, theta1, theta2;
};

struct Isol;
struct Foil;
struct Param;
struct Post;
struct Oper;
struct Geom;



void panel_info(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo& info);

Real norm2(const Real* x);
Real norm2_3D(const Real* x);

void panel_linvortex_stream(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo& panelInfo, Real& a, Real& b);

void panel_constsource_stream(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo& panelInfo, Real& a);

void panel_linsource_stream(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo& panelInfo, Real& a,Real& b);

void panel_linvortex_velocity(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo& info, Real& a1, Real& b1, Real& a2, Real& b2);


void panel_constsource_velocity(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo& info,
    Real& a1,Real& a2);

void panel_linsource_velocity(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo& info, Real& a1, Real& b1, Real& a2, Real& b2);

void inviscid_velocity(const Foil& foil,const Real* gammas, const Real& Vinf,
    const Real& alpha,const Real& CPx,const Real&CPy, Real*velocity);

void dvelocity_dgamma(const Foil& foil,const Real& CPx,const Real&CPy, Real*V_G);


#endif