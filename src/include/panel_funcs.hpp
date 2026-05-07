#pragma once

#include <cmath>

// Structure to hold panel information (templatised on Real)
template<typename Real>
struct PanelInfo {
    Real t[2];
    Real n[2];
    Real x, z, d, r1, r2, theta1, theta2;
};

// norm2 / norm2_3D are defined in real_type.h and real_type.hpp respectively.

template<typename Real>
void panel_info(const Real& panelStartX, const Real& panelStartY,
                const Real& panelEndX,   const Real& panelEndY,
                const Real& controlPointX, const Real& controlPointY, PanelInfo<Real>& info) {
    
    /* Finds commen panel quantities, such as tangent/normal vector, and
        distances/angles to some contorl point

        The struct PanelInfo<Real> is passed as reference (passed in input) to be re-used
        within loops for efficiency in terms of memory (i think).
    */

    // Panel-aligned tangent and normal vectors
    Real t_init[2] = {panelEndX - panelStartX, panelEndY - panelStartY};
    Real norm_t_init = norm2(t_init);
    
    info.t[0] = t_init[0] / norm_t_init;
    info.t[1] = t_init[1] / norm_t_init;
    
    info.n[0] = -info.t[1];
    info.n[1] = info.t[0];

    // Control point relative to panel start node
    Real xz[2] = {controlPointX - panelStartX, controlPointY - panelStartY};
    info.x = xz[0]*info.t[0] + xz[1]*info.t[1];  // Dot product in panel-aligned coord system
    info.z = xz[0]*info.n[0] + xz[1]*info.n[1];  // Dot product in panel-aligned coord system

    // Panel length
    info.d = std::sqrt(std::pow(t_init[0],2) + std::pow(t_init[1],2)) ;

    // Distances and angles
    info.r1 = std::sqrt(info.x*info.x + info.z*info.z);
    info.r2 = std::sqrt((info.x - info.d) * (info.x - info.d) + info.z*info.z);

    info.theta1 = std::atan2(info.z, info.x);
    info.theta2 = std::atan2(info.z, info.x - info.d);

}


template<typename Real>
void panel_linvortex_stream(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo<Real>& panelInfo, Real& a, Real& b) {
    
    // Get panel information by calling panel_info
    panel_info(panelStartX,panelStartY,panelEndX,panelEndY,controlPointX,controlPointY,panelInfo);
    
    // Check for r1, r2 zero
    const Real ep = 1e-9;
    Real logr1 = (panelInfo.r1 < ep) ? Real(0.0) : std::log(panelInfo.r1);
    Real logr2 = (panelInfo.r2 < ep) ? Real(0.0) : std::log(panelInfo.r2);

    // Streamfunction components
    Real P1 = (0.5/M_PI)*(panelInfo.z* (panelInfo.theta2 - panelInfo.theta1) - panelInfo.d + panelInfo.x*logr1 - (panelInfo.x - panelInfo.d)*logr2);
    Real P2 = panelInfo.x*P1 + (0.5/M_PI)*(0.5*panelInfo.r2*panelInfo.r2*logr2 - 0.5*panelInfo.r1*panelInfo.r1*logr1 - (panelInfo.r2*panelInfo.r2)/4.0 + (panelInfo.r1*panelInfo.r1)/4.0);

    // Influence coefficients
    a = P1 - P2 / panelInfo.d;
    b = P2 / panelInfo.d;
}

template<typename Real>
void panel_constsource_stream(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo<Real>& panelInfo, Real& a) {
    
    
    // Get panel information by calling panel_info
    panel_info(panelStartX,panelStartY,panelEndX,panelEndY,controlPointX,controlPointY,panelInfo);

    // Streamfunction
    Real ep = 1e-9;
    Real logr1, logr2;

    if (panelInfo.r1 < ep) {
        logr1 = 0.0;
        panelInfo.theta1 = M_PI;
        panelInfo.theta2 = M_PI;
    } 
    else {
        logr1 = std::log(panelInfo.r1);
    }

    if (panelInfo.r2 < ep) {
    logr2 = 0.0;
    panelInfo.theta1 = 0.0;
    panelInfo.theta2 = 0.0;
    } else {
    logr2 = std::log(panelInfo.r2);
    }

    Real P = (panelInfo.x * (panelInfo.theta1 - panelInfo.theta2) + panelInfo.d * panelInfo.theta2 + panelInfo.z * logr1 - panelInfo.z * logr2) / (2 * M_PI);
    if ((panelInfo.theta1 + panelInfo.theta2) > M_PI) {
        P = P - 0.25*panelInfo.d;
    } 
    else {
        P = P + 0.75*panelInfo.d;
    }

    // Influence coefficient
    a = P;
}

template<typename Real>
void panel_linsource_stream(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo<Real>& panelInfo, Real& a,Real& b) {
    
    
    // Get panel information by calling panel_info
    panel_info(panelStartX,panelStartY,panelEndX,panelEndY,controlPointX,controlPointY,panelInfo);

    // Streamfunction
    Real ep = 1e-9;
    Real logr1, logr2;


    if (panelInfo.theta1<0){ panelInfo.theta1 += 2*M_PI;}
    if (panelInfo.theta2<0){ panelInfo.theta2 += 2*M_PI;}

    if (panelInfo.r1 < ep) {
        logr1 = 0.0;
        panelInfo.theta1 = M_PI;
        panelInfo.theta2 = M_PI;
    } 
    else {
        logr1 = std::log(panelInfo.r1);
    }

    if (panelInfo.r2 < ep) {
    logr2 = 0.0;
    panelInfo.theta1 = 0.0;
    panelInfo.theta2 = 0.0;
    } else {
    logr2 = std::log(panelInfo.r2);
    }

   
    Real P1 = (0.5*M_1_PI)*(panelInfo.x*(panelInfo.theta1-panelInfo.theta2)+panelInfo.theta2*panelInfo.d + panelInfo.z*logr1 - panelInfo.z*logr2);
    Real P2 = panelInfo.x*P1 + (0.5/M_PI)*(0.5*panelInfo.r2*panelInfo.r2*panelInfo.theta2 - 0.5*panelInfo.r1*panelInfo.r1*panelInfo.theta1 - 0.5*panelInfo.z*panelInfo.d);

    a = P1 - P2/panelInfo.d;
    b = P2/panelInfo.d;

}



template<typename Real>
void panel_linvortex_velocity(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo<Real>& info, Real& a1, Real& b1, Real& a2, Real& b2){

    panel_info(panelStartX,panelStartY,panelEndX,panelEndY,controlPointX,controlPointY,info);
    
    // Compute velocity components
    Real temp1 = (info.theta2 - info.theta1) / (2.0*M_PI);
    Real temp2 = (2.0 * info.z * std::log(info.r1 / info.r2) - 2.0 * info.x * (info.theta2 - info.theta1)) / (4.0 * M_PI * info.d);
    Real ug1 = temp1 + temp2;
    Real ug2 = -temp2;

    temp1 = std::log(info.r2 / info.r1) / (2.0 * M_PI);
    temp2 = (info.x * std::log(info.r1 / info.r2) - info.d + info.z * (info.theta2 - info.theta1)) / (2.0 * M_PI * info.d);
    Real wg1 = temp1 + temp2;
    Real wg2 = -temp2;

    // Compute influence coefficients in original coordinate system
    a1 = ug1 * info.t[0] + wg1 * info.n[0];
    a2 = ug1 * info.t[1] + wg1 * info.n[1];
    b1 = ug2 * info.t[0] + wg2 * info.n[0];
    b2 = ug2 * info.t[1] + wg2 * info.n[1];

}


template<typename Real>
void panel_constsource_velocity(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo<Real>& info,
    Real& a1,Real& a2) {

    panel_info(panelStartX, panelStartY, panelEndX, panelEndY, controlPointX, controlPointY, info);


    // Handle small radius cases to avoid log(0)
    Real ep = 1e-9;

    Real logr1 = (info.r1 < ep) ? Real(0.0) : std::log(info.r1);
    Real logr2 = (info.r2 < ep) ? Real(0.0) : std::log(info.r2);

    if (info.r1 < ep) { info.theta1 = M_PI; info.theta2 = M_PI; }
    if (info.r2 < ep) { info.theta1 = 0; info.theta2 = 0; }

    // Velocity in panel-aligned coordinate system
    Real u = (0.5/M_PI)*(logr1 - logr2);
    Real w = (0.5/M_PI)*(info.theta2 - info.theta1);

    // Transform velocity to original coordinates
    a1 = u * info.t[0] + w * info.n[0];
    a2 = u * info.t[1] + w * info.n[1];

};

template<typename Real>
void panel_linsource_velocity(const Real& panelStartX, const Real& panelStartY,
    const Real& panelEndX,   const Real& panelEndY,
    const Real& controlPointX, const Real& controlPointY, PanelInfo<Real>& info, Real& a1, Real& b1, Real& a2, Real& b2){

    panel_info(panelStartX,panelStartY,panelEndX,panelEndY,controlPointX,controlPointY,info);
    
    // Compute velocity components
    Real temp1 = std::log(info.r1/info.r2) / (2*M_PI);
    Real temp2 = ((info.x * std::log(info.r1 / info.r2)) - info.d + info.z*(info.theta2 - info.theta1)) / (2.0 * M_PI * info.d);
    Real ug1 = temp1 - temp2;
    Real ug2 = temp2;

    temp1 = (info.theta2-info.theta1) / (2.0 * M_PI);
    temp2 = (-info.z * std::log(info.r1 / info.r2) + info.x*(info.theta2 - info.theta1)) / (2.0 * M_PI * info.d);
    Real wg1 = temp1 - temp2;
    Real wg2 = temp2;

    // Compute influence coefficients in original coordinate system
    a1 = ug1 * info.t[0] + wg1 * info.n[0];
    a2 = ug1 * info.t[1] + wg1 * info.n[1];
    b1 = ug2 * info.t[0] + wg2 * info.n[0];
    b2 = ug2 * info.t[1] + wg2 * info.n[1];

}
