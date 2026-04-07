#pragma once

#include "real_type.hpp"
#include "data_structs.hpp"
#include <iostream>
#include <cmath>

#include "WPSmodels.hpp"
#include "newAmiet.hpp"

#include <fstream>
#include <sstream>

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;


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

template<typename Real>
void calc_WPS(const std::string& model, const Real theta,const Real deltaStar,const Real delta,const Real tauW,
                    const Real tauMax,const Real edgeVel,const Real dpdx, const Real (&omega)[Nsound], Real nu,
                    const Real Uinf,
                    const Real X,const Real Y,const Real Z,const  Real S,const Real rho,const int isSuction,
                    Real (&WPS)[Nsound]){
 
    
    Real useTauW = tauW;
    if (tauW > tauMax){
        useTauW = tauMax;
    }

    
    if (model == "roz"){
            calc_WPS_Rozenburg(theta,deltaStar,delta,useTauW,tauMax,edgeVel,dpdx,omega,rho,nu,WPS);
        }
    else if (model == "goo")
    {
        calc_WPS_Goody(theta,deltaStar,delta,useTauW,tauMax,edgeVel,dpdx,omega,rho,nu,Uinf,WPS);
    }
    else if (model == "kam")
    {
        calc_WPS_Kamruzzaman(theta,deltaStar,useTauW,edgeVel,dpdx,omega,rho,nu,WPS);
    }
    else if (model == "tno")
    {
        calc_WPS_TNO(delta,useTauW,edgeVel,omega,rho,nu,isSuction,WPS);
    }

}

template<typename Real>
Real calc_OASPL_AD(const Real* botStates, const Real* topStates, const Real chordScale, const Real Uinf,
    const Real X,const Real Y,const Real Z, const Real S, const Real nu, const Real rho,
    const std::string& model){

    const Real f_min = 200.0;
    const Real df    = 100.0;

    Real omega[Nsound];
    Real Freq[Nsound];

    for (int i = 0; i < Nsound; ++i) {
        Freq[i]  = f_min + i * df;
        omega[i] = Freq[i] * 2.0 * M_PI;
    }

    //Real SppUpper[Nsound]={0}, SppLower[Nsound]={0};
    Real WPSUpper[Nsound]={0},WPSLower[Nsound]={0};
    
    Real theta = topStates[0];
    Real deltaS = topStates[1];
    Real tauMax = topStates[2];
    Real edgeVel = topStates[3];
    Real dpdx = topStates[4];
    Real tauWall = topStates[5];
    Real delta = topStates[6];    

    if (tauWall < 0.0) {
        tauWall *= -1.0;
    }
    
    if (tauMax > 0.0){ 
        
        calc_WPS(model,theta,deltaS,delta,tauWall,tauMax,edgeVel,dpdx,
                    omega,nu,Uinf,X,Y,Z,S,rho,1,WPSUpper);
    }
    
    
    theta = botStates[0];
    deltaS = botStates[1];
    tauMax = botStates[2];
    edgeVel = botStates[3];
    dpdx = botStates[4];
    tauWall = botStates[5];
    delta = botStates[6];  

    if (tauWall < 0.0) {
        tauWall *= -1.0;
    }

    if (tauMax > 0.0){ 
        calc_WPS(model,theta,deltaS,delta,tauWall,tauMax,edgeVel,dpdx,
                    omega,nu,Uinf,X,Y,Z,S,rho,0,WPSLower);
    }

    Real farfieldSpectra[Nsound] ;

    Real c = Uinf/340.0;
    TE_noise_outer<Real>(c,Uinf,X,Y,Z,chordScale/2,0.0,chordScale,S,340.0,rho,nu,
                omega,WPSLower,WPSUpper,farfieldSpectra);

    // integrate S_pp over frequency:
    Real integral = 0.0;
    for (int i = 0; i < Nsound - 1; ++i) {
        Real df = Freq[i+1] - Freq[i];
        integral += 0.5 * (farfieldSpectra[i] + farfieldSpectra[i+1]) * 2*M_PI* df;
    }

    // Now convert to dB re 20 µPa:
    Real pref2 = (20e-6)*(20e-6); // reference pressure squared
    Real OASPL = 10.0 * std::log10(integral / pref2);

    return OASPL;
}