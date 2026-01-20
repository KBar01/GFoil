#pragma once

#include <cmath>
#include "real_type.hpp"

template<typename Real>
void calc_WPS_Goody(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real edgeVel,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real Uinf,
                    Real (&phiqq)[Nsound]){
    
    Real a = 3;
    Real b = 2;
    Real c = 0.75;
    Real d = 0.5;
    Real e = 3.7;
    Real f = 1.1;
    Real g = -0.57;
    Real h = 7;
    Real i = 1;
    Real Ue = edgeVel;
    Real u_t = std::sqrt(tauWall/rho);
    Real Rt= (delta/Ue)/(nu/(u_t*u_t));
    Real SS   = Ue / (tauWall*tauWall*delta);
    Real FS   = delta/Ue ;

    for (int n=0;n<Nsound;++n){
        Real omegaBar= omega[n]*FS ;
        phiqq[n] = ((a*std::pow(omegaBar,b))/(std::pow(i*std::pow(omegaBar, c) + d, e) + std::pow((f*std::pow(Rt, g)*omegaBar), h))) / SS;
    }

}

template<typename Real>
void calc_WPS_Kamruzzaman(Real theta,
                    Real deltaS,
                    Real tauWall,
                    Real Ue,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real (&phiqq)[Nsound])
    {


    /*
    A semi-empirical surface pressure
    spectrum model for airfoil trailing-edge
    noise prediction

    M. Kamruzzaman, 2017 
    */

    Real H = deltaS/theta ;
    Real beta_c = std::max((theta/tauWall)*(dpdx),-0.5);
    Real Pi = 0.227;
    if (beta_c > -0.5){
        Pi = 0.8*std::pow(beta_c+0.5, 0.75);
    }
    
    // Table 1 in Lee : Comparison and Assessment of RecentEmpirical Models for Turbulent BoundaryLayer Wall Pressure Spectrum
    Real m = 0.5*std::pow(H/1.31, 0.3);
    Real a = 0.45*(1.75*std::pow(Pi*Pi*beta_c*beta_c, m) + 15);
    Real u_t = std::sqrt(tauWall/rho);
    Real Rt = (deltaS/Ue) / (nu/(u_t*u_t));
    Real B3 = std::pow(1.15*Rt, -2.0/7.0);
   
    Real SS   = Ue / (tauWall*tauWall*deltaS);
    Real FS   = deltaS/Ue ;

    for (int n=0;n<Nsound;++n){
        Real omegaBar= omega[n]*FS ;
        phiqq[n] = ((a*std::pow(omegaBar,2.0)) / (std::pow(std::pow(omegaBar, 1.637) + 0.27, 2.47) + std::pow((B3*omegaBar), 7.0)))/SS;
    }

}

template<typename Real>
void calc_WPS_Rozenburg(Real theta,
                    Real deltaS,
                    Real delta,
                    Real tauWall,
                    Real tauMax,
                    Real Ue,
                    Real dpdx,
                    const Real (&omega)[Nsound],
                    Real rho,
                    Real nu,
                    Real (&phiqq)[Nsound]){

    
    
    /*
    Wall-Pressure Spectral Model Including the Adverse
    Pressure Gradient Effects
    Yannick Rozenberg, 2012
    */

    Real Delta = delta/deltaS  ;
    
    Real beta_c = std::max((theta/tauWall)*(dpdx),-0.5);
    Real Pi = 0.227;
    if (beta_c > -0.5){
        Pi = 0.8*std::pow(beta_c+0.5, 0.75);
    }
    Real u_t = std::sqrt(tauWall/rho);
    Real Rt = (delta/Ue)/(nu/(u_t*u_t));


    Real b = 2; // Done
    Real c = 0.75; // Done
    Real A1 = 3.7 + 1.5*beta_c ; //Done - A1
    Real F1 = 4.76*std::pow((1.4/Delta), 0.75) * (0.375*A1 -1) ; // F1
    
    Real a = (2.82*Delta*Delta*std::pow((6.13*std::pow(Delta,-0.75) + F1), A1))  *  (4.2*(Pi/Delta) + 1); //Done
    Real f = 8.8; //done
    Real g = -0.57; //done
    Real F2 = std::min(3.0,19.0/ std::sqrt(Rt)); // done
    Real i = 4.76; //done

    Real SS   = Ue / (tauMax*tauMax*deltaS);
    Real FS   = deltaS/Ue ;

    Real C3prime = 8.8*std::pow(Rt, -0.57);

    for (int n=0;n<Nsound;++n){
        Real omegaBar= omega[n]*FS ;
        Real top = (a*std::pow(omegaBar, 2));
        Real bot = std::pow(4.76*std::pow(omegaBar, 0.75) + F1, A1)   +   std::pow((C3prime*omegaBar), F2);
        phiqq[n] = ( top/bot )/SS;
    }

}


/////////////////////////////////// All TNO Funcs /////////////////////////////////////
template<typename Real>
void mean_velocity_profile(const Real (&y)[NblPoints],
                           Real delta,
                           Real u_t,
                           Real nu,          // fluid kinematic viscosity
                           Real Ue,
                           Real (&U)[NblPoints],
                           Real (&dUdy)[NblPoints])
{

    // Following Eq 19

    const Real kappa = 0.41;
    const Real B = 5.5;
    for (int i = 0; i < NblPoints; ++i) {
        
        Real y_plus = y[i] * u_t / nu;

        if (y_plus <= 5.0){
            
            //u_plus = y_plus in viscous sub layer 
            U[i]  = u_t * y_plus ;
            // du+/dy = u_t/nu  so  dU/dy = u_t * (u_t/nu) = u_t^2 / nu
            dUdy[i] = (u_t * u_t) / nu;
        }
        else{

        Real W = 1-std::cos(M_PI*y[i] / delta); // Eq 20
        
        Real u_plus = (1.0 / kappa) * std::log(y_plus) + B +
                        0.5*W*((Ue/u_t) - (1.0/kappa)*std::log((u_t*delta)/nu) - B);
        
        // streamwise velocity
        U[i] = u_plus * u_t;
        // Derivative dU/dy
        // du+/dy:
        Real duplus_dy = (1.0 / (kappa * y[i])) + 0.5*(M_PI/delta)*std::sin(M_PI*y[i]/delta)*
                            ((Ue/u_t) - (1/(kappa))*std::log((u_t*delta)/nu) - B );

        dUdy[i] = u_t * duplus_dy;
        }
    }
}

template<typename Real>
void Turb_shear_stress(
    

    // calculate turbulent shear stress term u2^2 bar
    const Real (&dUdy)[NblPoints],
    const Real (&l_mix)[NblPoints],
    const int isSuction,
    Real (&u22)[NblPoints])
{
   
    // Using Eq 25-27
    for (int i = 0; i < NblPoints; ++i)
    {
       Real nu_t = l_mix[i]*l_mix[i]*std::sqrt(dUdy[i]*dUdy[i]);
       Real kt = std::sqrt((nu_t*nu_t*dUdy[i]*dUdy[i]) / 0.09);
        
       if (isSuction){
        u22[i] = 0.45*kt;
       }
       else{
        u22[i] = 0.3*kt;
       }
    }
}

template<typename Real>
void Integral_Length_scale(

    // calc vertical integral length scale L_2
    const Real delta,
    const Real (&y)[NblPoints],
    Real (&L2)[NblPoints],
    Real (&l_mix)[NblPoints])
{

    // Using Eq 22-24
    const Real k = 0.41;
    const Real B = 5.5 ;

    for (int i = 0; i < NblPoints; ++i)
    {
        l_mix[i] = (0.085 * delta * std::tanh( (k*y[i]) / (0.085*delta))) /
                        std::sqrt( std::pow(1+B*(y[i]/delta), 6.0) );
        
        L2[i] = l_mix[i] / 0.41 ;
    }
}

template<typename Real>
void Energy_density_spectrum(
    
    // calculating phi22 for midspan observer in far-field

    const Real k1,
    const Real (&L2)[NblPoints],
    Real (&phi22)[NblPoints])
{

    // Using Eq 28-32
    Real beta1 = 1.0;
    Real beta3 = 0.75;

    for (int i=0;i<NblPoints;++i){

        Real ke = 0.7468 / (2*L2[i]);
        Real term = ((beta1*k1)/ke) * ((beta1*k1)/ke);
        phi22[i] = (4.0/(9.0*M_PI)) * ((beta1*beta3)/(ke*ke)) * (term / std::pow(1+term, (7.0/3.0))) ;
    }
}

template<typename Real>
void calc_WPS_TNO(
    const Real delta,
    Real tauWall,
    const Real edgeVel,
    const Real (&omega)[Nsound],
    const Real rho,
    const Real nu,
    const int isSuction,
    Real (&phiqq)[Nsound])
{
    
    
    /*
    Source Characterization of Turbulent
    Boundary Layer Trailing Edge Noise Using an
    Improved TNO Model

    Seongkyu Lee, 2016
    */
    
    
    
    // Compute shear velocity and min y from target y+
    Real u_t = std::sqrt(tauWall / rho);

    const Real yplus_target = 1.0;  
    Real y_min = yplus_target * nu / u_t;

    
    //Build wall-normal grid with either cosine or linear spacing
    Real y[NblPoints];
    const Real y_max = delta;

    // flag: true = cosine stretching, false = linear spacing
    bool useCosineStretch = false;

    for (int i = 0; i < NblPoints; ++i)
    {
        Real eta = static_cast<Real>(i) / static_cast<Real>(NblPoints - 1);

        if (useCosineStretch)
        {
            // Cosine stretching
            Real y_stretch = 0.5 * (1.0 - std::cos(M_PI * eta));
            y[i] = y_min + (y_max - y_min) * y_stretch;
        }
        else
        {
            // Linear spacing
            Real y_stretch = eta;  
            y[i] = y_min + (y_max - y_min) * y_stretch;
        }
    }

    Real Uc = 0.65 * edgeVel;

    Real U[NblPoints], dUdy[NblPoints];
    mean_velocity_profile(y, delta, u_t, nu, edgeVel, U, dUdy);

    Real L2[NblPoints], l_mix[NblPoints];
    Integral_Length_scale(delta, y, L2, l_mix);

    Real u22[NblPoints];
    Turb_shear_stress(dUdy, l_mix, isSuction, u22);

    Real phi22[NblPoints];

    Real F[NblPoints];
    for (int i = 0; i < NblPoints; ++i){
        F[i] =  L2[i] * (1/Uc) * (dUdy[i] * dUdy[i]) * u22[i];
    }

    for (int w = 0; w < Nsound; ++w)
    {
        
        // Using Eq 37
        Real k1 = omega[w] / Uc;
        Real k = std::abs(k1);

        Energy_density_spectrum(k1, L2, phi22);

        Real integrand[NblPoints];
        for (int i = 0; i < NblPoints; ++i)
        {
            integrand[i] = F[i]*phi22[i]*std::exp(-2.0 * y[i] * k);
        }

        // trapezoidal integration
        Real integral = 0.0;
        for (int i = 1; i < NblPoints; ++i)
        {
            Real dy_local = y[i] - y[i - 1];
            integral += 0.5 * (integrand[i] + integrand[i - 1]) * dy_local;
        }
        

        // pretty sure this code gives you the wavenumber frequency PSD, 
        // so convert using same method in R&M with span corr length: 

        Real ly = 1.4*Uc / omega[w] ;

        phiqq[w] = (4.0 * rho*rho * integral) * M_PI * (1.0/ly) *  2.0;
    }
}