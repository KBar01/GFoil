#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.h"
#include "data_structs.h"
#include "main_func.h"
#include "get_funcs.h"
#include "panel_funcs.h"
#include "amiet.h"
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

#include "nlohmann/json.hpp"  // nlohmann/json

using json = nlohmann::json;

bool runCode(
    
    // restart logic for optimisations
    bool fromRestart,
    // geometry parameters
    const Real nCrit,
    const Real Ufac, 
    const Real TEfac,
    const Real chordScaling,
    const Real (&inXcoords)[Nin], 
    Real (&inYcoords)[Nin],
    // actual aero parameters
    Real alphad,
    Real Re, 
    Real Ma,
    Real rhoInf,
    Real kinViscInf,
    const Real &topTransPos,
    const Real &botTransPos,
    const bool force,
    // acoustic model parameters
    const std::string model,
    const Real sampleTE,
    const Real X,
    const Real Y,
    const Real Z,
    const Real S,
    // output parameters
    const int doCps

    ){
    Real alpha = (alphad/180)*M_PI;
    Oper oper(alpha,Re,Ma);
    oper.rho = rhoInf;

    Geom geom;
    
    Real flattenedCoords[2*Ncoords]={0};
    Real inCoords[2*Nin]={0};
    for (int i=0;i<Nin;++i){
        inCoords[colMajorIndex(0,i,2)] = inXcoords[i];
        inCoords[colMajorIndex(1,i,2)] = inYcoords[i];
    }
    make_panels(inCoords,flattenedCoords,Ufac,TEfac); // does spline to redist nodes over aerofoil for fixed number of 200 nodes
    
    // finding node positions to force transition ------------------------
    Real xTransTop = topTransPos * geom.chord ;
    Real xTransBot = botTransPos * geom.chord ;
    
    int idx_closest_bot = 0;
    int idx_closest_top = Ncoords - 1;  // top TE node
    
    if (force) {
        
        // Bottom surface: from bottom TE forward
        for (int i = 1; i < Ncoords; ++i) {
            Real x = flattenedCoords[colMajorIndex(0, i, 2)];
            Real dist = x - xTransBot;
            if (dist < 0.0) {
                idx_closest_bot = i;
                break;
            }
        }

        // Top surface: from top TE backward
        for (int i = Ncoords - 2; i >= 0; --i) {
            Real x = flattenedCoords[colMajorIndex(0, i, 2)];
            Real dist = x - xTransTop;
            if (dist < 0.0) {
                idx_closest_top = i;
                break;
            }
        }
    }

    Trans tdata;

    tdata.transNode[0] = idx_closest_bot;
    tdata.transNode[1] = idx_closest_top;
    tdata.transPos[0] = xTransBot ;
    tdata.transPos[1] = xTransTop ;

    Foil foil(flattenedCoords);
    Isol isol;
    Param param;
    param.ncrit = nCrit;
    Wake wake;
    Vsol vsol;
    Glob glob;


    build_gamma_codi(isol,foil,oper);
    init_thermo(oper,param,geom);
    build_wake(foil,geom,oper,isol,wake);
    stagpoint_find(isol,foil,wake);
    identify_surfaces(isol,vsol);
    set_wake_gap(foil,isol,vsol);
    calc_ue_m(foil,wake,isol,vsol);
    rebuild_ue_m(foil,wake,isol,vsol,false);

    if (fromRestart){

        std::ifstream prevfile("restart.json");

        json j;
        prevfile >> j;

        for (int i = 0; i < RVdimension; ++i) {
        double val = j["states"][i].get<double>();
        glob.U[i] = val;  // assigns numeric value to RealReverse
        }

        for (int i = 0; i < (Ncoords + Nwake); ++i) {
            bool val = j["turb"][i].get<bool>();
            vsol.turb[i] = val;
        }

        if (force){
            tdata.isForced[0]  = 1;
            tdata.isForced[1]  = 1;
        }
    }
    else {
        init_boundary_layer(oper,foil,param,isol,vsol,glob,tdata,force);
    }
    
    stagpoint_move(isol,glob,foil,wake,vsol);
    bool converged = solve_coupled(oper,foil,wake,param,vsol,isol,glob,tdata,force);
    Post post;
    calc_force(oper,geom,param,isol,foil,glob,post);
    
    ////////////////////////////// Aero done, now acoustics ///////////////////////////////////////////////////////////////////

    Real Uinf = (Re*kinViscInf)/(chordScaling) ;
    Real tauWall[Ncoords];
    if (doCps){
        Real cf_U[4]={0};
        for (int i=0;i<Ncoords;++i){
            tauWall[i] = get_cf(
                glob.U[colMajorIndex(0,i,4)],
                glob.U[colMajorIndex(1,i,4)],
                glob.U[colMajorIndex(2,i,4)],
                glob.U[colMajorIndex(3,i,4)],
                vsol.turb[i],
                false,
                param,
                cf_U
            );
            tauWall[i] *=  (oper.rho*(glob.U[colMajorIndex(3,i,4)]*Uinf * glob.U[colMajorIndex(3,i,4)]*Uinf))/2;
        }
    }
    
    Real topsurf[7],botsurf[7];
    Real xcoords[Ncoords]={0};
    Real ycoords[Ncoords]={0};

    for (int i=0;i<Ncoords;++i){
        xcoords[i] = flattenedCoords[colMajorIndex(0,i,2)];
        ycoords[i] = flattenedCoords[colMajorIndex(1,i,2)];
    }

    interpolate_at_95_both_surfaces(xcoords,glob.U,post.cp,oper,vsol,param,topsurf,botsurf,Uinf,sampleTE,chordScaling);
    Real OASPL = calc_OASPL(botsurf,topsurf,chordScaling,Uinf,X,Y,Z,S,kinViscInf,rhoInf,doCps,model);

    // check OASPL validity
    if (std::isnan(OASPL) || std::isinf(OASPL)) {
        converged = false;
    }
    
    if (converged){

        json out;
        
        out["conv"] = 1;
        out["aerofoilChord"] = chordScaling.getValue();
        out["freestreamVelocity"] = Uinf.getValue();
        out["samplingLoc"] = sampleTE.getValue();
        out["CL"]  = post.cl.getValue();
        out["CD"]  = post.cd.getValue();
        out["CM"]  = post.cm.getValue();
        out["OASPL"] = OASPL.getValue();

        if (doCps){
            
            json restart;
            //std::vector<double> states_d(RVdimension);
            // Extract numeric values from CoDiPack types
            //for (size_t i = 0; i < RVdimension; ++i){
            //    states_d[i] = glob.U[i].getValue();
            //} 
            //restart["states"] = states_d;
            //restart["turb"]   = vsol.turb;
            //std::ofstream restartFile("restart.json");
            //restartFile << restart.dump(4);  // pretty print with 4 spaces indentation
            //restartFile.close();

            //  calc transition point
            Real botTransX = geom.chord;
            for(int i=0;i<isol.stagIndex[0];++i){
                int isTurb = vsol.turb[isol.stagIndex[0] - i];
                if (isTurb){
                    botTransX = foil.x[colMajorIndex(0,isol.stagIndex[0]-i,2)];
                    break;
                } 
            }
            Real topTransX = geom.chord;
            for(int i=0;i<200-isol.stagIndex[1];++i){
                int isTurb = vsol.turb[isol.stagIndex[1] + i];
                if (isTurb){
                    topTransX = foil.x[colMajorIndex(0,isol.stagIndex[1]+i,2)];
                    break;
                } 
            }

            double inner[2*Ncoords] ;
            double cps[Ncoords];
            double tauWallOut[Ncoords];
            for (int i=0;i<(2*Ncoords);++i){
                inner[i] = foil.x[i].getValue() ;
            }

            for (int i=0;i<(Ncoords);++i){
                cps[i] = post.cp[i].getValue();
                tauWallOut[i] = tauWall[i].getValue();
            }
            out["innerFoil"] = inner;
            out["Cp"] = cps;
            out["tauWall"] = tauWallOut;
            
            out["stagnation"] = isol.stagIndex;
            out["topTransX"]  = topTransX.getValue();
            out["botTransX"]  = botTransX.getValue();
            
            
            //out["tauWall"] = tauWall;
            std::vector<std::string> BLoutputNames = {"CL", "CD",
                "thetaUpper", "deltaStarUpper", "tauMaxUpper","edgeVelocityUpper", "dpdxUpper", "tauWallUpper", "delta99Upper",
                "thetaLower", "deltaStarLower", "tauMaxLower","edgeVelocityLower", "dpdxLower", "tauWallLower", "delta99Lower"
            };
            out[BLoutputNames[2]] = topsurf[0].getValue();
            out[BLoutputNames[3]] = topsurf[1].getValue();
            out[BLoutputNames[4]] = topsurf[2].getValue();
            out[BLoutputNames[5]] = topsurf[3].getValue();
            out[BLoutputNames[6]] = topsurf[4].getValue();
            out[BLoutputNames[7]] = topsurf[5].getValue();
            out[BLoutputNames[8]] = topsurf[6].getValue();

            out[BLoutputNames[9]] = botsurf[0].getValue();
            out[BLoutputNames[10]] = botsurf[1].getValue();
            out[BLoutputNames[11]] = botsurf[2].getValue();
            out[BLoutputNames[12]] = botsurf[3].getValue();
            out[BLoutputNames[13]] = botsurf[4].getValue();
            out[BLoutputNames[14]] = botsurf[5].getValue();
            out[BLoutputNames[15]] = botsurf[6].getValue();
        }
        
        
        std::ofstream outFile("out.json");
        outFile << out.dump(4);  // pretty print with 4 spaces indentation
        outFile.close();
        
    }
    else{

        json out;
        out["conv"] = 0;
        std::ofstream outFile("out.json");
        outFile << out.dump(4);  // pretty print with 4 spaces indentation
        outFile.close();

    }

    return converged;
};


int main(){

    // Open JSON file
    std::ifstream infile("input.json");
    if (!infile) {
        std::cerr << "Failed to open input.json\n";
        return 1;
    }

    // Parse the JSON
    json j;
    infile >> j;

    const int doWPSonly = j["WPSonly"].get<int>();
    const int doCustWPS = j["custSpectra"].get<int>();

    if (doWPSonly){

        Real rhoInf = j["rho"].get<double>();
        Real kinViscInf = j["nu"].get<double>();
        Real chordScaling = j["chord"].get<double>();
        Real Re = j["Re"].get<double>();
        const Real X = j["X"].get<double>();
        const Real Y = j["Y"].get<double>();
        const Real Z = j["Z"].get<double>();
        const Real S = j["S"].get<double>();
        const Real toptheta   = j["toptheta"].get<double>();
        const Real topdstar   = j["topdstar"].get<double>();
        const Real topdelta   = j["topdelta"].get<double>();
        const Real toptauw    = j["toptauw"].get<double>();
        const Real toptaumax  = j["toptaumax"].get<double>();
        const Real topue      = j["topue"].get<double>();
        const Real topdpdx    = j["topdpdx"].get<double>();

        const Real bottheta   = j["bottheta"].get<double>();
        const Real botdstar   = j["botdstar"].get<double>();
        const Real botdelta   = j["botdelta"].get<double>();
        const Real bottauw    = j["bottauw"].get<double>();
        const Real bottaumax  = j["bottaumax"].get<double>();
        const Real botue      = j["botue"].get<double>();
        const Real botdpdx    = j["botdpdx"].get<double>();
        const std::string model = j["model"].get<std::string>();
        // HERE: optional assignment of BL properties
        Real topsurf[7],botsurf[7];
        topsurf[0] = toptheta;
        botsurf[0] = bottheta;
        topsurf[1] = topdstar;
        botsurf[1] = botdstar;
        topsurf[2] = toptaumax;
        botsurf[2] = bottaumax;
        topsurf[3] = topue;
        botsurf[3] = botue;
        topsurf[4] = topdpdx;
        botsurf[4] = botdpdx;
        topsurf[5] = toptauw;
        botsurf[5] = bottauw;
        topsurf[6] = topdelta;
        botsurf[6] = botdelta;

        Real Uinf = (Re*kinViscInf)/(chordScaling) ;
        Real OASPL = calc_OASPL(botsurf,topsurf,chordScaling,Uinf,X,Y,Z,S,kinViscInf,rhoInf,1,model);
        return 1 ;
    }

    if (doCustWPS){

        Real rhoInf = j["rho"].get<double>();
        Real kinViscInf = j["nu"].get<double>();
        Real chordScaling = j["chord"].get<double>();
        Real Re = j["Re"].get<double>();
        const Real X = j["X"].get<double>();
        const Real Y = j["Y"].get<double>();
        const Real Z = j["Z"].get<double>();
        const Real S = j["S"].get<double>();
        Real Uinf = (Re*kinViscInf)/(chordScaling) ;

        Real WPSUpper[Nsound];
        Real WPSLower[Nsound];
        Real Freq[Nsound];
        Real omega[Nsound];

        for (int i=0;i<Nsound;++i){
            WPSLower[i] = j["custWPSLower"][i].get<double>();
            WPSUpper[i] = j["custWPSUpper"][i].get<double>();
            Freq[i]    = j["custFreq"][i].get<double>();
            omega[i] = Freq[i] * 2.0 * M_PI;
        }
        
        Real farfieldSpectra[Nsound] ;
        Real FF_up[Nsound] ;
        Real FF_down[Nsound] ;
        TE_noise_outer((Uinf/340),Uinf,X,Y,Z,chordScaling/2,0.0,chordScaling,S,340,rhoInf,kinViscInf,
                omega,WPSLower,WPSUpper,farfieldSpectra,FF_up,FF_down);

        // integrate S_pp over frequency:
        Real integral = 0.0;
        for (int i = 0; i < Nsound - 1; ++i) {
            Real df = Freq[i+1] - Freq[i];
            integral += 0.5 * (farfieldSpectra[i] + farfieldSpectra[i+1]) * 8*M_PI* df;
        }

        // Now convert to dB re 20 µPa:
        Real pref2 = (20e-6)*(20e-6); // reference pressure squared
        Real OASPL = 10.0 * std::log10(integral / pref2);
        json j;
        double prefSqrd = pref2.getValue();
        std::vector<double> freq_d(Nsound);
        std::vector<double> wpsupper_d(Nsound);
        std::vector<double> wpslower_d(Nsound);
        std::vector<double> spectra_d(Nsound);

        std::vector<double> spectra_up(Nsound);
        std::vector<double> spectra_down(Nsound);
        
        for (int i=0; i<Nsound; ++i){
            freq_d[i]     = Freq[i].getValue();
            wpsupper_d[i] = (WPSUpper[i].getValue())/prefSqrd;
            wpslower_d[i] = (WPSLower[i].getValue())/prefSqrd;
            spectra_d[i]  = (farfieldSpectra[i].getValue())/prefSqrd;
            spectra_up[i]  = (FF_up[i].getValue())/prefSqrd;
            spectra_down[i]  = (FF_down[i].getValue())/prefSqrd;
        }
        j["frequency_Hz"]        = freq_d;
        j["WPS_upper/prefSqrd"]  = wpsupper_d;
        j["WPS_lower/prefSqrd"]  = wpslower_d;
        j["FF_spectra/prefSqrd"] = spectra_d;
        j["FF_spectra_up/prefSqrd"] = spectra_up;
        j["FF_spectra_down/prefSqrd"] = spectra_down;
        j["OASPL_dB"]     = OASPL.getValue();

        std::ofstream file("cutomSpectra.json");
        file << j.dump(4);     // pretty print, 4 spaces
        file.close();

        return 1 ;
    }
    else{

    Real inXcoords[Nin]={0}, inYcoords[Nin]={0};

    for (int i = 0; i < Nin; ++i) {
        inXcoords[i] = j["xcoords"][i];       // X[0][i] -> x[i] (x-coordinates)
        inYcoords[i] = j["ycoords"][i];   // X[1][i] -> x[N + i] (y-coordinates)
    }
   

    // Read input variables
    Real targetAlphaDeg = j["alpha_degrees"].get<double>();
    Real Re = j["Re"].get<double>();
    Real Ma = j["Ma"].get<double>();
    Real rhoInf = j["rho"].get<double>();
    Real nuInf = j["nu"].get<double>();
    Real custChord = j["chord"].get<double>();

    int doRestart = j["restart"].get<int>();

    Real sampleTE = j["sampleTE"].get<double>();
    const Real X = j["X"].get<double>();
    const Real Y = j["Y"].get<double>();
    const Real Z = j["Z"].get<double>();
    const Real S = j["S"].get<double>();
    const Real Ncrit = j["ncrit"].get<double>();
    const int doCps = j["returnData"].get<int>();
    const Real Ufac = j["Ufac"].get<double>();
    const Real TEfac = j["TEfac"].get<double>();
    const std::string model = j["model"].get<std::string>();
    // forcing transition variables
    const bool force = j["forcetrans"].get<int>();
    const Real topTransPos = j["toptrans"].get<double>();
    const Real botTransPos = j["bottrans"].get<double>();

    bool converged = runCode(doRestart,
        Ncrit,Ufac,TEfac,custChord,inXcoords,inYcoords,
        targetAlphaDeg,Re,Ma,rhoInf,nuInf,
        topTransPos,botTransPos,force,
        model,sampleTE,X,Y,Z,S,doCps);
    
    return converged;
    }
    
};