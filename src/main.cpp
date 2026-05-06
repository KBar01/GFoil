#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>

// Unified physics headers
#include "codi.hpp"
using Real = codi::RealReverse;
using Tape = typename Real::Tape;

#include "real_type.hpp"
#include "data_structs.hpp"
#include "spline.hpp"
#include "main_func.hpp"
#include "get_funcs.hpp"
#include "panel_funcs.hpp"
#include "calc_ue_m.hpp"
#include "solve_inv.hpp"
#include "build_global_sys.hpp"
#include "solve_glob.hpp"
#include "clear_RV.hpp"
#include "stagmove.hpp"
#include "coupled.hpp"
#include "init_BL.hpp"
#include "extract_BL_TE.hpp"
#include "sound.hpp"
#include "nlohmann/json.hpp"

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
    Oper<Real> oper(alpha,Re,Ma);
    oper.rho = rhoInf;

    Geom<Real> geom;

    Real flattenedCoords[2*Ncoords]={0};
    Real inCoords[2*Nin]={0};
    for (int i=0;i<Nin;++i){
        inCoords[colMajorIndex(0,i,2)] = inXcoords[i];
        inCoords[colMajorIndex(1,i,2)] = inYcoords[i];
    }
    make_panels<Real>(inCoords,flattenedCoords,Ufac,TEfac);

    Foil<Real> foil(flattenedCoords);
    Isolc<Real> isolc;
    Isolv<Real> isol_v;
    Param<Real> param;
    param.ncrit = nCrit;
    Wake<Real> wake;
    Vsol<Real> vsol;
    Glob<Real> glob;

    build_gamma_codi<Real>(isolc,foil,oper);
    init_thermo<Real>(oper,param,geom);
    build_wake<Real>(foil,geom,oper,isolc,wake);
    stagpoint_find<Real>(isolc,isol_v,foil,wake);
    identify_surfaces<Real>(isol_v,vsol);
    set_wake_gap<Real>(foil,isol_v,vsol);

    // Find forced transition nodes after surface identification so searches use
    // the correct surface node lists. geom.chord is correct after init_thermo.
    Real xTransTop = topTransPos * geom.chord;
    Real xTransBot = botTransPos * geom.chord;

    int idx_closest_bot = vsol.Is[0].empty() ? 0 : vsol.Is[0].front();
    int idx_closest_top = vsol.Is[1].empty() ? Ncoords - 1 : vsol.Is[1].back();

    if (force) {
        for (int i : vsol.Is[0]) {
            if (flattenedCoords[colMajorIndex(0,i,2)] >= xTransBot) { idx_closest_bot=i; break; }
        }
        for (int i : vsol.Is[1]) {
            if (flattenedCoords[colMajorIndex(0,i,2)] >= xTransTop) { idx_closest_top=i; break; }
        }
    }

    Trans<Real> tdata;
    tdata.transNode[0] = idx_closest_bot;
    tdata.transNode[1] = idx_closest_top;
    tdata.transPos[0] = xTransBot;
    tdata.transPos[1] = xTransTop;

    calc_ue_m<Real>(foil,wake,isolc,vsol);
    rebuild_ue_m<Real>(foil,wake,isol_v,vsol,false);

    if (fromRestart){
        std::ifstream prevfile("restart.json");
        json j; prevfile >> j;
        const int nsys = (Ncoords+Nwake);
        for (int i=0;i<4*nsys;++i) glob.U[i] = j["states"][i].get<double>();
        for (int i=0;i<nsys;++i) vsol.turb[i] = j["turb"][i].get<int>();
        isol_v.stagIndex[0] = j["stag"][0]; isol_v.stagIndex[1] = j["stag"][1];
        if (force) { tdata.isForced[0]=1; tdata.isForced[1]=1; }
    } else {
        init_boundary_layer<Real>(oper,foil,param,isolc,isol_v,vsol,glob,tdata,force);
    }

    stagpoint_move<Real>(isol_v,glob,foil,wake,vsol);
    bool converged = solve_coupled<Real>(oper,foil,wake,param,vsol,isolc,isol_v,glob,tdata,force);
    Post<Real> post;
    calc_force<Real>(oper,geom,param,foil,glob,post);
    
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

    interpolate_at_95_both_surfaces<Real>(xcoords,glob.U.data(),post.cp.data(),oper,vsol.turb.data(),param,topsurf,botsurf,Uinf,sampleTE,chordScaling);
    Real OASPL = calc_OASPL_AD<Real>(botsurf,topsurf,chordScaling,Uinf,X,Y,Z,S,kinViscInf,rhoInf,model);

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

            //  calc transition point
            Real botTransX = geom.chord;
            for(int i=0;i<isol_v.stagIndex[0];++i){
                int isTurb = vsol.turb[isol_v.stagIndex[0] - i];
                if (isTurb){
                    botTransX = foil.x[colMajorIndex(0,isol_v.stagIndex[0]-i,2)];
                    break;
                } 
            }
            Real topTransX = geom.chord;
            for(int i=0;i<200-isol_v.stagIndex[1];++i){
                int isTurb = vsol.turb[isol_v.stagIndex[1] + i];
                if (isTurb){
                    topTransX = foil.x[colMajorIndex(0,isol_v.stagIndex[1]+i,2)];
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
            
            out["stagnation"] = isol_v.stagIndex;
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
        Real OASPL = calc_OASPL_AD<Real>(botsurf,topsurf,chordScaling,Uinf,X,Y,Z,S,kinViscInf,rhoInf,model);
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