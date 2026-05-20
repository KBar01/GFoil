#include <iostream>
#include <vector>
#include <cmath>
#include "codi.hpp"
#include "real_type.h"
#include "data_structs.h"
#include "run_forward.h"
#include "sound.hpp"
#include <fstream>
#include <string>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

int main() {

    std::ifstream infile("input.json");
    if (!infile) {
        std::cerr << "Failed to open input.json\n";
        return 1;
    }
    json j;
    infile >> j;

    const int doWPSonly = j["WPSonly"].get<int>();

    if (doWPSonly) {

        Real rhoInf      = j["rho"].get<double>();
        Real kinViscInf  = j["nu"].get<double>();
        Real chordScaling = j["chord"].get<double>();
        Real Re          = j["Re"].get<double>();
        const Real S     = j["S"].get<double>();

        std::vector<double> obsX_d, obsY_d, obsZ_d;
        if (j["X"].is_array()) {
            obsX_d = j["X"].get<std::vector<double>>();
            obsY_d = j["Y"].get<std::vector<double>>();
            obsZ_d = j["Z"].get<std::vector<double>>();
        } else {
            obsX_d = { j["X"].get<double>() };
            obsY_d = { j["Y"].get<double>() };
            obsZ_d = { j["Z"].get<double>() };
        }
        int nObs = static_cast<int>(obsX_d.size());
        std::vector<Real> obsX(nObs), obsY(nObs), obsZ(nObs);
        for (int i = 0; i < nObs; ++i) {
            obsX[i] = obsX_d[i];
            obsY[i] = obsY_d[i];
            obsZ[i] = obsZ_d[i];
        }

        const Real toptheta  = j["toptheta"].get<double>();
        const Real topdstar  = j["topdstar"].get<double>();
        const Real topdelta  = j["topdelta"].get<double>();
        const Real toptauw   = j["toptauw"].get<double>();
        const Real toptaumax = j["toptaumax"].get<double>();
        const Real topue     = j["topue"].get<double>();
        const Real topdpdx   = j["topdpdx"].get<double>();

        const Real bottheta  = j["bottheta"].get<double>();
        const Real botdstar  = j["botdstar"].get<double>();
        const Real botdelta  = j["botdelta"].get<double>();
        const Real bottauw   = j["bottauw"].get<double>();
        const Real bottaumax = j["bottaumax"].get<double>();
        const Real botue     = j["botue"].get<double>();
        const Real botdpdx   = j["botdpdx"].get<double>();
        const std::string model = j["model"].get<std::string>();
        const int aWeighting    = j.value("aWeighting", 0);

        Real topsurf[7], botsurf[7];
        topsurf[0] = toptheta;  botsurf[0] = bottheta;
        topsurf[1] = topdstar;  botsurf[1] = botdstar;
        topsurf[2] = toptaumax; botsurf[2] = bottaumax;
        topsurf[3] = topue;     botsurf[3] = botue;
        topsurf[4] = topdpdx;  botsurf[4] = botdpdx;
        topsurf[5] = toptauw;   botsurf[5] = bottauw;
        topsurf[6] = topdelta;  botsurf[6] = botdelta;

        Real Uinf = (Re * kinViscInf) / (chordScaling);
        Real OASPL = calc_OASPL<Real, true>(botsurf, topsurf, chordScaling, Uinf,
            obsX.data(), obsY.data(), obsZ.data(), nObs, S, kinViscInf, rhoInf, model, 1, aWeighting);
        return 1;

    } else {

        Real inXcoords[Nin] = {0}, inYcoords[Nin] = {0};
        for (int i = 0; i < Nin; ++i) {
            inXcoords[i] = j["xcoords"][i];
            inYcoords[i] = j["ycoords"][i];
        }

        Real targetAlphaDeg = j["alpha_degrees"].get<double>();
        Real Re             = j["Re"].get<double>();
        Real Ma             = j["Ma"].get<double>();
        Real rhoInf         = j["rho"].get<double>();
        Real nuInf          = j["nu"].get<double>();
        Real custChord      = j["chord"].get<double>();
        int  doRestart      = j["restart"].get<int>();
        Real sampleTE       = j["sampleTE"].get<double>();
        const Real S        = j["S"].get<double>();
        const Real Ncrit      = j["ncrit"].get<double>();
        const Real Ncrithyst  = j.value("ncrithyst", 0.2);
        const int  doCps    = j["returnData"].get<int>();
        const Real Ufac     = j["Ufac"].get<double>();
        const Real TEfac    = j["TEfac"].get<double>();
        const std::string model = j["model"].get<std::string>();
        const int aWeighting    = j.value("aWeighting", 0);

        std::vector<double> obsX_d, obsY_d, obsZ_d;
        if (j["X"].is_array()) {
            obsX_d = j["X"].get<std::vector<double>>();
            obsY_d = j["Y"].get<std::vector<double>>();
            obsZ_d = j["Z"].get<std::vector<double>>();
        } else {
            obsX_d = { j["X"].get<double>() };
            obsY_d = { j["Y"].get<double>() };
            obsZ_d = { j["Z"].get<double>() };
        }
        int nObs = static_cast<int>(obsX_d.size());
        std::vector<Real> obsX(nObs), obsY(nObs), obsZ(nObs);
        for (int i = 0; i < nObs; ++i) {
            obsX[i] = obsX_d[i];
            obsY[i] = obsY_d[i];
            obsZ[i] = obsZ_d[i];
        }

        bool converged = runCode(doRestart,
            Ncrit, Ufac, TEfac, custChord, inXcoords, inYcoords,
            targetAlphaDeg, Re, Ma, rhoInf, nuInf,
            model, sampleTE, obsX.data(), obsY.data(), obsZ.data(), nObs, S, doCps,
            nullptr, nullptr, nullptr, aWeighting, Ncrithyst);

        return converged;
    }
}
