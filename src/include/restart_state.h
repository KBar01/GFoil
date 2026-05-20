#pragma once
#include <vector>

struct RestartState {
    std::vector<double> states;
    std::vector<int>    turb;
    std::vector<int>    stag;
    std::vector<double> RVvals;
    std::vector<int>    RVrows;
    std::vector<int>    RVcols;
    int                 RVnz = 0;
};

struct ForwardResult {
    bool   converged = false;
    double CL = 0.0, CD = 0.0, CM = 0.0, OASPL = 0.0;
};
