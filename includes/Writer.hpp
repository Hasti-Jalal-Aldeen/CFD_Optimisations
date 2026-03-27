#pragma once
#include "Global.hpp"
#include <string>

class Writer {
    public:
    Writer(int nvars, const std::string fname, const std::string varnames[NVARS], int Nc, double* arrs[NVARS]);
    void write_soln();

    private:
    int runcount;
    const std::string fname;
    std::string fmt;

    int nchars;
    const int nvars, Nc;
    double* arrs[NVARS];
    std::string varnames[NVARS];
};
