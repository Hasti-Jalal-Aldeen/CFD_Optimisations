#include "Writer.hpp"
#include <cstdio>
#include <format>
#include <iostream>

Writer::Writer(int nvars, const std::string fname, const std::string varnames[NVARS], int Nc, double* arrs[NVARS]) 
: nvars(nvars), Nc(Nc), fname(fname)
{
    runcount = 0;
    
    for (int var = 0; var < nvars; var++) {
        this->arrs[var]     = arrs[var];
        this->varnames[var] = varnames[var];
    }

    nchars = 0;
    for (int var = 0; var < nvars; var++) {
        nchars += varnames[var].length()+1;
    }
}

void Writer::write_soln() {
    fmt = std::format("{}_{}.hsln", fname, runcount);

    FILE *file = fopen(fmt.c_str(), "wb");

    fwrite(    &Nc, sizeof(int), 1, file);
    fwrite( &nvars, sizeof(int), 1, file);
    fwrite(&nchars, sizeof(int), 1, file);
    for (int var = 0; var < nvars; var++) {
        fwrite(varnames[var].c_str(), sizeof(char), varnames[var].length()+1, file);
    }
    for (int var = 0; var < nvars; var++) {
        fwrite(arrs[var], sizeof(double), Nc, file);
    }

    fclose(file);
    runcount++;
}
