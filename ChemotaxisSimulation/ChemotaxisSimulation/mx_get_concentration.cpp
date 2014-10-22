//
//  mx_get_concentration.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 17.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#include <vector>
#include <stdio.h>
#include <mex.h>
#include <matrix.h>

#include "simulation_wrapper.h"

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    
    
    
    
    auto pointer_to_obj = mxGetPr(prhs[0]);
    simulation_wrapper* myObj = (simulation_wrapper*) (long)pointer_to_obj[0];
    
    auto dim = myObj->get_concentration_dimensions();
    auto conc = myObj->get_concentration();
    
    plhs[0] = mxCreateDoubleMatrix(dim.first, dim.second, mxREAL);
    auto pointer_to_conc = mxGetPr(plhs[0]);
    
    memcpy(pointer_to_conc, conc, sizeof(double) * dim.first * dim.second);
    
    delete conc;
    
}