//
//  mx_update_simulation.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 16.10.14.
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
    
    
    auto time_update = mxGetPr(prhs[1]);
    auto dt_concentration = mxGetPr(prhs[2]);
    auto dt_chemotaxis = mxGetPr(prhs[3]);
    
    myObj->update(time_update, dt_concentration, dt_chemotaxis);
}