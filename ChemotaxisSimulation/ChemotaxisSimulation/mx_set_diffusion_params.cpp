//
//  mx_set_diffusion_params.cpp
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
    
    
    auto D = mxGetPr(prhs[1]);
    auto a = mxGetPr(prhs[2]);
    auto p = mxGetPr(prhs[3]);

    
    myObj->setup_ode_params(D, a, p);
}