//
//  mx_alloc_simulation.cpp
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
    
        
        simulation_wrapper* myObj = new simulation_wrapper();
        plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
        auto pointer_to_obj = mxGetPr(plhs[0]);
        pointer_to_obj[0] = (long) myObj;
    
    
}