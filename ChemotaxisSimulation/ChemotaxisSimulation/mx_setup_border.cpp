//
//  mx_setup_border.cpp
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
    
    
    auto border = mxGetPr(prhs[1]);
    auto size1 = mxGetM(prhs[1]);
    auto size2 = mxGetN(prhs[1]);
    
    auto grid_width = mxGetPr(prhs[2]);
    

    myObj->setup_border(border, size1 * size2, grid_width);
    
}