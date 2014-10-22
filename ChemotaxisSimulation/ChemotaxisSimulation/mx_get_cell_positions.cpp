//
//  mx_get_cell_positions.cpp
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
    
    auto n_cells = myObj->get_n_of_cells();
    auto pos = myObj->get_cell_positions();
    
    plhs[0] = mxCreateDoubleMatrix(2,n_cells,mxREAL);
    auto pointer_to_cells = mxGetPr(plhs[0]);
    
    memcpy(pointer_to_cells, pos, sizeof(double) * n_cells * 2);
    
    delete pos;

}