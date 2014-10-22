//
//  mx_get_producer_conc.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 17.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#include <stdio.h>


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
    
    auto conc = myObj->get_prod_conc_at_colony();
    
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    auto pointer_to_conc = mxGetPr(plhs[0]);
    pointer_to_conc[0] = conc;
    
    
}