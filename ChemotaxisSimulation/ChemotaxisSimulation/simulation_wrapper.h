//
//  simulation_wrapper.h
//  ChemotaxisSimulation
//
//  Created by David Hörl on 15.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#ifndef __ChemotaxisSimulation__simulation_wrapper__
#define __ChemotaxisSimulation__simulation_wrapper__

#include <stdio.h>
#include <vector>

#include "concentration_grid.h"
#include "ecoli_bacterium.h"
#include "random.h"
#include "geometry.h"


class simulation_wrapper {
    double t_;
    rng_type rng_;
    
    polygon_type border_;
    
    concentration_grid conc_;
    std::vector<ecoli_bacterium> cells_;
    
public:
    
    simulation_wrapper();
    
    void insert_cells(double* cells, size_t n);
    double* get_cell_positions();
    
    void setup_border(double* points, size_t n, double* grid_width);
    
    void setup_producers(double* rect, double* conc);
    
    void setup_ode_params(double* D, double* a, double*p){
        conc_.set_ode_params(*D, *a, *p);
    }
    
    void setup_growth_params(double* mu, double* lambda, double* a){
        conc_.setup_producer_growth(*mu, *lambda, *a);
    }
    
    void update(double* total_time, double* dt_concentration, double* dt_chemotaxis);
    
    std::pair<int, int> get_concentration_dimensions();
    double* get_concentration();
    
    int get_n_of_cells();
    
    double get_prod_conc_at_colony()
    {
        return conc_.get_producer_conc_at_colony();
    }
    
    double* get_pathway(int idx){
        double* res = new double[10];
        
        auto pw = cells_[idx].get_pathway();
        
        for (int i = 0; i < 10; i++){
            res[i] = pw(i);
        }
        
        return res;           
        
    }
    
    
};

#endif /* defined(__ChemotaxisSimulation__simulation_wrapper__) */
