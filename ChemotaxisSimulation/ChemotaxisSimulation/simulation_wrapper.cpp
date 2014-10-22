//
//  simulation_wrapper.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 15.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#include "simulation_wrapper.h"


simulation_wrapper::simulation_wrapper()
{
    conc_ = concentration_grid(1);
    rng_ = get_rng_seeded();
    
    //conc_.set_ode_params(1, .3, .01);
    //conc_.setup_producer_growth(.01, 100, 100);
    
    t_ = 0;
}

void simulation_wrapper::insert_cells(double* cells, size_t n)
{
    for (int i = 0; i < n; i = i + 2)
    {
        ecoli_bacterium t_bact = ecoli_bacterium(point2d(cells[i], cells[i+1]), rng_);
        cells_.push_back(t_bact);
    }
}

double* simulation_wrapper::get_cell_positions()
{
    size_t res_n = cells_.size() * 2;
    double* res = new double [res_n];
    
    for (int i = 0; i < cells_.size(); i++){
        auto t_p = cells_[i].get_position();
        res[2*i] = t_p.x;
        res[2*i+1] = t_p.y;
    }
    
    return res;
}

int simulation_wrapper::get_n_of_cells()
{
    return cells_.size();
}

void simulation_wrapper::setup_border(double* points, size_t n, double* grid_width)
{
    for (int i = 0; i < n; i = i + 2)
    {
        point2d t_p(points[i], points[i+1]);
        border_.push_back(t_p);
    }
    
    conc_ = concentration_grid(*grid_width);
    conc_.setup_from_polygon(border_);
}

void simulation_wrapper::setup_producers(double* rect, double* conc)
{
    conc_.setup_producer_conc(*conc, point2d(rect[0], rect[1]), point2d(rect[2], rect[3]));
}

void simulation_wrapper::update(double* total_time, double* dt)
{
    double target_time = t_ + *total_time;
    while (t_ < target_time) {
        conc_.update_producer_conc(t_);
        conc_.update_conc(*dt);
        
        for (auto& b : cells_)
        {
            auto t_p = b.get_position();
            double t_c = conc_.get_conc_at_point(t_p);
            b.set_concentration(t_c);
            
            b.update_bacterium(*dt, rng_, border_);
        }
        
        t_ += *dt;
    }
    
}

std::pair<int, int> simulation_wrapper::get_concentration_dimensions()
{
    return conc_.get_dimensions();
}

double* simulation_wrapper::get_concentration(){
    return conc_.get_raw_conc();
}