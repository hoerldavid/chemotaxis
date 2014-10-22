//
//  concentration_grid.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 14.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#include <limits>
#include <boost/geometry.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "concentration_grid.h"

typedef boost::geometry::model::d2::point_xy<double> geom_point_type;
typedef boost::geometry::model::polygon<geom_point_type> geom_poly_type;


void concentration_grid::setup_from_polygon(polygon_type& poly)
{
    double min_x = std::numeric_limits<double>::max(), min_y = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::min(), max_y = std::numeric_limits<double>::min();
    
    geom_poly_type t_poly;
    
    for (auto& p : poly)
    {
        if (p.x < min_x) min_x = p.x;
        if (p.x > max_x) max_x = p.x;
        if (p.y < min_y) min_y = p.y;
        if (p.y > max_y) max_y = p.y;
        
        geom_point_type t_point(p.x, p.y);
        boost::geometry::append(t_poly, t_point);
    }
    
    min = point2d(min_x, min_y);
    max = point2d(max_x, max_y);
    
    int g_width = (int)floor((max_x - min_x) / grid_width);
    int g_height = (int)floor((max_y - min_y) / grid_width);
    
    conc_ = concentration_state_type(g_height, g_width);
    
    for (int y = 0; y < conc_.size1(); y++){
        for (int x = 0; x < conc_.size2(); x++)
        {
            point2d t_p = grid_to_space(x, y, grid_width, min, max);
            geom_point_type t_p_geo(t_p.x, t_p.y);
            
            if (boost::geometry::covered_by(t_p_geo, t_poly)) {
                conc_(y,x) = 0.0;
            } else {
                conc_(y,x) = -1.0;
            }
            
        }
    }
    
    std::cout << conc_ << std::endl;
    
}


concentration_state_type concentration_grid::get_conc()
{
    return conc_;
}

double concentration_grid::get_conc_at_point(point2d p)
{
    // outside of grid -> concentration = 0 (should not happen)
    if (p.x < min.x && p.x > max.x && p.y < min.y && p.y > max.y) return 0;
    
    double num = 0;
    double denom = 0;
    
    auto idx = space_to_grid(p, grid_width, min, max);
    
    num += conc_(idx.second, idx.first);
    denom += 1;
    
    if (idx.first < conc_.size2() - 1) {
        num += conc_(idx.second, idx.first + 1);
        denom += 1;
    }
    
    if (idx.second < conc_.size1() - 1) {
        num += conc_(idx.second + 1, idx.first);
        denom += 1;
    }
    
    if (idx.first < conc_.size2() - 1 && idx.second < conc_.size1() - 1) {
        num += conc_(idx.second + 1, idx.first + 1);
        denom += 1;
    }
    
    return num/denom;
    
}


void concentration_grid::setup_producer_conc(double initial_conc, point2d from, point2d to){
    
    prod_conc_ = concentration_state_type(conc_.size1(), conc_.size2(), 0.0);
    prod_from_idx = space_to_grid(from, grid_width, min, max);
    prod_to_idx = space_to_grid(to, grid_width, min, max);
    
    int points = (prod_to_idx.first - prod_from_idx.first) * (prod_to_idx.second - prod_from_idx.second);
    initial_prod_conc_per_point = initial_conc / points;
    
    using namespace boost::numeric::ublas;
    
    matrix_range<concentration_state_type > mr (prod_conc_ , range (prod_from_idx.second, prod_to_idx.second), range (prod_from_idx.first, prod_to_idx.first));
    for (int i = 0; i < mr.size1 (); ++ i)
        for (int j = 0; j < mr.size2 (); ++ j)
            mr (i, j) = initial_prod_conc_per_point;
    
    
}

concentration_state_type concentration_grid::get_producer_conc()
{
    return prod_conc_;
}

void concentration_grid::update_producer_conc(double t)
{
    using namespace boost::numeric::ublas;
    
    matrix_range<concentration_state_type > mr (prod_conc_ , range (prod_from_idx.second, prod_to_idx.second), range (prod_from_idx.first, prod_to_idx.first));
    for (int i = 0; i < mr.size1 (); ++ i)
        for (int j = 0; j < mr.size2 (); ++ j)
            mr (i, j) = prod_conc_updater.get_conc(t, initial_prod_conc_per_point);
}

void concentration_grid::update_conc(double dt){
    
    using namespace boost::numeric::odeint;
    
    //typedef runge_kutta_dopri5<concentration_state_type> stepper_t;
    
    typedef runge_kutta4<concentration_state_type> stepper_t;

    //auto stepper = make_controlled(1.0e-6, 1.0e-6, stepper_t());
    
    concentration_ode ode(&prod_conc_, D, a, p);
    size_t steps = integrate_const( stepper_t() , ode ,
                    conc_ , 0.0 , dt, 0.01 ); // TODO: implicit euler doesnt work! stepsize??
}