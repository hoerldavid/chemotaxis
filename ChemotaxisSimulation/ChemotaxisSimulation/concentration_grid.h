//
//  concentration_grid.h
//  ChemotaxisSimulation
//
//  Created by David Hörl on 14.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#ifndef __ChemotaxisSimulation__concentration_grid__
#define __ChemotaxisSimulation__concentration_grid__

#include <stdio.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <cmath>

#include "geometry.h"


typedef boost::numeric::ublas::matrix< double > concentration_state_type;

struct concentration_ode
{
    concentration_state_type* prod_conc;
    double D; // diffusion coefficient
    double p; // production rate
    double a; // degradation rate
    
    concentration_ode(concentration_state_type* _prod_conc, double _D, double _a, double _p) :
            prod_conc(_prod_conc), D(_D), a(_a), p(_p) {}
    
    void operator()( const concentration_state_type &x , concentration_state_type &dxdt , double /* t */ ) const
    {
        size_t size1 = x.size1() , size2 = x.size2();
        
        for( size_t i=0 ; i<size1 ; ++i )
        {
            for( size_t j=0 ; j<size2 ; ++j )
            {
                
                if (x(i, j) == -1){ // no change outside of arena
                    dxdt(i,j) = 0;
                    continue;
                }
                
                double diff = 0;
                
                if (i > 0 && x( i - 1 , j ) != -1) {diff += x( i - 1 , j ) - x( i , j );}
                if (i < size1 -1 && x( i + 1 , j ) != -1) {diff += x( i + 1 , j ) - x( i , j );}
                if (j > 0 && x( i , j - 1 ) != -1) {diff += x( i , j - 1 ) - x( i , j );}
                if (j < size2 -1 && x( i , j + 1 ) != -1) {diff += x( i , j + 1 ) - x( i , j );}
                
                dxdt( i , j ) = diff * D + p * (*prod_conc)(i, j) - a * x(i, j);
            }
        }
    }
};


struct logistic_updater{
    double mu_, lambda_, a_;
    logistic_updater() {};
    logistic_updater(double mu, double lambda, double a) : mu_(mu), lambda_(lambda), a_(a) {}
    double get_conc(double t, double initial_conc)
    {
        double y = a_ / (1 + exp(4 * mu_ / a_ * (lambda_ - t) + 2));
        return exp(y) * initial_conc;
    }
};

class concentration_grid {
    concentration_state_type conc_;
    concentration_state_type prod_conc_;
    point2d min;
    point2d max;
    double grid_width = 0.5; // TODO: adaptive
    double initial_prod_conc_per_point;
    
    logistic_updater prod_conc_updater;
    
    std::pair<int, int> prod_from_idx, prod_to_idx;
    
    double D; // diffusion coefficient
    double p; // production rate
    double a; // degradation rate
    
public:
    
    concentration_grid() {}
    concentration_grid(double _grid_width) : grid_width(_grid_width) {}
    void setup_from_polygon(polygon_type& poly);
    
    concentration_state_type get_conc();
    concentration_state_type get_producer_conc();
    double get_conc_at_point(point2d p);
    void setup_producer_conc(double initial_conc, point2d from, point2d to);
    void setup_producer_growth(double mu, double lambda, double a) {prod_conc_updater = logistic_updater(mu, lambda, a);}
    
    void update_producer_conc(double t);
    void update_conc(double dt);
    void set_ode_params(double _D, double _a, double _p) {
        D = _D; a = _a; p = _p;
    }
    
    std::pair<int, int> get_dimensions() {
        return std::make_pair(conc_.size1(), conc_.size2());
    }
    
    double* get_raw_conc(){
        
        double* res = new double [conc_.size1() * conc_.size2()];
        
        int size2 = conc_.size2();
        for (int i = 0; i < conc_.size1(); i ++) {
            for (int j = 0; j < conc_.size2(); j ++) {
                res[i*size2 + j] = conc_(i, j);
            }
        }
        
        return res;
    }
    
    double get_producer_conc_at_colony(){
        return prod_conc_(prod_from_idx.second, prod_from_idx.first);
    }
    
};




#endif /* defined(__ChemotaxisSimulation__concentration_grid__) */
