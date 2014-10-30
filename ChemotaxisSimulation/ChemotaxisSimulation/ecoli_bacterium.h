#pragma once

#include <boost/operators.hpp>
#include <utility>
#include <vector>
#include <boost/array.hpp>
#include <cmath>
#include <boost/numeric/ublas/vector.hpp>
#include "geometry.h"
#include "random.h"






class ecoli_bacterium {
    
    //constexpr static const double h = 0.01;
    double speed_;
    point2d position_, direction_;
    boost::numeric::ublas::vector< double > pathway_;
    
    
    
public:
    
    ecoli_bacterium(point2d pos, rng_type& rng);
    
    void update_pathway(double dt);
    void fluctuate_r(double dt, rng_type& rng);
    void update_bacterium(double dt, rng_type& rng, polygon_type& poly);
    void tumble(double dt, rng_type& rng);
    void run(double dt, polygon_type& poly);
    
    point2d get_position();
    
    boost::numeric::ublas::vector< double > get_pathway(){return pathway_;}
    void set_concentration(double conc);
    
};