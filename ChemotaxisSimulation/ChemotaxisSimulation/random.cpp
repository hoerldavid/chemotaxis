//
//  random.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 14.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#include "random.h"
#include <chrono>

double get_rand_uniform(rng_type &rng, double from, double to)
{
    std::uniform_real_distribution<double> distribution(from, to);
    return distribution(rng);
}
double get_rand_gamma(rng_type& rng, double alpha, double beta){
    std::gamma_distribution<double> distribution(alpha, beta);
    return distribution(rng);
}

rng_type get_rng_seeded(){
    rng_type res = rng_type();
    long seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    res.seed((int)seed);
    return res;
}