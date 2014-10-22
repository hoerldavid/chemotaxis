//
//  random.h
//  ChemotaxisSimulation
//
//  Created by David Hörl on 14.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#ifndef ChemotaxisSimulation_random_h
#define ChemotaxisSimulation_random_h

#include <random>

typedef std::default_random_engine rng_type;

double get_rand_uniform(rng_type& rng, double from, double to);
double get_rand_gamma(rng_type& rng, double alpha, double beta);
rng_type get_rng_seeded();
#endif
