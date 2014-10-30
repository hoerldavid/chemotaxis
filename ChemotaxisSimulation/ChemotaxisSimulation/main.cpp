//
//  main.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 13.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#include <iostream>
#include <boost/container/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "geometry.h"
#include "random.h"
#include "concentration_grid.h"
#include "simulation_wrapper.h"

using namespace std;
using namespace boost::numeric::ublas;

int main(int argc, const char * argv[]) {
    // insert code here...
    
    simulation_wrapper sim;
    double border[26] = { 0,0,0,50,20,50,20,30,30,30,30,50,50,50,50,0,30,0,30,20,20,20,20,0,0,0 };
    double border_width = 1;
    sim.setup_border(border, 26, &border_width);
    double prods[4] = {5, 20, 15, 30};
    double prod_conc = 20;
    sim.setup_producers(prods, &prod_conc);
    double cells[2] = {40, 25};
    sim.insert_cells(cells, 2);
    
    
    double update_t = 100, dt = 1;
    sim.update(&update_t, &dt);

    
    return 0;
}
