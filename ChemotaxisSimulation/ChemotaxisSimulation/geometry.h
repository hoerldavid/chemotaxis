//
//  geometry.h
//  ChemotaxisSimulation
//
//  Created by David Hörl on 14.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#ifndef ChemotaxisSimulation_geometry_h
#define ChemotaxisSimulation_geometry_h

#include <boost/operators.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <vector>
#include <limits>
#include <utility>




struct point2d
: boost::addable< point2d          // point + point
, boost::subtractable< point2d     // point - point
, boost::dividable2< point2d, double    // point / T
, boost::multipliable2< point2d, double // point * T, T * point
> > > >
{
    double x, y;
    point2d(double _x, double _y) : x(_x), y(_y) {}
    point2d() : x(0), y(0) {}
    
    point2d& operator+=(const point2d& rhs) {x += rhs.x; y += rhs.y; return *this;}
    
    point2d& operator-=(const point2d& rhs) {x -= rhs.x; y -= rhs.y; return *this;}
    point2d& operator*=(double rhs) {x *= rhs; y *= rhs; return *this;}

    point2d& operator/=(double rhs) {x /= rhs; y /= rhs; return *this;}
    
    double squared_length() {return x*x + y*y;}
    double length() {return sqrt(squared_length());}
    
    void rotate(double degs){
        
        double rads = (degs / 180) * boost::math::constants::pi<double>();
        
        double t_x = cos(rads) * x - sin (rads) * y;
        double t_y = sin(rads) * x + cos(rads) * y;
        
        x = t_x;
        y = t_y;
        
    }
    
};


double cross2d(point2d a, point2d b);
bool intersect_lines(point2d& a1, point2d& a2, point2d& b1, point2d& b2, point2d& res);
double dot2d(point2d v1, point2d v2);
double get_angle(point2d v1, point2d v2);
point2d reflect_vector(point2d vi, point2d vs);

bool get_nearest_intersection(point2d& a1, point2d& a2, std::vector<point2d>& poly, bool ignore_start,
                              point2d& ress1, point2d& ress2, point2d& intersection );

void move_and_reflect(point2d& from, point2d& to, std::vector<point2d>& poly, bool ignore_start,
                      point2d& res_point, point2d& res_direction);

point2d grid_to_space(int x, int y, double grid_size, point2d min, point2d max);
std::pair<int, int> space_to_grid(point2d p, double grid_size, point2d min, point2d max);

typedef std::vector<point2d> polygon_type;


#endif
