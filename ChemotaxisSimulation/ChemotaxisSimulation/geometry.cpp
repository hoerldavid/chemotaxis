//
//  geometry.cpp
//  ChemotaxisSimulation
//
//  Created by David Hörl on 14.10.14.
//  Copyright (c) 2014 David Hörl. All rights reserved.
//

#include "geometry.h"


double cross2d(point2d a, point2d b)
{
    return a.x * b.y - b.x * a.y;
}

bool intersect_lines(point2d& a1, point2d& a2, point2d& b1, point2d& b2, point2d& res)
{
    auto va = a2 - a1;
    auto vb = b2 - b1;
    
    double denom = cross2d(va, vb);
    
    if (denom == 0) return false; //parallel lines
    
    double t = cross2d(b1 - a1, vb) / denom;
    double u = cross2d(b1 - a1, va) / denom;
    
    if ((t >= 0 && t <= 1) && (u >= 0 && u <= 1))
    {
        res = a1 + t * va;
        return true;
    }
    
    return false;
    
}

double dot2d(point2d v1, point2d v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

double get_angle(point2d v1, point2d v2)
{
    double d = dot2d(v1, v2);
    double n1 = v1.length();
    double n2 = v2.length();
    
    return acos(d/(n1*n2));
}

point2d reflect_vector(point2d vi, point2d vs)
{
    auto di = vi / vi.length();
    auto ds = vs / vs.length();
    ds = point2d(- ds.y, ds.x);
    
    double cosA = dot2d(ds, di);
    if (cosA > 0){
        ds *= -1;
        cosA = dot2d(ds, di);
    }
    point2d res = (di + ds * (-2 * cosA)) * vi.length();
    return res;
}

bool get_nearest_intersection(point2d& a1, point2d& a2, std::vector<point2d>& poly, bool ignore_start,
                              point2d& ress1, point2d& ress2, point2d& intersection )
{
    
    double dist = std::numeric_limits<double>::max();
    bool found_intersection = false;
    
    for (int i = 0; i < poly.size() - 1; i++ )
    {
        point2d t_intersection;
        bool does_intersect = intersect_lines(a1, a2, poly[i], poly[i+1], t_intersection);
        
        if(!does_intersect) continue;
        
        if (ignore_start && (t_intersection-a1).length() <= std::numeric_limits<double>::epsilon()) continue;
        
        if ((t_intersection-a1).length() < dist){
            dist = (t_intersection-a1).length();
            found_intersection = true;
            intersection = t_intersection;
            ress1 = poly[i];
            ress2 = poly[i+1];
        }
        
    }
    return found_intersection;
    
}


void move_and_reflect(point2d& from, point2d& to, std::vector<point2d>& poly, bool ignore_start,
                      point2d& res_point, point2d& res_direction)
{
    point2d ress1;
    point2d ress2;
    point2d intersection;
    
    bool found_intersection = get_nearest_intersection(from, to, poly, ignore_start, ress1, ress2, intersection);
    
    if (!found_intersection){
        res_point = to;
        res_direction = (to - from) / (to-from).length();
        return;
    }
    
    point2d reflected = reflect_vector(to-intersection, ress2-ress1);
    
    auto t_from = intersection;
    auto t_to = intersection + reflected;
    
    move_and_reflect(t_from, t_to, poly, true, res_point, res_direction);
}

/* convert indices in grid to space coordinates */
point2d grid_to_space(int x, int y, double grid_size, point2d min, point2d max)
{
    return point2d(min.x + grid_size * x, min.y + grid_size * y);
}

/* convert space coordinates to indices in grid starting at min */
std::pair<int, int> space_to_grid(point2d p, double grid_size, point2d min, point2d max)
{
    int x_coord = (int) floor((p - min).x / grid_size);
    int y_coord = (int) floor((p - min).y / grid_size);
    
    return std::make_pair(x_coord, y_coord);
}