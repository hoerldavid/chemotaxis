#include <stdio.h>
/*
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4.hpp>
#include <boost/numeric/odeint/stepper/rosenbrock4_controller.hpp>
#include <boost/numeric/odeint/stepper/generation.hpp>
 */
#include <boost/numeric/odeint.hpp>
//#include <boost/numeric/ublas/matrix.hpp>
//#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <utility>

#include "ecoli_bacterium.h"



using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector< double > ecoli_chemotaxis_state_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;
//typedef runge_kutta_dopri5< ecoli_chemotaxis_state_type > ecoli_chemotaxis_stepper_type;
typedef rosenbrock4< double > ecoli_chemotaxis_stepper_type;


/* ODE function specifying the chemotaxis pathway */
struct ecoli_pathway_ode{
    void operator() ( const ecoli_chemotaxis_state_type &x , ecoli_chemotaxis_state_type &dxdt , const double /*t */ )
    {
        
        /* MEANING OF THE STATE VECTOR */
        /*
         x(0) - x(4) :=  T0 - T4
         % x(5) := Ap
         % x(6) := Yp
         % x(7) := Bp
         % x(8) := L (ligand concentration)
         % x(9) := R concentration
         */
        /* CONSTANTS */
        const double kR = 0.39;
        const double kB = 6.3;
        const double kpB = 3;
        const double kA = 50;
        const double kY = 100;
        
        const double K0 = 27e-04;
        const double K1 = 20e-03;
        const double K2 = 150e-03;
        const double K3 = 150e-02;
        const double K4 = 60e+00;
        
        const double V0 = 0;
        const double V1 = .25;
        const double V2 = .5;
        const double V3 = .75;
        const double V4 = 1;
        
        const double kZ = 30 / 3.8; // TODO: correct unit?
        const double KR = .099;
        const double KB = 2.5;
        const double gB = 1;
        const double gZ = .1;
        
        const double AT = 5.3;
        const double TT = 5.3; // same as AT
        const double BT = .28;
        const double R = .08;
        const double YT = 9.7;
        //Z == 3.8!!
        const double Z = 3.8;
        
        const double Hm = 1.2;
        
        // probabilities to be in active state %
        double p0L = V0 * (1-((pow(x[8],Hm))/(pow(x[8],Hm) + pow(K0,Hm))));
        double p1L = V1 * (1-((pow(x[8],Hm))/(pow(x[8],Hm) + pow(K1,Hm))));
        double p2L = V2 * (1-((pow(x[8],Hm))/(pow(x[8],Hm) + pow(K2,Hm))));
        double p3L = V3 * (1-((pow(x[8],Hm))/(pow(x[8],Hm) + pow(K3,Hm))));
        double p4L = V4 * (1-((pow(x[8],Hm))/(pow(x[8],Hm) + pow(K4,Hm))));
        
        //active receptors %
        double T0A = p0L * x[0];
        double T1A = p1L * x[1];
        double T2A = p2L * x[2];
        double T3A = p3L * x[3];
        double T4A = p4L * x[4];
        double TA = T0A + T1A + T2A + T3A + T4A;
        
        // dTM/dt
        dxdt[0] = kB * x[7] * ( T1A /(KB + TA)) - kR * x[9] * (x[0]/(KR + TT));
        dxdt[1] = kR * x[9] * (x[0] / (KR + TT)) + kB * x[7] * ( T2A /(KB + TA))
        - kR * x[9] * (x[1]/(KR + TT)) - kB * x[7] * (T1A / (KB + TA));
        dxdt[2] = kR * x[9] * (x[1] / (KR + TT)) + kB * x[7] * ( T3A /(KB + TA))
        - kR * x[9] * (x[2]/(KR + TT)) - kB * x[7] * (T2A / (KB + TA));
        dxdt[3] = kR * x[9] * (x[2] / (KR + TT)) + kB * x[7] * ( T4A /(KB + TA))
        - kR * x[9] * (x[3]/(KR + TT)) - kB * x[7] * (T3A / (KB + TA));
        dxdt[4] = kR * x[9] * (x[3] / (KR + TT)) - kB * x[7] * (T4A / (KB + TA));
        
        
        //dAp/dt
        dxdt[5] = kA * (AT - x[5]) * TA - kY * x[5] * (YT - x[6]) - kpB * x[5] * (BT - x[7]);
        // dYp/dt
        dxdt[6] = kY * x[5] * (YT - x[6]) - kZ * x[6] * Z - gZ * x[6];
        // dBp/dt
        dxdt[7] = kpB * x[5] * (BT - x[7]) - gB * x[7];
        
        // CONCENTRATION unaffected
        dxdt[8] = 0;
        // R doesnt change by ODE
        dxdt[9] = 0;
    }
};

/* Approximate Jacobian by forward differences */
/* adapted from http://www.maths.lth.se/na/courses/FMN081/FMN081-06/lecture7.pdf */
struct ecoli_pathway_jacobian{
    void operator() (const ecoli_chemotaxis_state_type &x , matrix_type &J , const double t, ecoli_chemotaxis_state_type &dfdt)
    {
        ecoli_pathway_ode f;
        const double eps = 1e-8;
        
        ecoli_chemotaxis_state_type fx(10);
        f(x, fx, t);
        
        auto x_perturb = x;
        
        for (int i = 0; i < x.size(); i++){
            x_perturb(i) += eps;
            
            ecoli_chemotaxis_state_type fx_perturb(10);
            f(x_perturb, fx_perturb, t);
            
            boost::numeric::ublas::column(J, i) = fx_perturb;
            x_perturb(i) = x(i);
            
        }
    }
};

ecoli_bacterium::ecoli_bacterium(point2d pos, rng_type& rng)
{
    //initialize pathway
    pathway_ = boost::numeric::ublas::vector<double> (10, 0.0);
    
    pathway_[0] = .15; // all receptors unmethylated
    pathway_[1] = 2;
    pathway_[2] = 2;
    pathway_[3] = 1;
    pathway_[4] = .15;
    
    pathway_[9] = .07; // initial R concentration
    
    /*
     T = randSum(5) * 5.3;
     
     b.T0 = T(1);
     b.T1 = T(2);
     b.T2 = T(3);
     b.T3 = T(4);
     b.T4 = T(5);
     
     % ?? okay to set to 0 ??
     b.Ap = 0;%rand() * 5.3;
     b.Yp = 0;%rand() * 9.7;
     b.Bp = 0;%rand() * .28;
     b.conc = 0;
     b.R = .07;
     */
    
    speed_ = .03; // mm/s
    position_ = pos;
    
    direction_ = point2d(get_rand_uniform(rng, -1, 1), get_rand_uniform(rng, -1, 1));
    direction_ /= direction_.length();
}

void ecoli_bacterium::update_pathway(double dt){
    
    // internal step size
    double h = (dt/10 > 0.01) ? 0.01 : dt/10;
    
    //typedef rosenbrock4_controller<ecoli_chemotaxis_stepper_type> controller;
    //controller my_controller(1.0e-6, 1.0e-6);
    
    auto my_controller = make_controlled(1.0e-6, 1.0e-6, runge_kutta_dopri5<ecoli_chemotaxis_state_type>());
    /* TODO: errors ok that way? */
    
    auto ode = ecoli_pathway_ode();
    //auto ode = std::make_pair(ecoli_pathway_ode(), ecoli_pathway_jacobian());
    
    size_t steps = integrate_adaptive( my_controller,
                                      ode, pathway_, static_cast<double>(0), dt, h);
}


void ecoli_bacterium::fluctuate_r(double dt, rng_type& rng)
{
    
    const double koff = .068;
    
    // TODO: 0 not handled
    double upd = get_rand_uniform(rng, -1, 1) * 2 * koff * pathway_[9] * dt;
    pathway_[9] += upd;
    
}


void ecoli_bacterium::update_bacterium(double dt, rng_type& rng, polygon_type& poly)
{
    update_pathway(dt);
    fluctuate_r(dt, rng);
    const double Hc = 10.3;
    const double Kc = 3.1;
    double tau = pow(pathway_[6], Hc)/(pow(pathway_[6], Hc) + pow(Kc, Hc));
    
    if (get_rand_uniform(rng, 0, 1) > tau){
        run(dt, poly);
    }else{
        tumble(dt, rng);
    }
    
    
}

void ecoli_bacterium::tumble(double dt, rng_type& rng)
{
    double degs =  get_rand_gamma(rng, 4, 18.32) - 4.6;
    direction_.rotate(degs);
}

void ecoli_bacterium::run(double dt, polygon_type& poly)
{
    auto t_pos = position_ + (direction_ * speed_ * dt);
    point2d r_pos, r_dir;
    move_and_reflect(position_, t_pos, poly, false, r_pos, r_dir);
    position_ = r_pos;
    direction_ = r_dir;
    
}

point2d ecoli_bacterium::get_position(){
    return position_;
}

void ecoli_bacterium::set_concentration(double conc)
{
    pathway_[8] = conc;
}