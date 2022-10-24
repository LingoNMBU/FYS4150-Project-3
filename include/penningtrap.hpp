  
#ifndef __penningtrap_hpp__  
#define __penningtrap_hpp__

//Includes
#include <armadillo>
#include <vector>
#include <math.h>
#include "particle.hpp"
#include <cmath>
#include <complex>

//Penning trap class
class PenningTrap
{
  private:
    double iter_;

  public:
    //Declaration
    double B0, V0, w_v, f, d;
    std::vector<Particle> particles;
    bool var_field;

    //constructor
    PenningTrap(double B0_in, double V0_in, double d_in, bool var_field_in = false, double w_v_in = 0, double f_in = 0);

    // Add a particle to the trap #single #multiple
    void add_particle(Particle p_in);

    //Adding n randomly placed particle within 0.1d from the middle of the trap
    void random_particles(int n, double q, double m);

    // External electric field at point r=(x,y,z) #single #multiple
    arma::vec external_E_field(arma::vec r, double t = 0); 

    // External magnetic field at point r=(x,y,z) #single #multiple
    arma::vec external_B_field(arma::vec r); 

    // The total force on particle_i from the external fields #single
    arma::vec total_force_external(int i, double t = 0);

    // Force on particle_i from particle_j #multiple
    arma::vec force_particle(int i, int j);

    // The total force on particle_i from the other particles #multiple
    arma::vec total_force_particles(int i);

    // The total force on particle_i from both external fields and other particles #multiple
    arma::vec total_force(int i, bool particle_interactions, double t = 0);

    // Evolve the system one time step (dt) using Forward Euler #single #multiple #test
    void evolve_forward_Euler(double dt, bool particle_interactions, double t = 0);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order #single #multiple
    void evolve_RK4(double dt, bool particle_interactions, double t = 0);

    // Evolve the system one time step (dt) using Runge-Kutta 4th order #single #multiple
    void evolve_RK4_2(double dt, bool particle_interactions, double t = 0);

    // Analytical solution for a single particle in the penningtrap
    arma::vec single_particle_analytical(double t, double x_0, double v_0, double z_0);



};

#endif  