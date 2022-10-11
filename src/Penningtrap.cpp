  
#include "Penningtrap.hpp"
#include "Particle.hpp"
//Penning trap 

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  B0 = B0_in; //Magnetic field strength
  V0 = V0_in; // Electric Field strength
  d = d_in; //Characteristic dimension
}



// Add a particle to the trap #single #multiple
void PenningTrap::add_particle(Particle p_in)
{
  particles.push_back(p_in);
}

//// External electric field at point r=(x,y,z) #single #multiple
//arma::vec external_E_field(arma::vec r);  
//
//// External magnetic field at point r=(x,y,z) #single #multiple
//arma::vec external_B_field(arma::vec r);  
//
//// Force on particle_i from particle_j #multiple
//arma::vec force_particle(int i, int j);
//
//// The total force on particle_i from the external fields #single
//arma::vec total_force_external(int i);
//
//// The total force on particle_i from the other particles #multiple
//arma::vec total_force_particles(int i);
//
//// The total force on particle_i from both external fields and other particles #multiple
//arma::vec total_force(int i);

// Evolve the system one time step (dt) using Runge-Kutta 4th order #single #multiple
//void evolve_RK4(double dt);
  // global error Oh^4 or Oh^n for n-RK method, 4 is a usual balanced order

  //solves dy/dt = f(t,y)

  //Improved version of Predictor Corrector, y_i+1 = yi + h* avg(fi+1 and fi) where fi+1 is found from euler

  //n-RK uses n forward estimates to estimate the slope for forward jumping

  //Ago: 5 steps, makes 4 slope shifts as k-values

  //Step1: 
  //k1 = hf(ti,yi)

  //Step2:
  //k2 = hf(t_i + 1/2 * h, ri + k1)

  //Step3:
  //k3 = hf(t_i + 1/2 * h, ri + 1/2 * k2)

  //Step4:
  //k4 = hf(t_i + h, ri + k3)

  //Step5:
  //ri+1 = ri + 1/6(k1 + 2*k2 + 2*k3 + k4)  //In a sense an integration of a timestep using simpsons rule

  //For coupled equations do every step in sync, do every step for velocity and position
  //in practice 6 different k values per step, three dimensions for positions and velocities




// Evolve the system one time step (dt) using Forward Euler #single #multiple #test
//void evolve_forward_Euler(double dt);

// will make two coupled equations to solve for the velocity and positon vectors of the particle by representing v as the derivative of r in the same equaiton

//end up with equations on the form
//dr/dt = v(t)

// dv/dt = 1/m F(t,x,v)

// y_i+1 = yi + h*fi, fi = dy/dt, trunc error Oh^2, global error Oh

//single step only need one evaluation of y, yi

//Euler for coupled:

//vi+1 = vi + hfi
//xi+1 = xi + hvi

//Euler cromer for coupled

//vi+1 = vi + hfi
//xi+1 = xi + hvi+1

