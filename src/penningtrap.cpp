
#include "penningtrap.hpp"
#include "particle.hpp"
// Penning trap

PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, bool var_field_in, double w_v_in, double f_in)
{
  B0 = B0_in;               // Magnetic field strength
  V0 = V0_in;               // Electric Field strength
  d = d_in;                 // Characteristic dimension
  var_field = var_field_in; // Is there a varying electromagnetic field
  w_v = w_v_in;
  f = f_in;
}

// Add a particle to the trap #single #multiple
void PenningTrap::add_particle(Particle p_in)
{
  particles.push_back(p_in);
}

void PenningTrap::random_particles(int n, double q, double m)
{
  for (int i = 0; i < n; i++)
  {
    arma::arma_rng::set_seed_random();

    arma::vec r = arma::vec(3).randn() * 0.1 * d; // random initial position
    arma::vec v = arma::vec(3).randn() * 0.1 * d; // random initial velocity
    add_particle(Particle(r, v, q, m));
  }
}

// External electric field at point r=(x,y,z) #single #multiple
arma::vec PenningTrap::external_E_field(arma::vec r, double t)
{
  // V = (V0/(2*d*d))*( 2*r[2]*r[2] - r[0]*r[0] - r[1]*r[1])
  // E = -(V0/(2*d*d))*(4*r[2] - 2*r[0], -2*r[1])
  double V1;
  if (var_field == true)
  {
    V1 = V0 + V0*f*cos(w_v * t);

    //std::cout << t << std::endl;
  }
  else
  {
    V1 = V0;
  }


  double k = (V1 / (2. * d * d));
  arma::vec E = arma::vec(3); // -gradient of V

  E(2) = -k * 4. * r(2);
  E(0) = +k * 2. * r(0);
  E(1) = +k * 2. * r(1);

  // std::cout << __LINE__;
  // std::cout << E;
  // std::cout << std::endl;
  return E;
}

// External magnetic field at point r=(x,y,z) #single #multiple
arma::vec PenningTrap::external_B_field(arma::vec r)
{
  arma::vec B = arma::vec(3); 
  B(0) = 0.;
  B(1) = 0.;
  B(2) = B0;
  return B;
}

// The total force on particle_i from the external fields #single
arma::vec PenningTrap::total_force_external(int i, double t)
{
  arma::vec r = particles[i].r;
  arma::vec v = particles[i].v;
  double q = particles[i].q;
  arma::vec F_b = arma::vec(3).fill(0.);
  arma::vec F_e = arma::vec(3).fill(0.);

  // arma::vec F = q *external_E_field(r) + arma::cross(q*v, external_B_field(r));

  // b_x/q = v_y * b_z - v_z*b_y = v_y*b_0
  // b_y/q = v_z * b_x - v_x*b_z = -v_x*b_0
  // b_z/q = v_x * b_y - v_y*b_x = 0

  double distance = sqrt((pow((r(0)), 2) + pow((r(1)), 2) + pow((r(2)), 2)));


  if (distance < d)
  {
    F_b(0) = (q * v(1)* B0);
    F_b(1) = (-q * v(0) * B0);
    F_b(2) = 0.;

    for (int j = 0; j < 3; j++)
    {
      F_e(j) = q * external_E_field(r, t)(j);
    }
    //F_b = {(q * v(1)* B0), (-q * v(0) * B0), 0.};
    
  }


  return F_b + F_e;
}

// Force on particle_i from particle_j #multiple
arma::vec PenningTrap::force_particle(int i, int j)
{
  double k_coulomb = 1.38935333e5; // [u*(micrometer)^3 / (microsecond)^2*e^2

  //std::cout << "WHY THE FUCK IS THIS PRINTING" << std::endl;

  arma::vec r_i = particles[i].r;
  arma::vec r_j = particles[j].r;
  double q_i = particles[i].q;
  double q_j = particles[j].q;

  double distance3 = pow((pow((r_i(0) - r_j(0)), 2) + pow((r_i(1) - r_j(1)), 2) + pow((r_i(2) - r_j(2)), 2)), (3. / 2.));

  arma::vec F_ij = q_i * (k_coulomb * q_j * (r_i - r_j) / (distance3)); // F_ij = E_ij * q_i
  return F_ij;
}

// The total force on particle_i from the other particles #multiple
arma::vec PenningTrap::total_force_particles(int i)
{
  arma::vec F_on_particle = arma::vec(3).fill(0.);
  int n_particles = particles.size();

  for (int j = 0; j < n_particles; j++)
  {
    if (j != i)
    {
      F_on_particle += force_particle(i, j);
    }
  }
  return F_on_particle;
}

// The total force on particle_i from both external fields and other particles #multiple
arma::vec PenningTrap::total_force(int i, bool particle_interactions, double t)
{
  arma::vec total_force = arma::vec(3);
  arma::vec F_e = arma::vec(3);
  arma::vec F_p = arma::vec(3);

  F_e = total_force_external(i, t);
  F_p = total_force_particles(i);

  for (int j = 0; j < 3 ; j++)
  {
    if (particle_interactions == true)
    {
      total_force(j) = F_e(j) + F_p(j);
    }
    else
    {
      total_force(j) = F_e(j);
    }

  }
  //std::cout << std::endl;
  //std::cout << total_force << std::endl;


  return total_force;
}

// Evolve the system one time step (dt) using Forward Euler #single #multiple #test
void PenningTrap::evolve_forward_Euler(double dt, bool particle_interactions, double t)
{
  int n = particles.size();
  double h = dt;

  for (int i = 0; i < n; i++)
  {
    Particle particle_i = particles[i];
    arma::vec r0 = arma::vec(3);
    arma::vec v0 = arma::vec(3);
    arma::vec v1 = arma::vec(3);
    arma::vec r1 = arma::vec(3);
    double m = particle_i.m;

    for (int j = 0; j < 3; j++)
    {
      r0(j) = particle_i.r(j);
      v0(j) = particle_i.v(j);

      v1(j) = v0(j) + h * (1. / m) * total_force(i, particle_interactions, t)(j);
      r1(j) = r0(j) + h * v0(j);
    }

    //std::cout << "r0" << std::endl;
    //std::cout << r0 << std::endl;


    for (int j = 0; j < 3; j++)
    {
      particles[i].r(j) = r1(j);
      particles[i].v(j) = v1(j);
    }

    // vi+1 = vi + hfi
    // xi+1 = xi + hvi

    // std::cout << std::endl<< "k1 euler ";
    // std::cout << h*v0;
    // std::cout << std::endl << "New r euler ";
    // std::cout << particles[i].r;
    // std::cout << std::endl;
  }
}

void PenningTrap::evolve_RK4(double dt, bool particle_interactions, double t)
{
  int n = particles.size();
  double h = dt;
  double h_2 = h / 2.;

  //std::cout << std::endl<< "t ";
  //std::cout << t;
  //std::cout << std::endl << "h ";
  //std::cout << h;
  //std::cout << std::endl << "dt ";
  //std::cout << dt;
  //std::cout << std::endl;

  //std::cout << __LINE__ << std::endl;

  std::vector<Particle> original_particles = particles;
  std::vector<Particle> new_particles = particles;

  //std::cout << __LINE__ << std::endl;

  for (int i = 0; i < n; i++)
  {
    double m = particles[i].m;
    double inv_m = 1./m;

    std::cout << __LINE__ << std::endl;

    // finding k1
    arma::vec r0 = original_particles[i].r;
    arma::vec v0 = original_particles[i].v;

    arma::vec k1_v = h * (inv_m) * total_force(i, particle_interactions, t+h);
    arma::vec k1_r = h * v0 ;

    std::cout << __LINE__ << std::endl;

    // updating particle positions and velocities a half step with k1
    particles[i].v = v0 + k1_v / 2.;
    particles[i].r = r0 + k1_r / 2.;

    // finding k2
    arma::vec r05 = particles[i].r;
    arma::vec v05 = particles[i].v;

    arma::vec k2_v = h * (inv_m) * total_force(i, particle_interactions, t+h_2);
    arma::vec k2_r = h * v05;

    // updating particle positions and velocities at the same point with k2
    particles[i].v = v0 + k2_v * h_2;
    particles[i].r = r0 + k2_r * h_2;

    // finding k3
    arma::vec r05_2 = particles[i].r;
    arma::vec v05_2 = particles[i].v;

    arma::vec k3_v = h * (inv_m) * total_force(i, particle_interactions, (t+h_2));
    arma::vec k3_r = h * v05;

    // updating particle positions and velocities at t0 + h with k3
    particles[i].v = v0 + k3_v * h;
    particles[i].r = r0 + k3_r * h;

    // finding k3
    arma::vec r1 = particles[i].r;
    arma::vec v1 = particles[i].v;

    arma::vec k4_v = h * (inv_m) * total_force(i, particle_interactions, t+h);
    arma::vec k4_r = h * v0;

    // updating particle positions and velocities at t0 + h with a kombination of the earlier slopes
    new_particles[i].v = v0 + ((1 / 6.) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v));
    new_particles[i].r = r0 + ((1 / 6.) * (k1_r + 2 * k2_r + 2 * k3_r + k4_r));

    particles[i] = original_particles[i]; 
  }
  particles = new_particles;
}

void PenningTrap::evolve_RK4_2(double dt, bool particle_interactions, double t)
{
  int n = particles.size();
  double h = dt;
  double h_2 = h / 2.;
  double m = particles[0].m;
  double inv_m = 1./m;

  //std::cout << std::endl<< "t ";
  //std::cout << t;
  //std::cout << std::endl << "h ";
  //std::cout << h;
  //std::cout << std::endl << "dt ";
  //std::cout << dt;
  //std::cout << std::endl;

  std::vector<Particle> original_particles = particles;


  arma::mat k1_r = arma::mat(3,n);
  arma::mat k1_v = arma::mat(3,n);
  for (int i = 0; i < n; i++)
  {// finding k1s
    k1_v.col(i) = h * (inv_m) * total_force(i, particle_interactions, t+h);
    k1_r.col(i) = h * particles[i].v;
  }

  //std::cout << std::endl<< "t ";
  //std::cout << k1_r.col(0);
  //std::cout << std::endl << "h ";
  //std::cout << k1_v.col(0);
  //std::cout << std::endl << "dt ";
  //std::cout << dt;
  //std::cout << std::endl;

  for (int i = 0; i < n; i++)
  {// updating particle positions and velocities a half step with k1
    particles[i].v = original_particles[i].v + k1_v.col(i) / 2.;
    particles[i].r = original_particles[i].r + k1_r.col(i) / 2.;
  }

  arma::mat k2_r = arma::mat(3,n).fill(0.);
  arma::mat k2_v = arma::mat(3,n).fill(0.);
  for (int i = 0; i < n; i++)
  {// finding k2s
      k2_v.col(i) = h * (inv_m) * total_force(i, particle_interactions, t+h_2);
      k2_r.col(i) = h *  particles[i].v;
  }

  for (int i = 0; i < n; i++)
  {// updating particle positions and velocities a half step from start with k2
    particles[i].v = original_particles[i].v + k2_v.col(i) / 2.;
    particles[i].r = original_particles[i].r + k2_r.col(i) / 2.;
  }

  arma::mat k3_r = arma::mat(3,n);
  arma::mat k3_v = arma::mat(3,n); 
  for (int i = 0; i < n; i++)
  {// finding k3s
    k3_v.col(i) = h * (inv_m) * total_force(i, particle_interactions, t+h_2);
    k3_r.col(i) = h * particles[i].v;

  }

  for (int i = 0; i < n; i++)
  {// updating particle positions and velocities a step from start with k3
    particles[i].v = original_particles[i].v + k3_v.col(i);
    particles[i].r = original_particles[i].r + k3_r.col(i);
  }

  arma::mat k4_r = arma::mat(3,n);
  arma::mat k4_v = arma::mat(3,n); 
  for (int i = 0; i < n; i++)
  {// finding k4s
    k4_v.col(i) = h * (inv_m) * total_force(i, particle_interactions, t+h);
    k4_r.col(i) = h * particles[i].v;

  }

  for (int i = 0; i < n; i++)
  {
    particles[i].v = original_particles[i].v + ((1. / 6.) * (k1_v.col(i) + 2. * k2_v.col(i) + 2. * k3_v.col(i) + k4_v.col(i)));
    particles[i].r = original_particles[i].r + ((1. / 6.) * (k1_r.col(i) + 2. * k2_r.col(i) + 2. * k3_r.col(i) + k4_r.col(i)));
  }

}



arma::vec PenningTrap::single_particle_analytical(double t, double x_0, double z_0, double v_0)
{
  double m = particles[0].m;
  double q = particles[0].q;

  arma::vec r1 = arma::vec(3).fill(0.);

  double w_z = std::sqrt((2 * q * V0) / (m * d * d));
  double w_0 = (q * B0) / m;
  //std::cout << w_z << std::endl;
  //std::cout << w_z << std::endl;
  //std::cout << w_z << std::endl;
  double w_p = (w_0 + sqrt((w_0 * w_0) - (2 * w_z * w_z))) / (2);
  double w_m = (w_0 - sqrt((w_0 * w_0) - (2 * w_z * w_z))) / (2);

  double A_p = (v_0 + (w_m * x_0)) / (w_m - w_p);
  double A_m = -((v_0 + (w_p * x_0)) / (w_m - w_p));

  double phi_p = 0;
  double phi_m = 0;

  // Decomposition using eulers rule

  //r*e(iphi) = r*cos(x) + r*isin(x)

  //f(t) = A_p*e^(-i*(w_p*t + phi_p)) + A_m*e^(-i(w_m*t + phi_m))

  //f(t) = A_p*cos(-wp*t - phi_p)+ A_p*i*sin( -wp*t - phi_p) + A_m*cos(-wm*t - phi_m)+ A_m*i*sin( -wm*t - phi_m) 

  //Re(f(t)) = A_p*cos(-wp*t - phi_p) + A_m*cos(-wm*t - phi_m)

  //Im(f(t)) = A_p*sin(-wp*t - phi_p) + A_m*sin(-wm*t - phi_m) 

  double f1_x = A_p * cos((-w_p * t) - phi_p);
  double f2_x = A_m * cos((-w_m * t) - phi_m);

  double f1_y = A_p * sin((-w_p * t) - phi_p);
  double f2_y = A_m * sin((-w_m * t) - phi_m);

  r1(0) = f1_x + f2_x;
  r1(1) = f1_y + f2_y;
  r1(2) = z_0 * cos(t * w_z);

  return r1;
}

// will make two coupled equations to solve for the velocity and positon vectors of the particle by representing v as the derivative of r in the same equaiton

// end up with equations on the form
// dr/dt = v(t)
//  dv/dt = 1/m F(t,x,v)

// y_i+1 = yi + h*fi, fi = dy/dt, trunc error Oh^2, global error Oh

// single step only need one evaluation of y, yi

// Euler for coupled:

// vi+1 = vi + hfi
// xi+1 = xi + hvi

// Euler cromer for coupled

// vi+1 = vi + hfi
// xi+1 = xi + hvi+1

// Evolve the system one time step (dt) using Runge-Kutta 4th order #single #multiple
// void evolve_RK4(double dt);
// global error Oh^4 or Oh^n for n-RK method, 4 is a usual balanced order

// solves dy/dt = f(t,y)

// Improved version of Predictor Corrector, y_i+1 = yi + h* avg(fi+1 and fi) where fi+1 is found from euler

// n-RK uses n forward estimates to estimate the slope for forward jumping

// Ago: 5 steps, makes 4 slope shifts as k-values

// Step1:
// k1 = hf(ti,yi)

// Step2:
// k2 = hf(t_i + 1/2 * h, ri + k1)

// Step3:
// k3 = hf(t_i + 1/2 * h, ri + 1/2 * k2)

// Step4:
// k4 = hf(t_i + h, ri + k3)

// Step5:
// ri+1 = ri + 1/6(k1 + 2*k2 + 2*k3 + k4)  //In a sense an integration of a timestep using simpsons rule

// For coupled equations do every step in sync, do every step for velocity and position
// in practice 6 different k values per step, three dimensions for positions and velocities
