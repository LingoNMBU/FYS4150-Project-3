
#include "particle.hpp"
#include "penningtrap.hpp"
#include <cassert>
#include <iostream>

int main()
{


    //Test for creating a particle with a set of values
    arma::vec v0, r0, v1, r1;
    double q0, m0, d0, V0, B0;    
    v0 = arma::vec(3).fill(0.);
    r0 = arma::vec(3).fill(0.);
    v1 = arma::vec(3).fill(1.);
    r1 = arma::vec(3).fill(1.);
    q0 = 2.; // [e]
    m0 = 40.078; // [u] calcium +2
    d0 = 500.; // [micrometer]
    V0 = 2.41e6; // [(u*(micrometer)^2 )/ ((microsecond)^2 * elementary charge)]
    B0 = 96.5; // [u*(micrometer)^2 / (microsecond)*(elementary charge)]

    Particle particle0 = Particle(r0, v0, q0, m0); // r= 0, 0, 0  v= 0, 0, 0
    Particle particle1 = Particle(r1, v1, q0, m0); // r= 1, 1, 1  v= 1,1,1
    Particle particle2 = Particle(r1, v0, q0, m0); // r= 1, 1, 1  v= 0, 0 ,0
    Particle particle3 = Particle(r0, v1, q0, m0); // r= 0, 0, 0  v= 1, 1, 1

    
    

    //Test for particle creation
    assert (arma::approx_equal(v0, particle0.v,"reldiff", 0.));
    assert (arma::approx_equal(r0, particle0.r,"reldiff", 0.));

    //Test for adding particles to a Penning trap
    PenningTrap penningtrap0 = PenningTrap(B0, V0, d0);
    assert (penningtrap0.particles.size() == 0);
    penningtrap0.add_particle(particle0);
    assert (penningtrap0.particles.size() == 1);

    //Testing magnetic and electric field at 0,0,0
    arma::vec E0 = arma::vec(3).fill(0.);
    arma::vec B000 = penningtrap0.external_E_field(r0);
    arma::vec E000 = penningtrap0.external_B_field(r0);
    //assert (arma::approx_equal(B0, B000(2),"reldiff", 0.));
    //assert (arma::approx_equal(E0, E000,"reldiff", 0.));
    //std::cout << B000;
    //std::cout << E000;

    //Testing particle forces;
    PenningTrap penningtrap1 = PenningTrap(B0, V0, d0);
    penningtrap1.add_particle(particle0);
    penningtrap1.add_particle(particle1);
    //std::cout << penningtrap1.particles.size();
    arma::vec force_ij = penningtrap1.force_particle(0,1);

    //std::cout << std::endl;
    //std::cout << force_ij;
    //std::cout << std::endl;


    //Test field forces    
    PenningTrap penningtrap2 = PenningTrap(B0, V0, d0);
    penningtrap2.add_particle(particle0);
    penningtrap2.add_particle(particle3);
    arma::vec F20 = penningtrap2.total_force_external(0);
    arma::vec F23 = penningtrap2.total_force_external(1);
    assert (arma::approx_equal(F20, r0, "reldiff", 0.));
    //std::cout << F23;


    //test forward euler, RK4 and analytical solution for 1 particle
    double t_fin = 50;
    double t_start = 0.;
    int npoints = 50000;
    double dt = t_fin/npoints;
    arma::vec r2 = arma::vec(3);
    arma::vec v2 = arma::vec(3);
    r2 = {20., 0. , 20.};
    v2 = {0. , 25., 0.};
    Particle particle4 = Particle(r2, v2, q0, m0); // r= 0.1, 0, 0.1      v= 0, 0.1 ,0  , for single particle analytical solution

    arma::mat positions_ana(3,npoints);
    arma::mat positions_eu(3,npoints);
    arma::mat positions_rk(3,npoints);

    PenningTrap penningtrap22 = PenningTrap(B0, V0, d0, false);
    penningtrap22.add_particle(particle4);

    PenningTrap penningtrap3 = PenningTrap(B0, V0, d0, false);
    penningtrap3.add_particle(particle4);

    PenningTrap penningtrap4 = PenningTrap(B0, V0, d0, false);
    penningtrap4.add_particle(particle4);

    for( int i= 0; i < npoints ; i++)

    {   double curr_time = t_start + dt*double(i);

        //velocities.push_back(penningtrap3.particles[0].v);
        //positions.push_back(penningtrap3.particles[0].r);        
        //positions_ana.push_back(penningtrap3.single_particle_analytical(curr_time, 20., 20., 25.));

        positions_ana.col(i) = penningtrap22.single_particle_analytical(curr_time, 20., 20., 25.);
        positions_eu.col(i) = penningtrap3.particles[0].r;
        positions_rk.col(i) = penningtrap4.particles[0].r;

        //positions_ana.col(i) = penningtrap3.single_particle_analytical(curr_time, 20., 20., 25.);
        //positions_eu.col(i) = penningtrap3.total_force(0,false);
        //positions_rk.col(i) = penningtrap4.total_force(0,false);

        penningtrap3.evolve_forward_Euler(dt, false, curr_time);
        penningtrap4.evolve_RK4_2(dt, false, curr_time);
        
    }  


    positions_eu.save("positions_euler", arma::csv_ascii);
    positions_ana.save("positions_ana", arma::csv_ascii);
    positions_rk.save("positions_rk4", arma::csv_ascii);


    //test cloning particle
    PenningTrap penningtrap5 = PenningTrap(B0, V0, d0);
    penningtrap5.add_particle(particle0);
    std::vector<Particle> particles_test = penningtrap5.particles;
    particles_test[0].m = 1;
    assert(penningtrap5.particles[0].m != particles_test[0].m);

    //std::cout << std::endl;
    //std::cout << penningtrap4.particles[0].m;
    //std::cout << std::endl;
    //std::cout << particles_test[0].m;
    //std::cout << std::endl;

    //std::cout << std::endl;
    //for (arma::vec i: velocities)
    //    std::cout << i << std::endl;
    //    std::cout << std::endl;
    //std::cout << std::endl;
    //for (arma::vec i: positions)
    //    std::cout << i << std::endl;
    //    std::cout << std::endl;
    //std::cout << std::endl;
    //std::cout << std::endl;
    //for (arma::vec i: positions_ana)
    //    std::cout << i << std::endl;
    //    std::cout << std::endl;

    //test distance

    





    return 0;
}

//void test_particle_creation(arma::vec r0, arma::vec v0, double q0, double m0){
//    
//}