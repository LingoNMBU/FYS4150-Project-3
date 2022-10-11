
#include "Particle.hpp"
#include "Penningtrap.hpp"
#include <cassert>

int main()
{


    //Test for creating a particle with a set of values
    arma::vec v0, r0;
    double q0, m0;    
    v0 = arma::vec(3).fill(0.);
    r0 = arma::vec(3).fill(0.);
    q0 = 1.;
    m0 = 1.;

    Particle particle0 = Particle(r0, v0, q0, m0);
    assert (arma::approx_equal(v0, particle0.v,"reldiff", 0.));
    assert (arma::approx_equal(r0, particle0.r,"reldiff", 0.));

    PenningTrap penningtrap0 = PenningTrap(1., 2., 3.);
    assert (penningtrap0.particles.size() == 0);
    penningtrap0.add_particle(particle0);
    assert (penningtrap0.particles.size() == 1);



    return 0;
}

//void test_particle_creation(arma::vec r0, arma::vec v0, double q0, double m0){
//    
//}