  
#include "particle.hpp"


//constructor
Particle::Particle(arma::vec r_in, arma::vec v_in, double q_in, double m_in)
{
    q = q_in; // e
    m = m_in; // u
    v = v_in; //micrometers/microsecond = m/s
    r = r_in; // micrometers
}

