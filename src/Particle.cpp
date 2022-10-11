  
#include "Particle.hpp"


//constructor
Particle::Particle(arma::vec r_in, arma::vec v_in, double q_in, double m_in)
{
    q_ = q_in;
    m_ = m_in;
    v = v_in;
    r = r_in;
}