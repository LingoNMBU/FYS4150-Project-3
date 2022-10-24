  
#ifndef __Particle_hpp__  
#define __Particle_hpp__

//includes
#include <armadillo>

//Particle class
class Particle
{
  private:
    friend class Penningtrap; // gives Penningtrap access to change variables in Particle    
    //Declaration    

  public:
    //Declaration
    arma::vec r,v;
    double q, m;
    //constructor
    Particle(arma::vec r_in, arma::vec v_in, double q_in, double m_in);
};
  

#endif  