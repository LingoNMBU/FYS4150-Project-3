  
#ifndef __Penningtrap_hpp__  
#define __Penningtrap_hpp__

//Includes
#include <armadillo>
#include <vector>

//Penning trap class
class PenningTrap
{
  private:
    double iter_, d;

  public:
    //Declaration
    double B0, V0;
    std::vector<arma::vec> particles;

    //constructor
    PenningTrap(double B0_in, double V0_in, double d_in);

    // Add a particle to the trap #single #multiple
    void PenningTrap::add_particle(Particle p_in);
};

#endif  