
#include "particle.hpp"
#include "penningtrap.hpp"

int main()
{

    double convergence_rate(std::vector<arma::mat> r_ex, std::vector<arma::mat> r_i, std::vector<int> ns);

    double q0, m0, d0, V0, B0; 
    q0 = 2.; // [e]
    m0 = 40.078; // [u] calcium +2
    d0 = 500.; // [micrometer]
    V0 = 2.41e6; // [(u*(micrometer)^2 )/ ((microsecond)^2 * elementary charge)]
    B0 = 96.5; // [u*(micrometer)^2 / (microsecond)*(elementary charge)]
    double t_fin = 50;
    int npoints = 30000;
    double dt = t_fin/npoints;

    //Problem 8
    arma::vec r1 = arma::vec(3);
    arma::vec v1 = arma::vec(3);
    r1 = {20., 0. , 20.};
    v1 = {0. , 25., 0.};

    arma::vec r2 = arma::vec(3);
    arma::vec v2 = arma::vec(3);
    r2 = {25., 25. , 0.};
    v2 = {0. , 40., 5.};

    Particle particle1 = Particle(r1, v1, q0, m0); 
    Particle particle2 = Particle(r2, v2, q0, m0); 

    std::vector<int> ns;
    ns = {4000, 8000, 16000, 32000};
    std::vector<arma::mat> rs_eu;
    std::vector<arma::mat> rs_ana;
    std::vector<arma::mat> rs_rk;

    arma::vec convergences(2);


    for( int k = 0; k < 4; k++)
    {

        
        rs_ana.push_back(arma::mat(3,ns[k]));
        rs_eu.push_back(arma::mat(3,ns[k]));
        rs_rk.push_back(arma::mat(3,ns[k]));

        arma::mat rs_ana2(3,ns[k]);
        arma::mat rs_eu2(3,ns[k]);
        arma::mat rs_rk2(3,ns[k]);

        PenningTrap penningtrap4 = PenningTrap(B0, V0, d0);
        penningtrap4.add_particle(particle1);

        PenningTrap penningtrap5 = PenningTrap(B0, V0, d0);
        penningtrap5.add_particle(particle1);

        double dt = 50./ns[k];

        //std::cout << std::endl;
        //std::cout << k;
        //std::cout << std::endl;
        //std::cout << ns[k];
        //std::cout << std::endl;
        //std::cout << dt;
        //std::cout << std::endl;
        

        for( int i= 0; i < ns[k] ; i++)
        {
            double curr_time = dt*(i);

            rs_ana2.col(i) = penningtrap3.single_particle_analytical(curr_time, 20., 20., 25.);

            rs_rk2.col(i) = penningtrap4.particles[0].r;
            penningtrap4.evolve_RK4_2(dt, false); 

            rs_eu2.col(i) = penningtrap5.particles[0].r;
            penningtrap5.evolve_forward_Euler(dt, false);  

            
        }

        //std::cout << std::endl;
        //std::cout << rs_ana2;
        //std::cout << std::endl;
            
        std::string filename_eu = "rs_eu2_n_" + std::to_string(ns[k]);
        std::string filename_rk = "rs_rk2_n_" + std::to_string(ns[k]);
        std::string filename_ana = "rs_ana2_n_" + std::to_string(ns[k]);

        rs_eu2.save(filename_eu, arma::csv_ascii);
        rs_rk2.save(filename_rk, arma::csv_ascii);
        rs_ana2.save(filename_ana, arma::csv_ascii);

        rs_eu[k] =rs_eu2;
        rs_rk[k] =rs_rk2;
        rs_ana[k] = rs_ana2;
        
//
        //rs_eu[k].save(filename_eu, arma::csv_ascii);
        //rs_rk[k].save(filename_rk, arma::csv_ascii);
        //rs_ana[k].save(filename_ana, arma::csv_ascii);


    }

    std::cout << std::endl;
    std::cout << "Forward Euler";
    std::cout << std::endl;

    convergences(0) = convergence_rate(rs_ana, rs_eu, ns);

    std::cout << std::endl;
    std::cout << "Runge Kutta" ;
    std::cout << std::endl;
    convergences(1) = convergence_rate(rs_ana, rs_rk, ns);


    std::cout << "Runge Kutta convergence  " << std::endl;
    std::cout << convergences(1) << std::endl;
    std::cout << "Euler convergence  " << std::endl;
    std::cout << convergences(0) << std::endl;
    std::cout << std::endl;


    convergences.save("convergence_rates", arma::csv_ascii);


    return 0;

}


double convergence_rate(std::vector<arma::mat> r_ex, std::vector<arma::mat> r_i, std::vector<int> ns)
{
    arma::vec dmax = {0.,0.,0.,0.};
    arma::vec h = {0.,0.,0.,0.};
    double convergence_rate = 0;
    int n = dmax.size();

    for (int k = 0; k < n;k++)
    {
        h(k) = 50./(ns[k] -1);
        for(int i = 0; i < ns[k]; i++)
        {
            arma::vec delta = r_ex[k].col(i) - r_i[k].col(i);
            arma::vec absdelta = abs(delta);
            double sumabsdelta = arma::accu(delta);

            if(dmax(k) < sumabsdelta)
            {
                dmax(k) = sumabsdelta;

            }
        }
    }

    for (int k = 1; k < n; k++)
    {
        double ddmax = dmax(k)/dmax(k-1);

        double ddh = h(k)/h(k-1);

        convergence_rate += (1./3.)* (log(ddmax)/(log(ddh)));
   
    }

    return convergence_rate;