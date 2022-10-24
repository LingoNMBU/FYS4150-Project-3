    
#include "particle.hpp"
#include "penningtrap.hpp"

int main()
{

    int particles_inside(std::vector<Particle> particles, double d);

    ////Problem 9a
    double q0, m0, d0, V0, B0; 
    q0 = 2.; // [e]
    m0 = 40.078; // [u] calcium +2
    d0 = 500.; // [micrometer]
    V0 = 2.41e6; // [(u*(micrometer)^2 )/ ((microsecond)^2 * elementary charge)]
    B0 = 96.5; // [u*(micrometer)^2 / (microsecond)*(elementary charge)]

    bool prob9a = false;
    if(prob9a == true)
    {
        arma::vec fs = arma::vec({0.1, 0.4, 0.7});
        double t_fin = 500.;
        double t_start = 0.;
        double w_v_min = 0.2; // MHz or 1/microsecond
        double w_v_max = 2.5; // MHz
        double w_v_step  = 0.02; // MHz
        int n_w_steps = int(round((w_v_max-w_v_min)/w_v_step));
        int n_t_steps = 30000;
        int n_f_steps = fs.size();
        double h = (t_fin-t_start)/n_t_steps; //step size

        arma::mat trapped = arma::mat(n_f_steps,n_w_steps);

        for(int i=0; i <  n_f_steps; i++)
        {// loop through amplitudes f    
            double f = fs(i);

            std::cout << "f " << std::endl;
            std::cout << f << std::endl;
            std::cout << std::endl;

            for (int j=0; j <  n_w_steps; j++)
            {//loop through frequencies w_v
                double w_v = w_v_min + j*w_v_step;

                //Initialization of Penning trap
                PenningTrap harmonic_trap = PenningTrap(B0, V0, d0, true, w_v, f);
                //PenningTrap harmonic_trap = PenningTrap(B0, V0, d0);
                harmonic_trap.random_particles(50, q0, m0);
                int remaining_particles_b = particles_inside(harmonic_trap.particles, harmonic_trap.d);
                std::cout << "remaining_particles before " << std::endl;
                std::cout << remaining_particles_b << std::endl;
                std::cout << std::endl;

                //bool any_remaining = true;

                for (int k=0; k < n_t_steps; k++)
                {//loop through timesteps t
                    //if (any_remaining == true)
                    //{
                    double t = t_start + h*k;
                    harmonic_trap.evolve_RK4_2(h, false, t);
                        //int remaining_particles_now = particles_inside(harmonic_trap.particles, harmonic_trap.d);       
                    
                        //if (remaining_particles_now == 0)
                        //{
                            //any_remaining = false;
                            ////int remaining_particles_m = particles_inside(harmonic_trap.particles, harmonic_trap.d);
//
                            //std::cout << "no remaining particles , time " << t << std::endl;
                            //std::cout << "k " << k << std::endl;
                            //std::cout << "h " << h << std::endl;
                            //std::cout << "h*t " << h*k << std::endl;
                            //std::cout << std::endl;
                            
                        //}
                    //}
                    

                }

                int remaining_particles = particles_inside(harmonic_trap.particles, harmonic_trap.d);
                std::cout << "remaining_particles after " << std::endl;
                std::cout << remaining_particles << std::endl;
                std::cout << std::endl;

                trapped(i,j) = remaining_particles;

            std::cout << "w_v " << std::endl;
            std::cout << w_v << std::endl;
            std::cout << std::endl;
            }

        trapped.save("trapped3", arma::csv_ascii);
        }
    }



    //problem 9b


    bool prob9b = true;
    if (prob9b==true)
    {
        double f = 0.7;
        double period = 500.; // micro second
        double w_v_min = 0.8; // MHz or 1/microsecond
        double w_v_max = 0.9; // MHz
        int n_w_v = 100;
        double w_v_step  = (w_v_max - w_v_min)/(n_w_v-1); // MHz
        int n_w_steps = int(round((w_v_max-w_v_min)/w_v_step));
        int n_t_steps = 30000;
        double h = period/n_t_steps; //step size
        std::vector<bool> interaction = {false, true};

        arma::mat trapped = arma::mat(2,n_w_steps);

        for(int i=0; i <  2; i++)
        {// loop through amplitudes f    
            bool inter = interaction[i];

            std::cout << "f " << std::endl;
            std::cout << f << std::endl;
            std::cout << std::endl;
            std::cout << "particle interactions" << std::endl;
            std::cout << inter << std::endl;
            std::cout << std::endl;

            for (int j=0; j <  n_w_steps; j++)
            {//loop through frequencies w_v
                double w_v = w_v_min + j*w_v_step;

                //Initialization of Penning trap
                PenningTrap harmonic_trap = PenningTrap(B0, V0, d0, true, w_v, f);
                //PenningTrap harmonic_trap = PenningTrap(B0, V0, d0);
                harmonic_trap.random_particles(50, q0, m0);
                

                for (int k=0; k < n_t_steps; k++)
                {//loop through timesteps t
                    double t = h*k;
                    harmonic_trap.evolve_RK4_2(h, inter, t);

                    
                }

                int remaining_particles = particles_inside(harmonic_trap.particles, harmonic_trap.d);
                std::cout << "remaining_particles " << std::endl;
                std::cout << remaining_particles << std::endl;
                std::cout << std::endl;

                trapped(i,j) = remaining_particles;

            std::cout << "w_v " << std::endl;
            std::cout << w_v << std::endl;
            std::cout << std::endl;
            }


        }
        trapped.save("trapped_particles_fine2", arma::csv_ascii);
    }

}

int particles_inside(std::vector<Particle> particles, double d)
{
    int n = 0;
    int n_particles = particles.size();
    for (int i = 0; i< n_particles; i++)
    {
        arma::vec r = particles[i].r;

         double distance = pow((pow((r(0)), 2) + pow((r(1)), 2) + pow((r(2)), 2)), (1./ 2.));

        if (distance < d)
        {
            n+=1;
        }
    }
    return n;
}