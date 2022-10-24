
#include "particle.hpp"
#include "penningtrap.hpp"

int main()
{
    int particles_inside(std::vector<Particle> particles, double d);
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

    //Penning trap with one particle
    PenningTrap penningtrap1 = PenningTrap(B0, V0, d0);
    penningtrap1.add_particle(particle1);


    //Initialization
    arma::mat r_rk_s(3, npoints);

    for( int i= 0; i < npoints ; i++)
    {  
        r_rk_s.col(i) = penningtrap1.particles[0].r;
        penningtrap1.evolve_RK4_2(dt, false);        
    } 
    r_rk_s.save("positions_single_rk4", arma::csv_ascii);


    //Penning trap with two particles, with and without interactions
    PenningTrap penningtrap2 = PenningTrap(B0, V0, d0);
    penningtrap2.add_particle(particle1);
    penningtrap2.add_particle(particle2);

    PenningTrap penningtrap3 = PenningTrap(B0, V0, d0);
    penningtrap3.add_particle(particle1);
    penningtrap3.add_particle(particle2);
    
    arma::mat r_rk_d1(3, npoints);
    arma::mat r_rk_d2(3, npoints);
    arma::mat v_rk_d1(3, npoints);
    arma::mat v_rk_d2(3, npoints);

    arma::mat r_rk_d1_int(3, npoints);
    arma::mat r_rk_d2_int(3, npoints);
    arma::mat v_rk_d1_int(3, npoints);
    arma::mat v_rk_d2_int(3, npoints);

    for( int i= 0; i < npoints ; i++)
    {  
        r_rk_d1.col(i) = penningtrap2.particles[0].r;
        r_rk_d2.col(i) = penningtrap2.particles[1].r;
        v_rk_d1.col(i) = penningtrap2.particles[0].v;
        v_rk_d2.col(i) = penningtrap2.particles[1].v;
        penningtrap2.evolve_RK4_2(dt, false); 

        r_rk_d1_int.col(i) = penningtrap3.particles[0].r;
        r_rk_d2_int.col(i) = penningtrap3.particles[1].r;
        v_rk_d1_int.col(i) = penningtrap3.particles[0].v;
        v_rk_d2_int.col(i) = penningtrap3.particles[1].v;
        penningtrap3.evolve_RK4_2(dt, true);       
    }  

    r_rk_d1.save("r_double_rk4_1", arma::csv_ascii);
    r_rk_d2.save("r_double_rk4_2", arma::csv_ascii);
    v_rk_d1.save("v_double_rk4_1", arma::csv_ascii);
    v_rk_d2.save("v_double_rk4_2", arma::csv_ascii);

    r_rk_d1_int.save("r_double_rk4_1_int", arma::csv_ascii);
    r_rk_d2_int.save("r_double_rk4_2_int", arma::csv_ascii);
    v_rk_d1_int.save("v_double_rk4_1_int", arma::csv_ascii);
    v_rk_d2_int.save("v_double_rk4_2_int", arma::csv_ascii);


    //Penning trap single particles, different ns
    

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


    //Problem 9
    arma::vec fs = arma::vec({0.1, 0.4, 0.7});
    double period = 500; // micro second
    double w_v_min = 0.2; // MHz
    double w_v_max = 2.5; // MHz
    double w_v_step  = 0.1; // MHz
    int n_w_steps = int(round((w_v_max-w_v_min)/w_v_step));
    int n_t_steps = 10000;
    int n_f_steps = fs.size();
    double h = period/n_t_steps; //step size
    
    arma::mat trapped = arma::mat(n_f_steps,n_w_steps);

//    for(int i=0; i <  n_f_steps; i++)
//    {// loop through amplitudes f    
//        double f = fs(i);
//
//        std::cout << "f " << std::endl;
//        std::cout << f << std::endl;
//        std::cout << std::endl;
//
//        for (int j=0; j <  n_w_steps; j++)
//        {//loop through frequencies w_v
//            double w_v = w_v_min + j*w_v_step;
//
//            //Initialization of Penning trap
//            PenningTrap harmonic_trap = PenningTrap(B0, V0, d0, true, w_v, f);
//            harmonic_trap.random_particles(100, q0, m0);
//
//            
//
//            for (int k=0; k < n_t_steps; k++)
//            {//loop through timesteps t
//                double t = h*k;
//                harmonic_trap.evolve_RK4(h, false, t);
//
//                
//            }
//
//            int remaining_particles = particles_inside(harmonic_trap.particles, harmonic_trap.d);
//            std::cout << "remaining_particles " << std::endl;
//            std::cout << remaining_particles << std::endl;
//            std::cout << std::endl;
//
//            trapped(i,j) = remaining_particles;
//
//        std::cout << "w_v " << std::endl;
//        std::cout << w_v << std::endl;
//        std::cout << std::endl;
//        }
//    }
//
    //trapped.save("trapped_particles", arma::csv_ascii);
    
   //Test
   //double w_v_test = 0.2;
   //double f_test = 0.4;
   //int n_test = 100;
   //PenningTrap test_trap = PenningTrap(B0, V0, d0, true, w_v_test, f_test);
   ////PenningTrap test_trap = PenningTrap(B0, V0, d0);

   ////for (int i = 0; i < n_test; i++)
   ////{
   ////    test_trap.add_particle(particle1);
   ////}

   //test_trap.random_particles(n_test, q0, m0);



  // for (int i = 0; i < n_test ; i++)
  // {
  //     std::cout << "test r " << i << std::endl;
  //     std::cout << test_trap.particles[i].r << std::endl;
  //     std::cout << std::endl;
  // }

  // double test_h = 0.05;
  // int remaining_particles = particles_inside(test_trap.particles, test_trap.d);

  // std::cout << "remaining_particles " << std::endl;
  // std::cout << remaining_particles << std::endl;
  // std::cout << std::endl;

  // for(int i = 0; i < 10000; i++)
  // {
  //      double test_t = i*test_h;
  //      test_trap.evolve_forward_Euler(test_h, false, test_t);
  //      remaining_particles = particles_inside(test_trap.particles, test_trap.d);

  // }
  //  std::cout << "remaining_particles " << std::endl;
  //  std::cout << remaining_particles << std::endl;
  //  std::cout << std::endl;

   //for (int i = 0; i < n_test ; i++)
   //{
   //    std::cout << "test r " << i << std::endl;
   //    std::cout << test_trap.particles[i].r << std::endl;
   //    std::cout << std::endl;
   //}
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
                //std::cout << std::endl;
                //std::cout << absdelta;
                //std::cout << std::endl;

            }
        }
    }
    //std::cout << "dmax" << std::endl;
    //std::cout << dmax << std::endl;
    //std::cout << std::endl;

    //std::cout << "h" << std::endl;
    //std::cout << h << std::endl;
    //std::cout << std::endl;

    for (int k = 1; k < n; k++)
    {
        double ddmax = dmax(k)/dmax(k-1);

        double ddh = h(k)/h(k-1);

        convergence_rate += (1./3.)* (log(ddmax)/(log(ddh)));

        //std::cout << "k" << std::endl;
        //std::cout << k ;
        //std::cout << std::endl;

        //std::cout << "ddmax " << std::endl;
        //std::cout << ddmax ;
        //std::cout << std::endl;

        //std::cout << "ddh " << std::endl;
        //std::cout << ddh ;
        //std::cout << std::endl;        
    }

    return convergence_rate;

}

int particles_inside(std::vector<Particle> particles, double d)
{
    int n = 0;
    int n_particles = particles.size();
    for (int i = 0; i< n_particles; i++)
    {
        arma::vec r = particles[i].r;

         double distance = pow((pow((r(0)), 2) + pow((r(1)), 2) + pow((r(2)), 2)), (1./ 2.));


        //std::cout << "r " << std::endl;
        //std::cout << r << std::endl;
        //std::cout << std::endl;

        //std::cout << "distance " << std::endl;
        //std::cout << distance << std::endl;
        //std::cout << std::endl;

        if (distance < d)
        {
            n+=1;
        }
    }
    return n;
}
