#include<iostream>
#include<stdlib.h>
#include<limits.h>

#include "EVRP.hpp"
#include "ACO.hpp"
#include "stats.hpp"

using namespace std;

/****************************************************************/
/*                Main Function                                 */
/****************************************************************/
int main(int argc, char *argv[]) {

  int run, i;

  problem_instance = argv[1];
  change_degree = atof(argv[2]);
  change_speed = atoi(argv[3]);
  total_changes = atoi(argv[4]);
  alg_mode = atoi(argv[5]);

  read_problem(problem_instance);   //Read TSP from file

  change_speed = problem_size * 25;
  total_changes = total_changes + 1;
  max_iterations = TERMINATION;
  //max_iterations = 10;
  //cout << "termination " << TERMINATION << endl;
  max_trials = 10;
  set_algorithm_parameters();  //from ACO.h
  
  open_stats(); //open text files to store statistics stats.h
  
  //compute distances and nearest neigbours
  compute_distances();
  compute_nn_lists();         

  for(run = 1; run <= max_trials; run++) {
      //reset parameters for statistics
      initialize_environment();  //from EVRP.h
      //ACO algorithms initialization
      init_try(run);             //from ACO.h

      while(current_iteration < max_iterations) {
        //ACO algorithm iterative methods
        ACO();

        get_observation(current_iteration-1);           //from stat.hpp

        if(env_index == 0 && (current_iteration%(problem_size*100) == 0))
          change_environment();
        else 
          if(env_index > 0 && (current_iteration%change_speed==0)) 
            change_environment(); //from DTSPwc.hpp
      }

    //cout << best_so_far_ant->tour_length << " " << best_so_far_ant->no_recharges << " " << best_so_far_ant->no_customers << " " <<  endl;
    //offline_performance = get_current_best();
    //cout<< "Offline performance: " << get_offline_performance() << endl;      //from stats.h
    //cout << "Best found: " << get_current_best() << endl;
    //cout<< "Offline error: " << get_offline_error() << endl;                  //from stats.h
    //cout<< "Time best found: " << get_best_time() << " Iteration best found: " << get_best_iteration() << endl;
    //cout <<"Total time elapsed: " << get_time() << endl;
    //cout <<"Feasibility: " << counter << endl;            
    //store totals for each run
    get_mean(run-1,get_offline_performance()); //from stats.h
  }

  close_stats(); //close text files with stored statistics from stats.h
  //free memory
  free_ACO();
  free_EVRP();

  return 0;
}
