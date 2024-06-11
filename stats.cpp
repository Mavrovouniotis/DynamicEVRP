/*******************************************************************/
/*                                                                 */
/* This is an implementation of the dynamic perfomance measurements*/
/*                                                                 */
/*                                                                 */
/*******************************************************************/

#include<cmath>
#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<cstring>
#include<math.h>
#include<fstream>
#include<time.h>
#include<limits.h>

#include "ACO.hpp"
#include "EVRP.hpp"
#include "stats.hpp"

using namespace std;

//Used to output offline performance and population diversity

FILE *log_performance;
FILE *log_diversity;

//output files
char *perf_filename;
char *div_filename;

double* perf_of_trials;
double* perf_of_iterations;

double* div_of_trials;
double* div_of_iterations;


int max_iterations;
int max_trials;

bool flag = false;
double help = 0.0;

void open_stats(void){
    int i;
    //initialize
    perf_of_trials = new double[max_trials];
    perf_of_iterations = new double[max_iterations];
    div_of_trials = new double[max_trials];
    div_of_iterations = new double[max_iterations];
    for(i = 0; i < max_iterations; i++) {
        perf_of_iterations[i] = 0.0;
        div_of_iterations[i] = 0.0;
   }

   for(i = 0; i < max_trials; i++){
        perf_of_trials[i] = 0.0;
        div_of_trials[i] = 0.0;
   }

  //initialize and open output files
  perf_filename = new char[CHAR_LEN];
  sprintf(perf_filename, "Alg_%i_Performance_D_%.2f_S_%i_%s.txt",
	  alg_mode,change_degree,change_speed,problem_instance);
  //for performance
  if ((log_performance = fopen(perf_filename,"a")) == NULL) { exit(2); }


  //initialize and open output files
  div_filename = new char[CHAR_LEN];
  sprintf(div_filename, "Alg_%i_Div_D_%.2f_S_%i_%s.txt",
	  alg_mode,change_degree,change_speed,problem_instance);
  //for performance
  if ((log_diversity = fopen(div_filename,"a")) == NULL) { exit(2); }

}

void get_observation(int t){
  double performance;
  double diversity;
  int index_period;
  performance = get_current_best();
  //diversity = calc_lambda(0.05,problem_size,nn,pheromone);
  /*diversity = calc_diversity(n_ants);*/
  /*diversity = calc_entropy(problem_size,depth,pheromone,heuristic,alpha,beta);*/

  perf_of_iterations[t] += performance;
  //div_of_iterations[t] += diversity;
  offline_performance += performance;
  //pop_diversity += diversity;

  /*if(flag == true){
    index_period = (int)ceil(current_iteration/change_speed)-1;
    robustness = help / (double)get_current_best();
    //cout << help << " " <<  get_current_best() << endl;
    if(robustness > 1.0) robustness = 1.0;
    rob_of_changes[index_period] += robustness;
    tot_robustness += robustness;
    flag = false;
  }*/

 /* if(current_iteration%change_speed == 0){
    before = get_current_best();
    index_period = (int)ceil(current_iteration/change_speed)-1;
    before_of_changes[index_period] += before;
    before_change += before;
    flag = true;
    help = before; // store previous to use for robustness
  }*/

  cout<< "iteration: " << current_iteration << " best_so_far: " << performance << endl;
}

void get_mean(int r, double value) {

  perf_of_trials[r] = value;
  //before_of_trials[r] = before;
  //div_of_trials[r] = div;
  //rob_of_trials[r] = rob;
  //reset for next trials
  flag = false;
  help = 0.0;
}

double mean(double* values, int size){
  int i;
  double m = 0.0;
  for (i = 0; i < size; i++){
      m += values[i];
  }
  m = m / (double)size;
  return m; //mean
}

double stdev(double* values, int size, double average){
  int i;
  double dev = 0.0;

  if (size <= 1)
    return 0.0;
  for (i = 0; i < size; i++){
    dev += ((double)values[i] - average) * ((double)values[i] - average);
  }
  return sqrt(dev / (double)(size - 1)); //standard deviation
}

double calc_diversity(int size){

  int i, j;
   double div = 0.0;
   for (i = 0; i < size; i++){
     for (j = 0; j < size; j++) {
       if (i != j) {
         //div += distance_between(&ant_population[i], &ant_population[j]);//common edges
       }
     }
   }
   return (1.0 / (size * (size - 1))) * div; //population diversity
}

double calc_lambda(double l, int size, int depth, double** pheromone){

  int  i, m;
  double min, max, cutoff;
  double avg;
  double *num_branches = new double[size];
  double total;
  for (i = 0; i < problem_size; i++)
    num_branches[i] = 0.0;

  for ( m = 0 ; m < size ; m++ ) {
    /* determine max, min to calculate the cutoff value */
    min = pheromone[m][nn_list[m][1]];
    max = pheromone[m][nn_list[m][1]];
    for ( i = 1 ; i < depth ; i++ ) {
      if ( pheromone[m][nn_list[m][i]] > max )
        max = pheromone[m][nn_list[m][i]];
      if ( pheromone[m][nn_list[m][i]] < min )
        min = pheromone[m][nn_list[m][i]];
    }
    cutoff = min + l * (max - min);

    for ( i = 0 ; i < depth ; i++ ) {
      if (pheromone[m][nn_list[m][i]] > cutoff )
        num_branches[m] += 1.0;
    }
  }
  avg = 0.0;
  for ( m = 0 ; m < size ; m++ ) {
    avg += num_branches[m];
  }
  total =  avg / (double) (size);
  if (total < 2.0) total = 2.0;


  delete[] num_branches;

  return total;

}

double calc_entropy(int size, int depth, double** pheromone, double** heuristic, double a, double b){
	/*this behaviour measurement has not been tested*/
	int i, j;
	double sum_of_i = 0.0;
	double avg = 0.0;
        double prob;
	for(i = 0; i < size; i++){
       sum_of_i = 0.0;
	   for(j = 0; j < depth; j++){
         prob = pow(pheromone[i][nn_list[i][j]],a) * pow(heuristic[i][nn_list[i][j]],b);
		 sum_of_i -= ((prob) * log(prob));
	   }
	  avg += sum_of_i/(double)size;
	}
	return avg;
}

void close_stats(void){
  int i,j;
  double perf_mean_value, perf_stdev_value;
  double div_mean_value, div_stdev_value;

  //For graph plots
  for(i = 0; i < max_iterations; i++){
    perf_of_iterations[i] /= ((double)max_trials);
    div_of_iterations[i] /= ((double)max_trials);
    fprintf(log_performance, "%.2f", perf_of_iterations[i]);
    fprintf(log_performance,"\n");
    fprintf(log_diversity, "%.2f", div_of_iterations[i]);
    fprintf(log_diversity,"\n");
  }
  fprintf(log_performance,"\n");
  fprintf(log_performance,"Statistical results\n");
  fprintf(log_diversity,"\n");
  fprintf(log_diversity,"Statistical results\n");

  //For statistics
  for(i = 0; i < max_trials; i++){
    fprintf(log_performance, "%.2f", perf_of_trials[i]);
    fprintf(log_performance,"\n");
    fprintf(log_diversity, "%.2f", div_of_trials[i]);
    fprintf(log_diversity,"\n");
  }

  perf_mean_value = mean(perf_of_trials, max_trials);
  perf_stdev_value = stdev(perf_of_trials, max_trials, perf_mean_value);
  fprintf(log_performance,"Mean %f\t ", perf_mean_value);
  fprintf(log_performance,"\tStd Dev %f\t ", perf_stdev_value);
  cout << "offline performance: " << "<< mean " << perf_mean_value << " << stdev " << perf_stdev_value << endl;

  div_mean_value = mean(div_of_trials,max_trials);
  div_stdev_value = stdev(div_of_trials,max_trials,div_mean_value);
  fprintf(log_diversity, "Mean %f\t", div_mean_value);
  fprintf(log_diversity, "\tStd Dev %f\t ",div_stdev_value);
  cout << "population diversity: " << "<< mean " << div_mean_value << " << stdev " << div_stdev_value << endl;


  fclose(log_performance);
  fclose(log_diversity);


  //free memory
  delete[] perf_of_trials;
  delete[] perf_of_iterations;
  delete[] perf_filename;
  delete[] div_filename;
  delete[] div_of_trials;
  delete[] div_of_iterations;

}


