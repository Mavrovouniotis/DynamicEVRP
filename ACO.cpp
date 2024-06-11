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


int alg_mode; 

struct ant *ant_population;   //Population of ants or individuals


ant *best_so_far_ant;         //current best so far tour
ant *restart_best_ant;        //current best so far tour

bool flag_change_detected;

double **pheromone;          //Pheromone matrix
double **heuristic;          //Heuristic information matrix
double **total;              //Pheromone + Heuristic information matrix

double *prob_of_selection;   //Selection probabilities

//General ACO parameters
double alpha;
double be;
double q_0;
double trail0;
double rho;

double trail_max;       /* maximum pheromone trail in MMAS */
double trail_min;       /* minimum pheromone trail in MMAS */
int u_gb = INFTY;       /* every u_gb iterations update with best-so-far ant */
int restart_found_best;
int restart_iteration;

int iteration_best;

double found_branching;
double lambda;
bool new_best_found;

double branching_factor;
double branch_fac = 1.00001;

int n_ants;                   //Population size

int depth;                    //Candidate list size (nearest neighbour)
int **nn_list;                //Candidate lists
int nn;

bool *complete;

int current_iteration = 0;     //Used to calculate the period of change

bool ls_flag = false;          //indicates whether and which local search is used

int seed;                      //changing seed for the algorithms

int pseudonodes;


/****************************************************************/
/*                     Initialization                           */
/****************************************************************/

void allocate_ants(void){
  int i;
  int count = 0;
  int size_of_colonies;

   cout << problem_size << endl;
  //n_ants = problem_size+number_of_charging_stations;
  //n_ants = n_ants * num_of_colonies;
  n_ants = problem_size;

  pseudonodes = number_of_charging_stations * (problem_size_no_depot*2) + 250;
  ant_population = new ant[n_ants];

  complete = new bool[n_ants];
  //cout << n_ants << endl;
  for(i = 0; i < n_ants; i++) {
    ant_population[i].tour = new int[problem_size_no_depot+pseudonodes];
    ant_population[i].visited = new bool[problem_size];
    ant_population[i].routes = new int[problem_size_no_depot+pseudonodes];
    ant_population[i].acc_space = new double[problem_size_no_depot+pseudonodes];
    ant_population[i].space_available = new double[500];
    ant_population[i].id = i+1;
    ant_population[i].feasible = true;
  }

  best_so_far_ant = new ant;
  best_so_far_ant->tour = new int[problem_size_no_depot+pseudonodes];
  best_so_far_ant->visited = new bool[problem_size];
  best_so_far_ant->routes = new int[problem_size_no_depot+pseudonodes];
  best_so_far_ant->acc_space = new double[problem_size_no_depot+pseudonodes];
  best_so_far_ant->space_available = new double[500];
  best_so_far_ant->feasible = true;

  restart_best_ant = new ant;
  restart_best_ant->tour = new int[problem_size_no_depot+pseudonodes];
  restart_best_ant->visited = new bool[problem_size];
  restart_best_ant->routes = new int[problem_size_no_depot+pseudonodes];
  restart_best_ant->acc_space = new double[problem_size_no_depot+pseudonodes];
  restart_best_ant->space_available = new double[500];
  restart_best_ant->feasible = true;

  prob_of_selection = new double[depth+1];

  prob_of_selection[depth] = HUGE_VAL;

}

void allocate_structures(void){
  int i;
  heuristic = generate_2D_matrix_double(problem_size+number_of_charging_stations,problem_size+number_of_charging_stations);
  nn_list = generate_2D_matrix_int(problem_size+number_of_charging_stations,nn);


  pheromone = generate_2D_matrix_double(problem_size+number_of_charging_stations,problem_size+number_of_charging_stations);
  total = generate_2D_matrix_double(problem_size+number_of_charging_stations,problem_size+number_of_charging_stations);
  
  restart_best_ant = new ant;
  restart_best_ant->tour = new int[problem_size_no_depot+pseudonodes];
  restart_best_ant->visited = new bool[problem_size];
  restart_best_ant->routes = new int[problem_size_no_depot+pseudonodes];
  restart_best_ant->acc_space = new double[problem_size_no_depot+pseudonodes];
  restart_best_ant->space_available = new double[500];
  restart_best_ant->feasible = true;
  best_so_far_ant = new ant;
  best_so_far_ant->tour = new int[problem_size_no_depot+pseudonodes];
  best_so_far_ant->visited = new bool[problem_size];
  best_so_far_ant->routes = new int[problem_size_no_depot+pseudonodes];
  best_so_far_ant->acc_space = new double[problem_size_no_depot+pseudonodes];
  best_so_far_ant->space_available = new double[500];
  best_so_far_ant->feasible = true;
}



void set_algorithm_parameters(void){
  int i;
  //ACO common parameters
  alpha = 1;
  be = 5;
  depth = 20;
  //n_ants = 25;

  lambda = 0.05;
  rho = 0.8;
  q_0 = 0.0;
  u_gb = INFTY;
  


  if(depth >= problem_size) depth = problem_size - 1;
  nn = max(0,depth);
  if(nn >= problem_size) nn = problem_size-1;

  allocate_ants();
  allocate_structures();
}


void swap(int v[], int v2[],  int i, int j){
  int tmp;
  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
  tmp = v2[i];
  v2[i] = v2[j];
  v2[j] = tmp;
}


void sort(int v[], int v2[], int left, int right) {
  int k, last;

  if (left >= right)
    return;
  swap(v, v2, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      swap(v, v2, ++last, k);
  swap(v, v2, left, last);
  sort(v, v2, left, last);
  sort(v, v2, last+1, right);
}

double alg_random_number(int *idum){
  int k;
  double ans;
  //uniformly distributed random number [0,1]
  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}

void compute_nn_lists( void ){
  int i,j;
  int *distance_vector;
  int *help_vector;
  
  distance_vector = new int[problem_size];
  help_vector = new int[problem_size];
  //compute the nearest neigbhours of the objects
  for (j  = 0 ; j < problem_size; j++ ) {
    for (i = 0 ; i < problem_size; i++ ) {
      distance_vector[i] = distances[j][i];
      help_vector[i] = i;
    }
    distance_vector[j] = INFTY;
    sort(distance_vector, help_vector, 0, problem_size-1);
    for (i = 0 ; i < nn ; i++) {
      nn_list[j][i] = help_vector[i];
    }
  }

  //free memory
  delete[] distance_vector;
  delete[] help_vector;
}


void init_pheromone_trails(double init){
  int i,j;

  for(i = 0; i < problem_size+number_of_charging_stations; i++){
    for(j = 0; j <= i; j++){
      pheromone[i][j] = init;
      pheromone[j][i] = init;
      total[i][j] = init;
      total[j][i] = init;
    }
  }
}

void init_heuristic_info(void){
  int i,j;

  for(i = 0; i < problem_size+number_of_charging_stations; i++){
    for(j = 0; j <= i; j++){
     heuristic[i][j] = 1.0/((double) distances[i][j] + 0.1); //small value to avoid 1 div 0
      //heuristic[i][j] = distances[i][depot] + distances[depot][j] - distances[i][j] +EPSILON;
      heuristic[j][i] = heuristic[i][j];
    }
  }
}

void compute_total_info(void){
  int i,j;

  for(i = 0; i < problem_size+number_of_charging_stations; i++){
    for(j = 0; j < i; j++){
      total[i][j] = pow(pheromone[i][j],alpha) * pow(heuristic[i][j],be);
      total[j][i] = total[i][j];
    }
  }
}



void insert_station(ant *a, int node, int pos_c1){

  int i, n; 
  int id = 0;
  n = a->step + 1;

  for(i = n; i >= pos_c1; i--)
    a->tour[i] = a->tour[i-1];

  a->tour[pos_c1 - 1] = node;


  /*for(i = 0; i < n-1; i++){
    if(a->tour[i]==depot){
      id++;
    }
    a->routes[i] = id;
  }
  a->routes[n-1] = id;*/

  a->step++;
}

void remove_station(ant *a, int pos_c1){

  int i, n; 
  int id = 0;
  n = a->step + 1;

  i = pos_c1 + 1;
  while(i != n) {
    a->tour[i -1] = a->tour[i];
    i++;
  }


 /* for(i = 0; i < n-1; i++){
    if(a->tour[i]==depot){
      id++;
    }
    a->routes[i] = id;
  }
  a->routes[n-1] = id;*/

  a->step--;

}

void remove_stations(ant *a) {
  int i; 
  int from, to;
  i = 0;
  /*remove stations*/
  while(i!= a->step-1) {
    from  = a->tour[i];
    to = a->tour[i+1];
    if(is_charging_station(from) == true) {
      a->tour_length = a->tour_length - distances[a->tour[i-1]][from] - distances[from][to] + distances[a->tour[i-1]][to];
      remove_station(a,i);
      a->no_recharges--;
      i--;
    } else 
      i++;
  }
}


void satisfy_energy_constraint(ant *a){

  int i, j, h; 
  int from, to;
  double energy = battery_capacity;
  double capacity = max_capacity;
  double distance = 0.0;
  int station = 0;
  int offset = 0;
  double value_best = INFTY;
  double est_energy;
  double help;

  i = 0;
  /*remove stations*/

  /*if(alg_mode == 3) {
    while(i!= a->step-1) {
      from  = a->tour[i];
      to = a->tour[i+1];

      if(charging_station[from] == true) {
        //cout << "station" << from << endl;
        a->tour_length = a->tour_length - distances[a->tour[i-1]][from] - distances[from][to] + distances[a->tour[i-1]][to];
        remove_station(a,i);
        a->no_recharges--;
        i--;
      }

      i++;
    }
  }*/

  i = 0;
  while(i!=a->step-1) {
    from = a->tour[i];
    to = a->tour[i+1];

    if(from == depot) {
      energy = battery_capacity; 
      capacity = max_capacity;
    }
 
    help = energy - get_energy_consumption(from,to);
    //cout << help << endl;
    if(help < 0.0) {
      value_best = INFTY;
      station = 0;
      for(j = problem_size; j < (problem_size + number_of_charging_stations); j++){
      
          /*j != tour[i-1] avoids infite loop of adding the same station continously*/
         if((distances[a->tour[i-1]][j]+distances[j][a->tour[i]]) < value_best && j != a->tour[i-1] && (energy - get_energy_consumption(a->tour[i-1],j) < 0.0)){
            station = j;
            value_best = (distances[a->tour[i-1]][j]+distances[j][a->tour[i]]);
        }
      }
      /*int max = problem_size+number_of_charging_stations-1;
      int min = problem_size;
      station = floor(alg_random_number(&seed) * (max-min+1)+min);*/
     
      //no feasible station!!!
      //cannot add charging station next to another continously
      if(station == 0 || charging_station[a->tour[i-1]] == true) {
        a->feasible = false;
        break;
      } 

     insert_station(a, station, i+1);
     // a->no_recharges++;
     a->tour_length = a->tour_length -distances[a->tour[i-1]][a->tour[i+1]] + distances[a->tour[i-1]][a->tour[i]] + distances[a->tour[i]][a->tour[i+1]];
      //cout << "insert " << station << " between " << a->tour[i-1] << " " << a->tour[i] <<" " << a->tour[i+1]<< endl;
      // printTour(a);
      //cout << a->step << endl;
      a->no_recharges++;
      //remove stations inserted after station i
     /* h = i;
      while(a->tour[h]!=depot){
        h++;
        if(charging_station[a->tour[h]]==true) {
          //  printTour(a);
          // cout << "remove " << a->tour[h] << " between " << a->tour[h-1] << " " << a->tour[h] <<" " << a->tour[h+1]<< endl;
          a->tour_length =  a->tour_length -distances[a->tour[h-1]][a->tour[h]] -distances[a->tour[h]][a->tour[h+1]] + distances[a->tour[h-1]][a->tour[h+1]];
          remove_station(a,h);
          //cout << a->step << endl;
        }
      } //end while
      //start scanning from the begging of the tour
      i=0;*/

    } else {
      //update calculated energy
     
      energy -= get_energy_consumption(from,to);
      capacity -= init_objects[to].cust_demand;
    
      i++;
    }
    if(charging_station[to]==true) {
      energy = battery_capacity; 
    }

    if(i == a->step+pseudonodes) {
      a->feasible = false;
      break;
    }
  } //end of loop
  int id = 0;
  for(i = 0; i < a->step; i++) {
    if(a->tour[i]==depot){
      id++;
    }
    a->routes[i] = id;
  }
  a->routes[a->step-1] = id;
}


void checkTour(ant *a){
  int i, j ;
  int from, to;
  double energy = battery_capacity; 
  double capacity = max_capacity;
  double distance = 0.0;
  bool flag = true;
  double variable_rate;



  for(i = 0; i < a->step-1; i++){
    from = a->tour[i];
    to = a->tour[i+1];
    if(from > (problem_size_no_depot+number_of_charging_stations) || to > (problem_size_no_depot+number_of_charging_stations)) {
      flag = false;
      break;
    }
   
    if(from == depot){
      //if((i > 0) && (a->routes[i]!=a->routes[i-1])) {
        capacity = max_capacity;
        energy = battery_capacity;
      //}
    } 

    energy -= get_energy_consumption(from, to);
    if(energy < 0.0) {
     // cout << "error: energy below 0: " << " from " <<from << " to " << to << ", cap: " << capacity << " ene: "  << energy << endl;
         //printTour(a);
        // exit(1);
       flag = false;
       // exit(1);
       break;
    }

    if(charging_station[to]==true){
      energy = battery_capacity;
    }

    capacity -= init_objects[to].cust_demand;
    //cout << from << " " << to << " " << init_objects[to].cust_demand << " = " << capacity << endl;
    if(capacity < 0.0) {
         cout << "error: capacity below 0: " << from << " " << a->routes[i] << ", cap: " << capacity <<  endl;
         // printTour(a);
         // for(j = 0; j <= a->dummy_depots; j++){
          //  cout << a->space_available[j] << endl;
         // }
         // exit(1);
          flag = false;
           break;
    }

    distance += distances[from][to];

    //cout << time << endl;

   
    //cout << from << " " << to << " dis: " << distances[from][to] <<" : Cap: " << capacity << " Ene: " << energy << " Serv: " <<  time << endl;
   
  }
  //if(distance != a->tour_length) {
    //cout << "invalid solution quality: " <<  distance << endl;
    //printTour(a);
    //exit(1);
    //flag = false;
 // }
  if(a->no_customers != problem_size_no_depot) {
    //cout << "customers" << endl;
    //printTour(a);
    //exit(1);
    flag = false;
  }
  
  a->feasible = flag;
  //cout << a->feasible << endl;

}


/****************************************************************/
/*                    Construct Solutions                       */
/****************************************************************/
void ant_empty_memory(ant *a){
  
  int i;

  //clear previous ant solution
  for(i = 0; i < problem_size; i++) {

    a->visited[i] = false;
    a->routes[i] = 0;
  }

  for(i = 0; i < 100; i++) {
    a->space_available[i] = 0;
  }
  a->capacity = max_capacity; //initialize capacity
  a->energy_level = battery_capacity; //initialize energy level
  a->dummy_depots = 0;
  a->no_recharges = 0;
  a->no_customers = 0;
  a->step = 0;
  a->feasible = true;
}

void place_ant(ant *a){
  //place ants to the depot
  //a->start_new = true;
 
  a->tour[0] = depot;
  a->visited[depot] = true;
  a->dummy_depots+=1; //start first vehicle route
  a->routes[0]= a->dummy_depots;
  a->acc_space[0] = max_capacity;
  a->step = 1;

}

void update_energy(ant *a, int from, int to) {
  double variable_rate;

  a->energy_level -= get_energy_consumption(from,to);
  if(a->energy_level < 0) {
    cout << "energy level below 0" << endl;
    //exit(1);
  }
}


void update_capacity(ant *a, int to) {

  a->capacity -= init_objects[to].cust_demand;
  if(a->capacity < 0) {
    cout << "capacity below 0" << endl;
    //exit(1);
  }
}


bool check_energy(ant *a, int from, int to){
  if(a->energy_level - get_energy_consumption(from,to) >= 0.0) 
    return true;
  else 
    return false;

}

bool check_capacity(ant *a, int to) {

  if((a->capacity - init_objects[to].cust_demand) >= 0) 
    return true;
  else 
    return false;

}

//lookahead policy
bool check_range(ant *a, int from, int to){
  int i;
  double est_energy;
  double temp_cap;
  bool flag = false;
  est_energy = a->energy_level - get_energy_consumption(from,to);
  temp_cap = a->capacity - init_objects[to].cust_demand; 


  //cout << a->no_customers << " " <<problem_size_no_depot << endl;
  //check depot first
 
  if(get_energy_consumption(to,depot) <= est_energy) {
    flag = true;
  } else {
    for(i = problem_size; i < problem_size+number_of_charging_stations; i++){
      if(get_energy_consumption(to,i) <= est_energy) {
        flag = true;
        break;
      }
    }


  }

  return flag;
}



int find_closest_station(ant *a, int from, int to){

  int station_i = -1;
  int i; 
  double value_best = INFTY;
  //find the best station
  int previous = a->tour[a->step-1];
  //int previous2 = a->tour[a->step-3];
  //int previous3 = a->tour[a->step-4];
  //cout << previous << " " << previous2 << endl;
  //if(previous == depot) previous = -1;
 // if(previous2 == depot) previous2 = -1;
  //cout << previous << " " << previous2 << endl;
  //cout << previous << " " << previous2 << endl;
  //int previous3 = a->tour[a->step-4];
  for(i = problem_size; i < (problem_size + number_of_charging_stations); i++){

     if((distances[from][i]+distances[i][to]) < value_best && check_energy(a,from,i) == true && from !=i && i!=previous){
      station_i = i;
      value_best = (distances[from][i]+distances[i][to]);
    }
  }
  //check if the best station is the depot
 /* if((distances[from][depot]+distances[depot][to]) < value_best && check_energy(a,from,depot) == true &&from!=depot && depot!=previous){
    station_i = depot;
    value_best = (distances[from][depot]+distances[depot][to]);
  }*/
  /*if(station_i == -1){
      cout << "cannot find station from to"<< endl;
      for(i = problem_size; i < (problem_size + number_of_charging_stations); i++){
        cout << check_energy(a,from,i) << endl;
      }
     // cout << check_energy(a,from,depot) << endl;
      //exit(1);
  } else if (station_i == from) {
     cout << "Strangle!!" << endl;
    //exit(1);
  }*/

  return station_i;
}


void add_customer(ant *a, int to, int phase){

    a->visited[to] = true;
    a->tour[phase] = to;
    a->routes[phase] = a->dummy_depots;
    a->no_customers+=1;
    a->acc_space[phase] = a->capacity;
    //cout << "add customer: " << to << endl;

}

void add_depot(ant *a, int phase){
    a->space_available[a->dummy_depots] = a->capacity;
    a->dummy_depots+=1; // increase the route id
    a->capacity = max_capacity; //set capacity to the maximum
    a->energy_level = battery_capacity; //set energy level to full
    a->tour[phase] = depot; 
    a->routes[phase] = a->dummy_depots;
    a->acc_space[phase] = a->capacity;
    //cout << "return to depot " << depot << endl;
}

void add_station(ant *a, int station, int phase){
 
    a->tour[phase] = station;
    a->routes[phase] = a->dummy_depots;
    a->no_recharges+=1;
    a->acc_space[phase] = a->capacity;
    a->energy_level = battery_capacity;
    //cout << "station visited " << station <<endl;
}

void close_tour(ant *a, int phase) {
      
      int station;
      int from = a->tour[a->step-1];
      if(from!=depot){
        if(check_energy(a,from,depot)==true) {
          update_energy(a,from,depot);
           a->tour[a->step] = depot;
           a->routes[a->step] = a->dummy_depots;
          
           //add_depot(a,a->step); 
           a->step++;
        } else {
          station = find_closest_station(a,from,depot);
          //cout << check_energy(a,from,station) << endl;
         
          if(station == -1) {
            cout << "strange" << endl;
          } else {
            update_energy(a,from,station);
            add_station(a,station,phase);
            a->step++; 
          }
         // a->tour[a->step] = station;
          //a->routes[a->step] = a->dummy_depots;
       
          if(station!=depot && check_energy(a,station,depot)==true){
            update_energy(a,station,depot);
            a->tour[a->step] = depot;
            a->routes[a->step] = a->dummy_depots;
            a->step++; 
          } else {
            a->feasible = false;
            cout << "Cannot return to depot" << endl;
            //exit(1);
          }
        }
      }
   //update capacity for the last route
   a->space_available[a->dummy_depots] = a->capacity;
}




void add_component(ant *a, int from, int to, int phase) {

  int station;
  int help;
  //cout << " " << from << " " << to << endl;
  //cout << check_capacity(a,to) << check_energy(a,from,to) << check_range(a,from,to) <<a->visited[to] << endl;
  //cout << a->id << " " << a->step << " from " << from << " to " << to << " cust " << a->no_customers << endl;
  if(check_capacity(a,to) == true && 
    check_energy(a,from,to) == true && 
    check_range(a,from,to)==true && 
    from!=to && a->visited[to]==false ) {
    //add the customer and update energy and capacity
    
    //cout << a->no_customers << endl;
    // if(a->no_customers != problem_size_no_depot-1) {
     
    update_energy(a,from,to);
    update_capacity(a,to);
    add_customer(a,to,phase);
    
    /*} else { 
     // printTour(a);
      //cout << from << " " << to << endl; 
      if(check_energy(a,to,depot) == true) {
              update_energy(a,from,to);
              update_capacity(a,to);
              add_customer(a,to,phase);

      } else {
          if(check_energy(a,from,depot) == true) { 
              update_energy(a,from,depot);
              add_depot(a,phase);
          } else {
            station = find_closest_station(a,from,depot);
            if(station!=-1){
              update_energy(a,from,station);
              add_station(a,station,phase);
             } else {
              a->feasible = false;
              cout << "Strange!!!!" << endl;
            }
          }
        }
     } */
    
  } else if (check_capacity(a,to) == true && check_energy(a,from,to) == false) {
      //find a charging station to add instead of the "to" customer
     //cout << "energy" << endl;
      station = find_closest_station(a,from,to);
      if(station != -1) {
        update_energy(a,from,station);
        add_station(a,station,phase);
      } else if (station == -1 && check_energy(a,from,depot)==true) {
      //in case check_range returns only the depot within the range
        update_energy(a,from,depot);
        add_depot(a,phase);
      } else {
        a->feasible = false;
        add_depot(a,phase);
      }
       // cout << "Strange" << endl;
    

  } else if (check_capacity(a,to) == false && check_energy(a,from,depot)==true) {
     //EV cannot serve another customer
     //return to the depot when enough energy to close the route
     //out << capacity_left(a,from) << endl;
    //help = capacity_left(a,from);
    //if(help==-1){
     //all unvisited customer exceed capacity
     update_energy(a,from,depot);
     add_depot(a,phase);
     //} else { 
      //add the customer that leaves the least capacity in the vehicle
     // update_energy(a,from,help);
      //update_capacity(a,help);
      //add_customer(a,help,phase);
    //}
  } else if (check_capacity(a,to) == false && check_energy(a,from,depot)==false) {
      //EV cannot serve another customer
      //return to the depot to close the route but visit a charging station first
     //cout << "energy" << endl;
      station = find_closest_station(a,from,depot);
      update_energy(a,from,station);
      
     // if(station==depot) 
        //add_depot(a,phase); //in case the depot is in range for EV to return
     // else 
      add_station(a,station,phase);
 /*}else if (check_capacity(a,to) == true && check_energy(a,from,to)==true && check_range(a,from,to==false)) {
       station = find_closest_station(a,from);
       if (check_range(a,station,to) == true) {
        update_energy(a,from,station);
        add_station(a,station,phase);
       } else 
        cout << "false" << endl;*/
  } else { 
     //if (from!=to) if visited[to]==true or stucks in a cycle

     station = find_closest_station(a,from,to);
    // cout << "station " << station << endl;
     if(check_energy(a,from,depot)==true && station!=from && station == -1){
        update_energy(a,from,depot);
        add_depot(a,phase);
     } else if(check_energy(a,from,station) == true && station!=from && station != -1){
        update_energy(a,from,station);
        //if(station==depot) {
        //  add_depot(a,phase); //in case the depot is in range for EV to return
        //  cout << "depot " << station << endl;
      // } else {
          add_station(a,station,phase);
    
         // }
      }  else {
          cout << "Strange!!!" << endl;
          exit(1);
        }
  }
  a->step++; 
}


void choose_best_next(ant *a, int phase){
  int i,current, next;
  double value_best;
  next = problem_size_no_depot;
  current = a->tour[phase-1];

  value_best = -1.0;  //values in the list are always >=0.0
  //choose the next object with maximal (pheromone+heuristic) value
  //among all objects of the problem
  for(i = 0; i < problem_size; i++){
     if(a->visited[i])
    //if(a->visited[i] || check_capacity(a,i)==false)
      ;//if object visited
    else {
      if(total[current][i] > value_best){
        next = i;
        value_best = total[current][i];
      }
    }
  }

  add_component(a,current,next,phase);
}


void neighbour_choose_best_next(ant *a, int phase){
  int i,current,next,temp;
  double value_best, help;

    //next = depot;
  next = problem_size_no_depot;
  current = a->tour[phase-1];

  value_best = -1.0; //values in the list are always >=0.0

  //choose the next object with maximal (pheromone+heuristic) value
  //among all the nearest neighbour objects
  for(i =0; i < depth; i++){
    temp = nn_list[current][i];
    //if(a->visited[temp] || check_capacity(a,temp) == false)
    ////if(a->visited[temp] || charging_station[temp] == true  || check_range(a,current,temp)==false)
    if(a->visited[temp] || charging_station[temp] == true) 
      ;//if object visited
    else {
      help = total[current][temp];
      if(help > value_best){
        value_best = help;
        next = temp;
      }
    }
  }
  if(next == depot){
    //if all nearest neighnour objects are already visited
    choose_best_next(a,phase);
  } else {
    add_component(a,current,next,phase);
  }
}


void neighbour_choose_and_move_to_next(ant *a, int phase){
  int i,help, current, select;
  double rnd;
  double partial_sum = 0.0;
  double sum_prob = 0.0;
  double *prob_ptr;

  if((q_0 > 0.0) && (alg_random_number(&seed) < q_0)) {
    //with probability q_0 make the best possible choice
    neighbour_choose_best_next(a,phase);
    return;
  }
  prob_ptr = prob_of_selection; //selection probabilities of the nearest neigbhour objects
  current = a->tour[phase-1];

  //compute selection probabilities of nearest neigbhour objects
  for(i = 0; i < nn; i++){
    //if(a->visited[nn_list[current][i]] || check_capacity(a,nn_list[current][i]) == false){
    //if(a->visited[nn_list[current][i]] || charging_station[nn_list[current][i]] == true || check_range(a,current,nn_list[current][i])==false){ 
    if(a->visited[nn_list[current][i]] || charging_station[nn_list[current][i]] == true){ 
      prob_ptr[i] = 0.0;
    } else {
      prob_ptr[i] = total[current][nn_list[current][i]];
      sum_prob +=prob_ptr[i];
    }
  }
  if(sum_prob <= 0.0){
    //in case all neigbhbour objects are visited
    choose_best_next(a,phase);
  } else{
    //proabilistic selection (roullete wheel)
    rnd = alg_random_number(&seed);
    rnd *= sum_prob; //normalize
    select = 0;
    partial_sum = prob_ptr[select];
    //This loop always stops because prob_ptr[nn_ants] == HUGE_VAL
    while(partial_sum<=rnd){
      select++;
      partial_sum+=prob_ptr[select];
    }
    //this may very rarely happen because of rounding if rnd is close to 1
    if(select==depth){
      neighbour_choose_best_next(a,phase);
      return;
    }
    help = nn_list[current][select];
    add_component(a,current,help,phase);
  }
}


void copy_from_to(ant *a1, ant *a2){

  int i;
  //ant2 is a copy of ant1
  for(i = 0; i < a1->step; i++){
    a2->tour[i] = a1->tour[i];
    a2->routes[i] = a1->routes[i];
    a2->acc_space[i] = a1->acc_space[i];
  }
  a2->id = a1->id;
  a2->tour_length = a1->tour_length;
  a2->dummy_depots = a1->dummy_depots;
  a2->capacity = a1->capacity;
  a2->no_customers = a1->no_customers;
  a2->step = a1->step;
  a2->no_recharges = a1->no_recharges;
  a2->feasible = a1->feasible;
  for(i = 0; i < a1->dummy_depots+1; i++) {
    a2->space_available[i] = a1->space_available[i];
  }

}




void choose_closest_next(ant *a, int phase){

  int i,current,next,min;

  next = depot;
  current = a->tour[phase-1];
  min = INFTY;
  //choose closest object used in the nn_tour()
  for(i = 0; i < problem_size; i++){
    if(a->visited[i] || charging_station[i] == true )
   // if(a->visited[i] || charging_station[i] == true || check_capacity(a,i) == false)
      ; //if object not visited
     else {
      if(distances[current][i] < min){
        next = i;
        min = distances[current][i];
      }
    }
  }

  add_component(a,current,next,phase);
  
}

int nn_tour(void){

  int phase, help;
  phase=help=0;
  ant_empty_memory(&ant_population[0]);
  place_ant(&ant_population[0]);
  //compute the tour length of the nearest neigbour heuristic
  //used to initialize the pheromone trails
  while(ant_population[0].no_customers < problem_size_no_depot){
    choose_closest_next(&ant_population[0],ant_population[0].step);
  }

  ant_population[0].tour[ant_population[0].step] = ant_population[0].tour[0];
  ant_population[0].step++;

  //if(ls_flag == true) two_opt_first(ant_population[0].tour, ant_population[0].routes,ant_population[0].dummy_depots);

  ant_population[0].tour_length = dummy_fitness_evaluation(ant_population[0].tour,ant_population[0].step);
  help = ant_population[0].tour_length;

  //printTour(&ant_population[0]);
  ant_empty_memory(&ant_population[0]);

  return help;
}



void constant_pheromone_deposit(ant *a, double dT){
  int i,j,h;

  for(i = 0; i < a->step-1; i++){
    j = a->tour[i];
    h = a->tour[i+1];
    pheromone[j][h]+= dT;
    pheromone[h][j] = pheromone[j][h];
    /*compute total*/
    total[j][h] = pow(pheromone[j][h],alpha) * pow(heuristic[j][h],be);
    total[h][j] = total[j][h];

  }
}

void constant_pheromone_removal(ant *a, double dT){
  int i,j,h;

  for(i = 0; i < a->step-1; i++){
    j = a->tour[i];
    h = a->tour[i+1];
    pheromone[j][h]-= dT;
    pheromone[h][j] = pheromone[j][h];
     /*compute total*/
    total[j][h] = pow(pheromone[j][h],alpha) * pow(heuristic[j][h],be);
    total[h][j] = total[j][h];
  }
}

void action_when_change_detected(void) {

  int i, j;

  if(flag_change_detected == true) {
    flag_change_detected = false;

    /*restart always R1*/ 
    if(alg_mode == 1) {
        restart_best_ant->tour_length = INFTY;
        best_so_far_ant->tour_length = INFTY; //reset best solution found
        restart_iteration = 1;
        restart_found_best = 0;
        trail_max = 1.0 / ( (rho) * nn_tour() );
        trail_min = trail_max / ( 2.0 * problem_size );
        trail0 = trail_max;
        init_pheromone_trails(trail0);  //initialize pheromone trails
        compute_total_info();           //combine heuristic+pheromone
    }
    /*detect infeasibility and restart algorithm R2*/
    if(alg_mode == 2) {
      if(check_solution(best_so_far_ant->tour,best_so_far_ant->step,best_so_far_ant->feasible) == true) {
        cout << "valid" << endl;
      } else { 
        cout << "invalid" << endl;
        restart_best_ant->tour_length = INFTY;
        best_so_far_ant->tour_length = INFTY; //reset best solution found
        restart_iteration = 1;
        restart_found_best = 0;
        trail_max = 1.0 / ( (rho) * nn_tour() );
        trail_min = trail_max / ( 2.0 * problem_size );
        trail0 = trail_max;
        init_pheromone_trails(trail0);  //initialize pheromone trails
        compute_total_info();           //combine heuristic+pheromone
      }
    }
    /*detect infeasibility, repair action and then restart R3*/
    if(alg_mode == 3) {
        if(check_solution(best_so_far_ant->tour,best_so_far_ant->step,best_so_far_ant->feasible) == true) {
        cout << "valid" << endl;
      } else { 
        cout << "invalid" << endl;
        satisfy_energy_constraint(best_so_far_ant);
        if(best_so_far_ant->feasible == true) {
          //cout << check_solution(best_so_far_ant->tour,best_so_far_ant->step,best_so_far_ant->feasible) << endl;
          best_so_far_ant->tour_length = fitness_evaluation(best_so_far_ant->tour,best_so_far_ant->step,best_so_far_ant->feasible);
        } else {
           best_so_far_ant->tour_length = INFTY;
        }
        satisfy_energy_constraint(restart_best_ant);
        if(restart_best_ant->feasible == true) { 
          restart_best_ant->tour_length = fitness_evaluation(restart_best_ant->tour, restart_best_ant->step, restart_best_ant->feasible);
        } else {
          restart_best_ant->tour_length = INFTY;
        }
        restart_iteration = 1;
        restart_found_best = 0;
        trail_max = 1.0 / ( (rho) * nn_tour() );
        trail_min = trail_max / ( 2.0 * problem_size );
        trail0 = trail_max;
        init_pheromone_trails(trail0);  //initialize pheromone trails
        compute_total_info();           //combine heuristic+pheromone
      }
    }

    /*repair the solution and adapt MMAS_A*/
    if(alg_mode == 4) {
      if(check_solution(best_so_far_ant->tour,best_so_far_ant->step,best_so_far_ant->feasible) == true) {
        cout << "valid" << endl;
      } else { 
        cout << "invalid" << endl;
        satisfy_energy_constraint(best_so_far_ant);
        if(best_so_far_ant->feasible == true) {
          //cout << check_solution(best_so_far_ant->tour,best_so_far_ant->step,best_so_far_ant->feasible) << endl;
          best_so_far_ant->tour_length = fitness_evaluation(best_so_far_ant->tour,best_so_far_ant->step,best_so_far_ant->feasible);
        } else {
           best_so_far_ant->tour_length = INFTY;
        }
      }
    }
  }
}


void construct_solutions(void){

  int k, all_ants;
  action_when_change_detected();
  //clear memory of ants
  for(k = 0; k < n_ants; k++) {
    ant_empty_memory(&ant_population[k]);
    complete[k] = false;
  }

  //place ants on a random object
  for(k = 0; k < n_ants; k++)
    place_ant(&ant_population[k]);

  all_ants = 0;

  //select object until all objects are visited
  while(all_ants < n_ants) {
    for(k = 0; k < n_ants; k++){
      if(ant_population[k].no_customers < problem_size_no_depot){
        neighbour_choose_and_move_to_next(&ant_population[k],ant_population[k].step);
        //cout << ant_population[k].step << " " << ant_population[k].no_customers << endl;
      } else {
        if(complete[k] == false) {
           complete[k] = true;
           all_ants++; //until all ants satisfied all customer demands
        }
      }
    }
  }

  for(k = 0; k < n_ants; k++) {
    //add depot at the end to close the tour
    close_tour(&ant_population[k],ant_population[k].step);  
    check_solution(ant_population[k].tour, ant_population[k].step, ant_population[k].feasible);
    if(ant_population[k].feasible == true) {
      ant_population[k].tour_length = fitness_evaluation(ant_population[k].tour,ant_population[k].step,ant_population[k].feasible);
      //cout << k << " " << ant_population[k].tour_length << endl;
    } 
  }
}



/****************************************************************/
/*                   MMAS Pheromone Update                      */
/****************************************************************/

void global_pheromone_deposit(ant *a){

  int i,j,h;
  double d_tau;
  d_tau = 1.0 / (double)a->tour_length;

  for(i = 0; i < a->step-1; i++) {
      j = a->tour[i];
      h = a->tour[i+1];
      pheromone[j][h] += d_tau;
      pheromone[h][j] = pheromone[j][h];
  }   
}
   

void check_pheromone_trail_limits( void ){
    int i, j;

    for ( i = 0 ; i < problem_size+number_of_charging_stations; i++ ) {
	    for ( j = 0 ; j < i ; j++ ) {
	       if ( pheromone[i][j] < trail_min ) {
		         pheromone[i][j] = trail_min;
		         pheromone[j][i] = trail_min;
	       } else if ( pheromone[i][j] > trail_max ) {
		         pheromone[i][j] = trail_max;
	  	       pheromone[j][i] = trail_max;
	       }
      }
    }
}


int find_best(){

  int k,min,k_min;
  min = ant_population[0].tour_length;
  k_min = 0;
  for(k = 1; k < n_ants; k++) {
    // cout << ant_population[k].id << " " << id << endl;
    //cout << ant_population[k].tour_length << " " << min << endl;
    if(ant_population[k].tour_length < min && ant_population[k].feasible == true) {
      min=ant_population[k].tour_length;
      k_min = k;
    }
  }
  //cout << " n_ants: " << n_ants << endl;
  return k_min; //population best ant index
}


void evaporation(void){
 int i, j;

 for (i = 0 ; i < problem_size+number_of_charging_stations; i++) {
 	  for (j = 0 ; j <= i; j++) {
      pheromone[i][j] = (1 - rho) * pheromone[i][j];
	    pheromone[j][i] = pheromone[i][j];
	  }
  }
}

void mmas_pheromone_update(void){
    int i;

    evaporation();

    if (current_iteration%u_gb) {
	     //iteration_best = find_best();
	     global_pheromone_deposit(&ant_population[iteration_best]);
    } else {
       if(u_gb == 1 && (current_iteration - restart_found_best > 50)){
	        global_pheromone_deposit(best_so_far_ant);
       } else {
          global_pheromone_deposit(restart_best_ant);
       }
    }
    if(ls_flag == true) {
      if ( ( current_iteration - restart_iteration ) < 25 )
        u_gb = 25;
      else if ( (current_iteration - restart_iteration) < 75 )
        u_gb = 5;
      else if ( (current_iteration - restart_iteration) < 125 )
        u_gb = 3;
      else if ( (current_iteration - restart_iteration) < 250 )
        u_gb = 2;
      else
        u_gb = 1;
    } else
       u_gb = 25;

    check_pheromone_trail_limits();
}

void pheromone_update(){

  mmas_pheromone_update();
  compute_total_info();
}


/****************************************************************/
/*                    Update Best Ants                          */
/****************************************************************/

void update_best(void){
  double min;
  double p_x;
  int i_min;
  min = INFTY;
  i_min = -1;
  int i;
    iteration_best = find_best();
    if(ant_population[iteration_best].tour_length < best_so_far_ant->tour_length) {
      copy_from_to(&ant_population[iteration_best], best_so_far_ant);
      //update_tour_length(ant_population[iteration_best].tour_length); 
      copy_from_to(&ant_population[iteration_best],restart_best_ant);
      restart_found_best = current_iteration;
      if(ls_flag == false) {
        p_x = exp(log(0.05)/problem_size);
        trail_min = 1.0 * (1.0 - p_x) / (p_x * (double)((depth + 1) / 2.0));
        trail_max = 1.0 / ( (rho) * best_so_far_ant->tour_length );
        trail0 = trail_max;
        trail_min = trail_max * trail_min;
      } else {
        trail_max = 1. / ( (rho) * best_so_far_ant->tour_length );
        trail_min = trail_max / ( 2. * problem_size );
        trail0 = trail_max;
      }
    }
    if(ant_population[iteration_best].tour_length < restart_best_ant->tour_length) {
      copy_from_to(&ant_population[iteration_best], restart_best_ant);
      restart_found_best = current_iteration;
      cout << " restart best: " << restart_best_ant->tour_length << " restart_found_best: " << restart_found_best << endl;
    }
}

/****************************************************************/
/*                    Update Statistics                         */
/****************************************************************/

double node_branching(double l) {
  int  i, m;
  double min, max, cutoff;
  double avg;
  double *num_branches = new double[problem_size];
  for (i = 0; i < problem_size; i++){
    num_branches[i] = 0.0;
  }

  for ( m = 0 ; m < problem_size; m++ ) {
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
      if (pheromone[m][nn_list[m][i]] > cutoff)
        num_branches[m] += 1.0;
    }
  }
  avg = 0.0;
  for ( m = 0 ; m < problem_size; m++ ) {
    avg += num_branches[m];
  }
  delete[] num_branches;

  /* Norm branching factor to minimal value 1 */
  return ( avg / (double) ((problem_size) * 2)  );
}


void statistics_and_output(){

  int i;
  current_iteration++;
  if (!(current_iteration%10)) {
	 branching_factor = node_branching(lambda);
   if ((branching_factor < branch_fac) && (current_iteration - restart_found_best > 250) ) {
	 //if ( (current_iteration - multi_colonies[i].restart_found_best > 250) ) {
      cout << "INIT TRAILS!!!" << endl;
      restart_best_ant->tour_length = INFTY;
      init_pheromone_trails(trail_max);
      compute_total_info();
      restart_iteration = current_iteration;
    }
  }
  //cout << get_current_best() << endl;
  //get_observation(current_iteration-1,get_current_best());
  //output results
  //cout<< "iteration: " << current_iteration << " evals " << get_evals() <<  " best_so_far: " << get_current_best() << " %: " << get_current_error() << endl;

}

/****************************************************************/
/*                    Initialization                            */
/****************************************************************/

void init_try(int r){
    int i;
  
    current_iteration = 0;
    seed = r; // change the seed for another run

    init_heuristic_info();          //initialize heuristic info
    compute_nn_lists();

    best_so_far_ant->tour_length = INFTY; //reset best solution found
    restart_best_ant->tour_length = INFTY;
    restart_iteration = 1;
    restart_found_best = 0;
    branching_factor = node_branching(lambda);
    trail_max = 1.0 / ( (rho) * nn_tour() );
    trail_min = trail_max / ( 2.0 * problem_size );
    trail0 = trail_max;
    
    init_pheromone_trails(trail0);  //initialize pheromone trails
    compute_total_info();           //combine heuristic+pheromone

    flag_change_detected = false;
}


void ACO(){

      //cout << "pass 0" << endl;
      construct_solutions();   //from ACO.h
      //cout << "pass 1" << endl;
      update_best();           //from ACO.h
      //cout << "pass 2" << endl;
      pheromone_update();      //from ACO.h
      //cout << "pass 3" << endl;
      statistics_and_output(); //from ACO.h
      //cout << "pass 4" << endl;
}


void free_ACO() {

  int i,j;

  for(i = 0; i < n_ants; i++){
    delete[] ant_population[i].tour;
    delete[] ant_population[i].visited;
    delete[] ant_population[i].routes;
    delete[] ant_population[i].acc_space;
    delete[] ant_population[i].space_available;
  }
  delete[] ant_population;

  for(j = 0; j < problem_size+number_of_charging_stations; j++){
    delete[] pheromone[j];
    delete[] total[j];
  }

  delete[] pheromone;
  delete[] total;
  delete[] best_so_far_ant;
  delete[] restart_best_ant;
  
  delete[] best_so_far_ant->tour;
  delete[] best_so_far_ant->visited;
  delete[] best_so_far_ant->routes;
  delete[] best_so_far_ant->acc_space;
  delete[] best_so_far_ant->space_available;

  delete[] restart_best_ant->tour;
  delete[] restart_best_ant->visited;
  delete[] restart_best_ant->routes;
  delete[] restart_best_ant->acc_space;
  delete[] restart_best_ant->space_available;
 
  delete[] complete;
  delete[] heuristic;
  delete[] nn_list;
  delete[] prob_of_selection;

}

