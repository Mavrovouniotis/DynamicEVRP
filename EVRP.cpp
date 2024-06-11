#include<iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<cstring>
#include<math.h>
#include<fstream>
#include<time.h>
#include<limits.h>
//#include <sys/time.h>
//#include <sys/resource.h>

#include "EVRP.hpp"
#include "stats.hpp"
#include "ACO.hpp"

using namespace std;

struct object *init_objects;     //Actual objects of the instance

int depot;                       //depot id (usually 0)
int problem_size_no_depot;       //Size of the instance (excludes depot)
int problem_size_station;

double **distances;              //Distance matrix
char* problem_instance;          //Name of the instance
int problem_size;                //Number of cities

double offline_performance;
double offline_error;
double global_optimum_value;
int evals;
double current_best;
double current_error;
double best_found_iteration;
double best_found_time;

int number_of_charging_stations;
bool* charging_station;

int battery_capacity; //maximum energy of vehicles
int max_capacity;        //capacity of vehicles

double energy_consumption;
double **modify_factors;         /*modified factors matrix*/
double change_degree;            /*ratio of swapped objects*/
int change_speed;                /*period of changes in algorithmic iterations*/
int total_changes;               /*how many environmental changes*/

/*this parameters can be modified to increase or decrease the amount of change*/
double mu = 0.0;
double sigma;

double before_change;
double pop_diversity;
double tot_robustness;

int env_index;

/****************************************************************/
/*Random number generator where the seed is the same in all runs*/
/****************************************************************/
double env_random_number(double low, double high) {
  return ((double)(rand()%10000)/10000.0)*(high - low) + low;
}

/****************************************************************/
/*Random seed used from the env_random_number() method          */
/****************************************************************/
void reset_random_seed(){
  /*sreset seed to get the same consecutive dynamic changes for each run*/
  srand(1);
}

/****************************************************************/
/*Normally distributed number used to increase or decrease      */
/*the weights                                                   */
/****************************************************************/
double normal(double mu, double sigma, double lower, double higher){
  double p, p1, p2, v;
  do {
    p1 = env_random_number(-1.0,1.0);
    p2 = env_random_number(-1.0,1.0);

    p = (p1 * p1) + (p2 * p2);
  } while (p >= 1.0);
   // if (p == 0 || p > 1) return rand_num(-1.0,1.0);
    v = mu + (sigma * p1 * sqrt(-2.0 * (log(p)/p)));
    return v;
}

/****************************************************************/
/*           Initialize all modify factors equally              */
/****************************************************************/
void initialize_modify_factors(){
  int i,j;
  for(i = 0; i < (problem_size+number_of_charging_stations); i++){
   for(j = 0; j <= i; j++) {
     modify_factors[i][j] = 0.0;
     modify_factors[j][i] = modify_factors[i][j];
   }
  }
}


/****************************************************************/
/*Compute and return the euclidean distance of two objects      */
/****************************************************************/
double euclidean_distance(int i, int j) {
  double xd,yd;
  double r = 0.0;
  xd = init_objects[i].x - init_objects[j].x;
  yd = init_objects[i].y - init_objects[j].y;
  r  = sqrt(xd*xd + yd*yd) + 0.5;
 
  return r;
}

/****************************************************************/
/*Compute the distance matrix of the problem instance           */
/****************************************************************/
void compute_distances(void) {
int i, j;
int e2d;
  for(i = 0; i < (problem_size+number_of_charging_stations); i++){
    for(j = 0; j < (problem_size+number_of_charging_stations); j++){
      e2d = euclidean_distance(i,j);
      distances[i][j] = e2d + modify_factors[i][j];
      //in case a negative distance occurs
      if(distances[i][j] <= 0 && i!=j) distances[i][j] = e2d + abs(modify_factors[i][j]); 
    }
  }
}


/****************************************************************/
/*          Modify the arc from one node to another node        */
/****************************************************************/
void modify_weight(int from, int to){
  double value;
  //cout << value << endl;
  sigma = 0.2 * euclidean_distance(from,to);
  value = normal(mu,sigma,-10,10);
  modify_factors[from][to] = value;
  modify_factors[to][from] = modify_factors[from][to];
}


/****************************************************************/
/* Perform random changes to a predefined arcs according to the */
/* magnitude of change                                          */
/*                                                              */
/****************************************************************/
void add_random_change(void){
  int i,j,k;
  //calculate total number of arcs for symmetric cases
  int num_arcs = (problem_size * (problem_size - 1))/2.0;
  //int k = (problem_size * (problem_size - 1)); for asymmetric
  int changes = (int)abs(change_degree * num_arcs);
  int** changed = generate_2D_matrix_int(problem_size,problem_size);

  /*Implementing the clustered method in DTSP with node changes it is also 
  possible to restrict the selection of the arcs from a set of clustered 
  cities. In this implementation the selected arcs are NOT restricted*/

  //vary the weights of the arcs
  k = 0;
  for(i = 0; i < problem_size; i++){
    for(j = 0; j < i; j++){
      if(env_random_number(0.0,1.0) <= change_degree){
        if(k>=changes) break;
        modify_weight(i,j);
        changed[i][j] = 1;
        changed[j][i] = 1;
        k++;
      }
    }
  }

  //vary the weights of more arcs in case the total changes number does not match k
  int h1, h2;
  while(k < changes){
    h1 = env_random_number(0,problem_size);
    h2 = env_random_number(0,problem_size);
    if(changed[h1][h2] == 0){
      modify_weight(h1,h2);
      changed[h1][h2] = 1;
      changed[h2][h1] = 1;
      k++;
    }
  }

  //cout << "total arcs" << k << " " << changes <<  endl;

  //free changed memory
  for(int i = 0; i < problem_size; i++) {
    delete[] changed[i];
  }
  delete[] changed;
}

/****************************************************************/
/*Generate and return a two-dimension array of type int         */
/****************************************************************/
int ** generate_2D_matrix_int(int n, int m){
  int **matrix;
  matrix = new int*[n];
  for ( int i = 0 ; i < n ; i++ ) {
    matrix[i] = new int[m];
  }
  //initialize 2-d array
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) {
      matrix[i][j] = 0;
    }
  }
  return matrix;
}

/****************************************************************/
/*Generate and return a two-dimension array of type double      */
/****************************************************************/
double ** generate_2D_matrix_double(int n, int m){
  double **matrix;

  matrix = new double*[n];
  for ( int i = 0 ; i < n ; i++ ) {
    matrix[i] = new double[m];
  }
  //initialize the 2-d array
  for(int i = 0; i < n; i++){
    for(int j = 0; j < m; j++) {
      matrix[i][j] = 0.0;
    }
  }
  return matrix;
}


/****************************************************************/
/* Read the problem instance and generate the initial object    */
/* vector.                                                      */
/****************************************************************/
void read_problem(char* filename){
    int i;
  char line[CHAR_LEN];
  char * keywords;
  char Delimiters[] = " :=\n\t\r\f\v";
  ifstream fin(filename);
  while((fin.getline(line, CHAR_LEN-1))){

    if(!(keywords = strtok(line, Delimiters)))
      continue;
    if(!strcmp(keywords, "DIMENSION")){
      if(!sscanf(strtok(NULL, Delimiters), "%d", &problem_size)){
	     cout<<"DIMENSION error"<<endl;
	     exit(0);
      }
    }
    else if(!strcmp(keywords, "EDGE_WEIGHT_TYPE")){
      char * tempChar;
      if(!(tempChar=strtok(NULL, Delimiters))){
	cout<<"EDGE_WEIGHT_TYPE error"<<endl;
	exit(0);
      }
      if(strcmp(tempChar, "EUC_2D")){
	cout<<"not EUC_2D"<<endl;
	exit(0);
      }
    }
    else if (!strcmp(keywords, "CAPACITY")){
       if(!sscanf(strtok(NULL,Delimiters), "%d", &max_capacity)){
          cout << "CAPACITY error" << endl;
          exit(0);
       }
    }
    else if (!strcmp(keywords, "ENERGY_CAPACITY")){
       if(!sscanf(strtok(NULL,Delimiters), "%d", &battery_capacity)){
          cout << "ENERGY_CAPACITY error" << endl;
          exit(0);
       }
    }
    else if (!strcmp(keywords, "ENERGY_CONSUMPTION")){
       if(!sscanf(strtok(NULL,Delimiters), "%lf", &energy_consumption)){
          cout << "ENERGY_CONSUMPTION error" << endl;
          exit(0);
       }
    }
    else if (!strcmp(keywords, "STATIONS")){
       if(!sscanf(strtok(NULL,Delimiters), "%d", &number_of_charging_stations)){
          cout << "STATIONS error" << endl;
          exit(0);
       }
    }
    else if (!strcmp(keywords, "OPTIMAL_VALUE")){
       if(!sscanf(strtok(NULL,Delimiters), "%lf", &global_optimum_value)){
          cout << "OPTIMAL_VALUE error" << endl;
          exit(0);
       }
    }
    else if(!strcmp(keywords, "NODE_COORD_SECTION")){
      if(problem_size!=0){
        problem_size = problem_size + number_of_charging_stations;
        problem_size = problem_size - number_of_charging_stations;
        problem_size_no_depot = problem_size - 1;
         //number_of_charging_stations = 0;
        cout << problem_size << " " << number_of_charging_stations << endl;
        init_objects = new object[problem_size + number_of_charging_stations];
      
        for (i=0; i < problem_size+number_of_charging_stations; i++){
	         //store initial objects
           fin>>init_objects[i].id;
	         fin>>init_objects[i].x>>init_objects[i].y;
           init_objects[i].id -=1;
	         //initialize masked_objects with initial objects
        }
        //compute the distances using initial objects
        distances = generate_2D_matrix_double(problem_size+number_of_charging_stations, problem_size+number_of_charging_stations);
        modify_factors = generate_2D_matrix_double(problem_size+number_of_charging_stations,problem_size+number_of_charging_stations);
        //initialize_modify_factors();
        //compute_distances();
     } else {
           cout << "wrong problem instance file" << endl;
            exit(1);
     }
    }
    else if(!strcmp(keywords, "DEMAND_SECTION")){
     if(problem_size!=0){

       int temp;
       //masked_demand = new int[problem_size];
       //cust_demand = new int[problem_size+number_of_charging_stations];
       charging_station = new bool[problem_size+number_of_charging_stations];
       for(i =0; i < problem_size; i++){
        fin >> temp;
        fin >> init_objects[temp-1].cust_demand;
       }

       //initialize the charging stations. 
       //the depot is set to a charging station in the DEPOT_SECTION
       for(i = 0; i < problem_size+number_of_charging_stations; i++){
          if(i < problem_size) {
            charging_station[i] = false;
          } else {
            charging_station[i] = true;
            init_objects[i].cust_demand = 0;
          }
          //cout << init_objects[i].id << " " << init_objects[i].x << " "  << init_objects[i].y << " " << charging_station[i] << " " << cust_demand[i] << endl;
       } 
       //compute_distances();
      }
     }
     else if(!strcmp(keywords, "DEPOT_SECTION")){
      fin >> depot;
      depot-=1;
      charging_station[depot] = false;
     }
   
  }
  fin.close();
  if(problem_size == 0) {
           cout << "wrong problem instance file" << endl;
            exit(1);
  }

  //cout << max_capacity << " " << battery_capacity << " " << energy_consumption << endl;
  //cout << problem_size << " " << problem_size_no_depot << endl;

}

/****************************************************************/
/* Initialize the environment with the initial objects and      */
/* perform the initial dynamic change/generate base states      */
/****************************************************************/
void initialize_environment(){
  //Set the random seed
  //reset the metrics
  reset_random_seed();
  env_index = 0;

  offline_performance = 0.0;
  evals = 0;
  current_best = INFTY;

  initialize_modify_factors();

  //add_random_change();
  compute_distances();
  current_best = INFTY;

}

/****************************************************************/
/* Perform the dynamic change every "period" iteration          */
/*                                                              */
/****************************************************************/
void change_environment(){

    add_random_change();
    compute_distances();          /*update weight matrix*/
    current_best = INFTY;         /*reset current best*/
    flag_change_detected = true;
    env_index++;
}


double fitness_evaluation(int *t, int s, bool f) {
  int i;
  double tour_length;
  if(f == true) 
    tour_length = 0.0;
  else
    tour_length = 100000.0;
  evals++;
  //evaluation
  /*-----------------------------------------------------------*/
  //tour_length = distances[0][t[0]]; //depot - first customer
  for (i = 0; i < s-1; i++) {
    // if(id[i] == id[i+1]){
    //tour_length += distances[t[i]][t[(i+1)%s]];
    tour_length += distances[t[i]][t[i+1]];
    // cout << "from " << t[i] << " to " << t[i+1] << endl;
  }
  /*-----------------------------------------------------------*/
  // if(s = problem_size_no_depot+pseudonodes){
  // tour_length += 100000;
  //}

  if(tour_length < current_best){
      current_best = tour_length;
  }
  //cout << tour_length << endl
  return tour_length;
}

void update_tour_length(int length) {

  if(length < current_best) {
    current_best = length;
  }

}

/****************************************************************/
/* Evaluate routes     and return the length without count being*/
/* counted   					                                      		*/
/****************************************************************/
double dummy_fitness_evaluation(int *t, int s) {
  int i;
  double tour_length = 0.0;
  //evaluation
  /*-----------------------------------------------------------*/
  //tour_length = distances[0][t[0]]; //depot - first customer
  for (i = 0; i < s-1; i++) {
   // if(id[i] == id[i+1]){
    //tour_length += distances[t[i]][t[(i+1)%s]];
    tour_length += distances[t[i]][t[i+1]];
    // cout << "from " << t[i] << " to " << t[i+1] << endl;

  }
  /*-----------------------------------------------------------*/

  return tour_length;
}


/*return performance and behaviour measurements */
double get_offline_performance(){
    return offline_performance/(double)max_iterations;
}


double get_pop_diversity(){
    return pop_diversity/(double)max_iterations;
}

double get_current_best(){
 return current_best;
}

int get_evals(){
 return evals;
}

int get_problem_size(){
   return problem_size;
}


/****************************************************************/
/* Returns the energy consumed when travelling between two      */
/* points: from and to.                                         */
/****************************************************************/
double get_energy_consumption(int from, int to) {

    /*DO NOT USE THIS FUNCTION MAKE ANY CALCULATIONS TO THE ROUTE COST*/
    return energy_consumption*distances[from][to];

}

double get_distance(int from, int to){
  //adds partial evaluation to the overall fitness evaluation count
  //It can be used when local search is used and a whole evaluation is not necessary

  return distances[from][to];

}


/****************************************************************/
/* Returns true when a specific node is a charging station;     */
/* and false otherwise                                          */
/****************************************************************/
bool is_charging_station(int node){

  bool flag = false; 
  if(charging_station[node] == true)
    flag = true; 
  else 
    flag = false;
  return flag;
  
}


/****************************************************************/
/* Returns the demand for a specific customer                   */
/* points: from and to.                                         */
/****************************************************************/
int get_customer_demand(int customer){

  return init_objects[customer].cust_demand;

}


/****************************************************************/
/* Outputs the routes of the solution. Taken as input           */
/* an array of node indeces and its length                      */
/****************************************************************/
void print_solution(int *routes, int size, bool f) {
  int i;

  for( i = 0 ; i < size; i++ ) {
    cout << routes[i] <<  " , ";
  }  

}


/****************************************************************/
/* Validates the routes of the solution. Taken as input         */
/* an array of node indeces and its length                      */
/****************************************************************/
/*void check_solution(int *t, int size, bool f){
  int i, from, to;
  double energy_temp = battery_capacity; 
  double capacity_temp = max_capacity;
  double distance_temp = 0.0;

  for(i = 0; i < size-1; i++){
    from = t[i];
    to = t[i+1];
    capacity_temp -= get_customer_demand(to);
    energy_temp -= get_energy_consumption(from,to);
    distance_temp += get_distance(from,to);

    if(capacity_temp < 0.0) {
      cout << "error: capacity below 0 at customer " << to <<  endl;
      print_solution(t,size,f);
      exit(1);
    }
    if(energy_temp < 0.0) {
       cout << "error: energy below 0 from " << from << " to " << to <<  endl;
       print_solution(t,size,f);
       exit(1);
    }
    if(to == depot) {
      capacity_temp = max_capacity;
    }
    if(is_charging_station(to)==true || to==depot){
      energy_temp = battery_capacity;
    }
  }
  if(distance_temp != fitness_evaluation(t,size,f)) {
    cout << "error: check fitness evaluation" << endl;
  }
}*/

bool check_solution(int *t, int size, bool f){
  int i, from, to;
  double energy_temp = battery_capacity; 
  double capacity_temp = max_capacity;
  double distance_temp = 0.0;

  bool flag = true;

  for(i = 0; i < size-1; i++){
    from = t[i];
    to = t[i+1];
    capacity_temp -= get_customer_demand(to);
    energy_temp -= get_energy_consumption(from,to);
    distance_temp += get_distance(from,to);

    if(capacity_temp < 0.0) {
        cout << "error: capacity below 0 at customer " << to <<  endl;
        flag = false;
        break;
    }
    if(energy_temp < 0.0) {
       cout << "error: energy below 0 from " << from << " to " << to <<  endl;
        flag = false;
        break;
    }
    if(to == depot) {
      capacity_temp = max_capacity;
    }
    if(is_charging_station(to)==true || to==depot){
      energy_temp = battery_capacity;
    }
  }
  if(distance_temp != fitness_evaluation(t,size,f)) {
    cout << "error: check fitness evaluation" << endl;
    flag = false;
  }

  return flag;
}


void free_EVRP(){

  int i;

  delete[] init_objects;
  delete[] charging_station;

  for(i = 0; i < problem_size+number_of_charging_stations; i++) {
    delete[] distances[i];
    delete[] modify_factors[i];
  }

  delete[] distances;
  delete[] modify_factors;


}
