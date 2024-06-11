#define IA       16807
#define IM       2147483647
#define AM       (1.0/IM)
#define IQ       127773
#define IR       2836
#define MASK     123459876
#define EPSILON  0.000000000000000000000001

//Ant or individual structure that represent the TSP tour and tour cost
struct ant{
  int *tour;
  int *routes;      //shows the vehicles that cover the customers
  bool *visited;
  int id;
  double tour_length;
  int capacity;     //vehicles current capacity
  double service;
  int dummy_depots; //No. of vehicles
  bool start_new;   //start new route
  int no_customers; //No. of customers
  int step;
  int no_recharges; //number of visits in a charging station
  double energy_level; //vehicles currrent energy
  double* space_available;
  double *acc_space; 
  bool feasible;
};

extern int seed;
extern int alg_mode;
//extern int seed;
extern int n_ants;
extern int current_iteration;
extern double rho;
extern double q_0;
extern int **nn_list; 
extern int pseudonodes;
extern int nn;
extern ant *best_so_far_ant;
extern bool flag_change_detected;

void set_algorithm_parameters(void);
void allocate_ants(void);
void allocate_structures(void);
bool termination_condition(void);
void ACO();
void construct_solutions(void);
void pheromone_update();
void update_best(void);
void statistics_and_output();
void init_try(int r);
void free_ACO();
void local_search();
void init_heuristic_info();
void compute_nn_lists();
void printTour(ant *a);
void checkTour(ant *a);
void printTest(ant *a);
double alg_random_number(int *idum);
void sort(int v[], int v2[], int left, int right);
void printTrail(void);
void print_nn_list(void);
