#define INFTY		INT_MAX
#define CHAR_LEN 	100
#define TERMINATION (total_changes*change_speed) + (problem_size*100)

struct object {
  int id;
  double x;
  double y;
  double cust_demand;
};


typedef enum type_timer {REAL, VIRTUAL} TIMER_TYPE;
//USER INPUT PARAMETERS
extern char* problem_instance;          //Name of the instance
extern double global_optimum_value;     //Global optimum of the problem instance (if known)
extern double energy_consumption;
extern int battery_capacity;
extern object *init_objects;


//METHODS AND PARAMETERS THAT CAN BE USED IN YOUR ALGORITHM IMPLEMENTATION
extern double **distances;              //Distance matrix
extern int problem_size;                //Size of the instance
extern int problem_size_no_depot;
extern int number_of_charging_stations;
extern bool* charging_station;
extern int max_capacity;
extern int depot;
extern int evals;

extern double change_degree;
extern int change_speed;
extern int total_changes;
extern int max_size;
extern int env_index;

extern double current_best;
extern double offline_performance;
extern double offline_error;

double** generate_2D_matrix_double(int n, int m);
int** generate_2D_matrix_int(int n, int m);
void start_timers(void);
double elapsed_time(TIMER_TYPE type);
void read_problem(char* filename);
void compute_distances(void);
void initialize_environment();
void change_environment();
double fitness_evaluation(int *t, int s, bool f);
double distance_output(int *t,int *id, int recharges, int customers);
double dummy_fitness_evaluation(int *t, int s);
double get_offline_performance();
double get_current_best();
double get_best_iteration();
double get_best_time();
double get_offline_error();
double get_current_error();
int get_evals();
double get_time();
void free_EVRP();
void update_tour_length(int length);
double get_energy_consumption(int from, int to);
bool is_charging_station(int node);
//void check_solution(int *t, int size, bool f);
bool check_solution(int *t, int size, bool f);

