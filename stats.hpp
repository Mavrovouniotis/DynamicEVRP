void open_stats(void);
void close_stats(void);
double mean(double* values, int size);
double stdev(double* values, int size, double average);

/*Performance measurements*/
void get_observation(int t);
/*Behaviour mesaurements */
double calc_diversity(int size);
double calc_lambda(double l, int size, int depth, double** pheromone);
double calc_entropy(int size, int depth, double** pheromone, double** heuristic, double a, double b);
/*Calculate averages of all trials*/
void get_mean(int r, double value);

extern int max_trials;
extern int max_iterations;
