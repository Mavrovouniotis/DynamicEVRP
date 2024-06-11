To generate dynamic cases of the EVRP the C++ implementation of the Dynamic Benchmark Framework is used published in the following paper:

M. Mavrovouniotis, S. Yang, M. Van, C. Li, and M. Polycarpou. Ant colony optimization algorithms for dynamic 
optimization: A case study of the dynamic travelling salesperson problem. ``IEEE Computational Intelligence 
Magazine'', 15(1), pp. 52-63, Feb. 2019


An example on how to execute the environment: 

./main X-n143-k7.evrp 0.05 14300 10 1

main = executable 
X-n143-k7 = benchmark EVRP instance
0.05 = magntitude of change
14300 = frequency of change (iterations)
10 = total number of dynamic changes
1 = ACO variant

CONTENTS
=======
Implementation of ACO frameworks and variations
ACO.cpp
ACO.hpp

The code of the ACO algorithms is based on the ACOTSP implementation 
of Thomas Stuetzle: ACO algorithms for the TSP, version 1.03 
www.aco-metaheuristic.org/aco-code


Implementation of the dynamic benchmark framework
EVRP.cpp 
EVRP.hpp

Implementation of performance and behaviour measurements
stats.cpp
stats.hpp

The main method of the implementation
main.cpp


Contact
=======
Michalis Mavrovouniotis
m.mavrovouniotis(at)hotmail.com
Please let me know for any bugs/improvements/additions.
