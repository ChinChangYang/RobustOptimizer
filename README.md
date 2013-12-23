RobustOptimizer
===============

The most robust algorithm for global optimization problems

Goal
----

* Single-objective real-parameter optimization
* Multi-objective real-parameter optimization
* Min-max real-parameter optimization
* Noisy optimization
* Mixed-integer optimization
* Linear/Nonlinear inequality/equality constraint handling
* The most robust optimizer which can satisfy the above goals

Description
-----------

Basically, the aim of this project is to find the best algorithm (solver) that costs the least function evaluations at finding a minimizer of any functions. Currently, the candidate solvers are 

* Covariance Matrix Adaptation Evolution Strategy and 
* Differential Evolution. 

A robustness measure is developed to determine the most best solver. Currently, these solvers are tested on a set of 28 and 24 benchmark functions in CEC'13 and BBOB'12, respectively, including 

* Unimodal functions,  
* Multimodal functions, and
* Composition functions. 

The solvers are tested over various maximal function evaluations to see their robustness on the balance between expoitation and exploration power. Furthermore, the solvers are compared under the empirical cumulative distribution function values, in which we can see their error differences to the best solver for each test case. To simplify the results, number of win cases of each solver is also adopted to evaluate its robustness.

Feature in 4th Milestone
------------------------

### Algorithm
* CMA-ES
* DE/best/1/bin
* DE/rand/1/bin
* DEGL
* SaDE
* JADE
* Drift-free DE (DFDE)
* D-CMA-EA
* Rank-based DE (RBDE)
* Novel eigenvector-based crossover operator
 
### Benchmark
* CEC 2005 benchmark
* CEC 2011 benchmark
* BBOB 2012 benchmark
* CEC 2013 benchmark
 
### Measure
* Solution error
* Wilcoxon's rank-sum test
* Empirical cumulative probability distribution

Development
-----------

* Matlab(R) with Statistics Toolbox(TM) and Global Optimization Toolbox(TM)

License
-------

* RobustOptimizer is available under BSD license.
