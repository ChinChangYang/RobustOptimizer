RobustOptimizer
===============

The most robust algorithm for global optimization problems

Goal
----

* Single-objective real-parameter optimization (main)
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

A robustness measure is developed to determine the most best solver. For single-objective real-parameter optimization, these solvers are currently tested on a set of 30 benchmark functions in CEC'14, including 

* Unimodal functions,  
* Multimodal functions, and
* Composition functions. 

The solvers are tested over various maximal function evaluations to see their robustness on the balance between expoitation and exploration power. Furthermore, the solvers are compared under the empirical cumulative distribution function values, in which we can see their error differences to the best solver for each test case. To simplify the results, number of win cases of each solver is also adopted to evaluate its robustness.

Feature
------------------------

### Algorithm
* DE/best/1/bin
* DE/rand/1/bin
* DEGL
* SaDE
* JADE
* Drift-free DE (DFDE)
* D-CMA-EA
* Rank-based DE (RBDE)
* SHADE
* L-SHADE
* DEsPA
* Eigenvector-based crossover operator
* Successful-parent-selecting framework
 
### Benchmark
* CEC 2005 benchmark
* CEC 2011 benchmark
* BBOB 2012 benchmark
* CEC 2013 benchmark
* CEC 2014 benchmark
* CEC 2015 benchmark
 
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
