# GALGO
Genetic Algorithm in C++ STL with lower and upper bounds for constrained problems optimization.

# Description
GALGO is a C++ template library, headers only, designed to solve a problem under constraints (or not) by maximizing or minimizing an objective function on given boundaries. GALGO can also achieve multi-objective optimization. It does not use any external C++ library, only the Standard Template Library. GALGO is fast and can use parallelism when required through OpenMP. GALGO is flexible and has been written in a way allowing the user to easily add new methods to the genetic algorithm. This library already contains some methods for selection, cross-over and mutation among the most widely used. The user can choose among these pre-existing methods or create new ones.

# Template class GeneticAlgorithm

This is the main class you need to instantiate to run a genetic algorithm.

- ## Constructor
   ```C++
   GeneticAlgorithm(Functor<T> Objective, int popsize, const std::vector<T>& lowerBound, const std::vector<T>& upperBound, int nbgen, bool output = false)
   ```
With:
   - *Objective* = objective function functor 
   - *popsize* = population size or number of chromosomes
   - *lowerBound* = vector containing the parameter(s) lower bound
   - *upperBound* = vector containing the parameter(s) upper bound
   - *nbgen* = number of generations to run
   - *output* = control for outputting results (set to false by default)

- ## Member variables (public)
   - *lowerBound* = vector containing the parameter(s) lower bound
   - *upperBound* = vector containing the parameter(s) upper bound
   - *initialSet* = vector containing the parameter(s) starting point (empty by default)
   - *covrate* = cross-over rate between 0 and 1 (set to 0.5 by default)
   - *mutrate* = mutation rate between 0 and 1 (set to 0.05 by default)
   - *SP* = selective pressure, for RSP selection method only, between 1 and 2 (set to 1.5 by default)
   - *tolerance* = terminal condition to stop the algorithm (inactive by default, set to 0)
   - *elitpop* = elit population size, must not be greater than population size (set to 1 by default)
   - *matsize* = mating population size (set to population size by default)
   - *tntsize* = tournament size, for TNT selection method only (set to 10 by default)
   - *genstep* = generation step for outputting results (set to 10 by default)
   - *precision* = number of decimals for outputting results (set to 5 by default)
   
- ## Member functions (public)
   - *run()* for running the genetic algorithm
   - *result()* for getting the population best chromosome

# Encoding and Decoding Chromosomes
GALGO is based on chromosomes represented as a binary string of 0 and 1 containing the encoded parameters to be estimated. The user is free to choose the number of bits N to encode the parameters (or genes) composing a chromosome within the interval [1,64]. When initializing a population of chromosomes, a random 64 bits unsigned integer, we will call it X, will be generated for each parameter to be estimated, X being inside the interval [0,MAXVAL] where MAXVAL is the greatest unsigned integer obtained for the chosen number of bits. If the chosen number of bits to represent a gene is N, we will have:
```
MAXVAL = 2^N - 1
```
X will be then converted into a binary string of 0 and 1, the binary string will be truncated to the desired number of bits and added to the new chromosome. Once selection, cross-over and mutation have been applied, the new binary string obtained is decoded to get the new estimated parameter value Y. To do so, the binary string is first converted back into an unsigned integer X and the following equation is applied, with [minY,maxY] representing the parameter boundaries: 
```
Y = minY + (X / MAXVAL) * (maxY - minY)
```
This method of generating a random ratio rather than a random real number allows to achieve faster convergence as only values inside the boundaries [minY,maxY] will be generated when initializing the chromosome population but also when recombining and mutating them.
The step between two successive parameter values within the boundaries will be:
```
step = (maxY - minY) / MAXVAL
```
Estimating this step value before running the genetic algorithm is of high importance as it will greatly influence its performance. If the chosen number of bits N is large the algorithm will struggle to achieve fast convergence as the number of possible solutions within the boundaries will be too great, on the contrary, if it is too small the algorithm will struggle to find the global extremum and risk to quickly stall on a local extremum due to the lack of diversity within the chromosome population.

The functions to encode and decode the chromosomes are in the header Converter.hpp.


# Evolution

The pre-existing methods to evolve a chromosome population contained in the header Evolution.hpp are:

## Selection methods
- proportional roulette wheel selection (RWS)
- stochastic universal sampling (SUS)
- classic linear rank-based selection (RNK)
- linear rank-based selection with selective pressure (RSP)
- tournament selection (TNT)
- transform ranking selection (TRS)

## Cross-Over methods
- one point cross-over (P1XO)
- two point cross-over (P2XO)
- uniform cross-over (UXO)

## Mutation methods
- boundary mutation (BDM)
- single point mutation (SPM)
- uniform mutation (UNM)

GALGO also contains a default method for adaptation to constraints (DAC).

By default GALGO is set to run without constraints and with RWS, P1XO and SPM.


# Example

```C++
#include "Galgo.hpp"

// objective class example
template <typename T>
class MyObjective
{
public:
   // objective function example : Rosenbrock function
   // minimizing f(x,y) = (1 - x)^2 + 100 * (y - x^2)^2
   static std::vector<T> Objective(const std::vector<T>& x)
   {
      T obj = -(pow(1-x[0],2)+100*pow(x[1]-x[0]*x[0],2));
      return {obj};
   }
   // NB: GALGO maximize by default so we will maximize -f(x,y)
};

// constraints example:
// 1) x * y + x - y + 1.5 <= 0
// 2) 10 - x * y <= 0
template <typename T>
std::vector<T> MyConstraint(const std::vector<T>& x)
{
   return {x[0]*x[1]+x[0]-x[1]+1.5, 10-x[0]*x[1]};
}


int main()
{
   // defining lower bound LB and upper bound UB
   std::vector<double> LB({0.0,0.0});
   std::vector<double> UB({1.0,13.0});

   // initiliazing genetic algorithm
   galgo::GeneticAlgorithm<double> ga(MyObjective<double>::Objective,100,LB,UB,50,true);
 
   // setting constraints
   ga.Constraint = MyConstraint;

   // running genetic algorithm
   ga.run();
}
```

You can run the example above contained in the source file example.cpp by first compiling with the following command:
```
$ g++ -std=c++11 -O3 -Wall example.cpp -o run
```
and then running:
```
$ ./run
```

In this example we have constructed a class called MyObjective containing the function to optimize, this does not have to be necessarily the case if you do not need a complex objective function needing more arguments than the vector of parameters only.
