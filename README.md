# GALGO-2.0
Genetic Algorithm in C++ with template metaprogramming and abstraction for constrained optimization.

# Description
GALGO is a C++ template library, header only, designed to solve a problem under constraints (or not) by maximizing or minimizing an objective function on given boundaries. GALGO can also achieve multi-objective optimization. It does not use any external C++ library, only the Standard Template Library. GALGO is fast and can use parallelism when required through OpenMP. GALGO is flexible and has been written in a way allowing the user to easily add new methods to the genetic algorithm (within header file Evolution.hpp). This library already contains some methods for selection, cross-over and mutation among the most widely used. The user can choose among these pre-existing methods or create new ones. The new version of this library (2.0) is using metaprogramming with variadic templates and abstraction to allow optimization on an arbitrary number of parameters encoded using different number of bits.

# Encoding and Decoding Chromosomes
GALGO is based on chromosomes represented as a binary string of 0 and 1 containing the encoded parameters to be estimated. The user is free to choose the number of bits N to encode each one of them within the interval [1,64]. In the previous version of this library, this number had to be the same for all parameters to be estimated, in the new version 2.0 they can be encoded using a different number of bits. When initializing a population of chromosomes, a random 64 bits unsigned integer, we will call it X, will be generated for each parameter to be estimated, X being inside the interval [0,MAXVAL] where MAXVAL is the greatest unsigned integer obtained for the chosen number of bits. If the chosen number of bits to represent a gene is N, we will have:
```
MAXVAL = 2^N - 1
```
X will be then converted into a binary string of 0 and 1, the binary string will be truncated to the desired number of bits and added to the new chromosome. Once selection, cross-over and mutation have been applied, the new binary string obtained is decoded to get the new estimated parameter value Y. To do so, the binary string is first converted back into an unsigned integer X and the following equation is applied, with [minY,maxY] representing the parameter boundaries:
```
Y = minY + (X / MAXVAL) * (maxY - minY)
```
This method of generating a random ratio rather than a random real number allows to achieve faster convergence as only values inside the boundaries [minY,maxY] will be generated when initializing the chromosome population but also when recombining and mutating them.
The step between two consecutive parameter values within the boundaries will be:
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

GALGO also contains a default method for adaptation to constraint(s) (DAC).

By default GALGO is set to run with no constraint and with RWS, P1XO and SPM.

# Design

## C++ template class *Parameter*

This class used to initialize the parameter(s) to be estimated by providing a lower bound, an upper bound and an initial value if required (optional), is declared and defined within the header file Parameter.hpp.

```C++
namespace galgo {
   template <typename T,int N = 16>
   class Parameter;
}
```
The template parameter T can be either float or double for the precision of the solution returned. The template parameter N corresponds to the number of bits to encode the parameter, it must be between 1 and 64. This class inherits from an abstract base class called *BaseParameter*, declared and defined within the same header file, for storing the parameter(s) inside the same container even if they have been declared using a different template parameter N.

### Constructor
```C++
template <typename T, int N>
Parameter(const std::vector<T>& data);
```
### Member variable (public)
   - *data* = std::vector containing the parameter lower and upper bounds and an initial value if required


## C++ template class *GeneticAlgorithm*

This is the class you need to instantiate to run a genetic algorithm, declared and defined within the namespace galgo in the header file GeneticAlgorithm.hpp.  

```C++
namespace galgo {
   template <typename T>
   class GeneticAlgorithm;
}
```
The template parameter T can be either float or double for the precision of the solution returned.

### Constructor
```C++
template <typename T> template<int...N>
GeneticAlgorithm(Functor<T> Objective,int popsize,int nbgen,bool output,const Parameter<T,N>&...args);
```
With:
   - *objective* = objective function (function to optimize)
   - *popsize* = population size or number of chromosomes
   - *nbgen* = number of generations to run
   - *output* = control for outputting results
   - *args* = parameter(s) to be estimated (the parameter pack allows to instantiate this class using an arbitrary number of objects of type *Parameter*)
   
### Member function pointers (public)
   - *Selection* = for setting the selection method (set to RWS by default)
   - *CrossOver* = for setting the cross-over method (set to P1XO by default)
   - *Mutation* = for setting the mutation method (set to SPM by default)
   - *Adaptation* = for setting the adaptation to constraint(s) (optional)
   - *Constraint* = for setting the constraint(s) function (optional)
  
### Member variables (public)
   - *covrate* = cross-over rate between 0 and 1 (set to 0.5 by default)
   - *mutrate* = mutation rate between 0 and 1 (set to 0.05 by default)
   - *SP* = selective pressure, for RSP selection method only, between 1 and 2 (set to 1.5 by default)
   - *tolerance* = terminal condition to stop the algorithm (inactive by default, set to 0)
   - *elitpop* = elite population size, must not be greater than population size (set to 1 by default)
   - *matsize* = mating population size (set to population size by default)
   - *tntsize* = tournament size, for TNT selection method only (set to 10 by default)
   - *genstep* = generation step for outputting results (set to 10 by default)
   - *precision* = number of decimals for outputting results (set to 5 by default)
   
### Member functions (public)
   - *run()* for running the genetic algorithm
   - *result()* for getting the population best chromosome

NB: in the previous version the user had access to the parameter(s) lower bound, upper bound and initial set, they are now private inside *GeneticAlgorithm* class as these values are set when declaring an object using *Parameter* class. They can be modified at any time once the parameter(s) declared as they are contained inside a public vector.

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
   return {x[0]*x[1]+x[0]-x[1]+1.5,10-x[0]*x[1]};
}
// NB: a penalty will be applied if one of the constraints is > 0
// using the default adaptation to constraint(s) method

int main()
{
   // initializing parameters lower and upper bounds
   galgo::Parameter<double> par1({0.0,1.0});
   galgo::Parameter<double> par2({0.0,13.0});
   // an initial value can be added inside the initializer list after the upper bound
   // both parameter will be encoded using 16 bits (default value)
   // this value can be modified but has to remain between 1 and 64
   // example using 20 and 32 bits and initial values of 0.5 and 1.0:
   // galgo::Parameter<double,20> par1({0.0,1.0,0.5});
   // galgo::Parameter<double,32> par2({0.0,13.0,1.0});

   // initiliazing genetic algorithm
   galgo::GeneticAlgorithm<double> ga(MyObjective<double>::Objective,100,50,true,par1,par2);

   // setting constraints
   ga.Constraint = MyConstraint;

   // running genetic algorithm
   ga.run();
}
```

## Compilation

You can run the example above contained in the source file example.cpp by first compiling with the following command:
```
$ g++ -std=c++11 -O3 -Wall example.cpp -o run
```
and then running:
```
$ ./run
```
In this example we have constructed a class called MyObjective containing the function to optimize, this does not have to be necessarily the case if you do not need a complex objective function needing more arguments than the vector of parameters only.

If the objective function is time consuming you can go parallel by compiling with OpenMP enabled:
```
$ g++ -fopenmp -std=c++11 -O3 -Wall example.cpp -o run
```
GALGO is set to use the maximum number of threads available by default when OpenMP is enabled, you can reduce this value inside the header Galgo.hpp.

## Ouput

```
 Running Genetic Algorithm...
 ----------------------------
 Generation =  0 | X1 =   0.22918 | X2 =  12.61854 | F(x) = -15791.07606
 Generation = 10 | X1 =   0.81093 | X2 =  12.39359 | F(x) = -13773.38448
 Generation = 20 | X1 =   0.81209 | X2 =  12.38605 | F(x) = -13751.28199
 Generation = 30 | X1 =   0.81201 | X2 =  12.33666 | F(x) = -13635.97220
 Generation = 40 | X1 =   0.81212 | X2 =  12.31742 | F(x) = -13590.66665
 Generation = 50 | X1 =   0.81224 | X2 =  12.31583 | F(x) = -13586.50453

 Constraint(s)
 -------------
 C1(x) = -0.00021
 C2(x) = -0.00338
```

NB: the random number generator being randomly seeded you will not get exactly the same results.

# References

For those who want to know more about the methods I have implemented for selection and adaptation to constraints, please find below the links to the documents I have used:
- [Transform Ranking: a New Method of Fitness Scaling in Genetic Algorithms](http://shura.shu.ac.uk/5638/1/AI2008.pdf)
- [Constraint handling strategies in Genetic Algorithms](https://core.ac.uk/download/pdf/12039642.pdf?repositoryId=437)
