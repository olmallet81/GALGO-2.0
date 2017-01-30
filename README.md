# GALGO
Genetic Algorithm in C++ STL with lower and upper bounds for constrained problems optimization.

# Description
GALGO is a C++ template library, headers only, designed to solve a problem under constraints (or not) by maximizing or minimizing an objective function on given boundaries. It does not use any external C++ library, only the Standard Template Library. GALGO is fast and can use parallelism when required through OpenMP. GALGO is flexible and has been written in a way to allow the user to easily add new functionalities to the genetic algorithm. This library already contains some methods for selection, cross-over and mutation among the most widely used. The user can choose among these pre-existing methods to build a genetic algorithm or create new ones.

## Encoding and Decoding Chromosomes
This genetic algorithm is using chromosomes represented as a binary string of 0 and 1. When initializing a population of chromosome, a random unsigned integer, we will call it X, will be generated for each parameter to be estimated, X being inside the interval [0,MAXVAL] where MAXVAL is the greatest unsigned integer obtained for the chosen number of bits. If the number of bits is N, we will have:
```
MAXVAL = 2^N - 1
```
X will be then converted into a binary string of 0 and 1 and added to the new chromosome. Once selection, cross-over and mutation have been applied, the new binary string obtained is decoded to get the new estimated parameter value Y. To do so, the binary string is first converted back into an unsigned integer X and the follozing equation is applied: 
```
Y = minY + (X / MAXVAL) * (maxY - minY)
```
Where minY is the parameter lower bound and maxY its upper bound.
This method allows to achieve faster convergence as only values inside the boundaries [minY,maxY] will be generated when initializing the chromosome population but also when recombining and mutating them.
The precision of the solution will be:
```
precision = (maxY - minY) / MAXVAL
```
If the number of bits is large the algorithm will struggle to achieve convergence as the number of possible solutions will be greater, if it is too small the algorithm will struggle to find the global extremum and risk to quickly stall on a local extremum as the number of possible solutions will be smaller.


The pre-existing methods are:

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


