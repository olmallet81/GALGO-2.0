# GALGO
Genetic Algorithm in C++ STL with lower and upper bounds for constrained optimization.


# Description
GALGO is a C++ template library, headers only, designed to solve a problem under constraints (or not) by maximizing or minimizing an objective function on given boundaries. GALGO can also achieve multi-objective optimization. It does not use any external C++ library, only the Standard Template Library. GALGO is fast and can use parallelism when required through OpenMP. GALGO is flexible and has been written in a way allowing the user to easily add new methods to the genetic algorithm. This library already contains some methods for selection, cross-over and mutation among the most widely used. The user can choose among these pre-existing methods or create new ones.


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
The step between two consecutive parameter values within the boundaries will be:
```
step = (maxY - minY) / MAXVAL
```
Estimating this step value before running the genetic algorithm is of high importance as it will greatly influence its performance. If the chosen number of bits N is large the algorithm will struggle to achieve fast convergence as the number of possible solutions within the boundaries will be too great, on the contrary, if it is too small the algorithm will struggle to find the global extremum and risk to quickly stall on a local extremum due to the lack of diversity within the chromosome population.

The functions to encode and decode the chromosomes are in the header Converter.hpp.


# Evolution

The pre-existing methods to evolve a chromosome population contained in the header Evolution.hpp are:

## Selection methods
