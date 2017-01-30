# GALGO
Genetic Algorithm in C++ STL with lower and upper bounds for constrained problems optimization.

# Description
GALGO is a C++ template library, headers only, designed to solve a problem under constraints (or not) by maximizing or minimizing an objective function on given boundaries. It does not use any external C++ library, only the Standard Template Library. GALGO is fast and can use parallelism when required through OpenMP. GALGO is flexible and has been written in a way to allow the user to easily add new functionalities to the genetic algorithm. This library already contains some methods for selection, cross-over and mutation among the most widely used. The user can choose among these pre-existing methods to build a genetic algorithm or create new ones.

## Encoding and Decoding Chromosomes
This genetic algorithm is using chromosomes containing the binary representation of an unsigned integer represented as a string of 0 and 1. In order to achieve faster convergence, only values inside the boundaries will be generated when initializing the chromosome population but also when recombining and mutating them. To do so, a random ratio between 0 and 1 is generated, the denominator being the greatest unsigned integer obtained for the chosen number of bits, we will call it MAXVAL, and the numerator a random unsigned integer between 0 and MAXVAL, we will call it X. Hence, to decode a chromosome, the binary string is converted back to an unsigned integer and the value of the parameter Y inside the boundaries [minY,maxY] will be obtained using the following equation:
```
Y = minY + (X / MAXVAL) * (maxY - minY)
```

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


