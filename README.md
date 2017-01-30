# GALGO
Genetic Algorithm in C++ STL with lower and upper bounds for constrained problems optimization.

# Description
GALGO is a C++ template library, headers only, designed to solve a problem under constraints (or not) by maximizing or minimizing an objective function on given boundaries. It does not use any external C++ library, only the Standard Template Library. GALGO is fast and can use parallelism when required through OpenMP. GALGO is flexible and has been written in a way to allow the user to easily add new functionalities to the genetic algorithm. This library already contains some methods for selection, cross-over and mutation among the most widely used. The user can choose among these pre-existing methods to build a genetic algorithm or create new ones.

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
