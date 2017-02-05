//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef GALGO_H
#define GALGO_H

//================================================================================================= 

#include <algorithm>
#include <bitset>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <climits>
#include <cmath>
#include <cstring>

/*-------------------------------------------------------------------------------------------------*/

namespace galgo {

// forward declarations
template <typename T>
class BaseParameter;

template <typename T, int N = 16>
class Parameter;

template <typename T>
class GeneticAlgorithm;

template <typename T>
class Population;

template <typename T>
class Chromosome;

// convenient typedefs
template <typename T>
using CHR = std::shared_ptr<Chromosome<T>>;

template <typename T>
using PAR = std::unique_ptr<BaseParameter<T>>;

template <typename T, int...N>
using TUP = std::tuple<const Parameter<T,N>&...>;

}

/*-------------------------------------------------------------------------------------------------*/

#ifdef _OPENMP 
  #include <omp.h>
  // getting maximum number of threads available
  static const int MAX_THREADS = omp_get_max_threads();
#endif

/*-------------------------------------------------------------------------------------------------*/

#include "Randomize.hpp"
#include "Converter.hpp"
#include "Parameter.hpp"
#include "Evolution.hpp"
#include "Chromosome.hpp"
#include "Population.hpp"
#include "GeneticAlgorithm.hpp"

//================================================================================================= 

#endif

