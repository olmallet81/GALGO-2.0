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
#include <type_traits>
#include <vector>

#include <climits>
#include <cmath>
#include <cstring>

/*-------------------------------------------------------------------------------------------------*/

namespace galgo {

template <typename T, int N>
class GeneticAlgorithm;

template <typename T, int N>
class Population;

template <typename T, int N>
class Chromosome;

// convenient typedef
template <typename T, int N>
using CHR = std::shared_ptr<Chromosome<T,N>>;

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
#include "Evolution.hpp"
#include "Chromosome.hpp"
#include "Population.hpp"
#include "GeneticAlgorithm.hpp"

//================================================================================================= 

#endif

