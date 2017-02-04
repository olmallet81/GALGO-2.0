//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef RANDOMIZE_H
#define RANDOMIZE_H

namespace galgo {

//=================================================================================================

// template metaprogramming for getting maximum unsigned integral value from N bits
template <unsigned int N>
struct MAXVALUE 
{
    enum : uint64_t{ value = 2 * MAXVALUE<N - 1>::value };
};

// template specialization for initial case N = 0
template <>
struct MAXVALUE<0> 
{
    enum { value = 1 };
}; 

/*-------------------------------------------------------------------------------------------------*/

// Mersenne Twister 19937 pseudo-random number generator
std::random_device rand_dev;
std::mt19937_64 rng(rand_dev());

// generate uniform random probability in range [0,1)
std::uniform_real_distribution<> proba(0, 1);

/*-------------------------------------------------------------------------------------------------*/

// generate a uniform random number within the interval [min,max)
template <typename T>
inline T uniform(T min, T max)
{   
    #ifndef NDEBUG
    if (min >= max) {
      throw std::invalid_argument("Error: in galgo::uniform(T, T), first argument must be < to second argument.");
    }
    #endif
 
    return min + proba(rng) * (max - min);
}

//=================================================================================================

}

#endif
