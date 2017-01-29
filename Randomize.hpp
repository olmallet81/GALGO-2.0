#ifndef RANDOMIZE_H
#define RANDOMIZE_H

namespace galgo {

//=================================================================================================

// template metaprogramming to get the number of different unsigned integers obtained from n + 1 bits
template <unsigned int n>
struct NVALUES 
{
    enum : uint64_t{value = 2 * NVALUES<n - 1>::value};
};

// template specialization for initial case n = 0 (1 bit leads to 2 possible values 0 and 1)
template <>
struct NVALUES<0> 
{
    enum {value = 2};
};

/*-------------------------------------------------------------------------------------------------*/

// Mersenne Twister 19937 pseudo-random number generator
std::random_device rand_dev;
std::mt19937_64 rng(rand_dev());
//std::mt19937_64 rng(0);

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
