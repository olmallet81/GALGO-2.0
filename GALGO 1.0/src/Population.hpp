//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef POPULATION_HPP
#define POPULATION_HPP

namespace galgo {

//=================================================================================================

template <typename T, int N = 16>
class Population
{
   static_assert(std::is_same<float,T>::value || std::is_same<double,T>::value, "variable type can only be float or double, please amend.");
   static_assert(N > 0 && N <= 64, "number of bits cannot be ouside interval [1,64], please choose an integer within this interval.");

public: 
   // nullary constructor
   Population() {}
   // constructor
   Population(const GeneticAlgorithm<T,N>& ga);
   // create a population of chromosomes
   void creation();
   // evolve population, get next generation
   void evolution();

   // access element in current population at position pos
   const CHR<T,N>& operator()(int pos) const;
   // access element in mating population at position pos
   const CHR<T,N>& operator[](int pos) const;  
   // return iterator to current population beginning
   typename std::vector<CHR<T,N>>::iterator begin();
   // return const iterator to current population beginning
   typename std::vector<CHR<T,N>>::const_iterator cbegin() const;
   // return iterator to current population ending  
   typename std::vector<CHR<T,N>>::iterator end();  
   // return const iterator to current population ending  
   typename std::vector<CHR<T,N>>::const_iterator cend() const;  
   // select element at position pos in current population and copy it into mating population
   void select(int pos);
   // set all fitness to positive values 
   void adjustFitness();
   // compute fitness sum of current population
   T getSumFitness() const;
   // get worst objective function total result from current population
   T getWorstTotal() const;
   // return population size
   int popsize() const;
   // return mating population size
   int matsize() const;
   // return tournament size
   int tntsize() const;
   // return numero of generation
   int nogen() const;
   // return number of generations
   int nbgen() const;
   // return selection pressure
   T SP() const; 

private:
   std::vector<CHR<T,N>> curpop;               // current population
   std::vector<CHR<T,N>> matpop;               // mating population
   std::vector<CHR<T,N>> newpop;               // new population

   const GeneticAlgorithm<T,N>* ptr = nullptr; // pointer to genetic algorithm              
   int nbrcrov;                                // number of cross-over
   int matidx;                                 // mating population index

   // elitism => saving best chromosomes in new population
   void elitism();
   // create new population from recombination of the old one
   void recombination();
   // complete new population randomly
   void completion();
   // update population (adapting, sorting)
   void updating();
};

/*-------------------------------------------------------------------------------------------------*/

// constructor
template <typename T, int N>
Population<T,N>::Population(const GeneticAlgorithm<T,N>& ga)
{
   ptr = &ga;   
   nbrcrov = floor(ga.covrate * (ga.popsize - ga.elitpop));
   // adjusting nbrcrov (must be an even number)
   if (nbrcrov % 2 != 0) nbrcrov -= 1;
   // for convenience, we add elitpop to nbrcrov
   nbrcrov += ga.elitpop;
   // allocating memory
   curpop.resize(ga.popsize);
   matpop.resize(ga.matsize);
}

/*-------------------------------------------------------------------------------------------------*/

// create a population of chromosomes
template <typename T, int N>
void Population<T,N>::creation()
{
   int start = 0;
   // initializing first chromosome
   if (!ptr->initialSet.empty()) {
      curpop[0] = std::make_shared<Chromosome<T,N>>(*ptr);
      curpop[0]->initialize();
      curpop[0]->evaluate();
      start++;
   }
   // getting the rest
   #ifdef _OPENMP 
   #pragma omp parallel for num_threads(MAX_THREADS)
   #endif
   for (int i = start; i < ptr->popsize; ++i) {
      curpop[i] = std::make_shared<Chromosome<T,N>>(*ptr);
      curpop[i]->create();
      curpop[i]->evaluate();
   }
   // updating population
   this->updating();
}

/*-------------------------------------------------------------------------------------------------*/

// population evolution (selection, recombination, completion, mutation), get next generation
template <typename T, int N>
void Population<T,N>::evolution()
{   
   // initializing mating population index
   matidx = 0;
   // selecting mating population
   ptr->Selection(*this);
   // applying elitism if required
   this->elitism(); 
   // crossing-over mating population
   this->recombination();
   // completing new population
   this->completion();
   // moving new population into current population for next generation
   curpop = std::move(newpop);
   // updating population
   this->updating(); 
}

/*-------------------------------------------------------------------------------------------------*/

// elitism => saving best chromosomes in new population, making a copy of each elit chromosome
template <typename T, int N>
void Population<T,N>::elitism()
{
   // (re)allocating new population
   newpop.resize(ptr->popsize);

   if (ptr->elitpop > 0) {
      // copying elit chromosomes into new population
      std::transform(curpop.cbegin(), curpop.cend(), newpop.begin(), [](const CHR<T,N>& chr)->CHR<T,N>{return std::make_shared<Chromosome<T,N>>(*chr);});
   }
}

/*-------------------------------------------------------------------------------------------------*/

// create new population from recombination of the old one
template <typename T, int N>
void Population<T,N>::recombination()
{
   // creating a new population by cross-over
   #ifdef _OPENMP 
   #pragma omp parallel for num_threads(MAX_THREADS)
   #endif
   for (int i = ptr->elitpop; i < nbrcrov; i = i + 2) {      
      // initializing 2 new chromosome
      newpop[i] = std::make_shared<Chromosome<T,N>>(*ptr);
      newpop[i+1] = std::make_shared<Chromosome<T,N>>(*ptr);
      // crossing-over mating population to create 2 new chromosomes
      ptr->CrossOver(*this, newpop[i], newpop[i+1]);
      // mutating new chromosomes
      ptr->Mutation(newpop[i]);   
      ptr->Mutation(newpop[i+1]);   
      // evaluating new chromosomes
      newpop[i]->evaluate();
      newpop[i+1]->evaluate();
   } 
}

/*-------------------------------------------------------------------------------------------------*/

// complete new population
template <typename T, int N>
void Population<T,N>::completion()
{
   #ifdef _OPENMP 
   #pragma omp parallel for num_threads(MAX_THREADS)
   #endif
   for (int i = nbrcrov; i < ptr->popsize; ++i) {
      // selecting chromosome randomly from mating population
      newpop[i] = std::make_shared<Chromosome<T,N>>(*matpop[uniform<int>(0, ptr->matsize)]);
      // mutating chromosome
      ptr->Mutation(newpop[i]);
      // evaluating chromosome
      newpop[i]->evaluate();
   }
}

/*-------------------------------------------------------------------------------------------------*/

// update population (adapting, sorting)
template <typename T, int N>
void Population<T,N>::updating()
{
   // adapting population to constraints
   if (ptr->Constraint != nullptr) {
      ptr->Adaptation(*this); 
   }
   // sorting chromosomes from best to worst fitness
   std::sort(curpop.begin(),curpop.end(),[](const CHR<T,N>& chr1,const CHR<T,N>& chr2)->bool{return chr1->fitness > chr2->fitness;});
}

/*-------------------------------------------------------------------------------------------------*/

// access element in current population at position pos
template <typename T, int N>
const CHR<T,N>& Population<T,N>::operator()(int pos) const
{
   #ifndef NDEBUG
   if (pos > ptr->popsize - 1) {
      throw std::invalid_argument("Error: in galgo::Population<T>::operator()(int), exceeding current population memory.");
   }
   #endif

   return curpop[pos];
}

/*-------------------------------------------------------------------------------------------------*/

// access element in mating population at position pos
template <typename T, int N>
const CHR<T,N>& Population<T,N>::operator[](int pos) const
{
   #ifndef NDEBUG
   if (pos > ptr->matsize - 1) {
      throw std::invalid_argument("Error: in galgo::Population<T>::operator[](int), exceeding mating population memory.");
   }
   #endif

   return matpop[pos];
}

/*-------------------------------------------------------------------------------------------------*/

// return iterator to current population beginning 
template <typename T, int N>
inline typename std::vector<CHR<T,N>>::iterator Population<T,N>::begin()
{
   return curpop.begin(); 
}

/*-------------------------------------------------------------------------------------------------*/

// return const iterator to current population beginning 
template <typename T, int N>
inline typename std::vector<CHR<T,N>>::const_iterator Population<T,N>::cbegin() const
{
   return curpop.cbegin(); 
}

/*-------------------------------------------------------------------------------------------------*/

// return iterator to current population ending
template <typename T, int N>
inline typename std::vector<CHR<T,N>>::iterator Population<T,N>::end()
{ 
   return curpop.end();
}

/*-------------------------------------------------------------------------------------------------*/

// return const iterator to current population ending
template <typename T, int N>
inline typename std::vector<CHR<T,N>>::const_iterator Population<T,N>::cend() const
{ 
   return curpop.cend();
}

/*-------------------------------------------------------------------------------------------------*/

// select element at position pos in current population and copy it into mating population
template <typename T, int N>
inline void Population<T,N>::select(int pos)
{
   #ifndef NDEBUG
   if (pos > ptr->popsize - 1) {
      throw std::invalid_argument("Error: in galgo::Population<T>::select(int), exceeding current population memory.");
   }
   if (matidx == ptr->matsize) {
      throw std::invalid_argument("Error: in galgo::Population<T>::select(int), exceeding mating population memory.");
   }
   #endif

   matpop[matidx] = curpop[pos];
   matidx++;
}

/*-------------------------------------------------------------------------------------------------*/
    
// set all fitness to positive values (used in RWS and SUS selection methods)
template <typename T, int N>
void Population<T,N>::adjustFitness()
{
   // getting worst population fitness
   T worstFitness = curpop.back()->fitness;

   if (worstFitness < 0) {
      // getting best fitness
      T bestFitness = curpop.front()->fitness;
      // case where all fitness are equal and negative
      if (worstFitness == bestFitness) {
         std::for_each(curpop.begin(), curpop.end(), [](CHR<T,N>& chr)->void{chr->fitness *= -1;});
       } else {
         std::for_each(curpop.begin(), curpop.end(), [worstFitness](CHR<T,N>& chr)->void{chr->fitness -= worstFitness;});
      }
   }
}

/*-------------------------------------------------------------------------------------------------*/

// compute population fitness sum (used in TRS, RWS and SUS selection methods)
template <typename T, int N>
inline T Population<T,N>::getSumFitness() const
{
   return std::accumulate(curpop.cbegin(), curpop.cend(), 0.0, [](T sum, const CHR<T,N>& chr)->T{return sum + T(chr->fitness);});
}

/*-------------------------------------------------------------------------------------------------*/

// get worst objective function total result from current population (used in constraint(s) adaptation)
template <typename T, int N>
inline T Population<T,N>::getWorstTotal() const
{
   auto it = std::min_element(curpop.begin(), curpop.end(), [](const CHR<T,N>& chr1, const CHR<T,N>& chr2)->bool{return chr1->getTotal() < chr2->getTotal();});

   return (*it)->getTotal();
}

/*-------------------------------------------------------------------------------------------------*/

// return population size
template <typename T, int N>
inline int Population<T,N>::popsize() const
{
   return ptr->popsize;
}

/*-------------------------------------------------------------------------------------------------*/

// return mating population size
template <typename T, int N>
inline int Population<T,N>::matsize() const
{
   return ptr->matsize;
}

/*-------------------------------------------------------------------------------------------------*/

// return tournament size
template <typename T, int N>
inline int Population<T,N>::tntsize() const
{
   return ptr->tntsize;
}

/*-------------------------------------------------------------------------------------------------*/

// return numero of generation
template <typename T, int N>
inline int Population<T,N>::nogen() const
{
   return ptr->nogen;
}

/*-------------------------------------------------------------------------------------------------*/

// return number of generations
template <typename T, int N>
inline int Population<T,N>::nbgen() const
{
   return ptr->nbgen;
}

/*-------------------------------------------------------------------------------------------------*/

// return selection pressure
template <typename T, int N>
inline T Population<T,N>::SP() const
{
   return ptr->SP;
}

//=================================================================================================

}

#endif


