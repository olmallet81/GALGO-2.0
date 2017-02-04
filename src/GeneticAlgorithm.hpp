//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef GENETICALGORITHM_HPP
#define GENETICALGORITHM_HPP

namespace galgo {

//=================================================================================================

template <typename T, int...N>
class GeneticAlgorithm
{
   static_assert(std::is_same<float,T>::value || std::is_same<double,T>::value, "variable type can only be float or double, please amend.");

   template <typename K, int...S>
   friend class Population;
   template <typename K, int...S>
   friend class Chromosome;

   template <typename K>
   using Functor = std::vector<K> (*)(const std::vector<K>&);

private:
   Population<T,N...> pop;                                                 // population of chromosomes
   std::tuple<typename std::remove_reference<Parameter<T,N>>::type...> tp; // tuple of parameters
   std::vector<T> lowerBound;                                              // parameter(s) lower bound
   std::vector<T> upperBound;                                              // parameter(s) upper bound
   std::vector<T> initialSet;                                              // initial set of parameter(s)

public: 
   // objective function functor
   Functor<T> Objective; 
   // selection method functor initialized to roulette wheel selection                                   
   void (*Selection)(Population<T,N...>&) = RWS;  
   // cross-over method functor initialized to 1-point cross-over                                
   void (*CrossOver)(const Population<T,N...>&, CHR<T,N...>&, CHR<T,N...>&) = P1XO;
   // mutation method functor initialized to single-point mutation 
   void (*Mutation)(CHR<T,N...>&) = SPM;  
   // adaptation to constraint(s) method functor                                        
   void (*Adaptation)(Population<T,N...>&) = nullptr; 
   // constraint(s) functor                               
   std::vector<T> (*Constraint)(const std::vector<T>&) = nullptr; 

   T covrate = .50;   // cross-over rate
   T mutrate = .05;   // mutation rate   
   T SP = 1.5;        // selective pressure for RSP selection method 
   T tolerance = 0.0; // terminal condition (inactive if equal to zero)
                 
   int elitpop = 1;   // elit population size
   int matsize;       // mating pool size, set to popsize by default
   int tntsize = 10;  // tournament size
   int genstep = 10;  // generation step for outputting results
   int precision = 5; // precision for outputting results

   // constructor
   GeneticAlgorithm(Functor<T> objective, int popsize, const Parameter<T,N>&...args, int nbgen, bool output = false);
   // run genetic algorithm
   void run();
   // return best chromosome 
   const CHR<T,N...>& result() const;

private:
   int nbbit;     // total number of bits per chromosome
   int nbgen;     // number of generations
   int nogen = 0; // numero of generation
   int nbparam;   // number of parameters to be estimated
   int popsize;   // population size
   bool output;   // control if results must be outputted

   // end of recursion for initializing parameter(s) data
   template <int I = 0>
   typename std::enable_if<I == sizeof...(N), void>::type init(); 
   // recursion for initializing parameter(s) data
   template <int I = 0>
   typename std::enable_if<I < sizeof...(N), void>::type init();

   // check inputs validity
   void check() const ;
   // print results for each new generation
   void print() const;
};

/*-------------------------------------------------------------------------------------------------*/
   
// constructor
template <typename T, int...N>
GeneticAlgorithm<T,N...>::GeneticAlgorithm(Functor<T> objective, int popsize, const Parameter<T,N>&...args, int nbgen, bool output)
{
   this->Objective = objective;
   this->popsize = popsize;
   this->matsize = popsize;
   this->tp = std::make_tuple(args...);
   this->nbbit = sum(N...);
   this->nbgen = nbgen;
   this->nbparam = sizeof...(N);
   this->output = output;
}

/*-------------------------------------------------------------------------------------------------*/

// end of recursion for initializing parameter(s) data
template <typename T, int...N> template <int I>
inline typename std::enable_if<I == sizeof...(N), void>::type 
GeneticAlgorithm<T,N...>::init() {}

// recursion for initializing parameter(s) data
template <typename T, int...N> template <int I>
inline typename std::enable_if<I < sizeof...(N), void>::type 
GeneticAlgorithm<T,N...>::init() 
{
   const auto& par = std::get<I>(tp);
   lowerBound.push_back(par.data[0]);
   upperBound.push_back(par.data[1]);
   if (par.data.size() > 2) {
      initialSet.push_back(par.data[2]);
   }
   init<I + 1>();
}

/*-------------------------------------------------------------------------------------------------*/

// check inputs validity
template <typename T, int...N>
void GeneticAlgorithm<T,N...>::check() const
{
   if (!initialSet.empty()) {
      for (int i = 0; i < nbparam; ++i) {
         if (initialSet[i] < lowerBound[i] || initialSet[i] > upperBound[i]) {
            throw std::invalid_argument(" Error: in class galgo::Parameter<T,N>, initial parameter value cannot be outside the parameter boundaries, please choose a value between its lower and upper bounds.");
            return false;
         }
      }
      if (initialSet.size() != (unsigned)nbparam) {
         throw std::invalid_argument(" Error: in class galgo::GeneticAlgorithm<T>, initial set of parameters does not have the same dimension than the number of parameters, please adjust.");
         return false;
      }
   }
   if (SP < 1.0 || SP > 2.0) {
      throw std::invalid_argument(" Error: in class galgo::GeneticAlgorithm<T>, selective pressure (SP) cannot be outside [1.0,2.0], please choose a real value within this interval.");
   }
   if (elitpop > popsize || elitpop < 0) {
      throw std::invalid_argument(" Error: in class galgo::GeneticAlgorithm<T>, elit population (elitpop) cannot outside [0,popsize], please choose an integral value within this interval.");
   }
   if (covrate < 0.0 || covrate > 1.0) {
      throw std::invalid_argument(" Error: in class galgo::GeneticAlgorithm<T>, cross-over rate (covrate) cannot outside [0.0,1.0], please choose a real value within this interval.");
      return false;
   }
   if (genstep <= 0) {
      throw std::invalid_argument(" Error: in class galgo::GeneticAlgorithm<T>, generation step (genstep) cannot be <= 0, please choose an integral value > 0.");
   }
}

/*-------------------------------------------------------------------------------------------------*/
   
// run genetic algorithm
template <typename T, int...N>
void GeneticAlgorithm<T,N...>::run()
{
   lowerBound.clear();
   upperBound.clear();
   initialSet.clear();

   // initializing parameter(s) data
   init();
   // checking inputs validity
   check();

   // setting adaptation method to default if needed
   if (Constraint != nullptr && Adaptation == nullptr) {
      Adaptation = DAC;
   }

   // initializing population
   pop = Population<T,N...>(*this);

   if (output) {
      std::cout << "\n Running Genetic Algorithm...\n";
      std::cout << " ----------------------------\n";
   }

   // creating population
   pop.creation();
   // initializing best result and previous best result
   T bestResult = pop(0)->getTotal();
   T prevBestResult = bestResult;
   // outputting results 
   if (output) print();
    
   // starting population evolution
   for (nogen = 1; nogen <= nbgen; ++nogen) {
      // evolving population
      pop.evolution();
      // getting best current result
      bestResult = pop(0)->getTotal();
      // outputting results
      if (output) print();
      // checking convergence
      if (tolerance != 0.0) {
         if (fabs(bestResult - prevBestResult) < fabs(tolerance)) {
            break;
         }
         prevBestResult = bestResult;
      }
   } 

   // outputting contraint value
   if (Constraint != nullptr) {
      // getting best parameter(s) constraint value(s)
      std::vector<T> cst = pop(0)->getConstraint(); 
      if (output) {
         std::cout << "\n Constraint(s)\n";
         std::cout << " -------------\n";
         for (unsigned i = 0; i < cst.size(); ++i) {
            std::cout << " C"; 
            if (nbparam > 1) {
               std::cout << std::to_string(i + 1);
            }
            std::cout << "(x) = " << std::setw(6) << std::fixed << std::setprecision(precision) << cst[i] << "\n"; 
         }
         std::cout << "\n"; 
      }
   }   
}

/*-------------------------------------------------------------------------------------------------*/

// return best chromosome
template <typename T, int...N>
inline const CHR<T,N...>& GeneticAlgorithm<T,N...>::result() const
{
   return pop(0);
}

/*-------------------------------------------------------------------------------------------------*/
   
// print results for each new generation
template <typename T, int...N>
void GeneticAlgorithm<T,N...>::print() const
{
   // getting best parameter(s) from best chromosome
   std::vector<T> bestParam = pop(0)->getParam();
   std::vector<T> bestResult = pop(0)->getResult();

   if (nogen % genstep == 0) {
      std::cout << " Generation = " << std::setw(std::to_string(nbgen).size()) << nogen << " |";
      for (int i = 0; i < nbparam; ++i) {
	      std::cout << " X";
         if (nbparam > 1) {
            std::cout << std::to_string(i + 1);
         }
         std::cout << " = "  << std::setw(9) << std::fixed << std::setprecision(precision) << bestParam[i] << " |";
	   }
      for (unsigned i = 0; i < bestResult.size(); ++i) {
	      std::cout << " F";
         if (bestResult.size() > 1) {
            std::cout << std::to_string(i + 1);
         }
         std::cout << "(x) = " << std::setw(12) << std::fixed << std::setprecision(precision) << bestResult[i];
         if (i < bestResult.size() - 1) {
            std::cout << " |";
         } else {
            std::cout << "\n";
         }
	   }
 
   }
}
   
//=================================================================================================

}

#endif
