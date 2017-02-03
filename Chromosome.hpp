//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef CHROMOSOME_HPP
#define CHROMOSOME_HPP

namespace galgo {

//=================================================================================================

template <typename T, int...N>
class Chromosome
{
   static_assert(std::is_same<float,T>::value || std::is_same<double,T>::value, "variable type can only be float or double, please amend.");

public:
   // constructor
   Chromosome(const GeneticAlgorithm<T,N...>& ga);
   // copy constructor
   Chromosome(const Chromosome<T,N...>& rhs);
   // end of recursion for creating new chromosome
   template <int I = 0>
   typename std::enable_if<I == sizeof...(N), void>::type create();
   // recursion for creating new chromosome
   template <int I = 0>
   typename std::enable_if<I < sizeof...(N), void>::type create();
   // end of recursion for initializing new chromosome
   template <int I = 0>
   typename std::enable_if<I == sizeof...(N), void>::type initialize();
   // recursion for initializing new chromosome
   template <int I = 0>
   typename std::enable_if<I < sizeof...(N), void>::type initialize();
   // end of recursion for converting gene to real value
   template <int I = 0>
   typename std::enable_if<I == sizeof...(N), void>::type convert();
   // recursion for converting gene to real value
   template <int I = 0>
   typename std::enable_if<I < sizeof...(N), void>::type convert();
   // evaluate chromosome 
   void evaluate(); 
   // reset chromosome
   void reset();
   // end of recursion for setting or replace kth gene by a new one 
   template <int I = 0>
   inline typename std::enable_if<I == sizeof...(N), void>::type setGene(int k);
   // recursion for setting or replace kth gene by a new one 
   template <int I = 0>
   inline typename std::enable_if<I < sizeof...(N), void>::type setGene(int k);
   // end of recursion for initializing or replacing kth gene by a know value
   template <int I = 0>
   typename std::enable_if<I == sizeof...(N), void>::type initGene(int k, T x);
   // recursion for initializing or replacing kth gene by a know value
   template <int I = 0>
   typename std::enable_if<I < sizeof...(N), void>::type initGene(int k, T x);
   // add bit to chromosome
   void addBit(char bit);
   // initialize or replace an existing chromosome bit   
   void setBit(char bit, int pos);
   // flip an existing chromosome bit
   void flipBit(int pos);
   // get chromosome bit
   char getBit(int pos) const;
   // initialize or replace a portion of bits with a portion of another chromosome
   void setPortion(const Chromosome<T,N...>& x, int start, int end);
   // initialize or replace a portion of bits with a portion of another chromosome
   void setPortion(const Chromosome<T,N...>& x, int start);
   // get parameter value(s) from chromosome
   const std::vector<T>& getParam() const;
   // get objective function result
   const std::vector<T>& getResult() const;
   // get the total sum of all objective function(s) result
   T getTotal() const;
   // get constraint value(s)
   const std::vector<T> getConstraint() const;
   // return chromosome size in number of bits
   int size() const;
   // return number of chromosome bits to mutate
   T mutrate() const;
   // return number of genes in chromosome
   int nbgene() const;
   // return numero of generation this chromosome belongs to
   int nogen() const;
   // return lower bound(s)
   const std::vector<T>& lowerBound() const;
   // return upper bound(s)
   const std::vector<T>& upperBound() const;

private:
   std::vector<T> param;                          // estimated parameters parameters
   std::vector<T> result;                         // chromosome objective function(s) result
   std::string chr;                               // string of bits representing chromosome
   const GeneticAlgorithm<T,N...>* ptr = nullptr; // pointer to genetic algorithm
public:
   T fitness;                                     // chromosome fitness, objective function(s) result that can be modified (adapted to constraint(s), set to positive values, etc...)
private:
   T total;                                       // total sum of objective function(s) result
   int chrsize;                                   // chromosome size (in number of bits)
   int numgen;                                    // numero of generation
   int idx;                                       // index for chromosome breakdown
};

/*-------------------------------------------------------------------------------------------------*/

// constructor
template <typename T, int...N>
Chromosome<T,N...>::Chromosome(const GeneticAlgorithm<T,N...>& ga) 
{
   param.resize(sizeof...(N));
   ptr = &ga;
   chrsize = ga.nbbit;
   numgen = ga.nogen;
}

/*-------------------------------------------------------------------------------------------------*/

// copy constructor
template <typename T, int...N>
Chromosome<T,N...>::Chromosome(const Chromosome<T,N...>& rhs) 
{
   param = rhs.param;
   result = rhs.result;
   chr = rhs.chr;
   ptr = rhs.ptr;
   // re-initializing fitness to its original value
   fitness = rhs.total;
   total = rhs.total;
   chrsize = rhs.chrsize;
   numgen = rhs.numgen;
}

/*-------------------------------------------------------------------------------------------------*/

// end of recursion for creating new chromosome
template <typename T, int...N> template <int I> 
inline typename std::enable_if<I == sizeof...(N), void>::type 
Chromosome<T,N...>::create() {}

// recursion for creating new chromosome
template <typename T, int...N> template <int I> 
inline typename std::enable_if<I < sizeof...(N), void>::type 
Chromosome<T,N...>::create() 
{
   std::string str = std::get<I>(ptr->tp).encode();
   chr.append(str);  
   create<I + 1>();
}

/*-------------------------------------------------------------------------------------------------*/

// intialize chromosome from a known value
template <typename T, int...N> template <int I> 
inline typename std::enable_if<I == sizeof...(N), void>::type 
Chromosome<T,N...>::initialize() {}

// recursion for initializing chromosome from a known value
template <typename T, int...N> template <int I> 
inline typename std::enable_if<I < sizeof...(N), void>::type 
Chromosome<T,N...>::initialize() 
{
   std::string str = std::get<I>(ptr->tp).encode(ptr->initialSet[I]);
   chr.append(str);  
   initialize<I + 1>();
}

/*-------------------------------------------------------------------------------------------------*/

// end of recursion for converting gene to real value
template <typename T, int...N> template <int I>
inline typename std::enable_if<I == sizeof...(N), void>::type 
Chromosome<T,N...>::convert() {}

// recursion for converting gene to real value
template <typename T, int...N> template <int I>
inline typename std::enable_if<I < sizeof...(N), void>::type 
Chromosome<T,N...>::convert() 
{
   const auto& par = std::get<I>(ptr->tp);
   param[I] = par.decode(chr.substr(idx, par.size()));
   idx += par.size();  
   convert<I + 1>();
}

/*-------------------------------------------------------------------------------------------------*/

// evaluate chromosome fitness
template <typename T, int...N>
inline void Chromosome<T,N...>::evaluate() {
   // setting index for chromosome breakdown to 0
   idx = 0;
   // converting chromosome to value(s)
   convert();
   // computing objective result(s) 
   result = ptr->Objective(param);
   // computing sum of all results (in case there is not only one objective functions)
   total = std::accumulate(result.begin(), result.end(), 0.0);
   // initializing fitness to this total
   fitness = total;
}

/*-------------------------------------------------------------------------------------------------*/

// reset chromosome
template <typename T, int...N>
inline void Chromosome<T,N...>::reset()
{
   chr.clear();
   result = 0.0;
   total = 0.0;
   fitness = 0.0;
}

/*-------------------------------------------------------------------------------------------------*/

// end of recursion for setting or replace kth gene by a new one 
template <typename T, int...N> template <int I>
inline typename std::enable_if<I == sizeof...(N), void>::type 
Chromosome<T,N...>::setGene(int k) {}

// recursion for setting or replace kth gene by a new one 
template <typename T, int...N> template <int I>
inline typename std::enable_if<I < sizeof...(N), void>::type 
Chromosome<T,N...>::setGene(int k) 
{
   #ifndef NDEBUG
   if (k < 0 || k >= ptr->nbparam) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::setGene(int), argument cannot be outside interval [0,nbparam-1], please amend.");
   }
   #endif

   // setting index for chromosome breakdown to 0
   if (I == 0) idx = 0;

   // getting Ith parameter in tuple
   const auto& par = std::get<I>(ptr->tp);
 
   // if gene found
   if (I == k) {
      // generating a new gene
      std::string s = par.encode();
      // adding or replacing gene in chromosome
      chr.replace(idx, par.size(), s, 0, par.size());
      // ending recursion
      setGene<sizeof...(N)>(k); 
   } else {
      idx += par.size();
      setGene<I + 1>(k);
   }
}

/*-------------------------------------------------------------------------------------------------*/

// end of recursion for initializing or replacing kth gene by a know value
template <typename T, int...N> template <int I>
inline typename std::enable_if<I == sizeof...(N), void>::type 
Chromosome<T,N...>::initGene(int k, T x) {}

// recursion for initializing or replacing kth gene by a know value
template <typename T, int...N> template <int I>
inline typename std::enable_if<I < sizeof...(N), void>::type 
Chromosome<T,N...>::initGene(int k, T x)
{
   #ifndef NDEBUG
   if (k < 0 || k >= ptr->nbparam) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::initGene(int), first argument cannot be outside interval [0,nbparam-1], please amend.");
   }
   #endif

   // setting index for chromosome breakdown to 0
   if (I == 0) idx = 0;

   // getting Ith parameter in tuple
   const auto& par = std::get<I>(ptr->tp);
 
   // if gene found
   if (I == k) {
      // generating a new gene
      std::string s = par.encode(x);
      // adding or replacing gene in chromosome
      chr.replace(idx, par.size(), s, 0, par.size());
      // ending recursion
      initGene<sizeof...(N)>(k, x); 
   } else {
      idx += par.size();
      initGene<I + 1>(k, x);
   }
}

/*-------------------------------------------------------------------------------------------------*/

// add chromosome bit to chromosome (when constructing a new one)
template <typename T, int...N>
inline void Chromosome<T,N...>::addBit(char bit)
{
   chr.push_back(bit);

   #ifndef NDEBUG
   if (chr.size() > chrsize) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::setBit(char), exceeding chromosome size.");
   }
   #endif
}

/*-------------------------------------------------------------------------------------------------*/

// initialize or replace an existing chromosome bit   
template <typename T, int...N>
inline void Chromosome<T,N...>::setBit(char bit, int pos)
{  
   #ifndef NDEBUG
   if (pos >= chrsize) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::replaceBit(char, int), second argument cannot be equal or greater than chromosome size.");
   }
   #endif

   std::stringstream ss;
   std::string str;
   ss << bit;
   ss >> str;
   chr.replace(pos, 1, str);
   std::cout << chr << "\n";
}

/*-------------------------------------------------------------------------------------------------*/
      
// flip an existing chromosome bit
template <typename T, int...N>
inline void Chromosome<T,N...>::flipBit(int pos)
{
   #ifndef NDEBUG
   if (pos >= chrsize) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::flipBit(int), argument cannot be equal or greater than chromosome size.");
   }
   #endif

   if (chr[pos] == '0') {
      chr.replace(pos, 1, "1");
   } else {
      chr.replace(pos, 1, "0");
   }
}

/*-------------------------------------------------------------------------------------------------*/

// get a chromosome bit
template <typename T, int...N>
inline char Chromosome<T,N...>::getBit(int pos) const
{
   #ifndef NDEBUG
   if (pos >= chrsize) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::getBit(int), argument cannot be equal or greater than chromosome size.");
   }
   #endif

   return chr[pos];
}

/*-------------------------------------------------------------------------------------------------*/

// initialize or replace a portion of bits with a portion of another chromosome (from position start to position end included)
template <typename T, int...N>
inline void Chromosome<T,N...>::setPortion(const Chromosome<T,N...>& x, int start, int end)
{
   #ifndef NDEBUG
   if (start > chrsize) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::setPortion(const Chromosome<T>&, int, int), second argument cannot be greater than chromosome size.");
   }
   #endif

   chr.replace(start, end - start + 1, x.chr, start, end - start + 1);
}

/*-------------------------------------------------------------------------------------------------*/

// initialize or replace a portion of bits with a portion of another chromosome (from position start to the end of he chromosome)
template <typename T, int...N>
inline void Chromosome<T,N...>::setPortion(const Chromosome<T,N...>& x, int start)
{
   #ifndef NDEBUG
   if (start > chrsize) {
      throw std::invalid_argument("Error: in galgo::Chromosome<T>::setPortion(const Chromosome<T>&, int), second argument cannot be greater than chromosome size.");
   }
   #endif

   chr.replace(start, chrsize, x.chr, start, x.chrsize);
}

/*-------------------------------------------------------------------------------------------------*/

// get parameter value(s) from chromosome
template <typename T, int...N>
inline const std::vector<T>& Chromosome<T,N...>::getParam() const
{
   return param;
}

/*-------------------------------------------------------------------------------------------------*/

// get objective function result
template <typename T, int...N>
inline const std::vector<T>& Chromosome<T,N...>::getResult() const
{
   return result;
}

/*-------------------------------------------------------------------------------------------------*/

// get the total sum of all objective function(s) result
template <typename T, int...N>
inline T Chromosome<T,N...>::getTotal() const
{
   return total;
}

/*-------------------------------------------------------------------------------------------------*/

// get constraint value(s) for this chromosome
template <typename T, int...N>
inline const std::vector<T> Chromosome<T,N...>::getConstraint() const
{
   return ptr->Constraint(param);
}

/*-------------------------------------------------------------------------------------------------*/

// return chromosome size in number of bits
template <typename T, int...N>
inline int Chromosome<T,N...>::size() const
{
   return chrsize;
}

/*-------------------------------------------------------------------------------------------------*/

// return mutation rate 
template <typename T, int...N>
inline T Chromosome<T,N...>::mutrate() const
{
   return ptr->mutrate;
}

/*-------------------------------------------------------------------------------------------------*/

// return number of genes in chromosome
template <typename T, int...N>
inline int Chromosome<T,N...>::nbgene() const
{
   return ptr->nbparam;
}

/*-------------------------------------------------------------------------------------------------*/

// return numero of generation this chromosome belongs to
template <typename T, int...N>
inline int Chromosome<T,N...>::nogen() const
{
   return numgen;
}

/*-------------------------------------------------------------------------------------------------*/

// return lower bound(s)
template <typename T, int...N>
inline const std::vector<T>& Chromosome<T,N...>::lowerBound() const
{
   return ptr->lowerBound;
}

/*-------------------------------------------------------------------------------------------------*/

// return upper bound(s)
template <typename T, int...N>
inline const std::vector<T>& Chromosome<T,N...>::upperBound() const
{
   return ptr->upperBound;
}

//=================================================================================================

}

#endif
