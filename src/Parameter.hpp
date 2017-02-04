//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef PARAMETER_H
#define PARAMETER_H

namespace galgo {

//=================================================================================================

// end of recursion for computing the sum of a parameter pack of integral numbers
int sum(int first) 
{
   return first;
}

// recursion for computing the sum of a parameter pack of integral numbers
template <typename...Args>
int sum(int first, Args...args) 
{
   return first + sum(args...);
}

/*-------------------------------------------------------------------------------------------------*/

template <typename T, int N = 16>
class Parameter
{
   template <typename K, int...S>
   friend class Chromosome;

public:
   std::vector<T> data; // contains lower bound, upper bound and initial value (optional)
   // nullary constructor
   Parameter() {}
   // constructor
   Parameter(const std::vector<T>& data) {
      this->data = data;
   }
private:
   // encoding random unsigned integer
   std::string encode() const {
      // generating and encoding a random unsigned integer
      std::string str = galgo::GetBinary(Randomize<N>::generate());
      // truncating string to desired number of bits N
      return str.substr(str.size() - N, N);;
   }
   // encoding known unsigned integer
   std::string encode(T z) const {
      // converting known value to unsigned integer
      uint64_t value = Randomize<N>::MAXVAL * (z - data[0]) / (data[1] - data[0]);
      // encoding it into a string
      std::string str = GetBinary(value);
      // truncating string to desired number of bits N
      return str.substr(str.size() - N, N);;
   }
   // decoding string to real value
   T decode(const std::string& str) const {
      return data[0] + (galgo::GetValue(str) / (double)Randomize<N>::MAXVAL) * (data[1] - data[0]);
   }
   // return encoded parameter size in number of bits 
   T size() const {
      return N;
   }
};

//=================================================================================================

}

#endif
