//=================================================================================================
//                    Copyright (C) 2017 Olivier Mallet - All Rights Reserved                      
//=================================================================================================

#ifndef CONVERTER_HPP
#define CONVERTER_HPP

namespace galgo {

//=================================================================================================

// convert unsigned long long integer to binary string
std::string GetBinary(uint64_t value)
{
   std::bitset<sizeof(uint64_t)*CHAR_BIT> bits(value);
   // NB: CHAR_BIT = number of bits in char usually 8 but not always on older machines
   return bits.to_string();
}

/*-------------------------------------------------------------------------------------------------*/

// convert binary string to unsigned long long integer
uint64_t GetValue(const std::string& s)
{
   uint64_t value, x = 0;
   for (std::string::const_iterator it = s.begin(), end = s.end(); it != end; ++it) {
      x = (x << 1) + (*it - '0');
   }
   memcpy(&value, &x, sizeof(uint64_t));

   return value;
}

//=================================================================================================

}

#endif
