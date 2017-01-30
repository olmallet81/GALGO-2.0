#include "Galgo.hpp"

// objective class example
template <typename T>
class MyObjective
{
public:
   // objective function example : Rosenbrock function
   // minimizing f(x,y) = (1 - x)^2 + 100 * (y - x^2)^2
   static std::vector<T> Objective(const std::vector<T>& x)
   {
      T obj = -(pow(1-x[0],2)+100*pow(x[1]-x[0]*x[0],2));
      return std::vector<T>({obj});
   }
   // NB: GALGO maximize by default so we will maximize -f(x,y)
};

// constraint example:
// 1) x * y + x - y + 1.5 <= 0
// 2) 10 - x * y <= 0
template <typename T>
std::vector<T> MyConstraint(const std::vector<T>& x)
{
   return std::vector<T>({x[0]*x[1]+x[0]-x[1]+1.5,10-x[0]*x[1]});
}

int main()
{
   // lower bounds LB and upper bounds UB
   std::vector<double> LB({0.0,0.0});
   std::vector<double> UB({1.0,13.0});

   // initiliazing genetic algorithm
   galgo::GeneticAlgorithm<double> ga(MyObjective<double>::Objective,200,LB,UB,50,true);
 
   // setting constraints
   ga.Constraint = MyConstraint;

   // running genetic algorithm
   ga.run();
}
