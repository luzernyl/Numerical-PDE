#include "Poisson.h"
#include <cmath>

double zero(double x, double y)
{
  return 0;
}

double u(double x, double y)
{
  return exp(y + sin(x));
}

double u_x(double x, double y)
{
  return cos(x)*exp(y + sin(x));
}

double u_y(double x, double y)
{
  return exp(y + sin(x));
}

double f(double x, double y)
{
  return exp(y + sin(x))*(sin(x) - 1)*(sin(x) + 2);
}

int main(int argc, char* argv[])
{
  Poisson prob(u, f, 8);
  prob.print_rhs();
  /*
  Poisson prob1(u, f, 8, u_x, u_y);
  Poisson prob2(u, f, 16, u_x, u_y);
  Poisson prob3(u, f, 32, u_x, u_y);
  Poisson prob4(u, f, 64, u_x, u_y);

  VectorXd sol1, sol2, sol3, sol4;
  std::cout << "n = 8 : " << std::endl;
  //std::cout << "u8 = ";
  sol1 = prob1.solve();
  //prob1.print_sol4matlab();
  std::cout << std::endl;
  
  std::cout << "n = 16 : " << std::endl;
  //std::cout << "u16 = ";
  sol2 = prob2.solve();
  //prob2.print_sol4matlab();
  std::cout << std::endl;
  
  std::cout << "n = 32 : " << std::endl;
  //std::cout << "u32 = ";
  sol3 = prob3.solve();
  //prob3.print_sol4matlab();
  std::cout << std::endl;
  
  std::cout << "n = 64 : " << std::endl;
  //std::cout << "u64 = ";
  sol4 = prob4.solve();
  //prob4.print_sol4matlab();
  std::cout << std::endl;
  */

  std::cout << "Program Finished!" << std::endl;
  
  return 0;
}
