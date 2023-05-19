#include "Multigrid.h"

Real xsquared(Real x)
{
  return x*x;
}

Real one(Real x)
{
  return 1;
}

Real u(Real x)
{
  return exp(sin(x));
}

Real u_x(Real x)
{
  return cos(x)*exp(sin(x));
}

Real f(Real x)
{
  return (-1)*(-sin(x) + pow(cos(x),2))*exp(sin(x));
}

double u2(double x, double y)
{
  return exp(y + sin(x));
}

double u_x2(double x, double y)
{
  return cos(x)*exp(y + sin(x));
}

double u_y2(double x, double y)
{
  return exp(y + sin(x));
}

double f2(double x, double y)
{
  return exp(y + sin(x))*(sin(x) - 1)*(sin(x) + 2);
}

int main(int argc, char* argv[])
{
  Multigrid<1> dir32(u, f, 32);
  Multigrid<1> dir64(u, f, 64);
  Multigrid<1> dir128(u, f, 128);
  Multigrid<1> dir256(u, f, 256);

  Multigrid<1> neu32(u, u_x, f, 32, "PureNeumann");
  Multigrid<1> neu64(u, u_x, f, 64, "PureNeumann");
  Multigrid<1> neu128(u, u_x, f, 128, "PureNeumann");
  Multigrid<1> neu256(u, u_x, f, 256, "PureNeumann");

  Multigrid<1> left32(u, u_x, f, 32, "LeftNeumann");
  Multigrid<1> left64(u, u_x, f, 64, "LeftNeumann");
  Multigrid<1> left128(u, u_x, f, 128, "LeftNeumann");
  Multigrid<1> left256(u, u_x, f, 256, "LeftNeumann");

  Multigrid<1> right32(u, u_x, f, 32, "RightNeumann");
  Multigrid<1> right64(u, u_x, f, 64, "RightNeumann");
  Multigrid<1> right128(u, u_x, f, 128, "RightNeumann");
  Multigrid<1> right256(u, u_x, f, 256, "RightNeumann");
  
  Real rtol = pow(10,-16);
  int it = 30;

  
  std::cout << "Pure Dirichlet Boundary Condition" << std::endl;
  std::cout << "n = 32 : " << std::endl;
  dir32.solve(rtol, it);
  std::cout << std::endl;
  //dir32.print_res();
  //std::cout << std::endl;
  //dir32.print_resrate();
  //std::cout << std::endl;
  
  std::cout << "n = 64 : " << std::endl;
  dir64.solve(rtol, it);
  std::cout << std::endl;
  //dir64.print_res();
  //std::cout << std::endl;
  //dir64.print_resrate();
  //std::cout << std::endl;
  
  std::cout << "n = 128 : " << std::endl;
  dir128.solve(rtol, it);
  std::cout << std::endl;
  //dir128.print_res();
  //std::cout << std::endl;
  //dir128.print_resrate();
  //std::cout << std::endl;
  
  std::cout << "n = 256 : " << std::endl;
  dir256.solve(rtol, it);
  std::cout << std::endl;
  //dir256.print_res();
  //std::cout << std::endl;
  //dir256.print_resrate();
  //std::cout << std::endl;
  

  /*
  std::cout << "Pure Neumann Boundary Condition" << std::endl;
  std::cout << "n = 32 : " << std::endl;
  neu32.solve(rtol, it);
  std::cout << std::endl;
  //neu32.print_res();
  //std::cout << std::endl;
  //neu32.print_resrate();
  //std::cout << std::endl;
  

  std::cout << "n = 64 : " << std::endl;
  neu64.solve(rtol, it);
  std::cout << std::endl;
  //neu64.print_res();
  //std::cout << std::endl;
  //neu64.print_resrate();
  //std::cout << std::endl;

  std::cout << "n = 128 : " << std::endl;
  neu128.solve(rtol, it);
  std::cout << std::endl;
  //neu128.print_res();
  //std::cout << std::endl;
  //neu128.print_resrate();
  //std::cout << std::endl;

  std::cout << "n = 256 : " << std::endl;
  neu256.solve(rtol, it);
  std::cout << std::endl;
  //neu256.print_res();
  //std::cout << std::endl;
  //neu256.print_resrate();
  //std::cout << std::endl;
  */

  /*
  std::cout << "Left Neumann Boundary Condition" << std::endl;
  std::cout << "n = 32 : " << std::endl;
  left32.solve(rtol, it);
  std::cout << std::endl;
  left32.print_res();
  std::cout << std::endl;
  left32.print_resrate();
  std::cout << std::endl;

  std::cout << "n = 64 : " << std::endl;
  left64.solve(rtol, it);
  std::cout << std::endl;
  left64.print_res();
  std::cout << std::endl;
  left64.print_resrate();
  std::cout << std::endl;

  std::cout << "n = 128 : " << std::endl;
  left128.solve(rtol, it);
  std::cout << std::endl;
  left128.print_res();
  std::cout << std::endl;
  left128.print_resrate();
  std::cout << std::endl;

  std::cout << "n = 256 : " << std::endl;
  left256.solve(rtol, it);
  std::cout << std::endl;
  left256.print_res();
  std::cout << std::endl;
  left256.print_resrate();
  std::cout << std::endl;
  */

  /*
  std::cout << "Right Neumann Boundary Condition" << std::endl;
  std::cout << "n = 32 : " << std::endl;
  right32.solve(rtol, it);
  std::cout << std::endl;
  right32.print_res();
  std::cout << std::endl;
  right32.print_resrate();
  std::cout << std::endl;

  std::cout << "n = 64 : " << std::endl;
  right64.solve(rtol, it);
  std::cout << std::endl;
  right64.print_res();
  std::cout << std::endl;
  right64.print_resrate();
  std::cout << std::endl;

  std::cout << "n = 128 : " << std::endl;
  right128.solve(rtol, it);
  std::cout << std::endl;
  right128.print_res();
  std::cout << std::endl;
  right128.print_resrate();
  std::cout << std::endl;

  std::cout << "n = 256 : " << std::endl;
  right256.solve(rtol, it);
  std::cout << std::endl;
  right256.print_res();
  std::cout << std::endl;
  right256.print_resrate();
  std::cout << std::endl;
  */
  
  //Multigrid<2> test2(u2, f2, 32);
  return 0;
}
