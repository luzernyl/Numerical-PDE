/**
 * @file   Poisson.h
 * @author Luzern Yuven Luis <student@student>
 * @date   Thu Mar 31 01:15:16 2022
 * 
 * @brief  A C++ package to solve two-dimensional Poisson equation
 *         in the domain (0,1)^2 by the finite difference method
 * 
 * 
 */

#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <cstring>

using namespace Eigen;
using SpMat = SparseMatrix<double>;
using T = Triplet<double>;

class Poisson {
 public:
  // Boundary Conditions
  // Dirichlet : u(0,y), u(1,y), u(x,0), u(x,1)
  // Neumann   : u'(0,y), u'(1,y), u'(x,0), u'(x,1)
  // M1 : u'(0,y), u(1,y), u(x,0), u(x,1) (Left Edge)
  // M2 : u(0,y), u'(1,y), u(x,0), u(x,1) (Right Edge)
  // M3 : u(0,y), u(1,y), u'(x,0), u(x,1) (Bottom Edge)
  // M4 : u(0,y), u(1,y), u(x,0), u'(x,1) (Top Edge)
  // M5 : u(0,y), u(1,y), u'(x,0), u'(x,1) (Neumann for Bottom and Top Edge)
  // M6 : u'(0,y), u'(1,y), u(x,0), u(x,1) (Neumann for Left and Right Edge)
  // M7 : u(0,y), u'(1,y), u(x,0), u'(x,1) (Neumann for Right and Top Edge)
  // M8 : u'(0,y), u(1,y), u'(x,0), u(x,1) (Neumann for Left and Bottom Edge)
  // M9 : u'(0,y), u(1,y), u(x,0), u'(x,1) (Neumann for Left and Top Edge)
  // M10 : u(0,y), u'(1,y), u'(x,0), u(x,1) (Neumann for Right and Bottom Edge)
  // M11 : u(0,y), u'(1,y), u'(x,0), u'(x,1) (Dirichlet for Left Edge)
  // M12 : u'(0,y), u(1,y), u'(x,0), u'(x,1) (Dirichlet for Right Edge)
  // M13 : u'(0,y), u'(1,y), u(x,0), u'(x,1) (Dirichlet for Bottom Edge)
  // M14 : u'(0,y), u'(1,y), u'(x,0), u(x,1) (Dirichlet for Top Edge)

  /** 
   * Constructor for Dirichlet boundary condition
   * 
   * @param u function u, used to get boundary conditions and real solution to
   *          calculate errors.
   * @param f RHS of Poisson equation
   * @param _m number of grids
   */
  Poisson(double (*u)(double, double), double (*f)(double, double), int _m);

  /** 
   * Constructor for pure Neumann boundary condition
   * 
   * @param u function u, used to get boundary conditions and real solution to
   *          calculate errors.
   * @param f RHS of Poisson equation
   * @param _m number of grids
   * @param u_x partial derivative of u with respect to x 
   * @param u_y partial derivative of u with respect to y
   */
  Poisson(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  
  /** 
   * Constructor for M1 - M6 boundary conditions
   * 
   * @param u function u, used to get boundary conditions and real solution to
   *          calculate errors
   * @param f RHS of Poisson equation
   * @param _m number of grids
   * @param BC type of mixed boundary condition
   * @param g corresponding boundary condition (either u_x or u_y)
   */
  Poisson(double (*u)(double, double), double (*f)(double, double), int _m, std::string BC, double (*g)(double, double));

  /** 
   * Constructor for M7 - M14 boundary conditions
   * 
   * @param u function u, used to get boundary conditions and real solution to
   *          calculate errors
   * @param f RHS of Poisson equation
   * @param _m number of grids
   * @param BC type of mixed boundary condition
   * @param u_x partial derivative of u with respect to x
   * @param u_y partial derivative of u with respect to y
   */
  Poisson(double (*u)(double, double), double (*f)(double, double), int _m, std::string BC, double (*u_x)(double, double), double (*u_y)(double, double));

  // printing functions

  /** 
   * Print matrix of system of linear equations
   * 
   */
  void print_matrix() const { std::cout << A << std::endl << std::endl; };

  /** 
   * Print RHS of system of linear equations
   * 
   */
  void print_rhs() const { std::cout << F << std::endl << std::endl; }

  /** 
   * Print real solution of system of linear equations
   * 
   */
  void print_real_sol() const { std::cout << X << std::endl << std::endl; }

  /** 
   * Print solutions to use in MATLAB
   * 
   */
  void print_sol4matlab() const
  {
    int size = 0;
    if(neumann)
      size = U.size() - 1;
    else
      size = U.size();
    
    std::cout << "[";
    for(int i = 0; i < size; ++i)
    {
      if((i+1)%M == 0)
	std::cout << U(i) << ";";
      else if(i%8 == 0 && i != 0)
	std::cout << U(i) << ",..." << std::endl;
      else
	std::cout << U(i) << ",";
    }
    std::cout << "];" << std::endl;
  }

  /** 
   * Solve system of linear equations, and calculate norms of error.
   * 
   * 
   * @return computed solution
   */
  VectorXd solve();
  
 private:
  int m; 
  int M;       // number of unknowns for x_i
  SpMat A;
  VectorXd U;  // to be solved
  VectorXd F;  // RHS
  VectorXd X;  // real solution
  bool neumann = false;  // true if BC is neumann

  void resize_system(int m);

  /************************************************************************************************************************************/

  // Functions for mixed boundary conditions

  void M1(double (*u)(double, double), double (*f)(double, double), int _m, double (*g)(double, double));
  void M2(double (*u)(double, double), double (*f)(double, double), int _m, double (*g)(double, double));
  void M3(double (*u)(double, double), double (*f)(double, double), int _m, double (*g)(double, double));
  void M4(double (*u)(double, double), double (*f)(double, double), int _m, double (*g)(double, double));
  void M5(double (*u)(double, double), double (*f)(double, double), int _m, double (*g)(double, double));
  void M6(double (*u)(double, double), double (*f)(double, double), int _m, double (*g)(double, double));
  
  void M7(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  void M8(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  void M9(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  void M10(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  void M11(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  void M12(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  void M13(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
  void M14(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double));
};
