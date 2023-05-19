/**
 * @file   Multigrid.h
 * @author Luzern Yuven Luis <student@student>
 * @date   Mon Apr 25 17:55:18 2022
 * 
 * @brief  A solver for Poisson equation using Multigrid Method.
 * 
 * 
 */


#ifndef __MULTIGRID_H__
#define __MULTIGRID_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <cstring>

using Real = double;
using namespace Eigen;
using SpMat = SparseMatrix<Real>;
using T = Triplet<Real>;
using Vec = Matrix<Real, Dynamic, 1>;

struct Grid
{
  Vec V;
  Vec F;
};

template<int dim>
class Multigrid
{
 public:
  /** 
   * Constructor for Dirichlet boundary condition
   * 
   * @param u analytic function
   * @param f RHS function
   * @param _n number of grids
   */
  Multigrid(Real (*u)(Real), Real (*f)(Real), int _n);

  /** 
   * Constructor for Neumann and mixed boundary conditions
   * 
   * @param u    analytic function
   * @param u_x  first derivative of u
   * @param f    RHS function
   * @param _n   number of grids
   * @param _BC  type of boundary condition ("PureNeumann", "LeftNeumann",
   *             "RightNeumann")
   */
  Multigrid(Real (*u)(Real), Real (*u_x)(Real), Real (*f)(Real), int _n, std::string _BC);

  // const functions
  /** 
   * Print matrix of the system
   * 
   */
  void print_matrix() const { std::cout << A << std::endl;}

  /** 
   * Print RHS of the system
   * 
   */
  void print_rhs() const { std::cout << F << std::endl; }

  /** 
   * Print analytic solution of the Poisson equation
   * 
   */
  void print_real_sol() const { std::cout << X << std::endl; }

  /** 
   * Print solution of system
   * 
   */
  void print_computed_rhs() const { std::cout << A*U << std::endl; }

  /** 
   * Print residuals of each V-cycle
   * 
   */
  void print_res() const
  {
    std::cout << "[";
    for(int i = 0; i < res.size(); ++i)
    {
      std::cout << res[i] << "; ";
      if(i%4 == 0 && i != 0)
	std::cout << "..." << std::endl;
    }
    std::cout << "];";
  }

  /** 
   * Print reduction rate of each V-cycle
   * 
   */
  void print_resrate() const
  {
    std::cout << "[";
    for(int i = 0; i < resrate.size(); ++i)
    {
      std::cout << resrate[i] << "; ";
      if(i%4 == 0 && i != 0)
	std::cout << "..." << std::endl;
    }
    std::cout << "];";
  }

  /** 
   * Solve the system of linear equations
   * 
   * @param rtol tolerance of residual
   * @param it   maximum number of iterations
   * 
   * @return solution of system
   */
  Vec solve(Real rtol, int it);

 private:
  int n;
  int COARSEST;                    // number of grids in the coarsest grid
  double h;
  SpMat A;
  Vec U;                           // computed solution
  Vec F;                           // RHS of system
  Vec X;                           // analytic solution
  int level;                       // level of V-cycle
  std::vector<Grid> Grids;         // vector to store U and F at each level
  std::string BC;                  // type of boundary condition
  std::vector<double> res;         // vector to store residual of each V-cycle
  std::vector<double> resrate;     // vector to store reduction rate of each V-cycle

  void resize_system(int m);
  /** 
   * Construct matrix
   * 
   * @param m    matrix size
   */
  void make_matrix(int m, std::string _BC);

  /** 
   * Relax _U with matrix
   * 
   * @param _U   initial U
   * @param _F   RHS 
   * @param v1   number of relaxations
   * 
   * @return relaxed U
   */
  Vec relax(Vec _U, Vec _F, int v1, std::string _BC);

  // restriction operators
  Vec injection(Vec& R);
  Vec full_weighting(Vec& R);

  // prolongation / interpolation operator
  Vec linear_interp(Vec& E);

  Vec VC(Vec _V, Vec _F, int v1, int v2);
  Vec FMG(Vec _F, int v1, int v2);
  
};

/*
template<>
class Multigrid<2>
{
 public:
  Multigrid(Real (*u)(Real, Real), Real(*f)(Real, Real), int _n);
}

inline Multigrid<2>::Multigrid(Real (*u)(Real, Real), Real (*f)(Real, Real), int _n)
{
  std::cout << "2D" << std::endl;
}
*/

#endif
