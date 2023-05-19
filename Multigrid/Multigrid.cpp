#ifndef __MULTIGRID_CPP__
#define __MULTIGRID_CPP__

#include "Multigrid.h"

template class Multigrid<1>;

template<int dim>
void Multigrid<dim>::resize_system(int m)
{
  A.resize(m, m);
  U.resize(m, 1);
  F.resize(m, 1);
  X.resize(m, 1);
}

template<int dim>
void Multigrid<dim>::make_matrix(int m, std::string _BC)
{
  std::vector<T> coef;
  //if(_BC == "Dirichlet")
  //{
    A.resize(m,m);
    h = 1.0/(m+1);
    for(int i = 0; i < m; ++i)
    {
      coef.push_back(T(i,i,2/(h*h)));
      if(i != m-1)
	{
          coef.push_back(T(i,i+1,-1/(h*h)));
	  coef.push_back(T(i+1,i,-1/(h*h)));
	}
    }
    //}

    /*
  else
  {
    A.resize(m,m);
    h = 1.0/(m-1);
    if(_BC == "PureNeumann" || _BC == "LeftNeumann")
    {
      coef.push_back(T(0, 0, -1/h));
      coef.push_back(T(0, 1, 1/h));
    }
    else
    {
      coef.push_back(T(0, 0, 1.0));
      coef.push_back(T(0, 1, 0.0));
    }
    coef.push_back(T(1, 0, 1/(h*h)));
    for(int i = 1; i < m-1; ++i)
    {
      coef.push_back(T(i,i,-2/(h*h)));
      if(i != m-2)
      {
        coef.push_back(T(i,i+1,1/(h*h)));
	coef.push_back(T(i+1,i,1/(h*h)));
      }
    }
    coef.push_back(T(m-2, m-1, 1/(h*h)));
    if(_BC == "LeftNeumann")
    {
      coef.push_back(T(m-1, m-2, 0.0));
      coef.push_back(T(m-1, m-1, 1.0));
    }
    else
    {
      coef.push_back(T(m-1, m-2, 1/h));
      coef.push_back(T(m-1, m-1, -1/h));
    }
  }
    */
  
  A.setFromTriplets(coef.begin(), coef.end());
}

template<int dim>
Multigrid<dim>::Multigrid(Real (*u)(Real), Real(*f)(Real), int _n)
{
  n = _n;
  resize_system(n-1);
  h = 1.0/n;
  level = 0;
  COARSEST = log2(n) - 2;
  BC = "Dirichlet";

  // compute vector of real solution
  int k = 0;
  for(int i = 1; i <= n-1; ++i)
  {
    X(k) = u(i*h);
    ++k;
  }

  k = 0;
  for(int i = 1; i <= n-1; ++i)
  {
    double sum = 0;
    if(i-1 == 0)
      sum += u(0);
    if(i+1 == n)
      sum += u(1);
    F(k) = f(i*h) + (sum/(h*h));
    ++k;
  }

  Grids.resize(COARSEST+1);
  for(int i = 0; i <= COARSEST; ++i)
  {
    Grids[i].V.resize(n/(i+1),1);
    Grids[i].F.resize(n/(i+1),1);
  }
  Grids[0].V = MatrixXd::Zero(n-1,1);
  
}

template<int dim>
Multigrid<dim>::Multigrid(Real (*u)(Real), Real (*u_x)(Real), Real (*f)(Real), int _n, std::string _BC)
{
    n = _n+1;
    resize_system(n);
    h = 1.0/_n;
    level = 0;
    COARSEST = log2(n) - 1;
    BC = _BC;

    // compute vector of real solution
    int k = 0;
    for(int i = 0; i <= n-1; ++i)
      {
	X(k) = u(i*h);
	++k;
      }

    k = 0;
    for(int i = 0; i <= n-1; ++i)
      {
	double sum = 0;
	if(i == 0)
        {
	  if(_BC == "LeftNeumann" || _BC == "PureNeumann")
	    F(k) = u_x(0) + h/2*f(0);
	  else if(_BC == "RightNeumann")
	    F(k) = u(0);
	}
	else if(i == n-1)
	{
	  if(_BC == "LeftNeumann")
	    F(k) = u(1);
	  else
	    F(k) = -u_x(1) + h/2*f(1);
	}
	else
	  F(k) = f(i*h);
	++k;
      }

    Grids.resize(COARSEST+1);
    for(int i = 0; i <= COARSEST; ++i)
      {
	Grids[i].V.resize(n/(i+1),1);
	Grids[i].F.resize(n/(i+1),1);
      }
    Grids[0].V = MatrixXd::Zero(n,1);
}

template<int dim>
Vec Multigrid<dim>::relax(Vec _U, Vec _F, int v1, std::string _BC)
{
  int m = _U.rows();
  make_matrix(m, _BC);
  MatrixXd T(m,m);
  T = MatrixXd::Identity(m,m) - ((2.0/3)*h*h/2) * A;
  MatrixXd Dinv = MatrixXd::Zero(m,m);
  for(int i = 0; i < m; ++i)
  {
    Dinv(i,i) = 1 / A.coeff(i,i);
  }
  
  // iterate
  for(int i = 1; i <= v1; ++i)
  {
    _U = T*_U + (2.0/3)*Dinv*_F;
  }

  return _U;
}

template<int dim>
Vec Multigrid<dim>::injection(Vec& R)
{
  Vec coarseR(R.rows());
  for(int i = 0; i < R.rows()/2 - 1; ++i)
  {
    coarseR(i) = R(2*i);
  }
  //level++;

  return coarseR;
}

template<int dim>
Vec Multigrid<dim>::full_weighting(Vec& R)
{
  SpMat I;
  std::vector<T> coef;
  if(BC == "Dirichlet")
  {
    I.resize((R.rows() + 1)/2 - 1, R.rows());
    for(int i = 0; i < (R.rows()+1)/2 - 1; ++i)
    {
      coef.push_back(T(i, 2*i, 1.0/4));
      coef.push_back(T(i, 2*i + 1, 2.0/4));
      coef.push_back(T(i, 2*i + 2, 1.0/4));
    }
  }
  else
  {
    I.resize(R.rows() / 2 + 1, R.rows());
    coef.push_back(T(0, 0, 1.0/2));
    coef.push_back(T(0, 1, 1.0/4));
    for(int i = 1; i < R.rows()/2; ++i)
    {
      coef.push_back(T(i, 2*i - 1, 1.0/4));
      coef.push_back(T(i, 2*i, 2.0/4));
      coef.push_back(T(i, 2*i + 1, 1.0/4));
    }
    coef.push_back(T(R.rows()/2, R.rows() - 2, 1.0/4));
    coef.push_back(T(R.rows()/2, R.rows() - 1, 1.0/2));
  }
  I.setFromTriplets(coef.begin(), coef.end());

  return I*R;
}

template<int dim>
Vec Multigrid<dim>::linear_interp(Vec& E)
{
  SpMat I;
  std::vector<T> coef;
  if(BC == "Dirichlet")
  {
    I.resize(2*(E.rows() + 1) - 1, E.rows());
    for(int i = 0; i < E.rows(); ++i)
      {
	coef.push_back(T(2*i, i, 1.0/2));
	coef.push_back(T(2*i + 1, i, 1.0));
	coef.push_back(T(2*i + 2, i, 1.0/2));
      }
  }
  else
  {
    I.resize(2*E.rows() - 1, E.rows());
    coef.push_back(T(0, 0, 1.0));
    coef.push_back(T(1, 0, 1.0/2));
    for(int i = 1; i < E.rows()-1; ++i)
    {
      coef.push_back(T(2*i - 1, i, 1.0/2));
      coef.push_back(T(2*i, i, 1.0));
      coef.push_back(T(2*i + 1, i, 1.0/2));
    }
    coef.push_back(T(2*E.rows() - 3, E.rows() - 1, 1.0/2));
    coef.push_back(T(2*E.rows() - 2, E.rows() - 1, 1.0));
  }
  I.setFromTriplets(coef.begin(), coef.end());

  return I*E;
}

template<int dim>
Vec Multigrid<dim>::VC(Vec _V, Vec _F, int v1, int v2)
{
  //if(level == 0)
  //Grids[level].V = relax(_V, _F, v1, BC);
  //else
    Grids[level].V = relax(_V, _F, v1, "Dirichlet");
  Grids[level].F = _F;
  
  if(level == COARSEST)
  {
    SparseLU<SpMat, COLAMDOrdering<int> > solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    Grids[level].V = solver.solve(Grids[level].F);
    return Grids[level].V;
  }
  else
  {
    level++;
    Vec R = Grids[level-1].F - A*Grids[level-1].V;
    Grids[level].F = full_weighting(R);
    Grids[level].V = MatrixXd::Zero(Grids[level].F.rows(), 1);
    Grids[level].V = VC(Grids[level].V, Grids[level].F, v1, v2);

    level--;
    Grids[level].V = Grids[level].V + linear_interp(Grids[level+1].V);
    //if(level == 0)
    //Grids[level].V = relax(Grids[level].V, Grids[level].F, v2, BC);
    //else
      Grids[level].V = relax(Grids[level].V, Grids[level].F, v2, "Dirichlet");
  }

  return Grids[level].V;
}

template<int dim>
Vec Multigrid<dim>::FMG(Vec _F, int v1, int v2)
{
  Grids[level].F = _F;
  if(level == COARSEST)
  {
    Grids[level].V = MatrixXd::Zero(Grids[level].F.rows(), 1);
    Grids[level].V = VC(Grids[level].V, Grids[level].F, v1, v2);
    return Grids[level].V;
  }
  else
  {
    level++;
    Grids[level].F = full_weighting(Grids[level-1].F);
    Grids[level].V = FMG(Grids[level].F, v1, v2);
    
    level--;
    Grids[level].V = linear_interp(Grids[level+1].V);
    Grids[level].V = VC(Grids[level].V, Grids[level].F, v1, v2);
  }
  return Grids[level].V;
}

template<int dim>
Vec Multigrid<dim>::solve(Real rtol, int it)
{
  for(int i = 0; i <= COARSEST; ++i)
     std::fill(Grids[level].V.begin(), Grids[level].V.end(), 0);
  
  U = FMG(F,3,3);
  Vec R = F - A*U;
  //U = MatrixXd::Zero(n,1);
  res.push_back(R.lpNorm<Infinity>());
  std::cout << "inf norm of Residual = " << res[0] << std::endl;
  
  
  int cycle = 1;
  while(R.lpNorm<Infinity>() >= rtol && cycle <= it)
  {
    std::cout << "V-cycle " << cycle << std::endl;
    U = VC(U,F,3,3);
    R = F - A*U;
    res.push_back(R.lpNorm<Infinity>());
    resrate.push_back((res[cycle-1] - res[cycle]) / res[cycle-1]);
    std::cout << "inf norm of Residual = " << res[cycle] << std::endl;
    std::cout << "Reduction rate = " << resrate[cycle-1] << std::endl;
    cycle++;
    }
  std::cout << "Maximum norm of error vector : " << (X-U).lpNorm<Infinity>() << std::endl;
  
  return U;
 
}

#endif
