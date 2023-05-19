#ifndef _POISSON_CPP_
#define _POISSON_CPP_

#include "Poisson.h"

void Poisson::resize_system(int m)
{
  A.resize(m*m, m*m);
  U.resize(m*m);
  F.resize(m*m);
  X.resize(m*m);
}

Poisson::Poisson(double (*u)(double, double), double (*f)(double, double), int _m)
{
  m = _m;
  M = m;
  resize_system(m);
  double h = 1.0/(_m+1);

  std::vector<T> coef;
  for(int i = 0; i < m*m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
      coef.push_back(T(j,j,4/(h*h)));
      // generate the first sub and superdiagonal
      if(j != i+m-1)
      {
	coef.push_back(T(j+1,j,-1/(h*h)));
	coef.push_back(T(j,j+1,-1/(h*h)));
      }
    }
  }
  // generate the m-th sub and superdiagonal
  for(int i = 0 ; i < m*m - m; ++i)
  {
    coef.push_back(T(i+m,i,-1/(h*h)));
    coef.push_back(T(i,i+m,-1/(h*h)));
  }
  
  A.setFromTriplets(coef.begin(), coef.end());
  
  // compute vector of real solution
  int k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      double sum = 0;
      if(i-1 == 0)
	sum += u(0,j*h);
      if(i+1 == m+1)
	sum += u(1,j*h);
      if(j-1 == 0)
	sum += u(i*h,0);
      if(j+1 == m+1)
	sum += u(i*h,1);
      F(k) = f(i*h, j*h) + (sum/(h*h));
      ++k;
    }
  }
}

Poisson::Poisson(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  neumann = true;
  m = _m + 2;
  M = m;
  A.resize(m*m + 1, m*m + 1);
  U.resize(m*m + 1);
  F.resize(m*m + 1);
  X.resize(m*m + 1);
  double h = 1.0/(_m+1);

  std::vector<T> coef;

  for(int i = 0; i < m*m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
        if(j == i)
	{
	  coef.push_back(T(j,j+1,-2/(h*h)));
	  coef.push_back(T(j+1,j,-1/(h*h)));
	}
	else if(j == i+m-2)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else if(j < i+m-1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+2)-th sub and superdiagonal
  for(int i = 0 ; i < m*m - m; ++i)
  {
    if(i >= 0 && i < m)
    {
      coef.push_back(T(i,i+m,-2/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
    else if(i >= m && i <  m*m - m*2)
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-2/(h*h)));
    }
  }

  int k = 0;
  double total = 0;
  for(int j = 0; j <= _m+1; ++j)
  {
    for(int i = 0; i <= _m+1; ++i)
    {
      X(k) = u(i*h,j*h);
      total += X(k);
      ++k;
    }
  }
  double mean = total / (m*m);
  X(m*m) = mean;

  for(int j = 0; j < m*m; ++j)
  {
    coef.push_back(T(j, m*m, 1));
    coef.push_back(T(m*m, j, 1));
  }
  coef.push_back(T(m*m,m*m, 0));
  
  A.setFromTriplets(coef.begin(), coef.end());
  
  k = 0;
  for(int j = 0; j <= _m+1; ++j)
  {
    for(int i = 0; i <= _m+1; ++i)
    {
      double sum = 0;
      if(i == 0 && j == 0)
	sum += -2*u_x(0,0)/h - 2*u_y(0,0)/h;
      else if(i == 0 && j == _m+1)
	sum += -2*u_y(0,1)/h + 2*u_x(0,1)/h;
      else if(i == _m+1 && j == 0)
	sum += -2*u_y(1,0)/h + 2*u_x(1,0)/h;
      else if(i == _m+1 && j == _m+1)
	sum += 2*u_x(1,1)/h + 2*u_y(1,1)/h;
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h;
      else if(i == _m+1)
	sum += 2*u_x(1,j*h)/h;
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h;
      else if(j == _m+1)
	sum += 2*u_y(i*h,1)/h;

      F(k) = f(i*h,j*h) + sum - mean;
      ++k;
    }
  }
  F(m*m) = total;  
}


Poisson::Poisson(double (*u)(double, double), double(*f)(double, double), int _m, std::string BC, double (*g)(double, double))
{
  if(BC == "M1") M1(u, f, _m, g);
  else if(BC == "M2") M2(u, f, _m, g);
  else if(BC == "M3") M3(u, f, _m, g);
  else if(BC == "M4") M4(u, f, _m, g);
  else if(BC == "M5") M5(u, f, _m, g);
  else if(BC == "M6") M6(u, f, _m, g);
}

Poisson::Poisson(double (*u)(double, double), double (*f)(double, double), int _m, std::string BC, double (*u_x)(double, double), double (*u_y)(double, double))
{
  if(BC == "M7") M7(u, f, _m, u_x, u_y);
  else if(BC == "M8") M8(u, f, _m, u_x, u_y);
  else if(BC == "M9") M9(u, f, _m, u_x, u_y);
  else if(BC == "M10") M10(u, f, _m, u_x, u_y);
  else if(BC == "M11") M11(u, f, _m, u_x, u_y);
  else if(BC == "M12") M12(u, f, _m, u_x, u_y);
  else if(BC == "M13") M13(u, f, _m, u_x, u_y);
  else if(BC == "M14") M14(u, f, _m, u_x, u_y);
}

void Poisson::M1(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double))
{
  m = _m;
  M = m + 1;
  A.resize(m*m + m, m*m +m);
  U.resize(m*m + m);
  F.resize(m*m + m);
  X.resize(m*m + m);
  double h = 1.0/(m+1);
 
  std::vector<T> coef;
  for(int i = 0; i < m*m; i+=(m+1))
  {
    for(int j = i; j < i+m+1; ++j)
    {
      coef.push_back(T(j,j,4/(h*h)));
      // generate the first sub and superdiagonal
      if(j != i+m)
      {
	if(j == i)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-2/(h*h)));
	}
	else
        {
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
      }
    }
  }
  // generate the m-th sub and superdiagonal
  for(int i = 0 ; i < m*m - 1; ++i)
  {
    coef.push_back(T(i+m+1,i,-1/(h*h)));
    coef.push_back(T(i,i+m+1,-1/(h*h)));
  }
  
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 0; i <= m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 0; i <= m; ++i)
    {
      double sum = 0;
      // RHS for left edge of (0,1)x(0,1) (Neumann)
      if(i == 0 && j-1 == 0)
	sum += -2*u_x(0,0)/h + u(0,0)/(h*h); // use Dirichlet
      else if(i == 0 && j+1 == m+1)
	sum += -2*u_x(0,1)/h + u(0,1)/(h*h); // use Dirichlet
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h;

      // RHS for the other edges (Dirichlet)
      if(i+1 == m+1)             // right edge
	sum += u(1,j*h)/(h*h); 
      if(j-1 == 0 && i != 0)   // bottom edge
	sum += u(i*h,0)/(h*h);
      if(j+1 == m+1 && i != 0) // top edge
	sum += u(i*h,1)/(h*h);
      
      F(k) = f(i*h,j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M2(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double))
{
  m = _m;
  M = m + 1;
  A.resize(m*m + m, m*m + m);
  U.resize(m*m + m);
  F.resize(m*m + m);
  X.resize(m*m + m);
  double h = 1.0/(m+1);

  std::vector<T> coef;
  for(int i = 0; i < m*m; i+=(m+1))
  {
    for(int j = i; j < i+m+1; ++j)
    {
      coef.push_back(T(j,j,4/(h*h)));
      // generate the first sub and superdiagonal
      if(j != i+m)
      {
	if(j == i+m-1)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else
        {
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
      }
    }
  }
  // generate the m-th sub and superdiagonal
  for(int i = 0 ; i < m*m - 1; ++i)
  {
    coef.push_back(T(i+m+1,i,-1/(h*h)));
    coef.push_back(T(i,i+m+1,-1/(h*h)));
  }
  
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 1; i <= m+1; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 1; i <= m+1; ++i)
    {
      double sum = 0;
      // RHS for right edge of (0,1)x(0,1) (Neumann)
      if(i == m+1 && j-1 == 0)
	sum += 2*u_x(1,0)/h + u(1,0)/(h*h); // point (1,0)
      else if(i == m+1 && j+1 == m+1)
	sum += 2*u_x(1,1)/h + u(1,1)/(h*h); // point (1,1)
      else if(i == m+1)
	sum += 2*u_x(1,j*h)/h;

      // RHS for the other edges (Dirichlet)
      if(i-1 == 0)                 // left edge
	sum += u(0,j*h)/(h*h);
      if(j-1 == 0 && i != m+1)   // bottom edge
	sum += u(i*h,0)/(h*h);
      if(j+1 == m+1 && i != m+1) // top edge
	sum += u(i*h,1)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M3(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_y)(double, double))
{
  m = _m;
  M = m;
  A.resize(m*m + m, m*m + m);
  U.resize(m*m + m);
  F.resize(m*m + m);
  X.resize(m*m + m);
  double h = 1.0/(m+1);
  
  std::vector<T> coef;
  for(int i = 0; i < m*m + m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
      coef.push_back(T(j,j,4/(h*h)));
      // generate the first sub and superdiagonal
      if(j != i+m-1)
      {
	coef.push_back(T(j+1,j,-1/(h*h)));
	coef.push_back(T(j,j+1,-1/(h*h)));
      }
    }
  }
  // generate the m-th sub and superdiagonal
  for(int i = 0 ; i < m*m; ++i)
  {
    coef.push_back(T(i+m,i,-1/(h*h)));
    if(i >= 0 && i < m)
      coef.push_back(T(i,i+m,-2/(h*h)));
    else
      coef.push_back(T(i,i+m,-1/(h*h)));
  }
  
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 0; j <= m; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 0; j <= m; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      double sum = 0;
      // RHS for bottom edge of (0,1)x(0,1) (Neumann)
      if(j == 0 && i-1 == 0)
	sum += -2*u_y(0,0)/h + u(0,0)/(h*h); // point (0,0), use Dirichlet
      else if(j == 0 && i+1 == m+1)
	sum += -2*u_y(1,0)/h + u(1,0)/(h*h); // point (1,0), use Dirichlet
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h;

      // RHS for the other edges (Dirichlet)
      if(i-1 == 0 && j != 0)    // left edge - (0,0)
	sum += u(0,j*h)/(h*h);
      if(i+1 == m+1 && j != 0)  // right edge - (1,0)
	sum += u(1,j*h)/(h*h);
      if(j+1 == m+1)            // top edge
	sum += u(i*h,1)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M4(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_y)(double, double))
{
  m = _m;
  M = m;
  A.resize(m*m + m, m*m + m);
  U.resize(m*m + m);
  F.resize(m*m + m);
  X.resize(m*m + m);
  double h = 1.0/(m+1);

  std::vector<T> coef;
  for(int i = 0; i < m*m + m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
      coef.push_back(T(j,j,4/(h*h)));
      // generate the first sub and superdiagonal
      if(j != i+m-1)
      {
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
      }
    }
  }
  // generate the m-th sub and superdiagonal
  for(int i = 0 ; i < m*m; ++i)
  {
    if(i >= m*m - m && i < m*m)
      coef.push_back(T(i+m,i,-2/(h*h)));
    else
      coef.push_back(T(i+m,i,-1/(h*h)));
    
    coef.push_back(T(i,i+m,-1/(h*h)));
  }
  
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 1; j <= m+1; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 1; j <= m+1; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      double sum = 0;
      // RHS for top edge of (0,1)x(0,1) (Neumann)
      if(j == m+1 && i-1 == 0)
	sum += 2*u_y(0,1)/h + u(0,1)/(h*h); // point (0,1), use Dirichlet
      else if(j == m+1 && i+1 == m+1)
	sum += 2*u_y(1,1)/h + u(1,1)/(h*h); // point (1,1), use Dirichlet
      else if(j == m+1)
	sum += 2*u_y(i*h,1)/h;

      // RHS for the other edges (Dirichlet)
      if(i-1 == 0 && j != m+1)   // left edge - (0,1)
	sum += u(0,j*h)/(h*h);
      if(i+1 == m+1 && j != m+1) // right edge - (1,1) 
	sum += u(1,j*h)/(h*h);
      if(j-1 == 0)               // bottom edge
	sum += u(i*h,0)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M5(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_y)(double, double))
{
  m = _m;
  double h = 1.0/(m+1);
  // U_ij for i = 1,...,m , j = 0,...,(m+1)
  M = m;
  A.resize(m*(m+2), m*(m+2));
  U.resize(m*(m+2));
  F.resize(m*(m+2));
  X.resize(m*(m+2));

  std::vector<T> coef;
  for(int i = 0; i < m*(m+2); i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
      coef.push_back(T(j,j,4/(h*h)));
      // generate the first sub and superdiagonal
      if(j != i+m-1)
      {
	coef.push_back(T(j+1,j,-1/(h*h)));
	coef.push_back(T(j,j+1,-1/(h*h)));
      }
    }
  }
  // generate the (m+2)-th sub and superdiagonal
  for(int i = 0 ; i < m*(m+2) - m; ++i)
  {
    if(i >= 0 && i < m)
    {
      coef.push_back(T(i+m,i,-1/(h*h)));
      coef.push_back(T(i,i+m,-2/(h*h)));
    }
    else if(i >= m*(m+2) - 2*m && i < m*(m+2) - m)
    {
      coef.push_back(T(i+m,i,-2/(h*h)));
      coef.push_back(T(i,i+m,-1/(h*h)));
    }
    else
    {
      coef.push_back(T(i+m,i,-1/(h*h)));
      coef.push_back(T(i,i+m,-1/(h*h)));
    }
  }
  
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 0; j <= m+1; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 0; j <= m+1; ++j)
  {
    for(int i = 1; i <= m; ++i)
    {
      double sum = 0;
      // RHS for top and bottom edge of (0,1)x(0,1) (Neumann)
      if(i-1 == 0 && j == 0)
	sum += -2*u_y(0,0)/h + u(0,0)/(h*h); // point (0,0), use Dirichlet
      else if(i-1 == 0 && j == m+1)
	sum += 2*u_y(0,1)/h + u(0,1)/(h*h); // point (0,1), use Dirichlet
      else if(i+1 == m+1 && j == 0)
	sum += -2*u_y(1,0)/h + u(1,0)/(h*h); // point (1,0), use Dirichlet
      else if(i+1 == m+1 && j == m+1)
	sum += 2*u_y(1,1)/h + u(1,1)/(h*h); // point (1,1), use Dirichlet
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h; // bottom edge
      else if(j == m+1)
	sum += 2*u_y(i*h,1)/h;  // top edge

      // RHS for the other edges (Dirichlet)
      if(i-1 == 0 && j != 0 && j != m+1)   // left edge
	sum += u(0,j*h)/(h*h);
      if(i+1 == m+1 && j != 0 && j != m+1) // right edge 
	sum += u(1,j*h)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M6(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double))
{
  m = _m;
  double h = 1.0/(m+1);
  // U_ij for i = 0,...,m+1 , j = 1,...,m
  M = m + 2;
  A.resize(m*(m+2), m*(m+2));
  U.resize(m*(m+2));
  F.resize(m*(m+2));
  X.resize(m*(m+2));

  std::vector<T> coef;

  for(int i = 0; i < m*(m+2); i+=(m+2))
  {
    for(int j = i; j < i+m+2; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
        if(j == i)
	{
	  coef.push_back(T(j,j+1,-2/(h*h)));
	  coef.push_back(T(j+1,j,-1/(h*h)));
	}
	else if(j == i+m)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else if(j < i+m+1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate the (m+2)-th sub and superdiagonal
  for(int i = 0 ; i < m*(m+2) - (m+2); ++i)
  {
    coef.push_back(T(i+m+2,i,-1/(h*h)));
    coef.push_back(T(i,i+m+2,-1/(h*h)));
  }
  
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 0; i <= m+1; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }
  
  k = 0;
  for(int j = 1; j <= m; ++j)
  {
    for(int i = 0; i <= m+1; ++i)
    {
      double sum = 0;
      // RHS for left and right edge of (0,1)x(0,1) (Neumann)
      if(i == 0 && j-1 == 0)
	sum += -2*u_x(0,0)/h + u(0,0)/(h*h); // point (0,0), use Dirichlet
      else if(i == 0 && j+1 == m+1)
	sum += -2*u_x(0,1)/h + u(0,1)/(h*h); // point (0,1), use Dirichlet
      else if(i == m+1 && j-1 == 0)
	sum += 2*u_x(1,0)/h + u(1,0)/(h*h); // point (1,0), use Dirichlet
      else if(i == m+1 && j+1 == m+1)
	sum += 2*u_x(1,1)/h + u(1,1)/(h*h); // point (1,1), use Dirichlet
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h; // left edge
      else if(i == m+1)
	sum += 2*u_x(1,j*h)/h;  // right edge

      // RHS for the other edges (Dirichlet)
      if(j-1 == 0 && i != 0 && i != m+1)   // bottom edge
	sum += u(i*h, 0)/(h*h);
      if(j+1 == m+1 && i != 0 && i != m+1) // top edge 
	sum += u(i*h,1)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M7(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double) , double (*u_y)(double, double))
{
  double h = 1.0/(_m+1);
  // unknowns are U_ij for i,j = 1,...,(m+1)
  m = _m + 1;
  M = m;
  resize_system(m);
  
  std::vector<T> coef;

  for(int i = 0; i < m*m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
	if(j == i+m-2)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else if(j < i+m-1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+1)-th sub and superdiagonal
  for(int i = 0 ; i < m*m - m; ++i)
  {
    if (i >= m*m - m*2 && i < m*m - m)
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-2/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 1; j <= _m+1; ++j)
  {
    for(int i = 1; i <= _m+1; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 1; j <= _m+1; ++j)
  {
    for(int i = 1; i <= _m+1; ++i)
    {
      double sum = 0;
      // RHS for right and top edge of (0,1)x(0,1) (Neumann)
      if(i-1 == 0 && j == _m+1)
	sum += 2*u_y(0,1)/h + u(0,1)/(h*h); // point (0,1), use Dirichlet
      else if(i == _m+1 && j-1 == 0)
	sum += 2*u_x(1,0)/h + u(1,0)/(h*h); // point (1,0), use Dirichlet
      else if(i == _m+1 && j == _m+1)
	sum += 2*u_x(1,1)/h + 2*u_y(1,1)/h;   // point (1,1), use ghost cell
      else if(i == _m+1)
	sum += 2*u_x(1,j*h)/h;              // right edge
      else if(j == _m+1)
	sum += 2*u_y(i*h,1)/h;              // top edge

      // RHS for the other edges (Dirichlet)     
      if(j-1 == 0 && i != _m+1)   // bottom edge - (1,0)
	sum += u(i*h,0)/(h*h);
      if(i-1 == 0 && j != _m+1)   // left edge - (0,1)
	sum += u(0,j*h)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M8(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  double h = 1.0/(_m+1);
  // U_ij for i = 0,...,m , j = 0,...,m
  m = _m + 1;
  M = m;
  resize_system(m);

  std::vector<T> coef;

  for(int i = 0; i < m*m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
        if(j == i)
	{
	  coef.push_back(T(j,j+1,-2/(h*h)));
	  coef.push_back(T(j+1,j,-1/(h*h)));
	}
	else if(j < i+m-1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+1)-th sub and superdiagonal
  for(int i = 0 ; i < m*m - m; ++i)
  {
    if(i >= 0 && i < m)
    {
      coef.push_back(T(i,i+m,-2/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 0; j <= _m; ++j)
  {
    for(int i = 0; i <= _m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 0; j <= _m; ++j)
  {
    for(int i = 0; i <= _m; ++i)
    {
      double sum = 0;
      // RHS for left and bottom edge of (0,1)x(0,1) (Neumann)
      if(i == 0 && j == 0)
	sum += -2*u_x(0,0)/h - 2*u_y(0,0)/h;    // point (0,0), use ghost cell
      else if(i+1 == _m+1 && j == 0)
	sum += -2*u_y(1,0)/h + u(1,0)/(h*h);   // point (1,0), use Dirichlet
      else if(i == 0 && j+1 == _m+1)
	sum += -2*u_x(0,1)/h + u(0,1)/(h*h);   // point (0,1), use Dirichlet
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h;                // left edge - (0,0) - (0,1)
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h;                // bottom edge - (0,0) - (1,0)

      // RHS for the other edges (Dirichlet)                  
      if(i+1 == _m+1 && j != 0) // right edge - (1,0)
	sum += u(1,j*h)/(h*h);
      if(j+1 == _m+1 && i != 0) // top edge - (0,1)
	sum += u(i*h,1)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M9(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  double h = 1.0/(_m+1);
  // U_ij for i = 0,...,m , j = 1,...,m+1
  m = _m + 1;
  M = m;
  resize_system(m);
  
  std::vector<T> coef;

  for(int i = 0; i < m*m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
        if(j == i)
	{
	  coef.push_back(T(j,j+1,-2/(h*h)));
	  coef.push_back(T(j+1,j,-1/(h*h)));
	}
	else if(j < i+m-1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+1)-th sub and superdiagonal
  for(int i = 0 ; i < m*m - m; ++i)
  {
    if(i >= m*m - m*2 && i <  m*m - m)
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-2/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());
  
  int k = 0;
  for(int j = 1; j <= _m+1; ++j)
  {
    for(int i = 0; i <= _m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }
  
  k = 0;
  for(int j = 1; j <= _m+1; ++j)
  {
    for(int i = 0; i <= _m; ++i)
    {
      double sum = 0;
      // RHS for left and top edge of (0,1)x(0,1) (Neumann)
      if(i == 0 && j == _m+1)
	sum += 2*u_x(0,1)/h - 2*u_y(0,1)/h;      // point (0,1), use ghost cell
      else if(i == 0 && j-1 == 0)
	sum += -2*u_x(0,0)/h + u(0,0)/(h*h);  // point (0,0), use Dirichlet
      else if(i+1 == _m+1 && j == _m+1)
	sum += 2*u_y(1,1)/h + u(1,1)/(h*h);  // point (1,1), use Dirichlet
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h;              // left edge - (0,0) - (0,1)
      else if(j == _m+1)
	sum += 2*u_y(i*h,1)/h;               // top edge - (0,1) - (1,1)

      // RHS for the other edges (Dirichlet)
      if(i+1 == _m+1 && j != _m+1)  // right edge - (1,1)
	sum += u(1,j*h)/(h*h);
      if(j-1 == 0 && i != 0)     // bottom edge - (0,0)
	sum += u(i*h,0)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M10(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  double h = 1.0/(_m+1);
  // U_ij for i = 1,...,m+1 , j = 0,...,m
  m = _m + 1;
  M = m;
  resize_system(m);
  
  std::vector<T> coef;

  for(int i = 0; i < m*m; i+=m)
  {
    for(int j = i; j < i+m; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
	if(j == i+m-2)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else if(j < i+m-1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+1)-th sub and superdiagonal
  for(int i = 0 ; i < m*m - m; ++i)
  {
    if(i >= 0 && i < m)
    {
      coef.push_back(T(i,i+m,-2/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m,-1/(h*h)));
      coef.push_back(T(i+m,i,-1/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 0; j <= _m; ++j)
  {
    for(int i = 1; i <= _m+1; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 0; j <= _m; ++j)
  {
    for(int i = 1; i <= _m+1; ++i)
    {
      double sum = 0;
      // RHS for right and bottom edge of (0,1)x(0,1) (Neumann)
      if(i == _m+1 && j == 0)
	sum += 2*u_x(1,0)/h - 2*u_y(1,0)/h;      // point (1,0), use ghost cell
      else if(i-1 == 0 && j == 0)
	sum += -2*u_y(0,0)/h + u(0,0)/(h*h);   // point (0,0), use Dirichlet
      else if(i == _m+1 && j+1 == _m+1)
	sum += 2*u_x(1,1)/h + u(1,1)/(h*h);    // point (1,1), use Dirichlet
      else if(i == _m+1)
	sum += 2*u_x(1,j*h)/h;                 // right edge - (1,0) - (1,1)
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h;                // bottom edge - (0,0) - (1,0)

      // RHS for the other edges (Dirichlet)
      if(i-1 == 0 && j != 0)     // left edge - (0,0)
	sum += u(0,j*h)/(h*h);
      if(j+1 == _m+1 && i != _m+1)   // top edge - (1,1)
	sum += u(i*h,1)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M11(double (*u)(double, double), double (*f)(double, double),int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  m = _m;
  M = m + 1;
  double h = 1.0/(m+1);
  // U_ij for i = 1,...,m+1 , j = 0,...,m+1
  A.resize((m+1)*(m+2), (m+1)*(m+2));
  U.resize((m+1)*(m+2));
  F.resize((m+1)*(m+2));
  X.resize((m+1)*(m+2));

  std::vector<T> coef;

  for(int i = 0; i < (m+1)*(m+2); i+=(m+1))
  {
    for(int j = i; j < i+m+1; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
	if(j == i+m-1)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else if(j < i+m)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+1)-th sub and superdiagonal
  for(int i = 0 ; i < (m+1)*(m+2) - (m+1); ++i)
  {
    if(i >= 0 && i < m+1)
    {
      coef.push_back(T(i,i+m+1,-2/(h*h)));
      coef.push_back(T(i+m+1,i,-1/(h*h)));
    }
    else if(i >= m+1 && i <  (m+1)*(m+2) - (m+1)*2)
    {
      coef.push_back(T(i,i+m+1,-1/(h*h)));
      coef.push_back(T(i+m+1,i,-1/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m+1,-1/(h*h)));
      coef.push_back(T(i+m+1,i,-2/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 0; j <= m+1; ++j)
  {
    for(int i = 1; i <= m+1; ++i)
    {
      X(k) = u(i*h, j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 0; j <= m+1; ++j)
  {
    for(int i = 1; i <= m+1; ++i)
    {
      double sum = 0;
      // RHS for right, bottom and top edge of (0,1)x(0,1) (Neumann)
      if(i == m+1 && j == 0)
	sum += 2*u_x(1,0)/h - 2*u_y(1,0)/h;      // point (1,0), use ghost cell
      else if(i == m+1 && j == m+1)
	sum += 2*u_x(1,1)/h + 2*u_y(1,1)/h;     // point (1,1), use ghost cell
      else if(i-1 == 0 && j == 0)
	sum += -2*u_y(0,0)/h + u(0,0)/(h*h);  // point (0,0), use Dirichlet
      else if(i-1 == 0 && j == m+1)
	sum += 2*u_y(0,1)/h + u(0,1)/(h*h);       // point (0,1), use Dirichlet
      else if(i == m+1)
	sum += 2*u_x(1,j*h)/h;          // right edge - (1,0) - (1,1)
      else if(j == m+1)
	sum += 2*u_y(i*h,1)/h;           // top edge - (0,1) - (1,1)
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h;          // bottom edge - (0,0) - (1,0)

      // RHS for the left edge (Dirichlet)
      if(i-1 == 0 && j != 0 && j != m+1)
	sum += u(0,j*h)/(h*h);
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M12(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  m = _m;
  M = m + 1;;
  double h = 1.0/(m+1);
  // U_ij for i = 0,...,m , j = 0,...,m+1
  A.resize((m+1)*(m+2), (m+1)*(m+2));
  U.resize((m+1)*(m+2));
  F.resize((m+1)*(m+2));
  X.resize((m+1)*(m+2));
  
  std::vector<T> coef;

  for(int i = 0; i < (m+1)*(m+2); i+=(m+1))
  {
    for(int j = i; j < i+m+1; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
        if(j == i)
	{
	  coef.push_back(T(j,j+1,-2/(h*h)));
	  coef.push_back(T(j+1,j,-1/(h*h)));
	}
	else if(j < i+m)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+1)-th sub and superdiagonal
  for(int i = 0 ; i < (m+1)*(m+2) - (m+1); ++i)
  {
    if(i >= 0 && i < m+1)
    {
      coef.push_back(T(i,i+m+1,-2/(h*h)));
      coef.push_back(T(i+m+1,i,-1/(h*h)));
    }
    else if(i >= m+1 && i <  (m+1)*(m+2) - (m+1)*2)
    {
      coef.push_back(T(i,i+m+1,-1/(h*h)));
      coef.push_back(T(i+m+1,i,-1/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m+1,-1/(h*h)));
      coef.push_back(T(i+m+1,i,-2/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 0; j <= m+1; ++j)
  {
    for(int i = 0; i <= m; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 0; j <= m+1; ++j)
  {
    for(int i = 0; i <= m; ++i)
    {
      double sum = 0;
      // RHS for left, bottom and top edge of (0,1)x(0,1) (Neumann)
      if(i == 0 && j == 0)
	sum += -2*u_x(0,0)/h - 2*u_y(0,0)/h;      // point (0,0), use ghost cell
      else if(i == 0 && j == m+1)
	sum += -2*u_x(0,1)/h + 2*u_y(0,1)/h;      // point (0,1), use ghost cell
      else if(i+1 == m+1 && j == 0)
	sum += -2*u_y(1,0)/h + u(1,0)/(h*h);  // point (1,0), use Dirichlet
      else if(i+1 == m+1 && j == m+1)
	sum += 2*u_y(1,1)/h + u(1,1)/(h*h);   // point (1,1), use Dirichlet
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h;          // left edge - (0,0) - (0,1)
      else if(j == m+1)
	sum += 2*u_y(i*h,1)/h;           // top edge - (0,1) - (1,1)
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h;          // bottom edge - (0,0) - (1,0)

      // RHS for the right edge (Dirichlet)
      if(i+1 == m+1 && j != 0 && j != m+1)
	sum += u(1,j*h)/(h*h);       // right edge - (1,0) - (1,1)
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M13(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  m = _m;
  M = m + 2;
  double h = 1.0/(m+1);
  // U_ij for i = 0,...,m+1 , j = 1,...,m+1
  A.resize((m+2)*(m+1), (m+2)*(m+1));
  U.resize((m+2)*(m+1));
  F.resize((m+2)*(m+1));
  X.resize((m+2)*(m+1));

  std::vector<T> coef;

  for(int i = 0; i < (m+2)*(m+1); i+=(m+2))
  {
    for(int j = i; j < i+m+2; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
        if(j == i)
	{
	  coef.push_back(T(j,j+1,-2/(h*h)));
	  coef.push_back(T(j+1,j,-1/(h*h)));
	}
	else if(j == i+m)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else if(j < i+m+1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+2)-th sub and superdiagonal
  for(int i = 0 ; i < (m+2)*(m+1) - (m+2); ++i)
  {
    if(i >= (m+1)*(m+2) - (m+2)*2)
    {
      coef.push_back(T(i,i+m+2,-1/(h*h)));
      coef.push_back(T(i+m+2,i,-2/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m+2,-1/(h*h)));
      coef.push_back(T(i+m+2,i,-1/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 1; j <= m+1; ++j)
  {
    for(int i = 0; i <= m+1; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }

  k = 0;
  for(int j = 1; j <= m+1; ++j)
  {
    for(int i = 0; i <= m+1; ++i)
    {
      double sum = 0;
      // RHS for left, right and top edge of (0,1)x(0,1) (Neumann)
      if(i == 0 && j == m+1)
	sum += -2*u_x(0,1)/h + 2*u_y(0,1)/h;      // point (0,1), use ghost cell
      else if(i == m+1 && j == m+1)
	sum += 2*u_x(1,1)/h + 2*u_y(1,1)/h;       // point (1,1), use ghost cell
      else if(i == 0 && j-1 == 0)
	sum += -2*u_x(0,0)/h + u(0,0)/(h*h);  // point (0,0), use Dirichlet
      else if(i == m+1 && j-1 == 0)
	sum += 2*u_x(1,0)/h + u(1,0)/(h*h);     // point (1,0), use Dirichlet
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h;          // left edge - (0,0) - (0,1)
      else if(j == m+1)
	sum += 2*u_y(i*h,1)/h;           // top edge - (0,1) - (1,1)
      else if(i == m+1)
	sum += 2*u_x(1,j*h)/h;          // right edge - (1,0) - (1,1)

      // RHS for the bottom edge (Dirichlet)
      if(j-1 == 0 && i != 0 && i != m+1)
	sum += u(i*h,0)/(h*h);           // bottom edge - (0,0) - (1,0) 
      
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

void Poisson::M14(double (*u)(double, double), double (*f)(double, double), int _m, double (*u_x)(double, double), double (*u_y)(double, double))
{
  m = _m;
  M = m + 2;
  double h = 1.0/(m+1);
  // U_ij for i = 0,...,m+1 , j = 0,...,m
  A.resize((m+1)*(m+2), (m+1)*(m+2));
  U.resize((m+1)*(m+2));
  F.resize((m+1)*(m+2));
  X.resize((m+1)*(m+2));

  std::vector<T> coef;

  for(int i = 0; i < (m+1)*(m+2); i+=(m+2))
  {
    for(int j = i; j < i+m+2; ++j)
    {
	coef.push_back(T(j,j,4/(h*h)));
	// generate first sub and superdiagonal
        if(j == i)
	{
	  coef.push_back(T(j,j+1,-2/(h*h)));
	  coef.push_back(T(j+1,j,-1/(h*h)));
	}
	else if(j == i+m)
	{
	  coef.push_back(T(j+1,j,-2/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
	else if(j < i+m+1)
	{
	  coef.push_back(T(j+1,j,-1/(h*h)));
	  coef.push_back(T(j,j+1,-1/(h*h)));
	}
    }
  }

  // generate (m+2)-th sub and superdiagonal
  for(int i = 0 ; i < (m+1)*(m+2) - (m+2); ++i)
  {
    if(i >= 0 && i < (m+2))
    {
      coef.push_back(T(i,i+m+2,-2/(h*h)));
      coef.push_back(T(i+m+2,i,-1/(h*h)));
    }
    else
    {
      coef.push_back(T(i,i+m+2,-1/(h*h)));
      coef.push_back(T(i+m+2,i,-1/(h*h)));
    }
  } 
  A.setFromTriplets(coef.begin(), coef.end());

  int k = 0;
  for(int j = 0; j <= m; ++j)
  {
    for(int i = 0; i <= m+1; ++i)
    {
      X(k) = u(i*h,j*h);
      ++k;
    }
  }
  
  k = 0;
  for(int j = 0; j <= m; ++j)
  {
    for(int i = 0; i <= m+1; ++i)
    {
      double sum = 0;
      // RHS for left, right and bottom edge of (0,1)x(0,1) (Neumann)
      if(i == 0 && j == 0)
	sum += -2*u_x(0,0)/h - 2*u_y(0,0)/h;     // point (0,0), use ghost cell
      else if(i == m+1 && j == 0)
	sum += 2*u_x(1,0)/h - 2*u_y(1,0)/h;      // point (1,0), use ghost cell
      else if(i == 0 && j+1 == m+1)
	sum += -2*u_x(0,1)/h + u(0,1)/(h*h);  // point (0,1), use Dirichlet
      else if(i == m+1 && j+1 == m+1)
	sum += 2*u_x(1,1)/h + u(1,1)/(h*h);    // point (1,1), use Dirichlet
      else if(i == 0)
	sum += -2*u_x(0,j*h)/h;          // left edge - (0,0) - (0,1)
      else if(j == 0)
	sum += -2*u_y(i*h,0)/h;           // bottom edge - (0,0) - (1,0)
      else if(i == m+1)
	sum += 2*u_x(1,j*h)/h;          // right edge - (1,0) - (1,1)

      // RHS for the top edge (Dirichlet)
      if(j+1 == m+1 && i != 0 && i != m+1)
	sum += u(i*h,1)/(h*h);         // top edge - (0,1) - (1,1)
       
      F(k) = f(i*h, j*h) + sum;
      ++k;
    }
  }
}

VectorXd Poisson::solve()
{
  SparseLU<SpMat, COLAMDOrdering<int> > solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  U = solver.solve(F);

  VectorXd E;
  if(neumann)
  {
    E.resize(M*M);
    for(int i = 0; i < m*m; ++i)
      E(i) = X(i) - U(i);
  }
  else
  {
    E = X - U;
  }
  std::cout << "1-norm = " << E.lpNorm<1>()*h << std::endl;
  std::cout << "2-norm = " << E.norm() << std::endl;
  std::cout << "inf-norm = " << E.lpNorm<Infinity>() << std::endl;
  
  return U;
}

#endif
