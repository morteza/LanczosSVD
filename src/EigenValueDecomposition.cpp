#include "EigenValueDecomposition.h"
#include <cmath>

using namespace std;

EigenValueDecomposition::EigenValueDecomposition(const Matrix &input)
{
    n = input.cols;
    d = new double[n];
    e = new double[n];
    V = new double*[n];
    for (int i = 0 ; i < n ; i++)
    {
	V[i] = new double[n];
    }

    //! Check if input matrix is symmetric
    bool isSymmetric = true;
    for (int i = 0; (i < n) & isSymmetric ; i++) {
       for (int j = 0; (j < n) & isSymmetric ; j++) {
	  isSymmetric = (input.values[i][j] == input.values[j][i]);
       }
    }

    if (isSymmetric) {
      for (int i = 0; i < n; i++)
      {
	for (int j = 0; j < n; j++)
	{
	  V[i][j] = input.values[i][j];
	}
      }

      //! TriDiagonalize matrix V (which is a copy of input matrix)
      triDiagonalize();

      //! QL Decomposition
      QLDecomposition();
    }
}

void EigenValueDecomposition::QLDecomposition() {

    //! shift e vector by 1, and remove the first entry
    for(int i = 0 ; i < n-1 ; i++)
    {
	e[i] = e[i+1];
    }
    e[n - 1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = 2.0e-52;

    for (int l = 0; l < n; l++) {

      // Find small subdiagonal element

      double tmp = abs(d[l]) + abs(e[l]);
      if (tmp > tst1)
	  tst1 = tmp;
      int m = l;
      while (m < n) {
	if (abs(e[m]) <= eps * tst1) {
	  break;
	}
	m++;
      }

      if (m > l) {
	int iter = 0;
	do {
	  iter++;

	  // Compute implicit shift

	  double g = d[l];
	  double p = (d[l + 1] - g) / (2.0 * e[l]);
	  double r = hypot(p, 1.0);
	  if (p < 0) {
	    r = -r;
	  }
	  d[l] = e[l] / (p + r);
	  d[l + 1] = e[l] * (p + r);
	  double dl1 = d[l + 1];
	  double h = g - d[l];
	  for (int i = l + 2; i < n; i++) {
	    d[i] -= h;
	  }
	  f += h;

	  // Implicit QL transformation.

	  p = d[m];
	  double c = 1.0;
	  double c2 = c;
	  double c3 = c;
	  double el1 = e[l + 1];
	  double s = 0.0;
	  double s2 = 0.0;
	  for (int i = m - 1; i >= l; i--) {
	    c3 = c2;
	    c2 = c;
	    s2 = s;
	    g = c * e[i];
	    h = c * p;
	    r = hypot(p, e[i]);
	    e[i + 1] = s * r;
	    s = e[i] / r;
	    c = p / r;
	    p = c * d[i] - s * g;
	    d[i + 1] = h + s * (c * g + s * d[i]);

	    // Accumulate transformation.

	    for (int k = 0; k < n; k++) {
	      h = V[k][i + 1];
	      V[k][i + 1] = s * V[k][i] + c * h;
	      V[k][i] = c * V[k][i] - s * h;
	    }
	  }
	  p = -s * s2 * c3 * el1 * e[l] / dl1;
	  e[l] = s * p;
	  d[l] = c * p;

	  // Check for convergence.

	} while (abs(e[l]) > eps * tst1);
      }
      d[l] += f;
      e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors

    for (int i = 0; i < n - 1; i++) {
      int k = i;
      double p = d[i];
      for (int j = i + 1; j < n; j++) {
	if (d[j] < p) {
	  k = j;
	  p = d[j];
	}
      }
      if (k != i) {
	d[k] = d[i];
	d[i] = p;
	for (int j = 0; j < n; j++) {
	  p = V[j][i];
	  V[j][i] = V[j][k];
	  V[j][k] = p;
	}
      }
    }
  }


//! Symmetric Householder reduction to tridiagonal form
void EigenValueDecomposition::triDiagonalize() {

    for (int j = 0; j < n; j++) {
       d[j] = V[n-1][j];
    }

    //! Householder reduction to tridiagonal form

    for (int i = n-1; i > 0; i--) {

       //! Scale to avoid under/overflow

       double scale = 0.0;
       double h = 0.0;
       for (int k = 0; k < i; k++) {
	  scale = scale + abs(d[k]);
       }
       if (scale == 0.0) {
	  e[i] = d[i-1];
	  for (int j = 0; j < i; j++) {
	     d[j] = V[i-1][j];
	     V[i][j] = 0.0;
	     V[j][i] = 0.0;
	  }
       } else {

	  // Generate Householder vector.

	  for (int k = 0; k < i; k++) {
	     d[k] /= scale;
	     h += d[k] * d[k];
	  }
	  double f = d[i-1];
	  double g = sqrt(h);
	  if (f > 0) {
	     g = -g;
	  }
	  e[i] = scale * g;
	  h = h - f * g;
	  d[i-1] = f - g;
	  for (int j = 0; j < i; j++) {
	     e[j] = 0.0;
	  }

	  // Apply similarity transformation to remaining columns.

	  for (int j = 0; j < i; j++) {
	     f = d[j];
	     V[j][i] = f;
	     g = e[j] + V[j][j] * f;
	     for (int k = j+1; k <= i-1; k++) {
		g += V[k][j] * d[k];
		e[k] += V[k][j] * f;
	     }
	     e[j] = g;
	  }
	  f = 0.0;
	  for (int j = 0; j < i; j++) {
	     e[j] /= h;
	     f += e[j] * d[j];
	  }
	  double hh = f / (h + h);
	  for (int j = 0; j < i; j++) {
	     e[j] -= hh * d[j];
	  }
	  for (int j = 0; j < i; j++) {
	     f = d[j];
	     g = e[j];
	     for (int k = j; k <= i-1; k++) {
		V[k][j] -= (f * e[k] + g * d[k]);
	     }
	     d[j] = V[i-1][j];
	     V[i][j] = 0.0;
	  }
       }
       d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n-1; i++) {
       V[n-1][i] = V[i][i];
       V[i][i] = 1.0;
       double h = d[i+1];
       if (h != 0.0) {
	  for (int k = 0; k <= i; k++) {
	     d[k] = V[k][i+1] / h;
	  }
	  for (int j = 0; j <= i; j++) {
	     double g = 0.0;
	     for (int k = 0; k <= i; k++) {
		g += V[k][i+1] * V[k][j];
	     }
	     for (int k = 0; k <= i; k++) {
		V[k][j] -= g * d[k];
	     }
	  }
       }
       for (int k = 0; k <= i; k++) {
	  V[k][i+1] = 0.0;
       }
    }
    for (int j = 0; j < n; j++) {
       d[j] = V[n-1][j];
       V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
 }

double EigenValueDecomposition::hypot(double a, double b) {
  double r;
  if (abs(a) > abs(b)) {
	r = b/a;
	r = abs(a)*sqrt(1+r*r);
  } else if (b != 0) {
	r = a/b;
	r = abs(b)*sqrt(1+r*r);
  } else {
	r = 0.0;
  }
  return r;
}

