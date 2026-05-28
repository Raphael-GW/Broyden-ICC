#ifndef __BROYDEN_H__
#define __BROYDEN_H__

void montaJacobiana(const double *x, double *a, double *b, double *c, int n);

void avaliaF(const double *x, double *f, int n);

void newton(double *a, double *b, double *c, double *f, double *s, double *x, double tol, int max_iter, int n);

double norma(const double *v, int n);

#endif