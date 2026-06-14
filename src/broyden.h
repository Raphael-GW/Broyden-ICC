#ifndef __BROYDEN_H__
#define __BROYDEN_H__

#include "utils.h"

// void montaJacobiana(const double *x, double **j, int n);

//Jacobiana Otimizada
void montaJacobiana(long double *x, long double *a,  long double *b, long double *c, int n);

void avaliaF(const long double *x, long double *f, int n);

//newton otimizado
void newton(long double *a, long double *b, long double *c,  long double *f, long double *s, long double *x, long double tol, 
int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao, FILE *out);
  
// void newton(double **j, double *f, double *s, double *x, double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao);

long double norma(long double *v, int n);

#endif