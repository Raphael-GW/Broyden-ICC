#ifndef __BROYDEN_H__
#define __BROYDEN_H__

#include "utils.h"

void montaJacobiana(const long double *x, long double **j, int n);

void avaliaF(const long double *x, long double *f, int n);

//newton otimizado
//void newton(long double *a, long double *b, long double *c, long double *f, long double *s, long double *x, long double tol, int max_iter, int n, rtime_t tempoJ, rtime_t tempo_resolucao);

void newton(long double **j, long double *f, long double *s, long double *x, long double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao);

long double norma(const long double *v, int n);

#endif