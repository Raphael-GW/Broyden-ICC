#ifndef __BROYDEN_H__
#define __BROYDEN_H__

#include "utils.h"

void montaJacobiana(const double *x, double **j, int n);

void avaliaF(const double *x, double *f, int n);

//newton otimizado
//void newton(double *a, double *b, double *c, double *f, double *s, double *x, double tol, int max_iter, int n, rtime_t tempoJ, rtime_t tempo_resolucao);

void newton(double **j, double *f, double *s, double *x, double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao);

double norma(const double *v, int n);

#endif