#ifndef __BROYDEN_H__
#define __BROYDEN_H__

#include "utils.h"

//Jacobiana Otimizada
void montaJacobiana(long double *x, long double *a,  long double *b, long double *c, int n);

void avaliaF(const long double *x, long double *f, int n);

//newton otimizado
void newton(long double *a, long double *b, long double *c,  long double *f, long double *s, long double *x, long double tol, 
int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao, FILE *out);

long double norma(long double *v, int n);

#endif