/*
 * Metodo de Newton aplicado ao Sistema Tridiagonal de Broyden
 * -----------------------------------------------------------
 * Resolve F(x) = 0, onde:
 *   f_1(x) = -2*x_1^2 + 3*x_1 - 2*x_2 + 1
 *   f_i(x) = -2*x_i^2 + 3*x_i - x_{i-1} - 2*x_{i+1} + 1,   2 <= i <= n-1
 *   f_n(x) = -2*x_n^2 + 3*x_n - x_{n-1}
 *
 * A Jacobiana e' tridiagonal:
 *   subdiagonal: a_i = -1                        (i = 2,...,n)
 *   diagonal  : b_i = -4*x_i + 3                 (i = 1,...,n)
 *   superdiag.: c_i = -2                         (i = 1,...,n-1)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <likwid.h>

#include "utils.h"

/* --------- F(x): calcula o vetor F no ponto x --------- */
void avaliaF(const double *x, double *f, int n) {
    f[0] = -2.0*x[0]*x[0] + 3.0*x[0] - 2.0*x[1] + 1.0;

    for (int i = 1; i < n - 1; i++) {
        f[i] = -2.0*x[i]*x[i] + 3.0*x[i] - x[i-1] - 2.0*x[i+1] + 1.0;
    }
    
    f[n-1] = -2.0*x[n-1]*x[n-1] + 3.0*x[n-1] - x[n-2];
}

// --------- Monta as 3 diagonais da Jacobiana J(x) ---------
//  * a: subdiagonal  (indices 1..n-1 usados; a[0] ignorado)
//  * b: diagonal     (indices 0..n-1)
//  * c: superdiag.   (indices 0..n-2 usados; c[n-1] ignorado)
 
void montaJacobiana(double *x, double *a, double *b, double *c, int n) {
    for (int i = 0; i < n; i++) {
        b[i] = -4.0 * x[i] + 3.0;     // diagonal principal varia com x 
    }
    for (int i = 0; i < n-1; i++) a[i] = -1.0;
    for (int i = 0; i < n-1; i++) c[i] = -2.0;
}

// void montaJacobiana(double *x, double **j, int n){
//     for (int i = 0; i < n; ++i){
//         for (int k = 0; k < n; ++k){
//             j[i][k] = 0.0; // preenche a matriz com zeros
//         }
//     }
    
//     j[0][0] = -4.0*x[0] + 3.0;
//     j[0][1] = -2.0;
//     for (int i = 1; i < n-1; ++i){ // monta as 3 diagonais
//         j[i][i-1] = -1.0;
//         j[i][i] = -4.0*x[i] + 3.0;
//         j[i][i+1] = -2.0;
//     }
//     j[n-1][n-2] = -1.0;
//     j[n-1][n-1] = -4.0*x[n-1] + 3.0;
// }


//Eliminação de Gauss otimizada
void eliminacaoGauss (double *a, double *b, double *c, double *s, double *f, int n){
    // triangularizacao
    for (int i = 0; i < n-1; ++i){
        double m = a[i] / b[i];
        a[i] = 0.0;
        b[i+1] -= c[i]*m;
        f[i+1] -= f[i]*m;
    }

    //retro-substituicao
    s[n-1] = f[n-1] / b[n-1];
    for (int i = n-2; i >= 0; --i){
        s[i] = (f[i] - c[i]*s[i+1]) / b[i];
    }
}

// void eliminacaoGauss (double **j, double *s, double *f, int n){
//     for (int i = 0; i < n; ++i){
//         for (int k = i+1; k < n; ++k){
//             double m = j[k][i] / j[i][i];
//             j[k][i] = 0;
//             for (int l = i+1; l < n; ++l)
//                 j[k][l] -= j[i][l]*m;
//             f[k] -= f[i]*m;
//         }
//     }

//     //retro-substituicao
//     s[n-1] = f[n-1] / j[n-1][n-1];
//     for (int i = n-2; i >= 0; --i){
//         s[i] = f[i];
//         for (int l = i+1; l < n; ++l)
//             s[i] -= j[i][l]*s[l];
//         s[i] /= j[i][i];
//     }
// }


/* --------- norma euclidiana --------- */
double norma(double *v, int n) {
    double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrt(s);
}

// void newton(double **j, double *f, double *s, double *x, double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao) {
//   double *menos_f = (double *) malloc(n * sizeof(double)); // -F(X)
//   rtime_t tempo;

//   LIKWID_MARKER_START("Monta_Jacobiana");
//   tempo = timestamp();
//   montaJacobiana(x, j, n);
//   *tempoJ += timestamp() - tempo;
//   LIKWID_MARKER_STOP("Monta_Jacobiana");


//   for (int k = 0; k < max_iter; k++) {
        
//     double nF = norma(f, n);

//     if (nF < tol) {
//       printf("%2d | %.6e |  (convergiu)\n", k, nF);
//       printf("\nConvergencia em %d iteracoes.\n", k);
//       break;
//     }
        

//     for (int i = 0; i < n; i++) menos_f[i] = -f[i];

	
//     LIKWID_MARKER_START("Resolucao_SL");
//     tempo = timestamp ();
//     eliminacaoGauss(j, s, menos_f, n); // resolve J(X) * s = -F(X)
// 	  *tempo_resolucao += timestamp () - tempo;
//     LIKWID_MARKER_STOP("Resolucao_SL");
	

//     for (int i = 0; i < n; i++) x[i] += s[i]; // atualiza X(i+1)

//     for (int i = 0; i < n; i++) {
//       printf ("x%d = %.6f \n", i+1, x[i]);
//     }
//     printf("#\n");

//     nF = norma(s, n);

//     if (nF < tol) {
//       printf("%2d | %.6e |  (convergiu)\n", k, nF);
//       printf("\nConvergencia em %d iteracoes.\n", k);
//       break;
//     }

//     avaliaF(x, f, n); //monta o vetor F(X)

//     tempo = timestamp ();
        
//     LIKWID_MARKER_START("Monta_Jacobiana");
//     montaJacobiana(x, j, n);
//     LIKWID_MARKER_STOP("Monta_Jacobiana");
        
// 	  *tempoJ += timestamp() - tempo;
//   }

//   free(menos_f);
// }

//Newton otimizado
void newton(double *a, double *b, double *c, double *f, double *s, double *x, double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao) {
  double *menos_f = (double *) malloc(n * sizeof(double));
  rtime_t tempo;

  LIKWID_MARKER_START("Monta_Jacobiana");
  tempo = timestamp();
  montaJacobiana(x, a, b, c, n); // Alterado aqui
  *tempoJ += timestamp() - tempo;
  LIKWID_MARKER_STOP("Monta_Jacobiana");

  for (int k = 0; k < max_iter; k++) {
    double nF = norma(f, n);
    if (nF < tol) { break; }

    for (int i = 0; i < n; i++) menos_f[i] = -f[i];

    LIKWID_MARKER_START("Resolucao_SL");
    tempo = timestamp ();
    eliminacaoGauss(a, b, c, s, menos_f, n); // Alterado aqui
	*tempo_resolucao += timestamp () - tempo;
    LIKWID_MARKER_STOP("Resolucao_SL");

    for (int i = 0; i < n; i++) x[i] += s[i];

    for (int i = 0; i < n; i++) {
      printf("x%d = %.6f \n", i+1, x[i]);
    }
    printf("#\n");

    nF = norma(s, n);
    if (nF < tol) { break; }

    avaliaF(x, f, n);

    tempo = timestamp ();
    LIKWID_MARKER_START("Monta_Jacobiana");
    montaJacobiana(x, a, b, c, n); // Alterado aqui
    LIKWID_MARKER_STOP("Monta_Jacobiana");
	*tempoJ += timestamp() - tempo;
  }
  free(menos_f);
}