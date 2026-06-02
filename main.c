#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "broyden.h"
#include "utils.h"

int main() {
    int n, max_iter; // define  o tamanho do sistema e o número máximo de iterações
    double x0, tol; // define o chute inicial e a tolerância para convergência
    rtime_t tempo, tempo_montaJ, tempo_newton, tempo_resolucao = 0.0;

    
    printf ("Digite a dimensao do sistema (N): ");
    scanf ("%d", &n);
    printf ("Digite o chute inicial (x0): ");
    scanf ("%lf", &x0);
    printf ("Digite a tolerancia (TOL): ");
    scanf ("%lf", &tol);
    printf ("Digite o numero maximo de iteracoes: ");
    scanf ("%d", &max_iter);

    double *x =  malloc(n * sizeof(double)); // vetor com as variáveis x
    double *f =  malloc(n * sizeof(double)); // vetor com cada resultado de f(x)
    double *s =  malloc(n * sizeof(double)); // delta usado para atualizar x (i+1)
    double *a =  malloc(n * sizeof(double)); // subdiagonal
    double *b =  malloc(n * sizeof(double)); // diagonal principal
    double *c =  malloc(n * sizeof(double)); // superdiagonal

    /* chute inicial para todos */
    for (int i = 0; i < n; i++) x[i] = x0;
    
    avaliaF(x, f, n); //monta o vetor F(X)

    tempo = timestamp();
    montaJacobiana(x, a, b, c, n);
    tempo_montaJ = timestamp() - tempo;

    tempo = timestamp ();
    newton(a, b, c, f, s, x, tol, max_iter, n, tempo_montaJ, tempo_resolucao);
    tempo_newton = timestamp() - tempo;

    
    printf ("##########\n");
    printf ("# Tempo Total: %f\n", tempo_newton);
    printf ("# Tempo Jacobiana: %f\n", tempo_montaJ);
    printf ("# Tampo SL: %f\n", tempo_resolucao);
    printf ("##########\n");

    free(x); 
    free(f); 
    free(s);
    free(a); 
    free(b); 
    free(c);
    return 0;
}