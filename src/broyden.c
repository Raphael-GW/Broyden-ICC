#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <likwid.h>

#include "utils.h"

/* --------- F(x): calcula o vetor F no ponto x --------- */
void avaliaF(const long double *x, long double *f, int n) {
    f[0] = -2.0L*x[0]*x[0] + 3.0L*x[0] - 2.0L*x[1] + 1.0L;

    int limit = 1 + (n - 2) / 8 * 8; // Limite para o loop unroll, garantindo que não ultrapasse n-1
    for (int i = 1; i < limit; i += 8) {
        f[i] = -2.0L*x[i]*x[i] + 3.0L*x[i] - x[i-1] - 2.0L*x[i+1] + 1.0L;
        f[i+1] = -2.0L*x[i+1]*x[i+1] + 3.0L*x[i+1] - x[i] - 2.0L*x[i+2] + 1.0L;
        f[i+2] = -2.0L*x[i+2]*x[i+2] + 3.0L*x[i+2] - x[i+1] - 2.0L*x[i+3] + 1.0L;
        f[i+3] = -2.0L*x[i+3]*x[i+3] + 3.0L*x[i+3] - x[i+2] - 2.0L*x[i+4] + 1.0L;
        f[i+4] = -2.0L*x[i+4]*x[i+4] + 3.0L*x[i+4] - x[i+3] - 2.0L*x[i+5] + 1.0L;
        f[i+5] = -2.0L*x[i+5]*x[i+5] + 3.0L*x[i+5] - x[i+4] - 2.0L*x[i+6] + 1.0L;
        f[i+6] = -2.0L*x[i+6]*x[i+6] + 3.0L*x[i+6] - x[i+5] - 2.0L*x[i+7] + 1.0L;
        f[i+7] = -2.0L*x[i+7]*x[i+7] + 3.0L*x[i+7] - x[i+6] - 2.0L*x[i+8] + 1.0L;
    }
    for (int i = limit; i < n - 1; i++) {
        f[i] = -2.0L*x[i]*x[i] + 3.0L*x[i] - x[i-1] - 2.0L*x[i+1] + 1.0L;
    }
    
    f[n-1] = -2.0L*x[n-1]*x[n-1] + 3.0L*x[n-1] - x[n-2];
}
 
void montaJacobiana(long double *x, long double *a,  long double *b, long double *c, int n) {
// Reconstrói a diagonal principal baseada no x atualizado
    for (int i = 0; i < n; i++) {
        b[i] = -4.0L * x[i] + 3.0L;
    }
    
    // Reseta completamente as subdiagonais e superdiagonais (limpando o lixo da iteração anterior)
    a[0] = 0.0L;
    c[0] = -2.0L;   // ← BUG FIX: c[0] nunca era setado!
    for (int i = 1; i < n - 1; i++) {
        a[i] = -1.0L;
        c[i] = -2.0L;
    }
    a[n-1] = -1.0L;
    c[n-1] = 0.0L;
    
}


//Eliminação de Gauss otimizada
void eliminacaoGauss (long double *a, long double *b, long double *c, long double *s, long double *f, int n){
    // Triangularização (Algoritmo de Thomas)
    int limit = 1 + ((n - 1) / 8) * 8; // Limite para o loop unroll, garantindo que não ultrapasse n-1
    
    for (int i = 1; i < limit; i += 8) {
        // Proteção essencial contra divisão por zero ou valores perigosamente próximos de zero
        
        
        long double m1 = a[i] / b[i-1];
        b[i] -= m1 * c[i-1];
        f[i] -= m1 * f[i-1];
        a[i] = 0.0;

        long double m2 = a[i+1] / b[i];
        b[i+1] -= m2 * c[i];
        f[i+1] -= m2 * f[i];
        a[i+1] = 0.0;

        long double m3 = a[i+2] / b[i+1];
        b[i+2] -= m3 * c[i+1];
        f[i+2] -= m3 * f[i+1];
        a[i+2] = 0.0;

        long double m4 = a[i+3] / b[i+2];
        b[i+3] -= m4 * c[i+2];
        f[i+3] -= m4 * f[i+2];
        a[i+3] = 0.0;

        long double m5 = a[i+4] / b[i+3];
        b[i+4] -= m5 * c[i+3];
        f[i+4] -= m5 * f[i+3];
        a[i+4] = 0.0;

        long double m6 = a[i+5] / b[i+4];
        b[i+5] -= m6 * c[i+4];
        f[i+5] -= m6 * f[i+4];
        a[i+5] = 0.0;

        long double m7 = a[i+6] / b[i+5];
        b[i+6] -= m7 * c[i+5];
        f[i+6] -= m7 * f[i+5];
        a[i+6] = 0.0;

        long double m8 = a[i+7] / b[i+6];
        b[i+7] -= m8 * c[i+6];
        f[i+7] -= m8 * f[i+6];
        a[i+7] = 0.0;
        
        
    }

    // Proteção para o último elemento da diagonal antes da retro-substituição
    if (fabsl(b[n-1]) < 1e-12) {
        fprintf(stderr, "Erro Numérico: b[%d] é zero na retro-substituição.\n", n-1);
        exit(EXIT_FAILURE);
    }

    // Retro-substituição
    s[n-1] = f[n-1] / b[n-1];
    int i = n - 2; // Declarado fora para manter o estado entre os loops

    // 1. Loop Desenrolado (Sempre que houver pelo menos 8 elementos restantes)
    for (; i >= 7; i -= 8) {
        s[i]   = (f[i]   - c[i]   * s[i+1]) / b[i];
        s[i-1] = (f[i-1] - c[i-1] * s[i])   / b[i-1];
        s[i-2] = (f[i-2] - c[i-2] * s[i-1]) / b[i-2];
        s[i-3] = (f[i-3] - c[i-3] * s[i-2]) / b[i-3];
        s[i-4] = (f[i-4] - c[i-4] * s[i-3]) / b[i-4];
        s[i-5] = (f[i-5] - c[i-5] * s[i-4]) / b[i-5];
        s[i-6] = (f[i-6] - c[i-6] * s[i-5]) / b[i-6];
        s[i-7] = (f[i-7] - c[i-7] * s[i-6]) / b[i-7];
    }

    // 2. Loop de Resto (Continua exatamente de onde o anterior parou até o 0)
    for (; i >= 0; i--) {
        s[i] = (f[i] - c[i] * s[i+1]) / b[i];
    }
}




/* --------- norma euclidiana --------- */
long double norma(long double *v, int n) {
    long double s = 0.0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return sqrtl(s);
}

//Newton otimizado
void newton(long double *a, long double *b, long double *c,  long double *f, long double *s, long double *x, long double tol, int max_iter, int n, rtime_t *tempoJ, rtime_t *tempo_resolucao) {
  long double *menos_f = (long double *) malloc(n * sizeof(long double));
  rtime_t tempo;
  

  
  //LIKWID_MARKER_START("Monta_Jacobiana");
  tempo = timestamp();
  montaJacobiana(x, a, b, c, n); // Alterado aqui
  *tempoJ += timestamp() - tempo;
  //LIKWID_MARKER_STOP("Monta_Jacobiana");

  
  for (int k = 0; k < max_iter; k++) {
    long double nF = norma(f, n);
    if (nF < tol) { break; }

    for (int i = 0; i < n; i++) menos_f[i] = -f[i];

    //LIKWID_MARKER_START("Resolucao_SL");
    tempo = timestamp ();
    eliminacaoGauss(a, b, c, s, menos_f, n); // Alterado aqui
	*tempo_resolucao += timestamp () - tempo;
    //LIKWID_MARKER_STOP("Resolucao_SL");

    for (int i = 0; i < n; i++) x[i] += s[i];

    for (int i = 0; i < n; i++) {
      printf("x%d = %.6Lf \n", i+1, x[i]);
    }
    printf("#\n");

    nF = norma(s, n);
    if (nF < tol) { break; }

    avaliaF(x, f, n);

    tempo = timestamp ();
    //LIKWID_MARKER_START("Monta_Jacobiana");
    montaJacobiana(x, a, b, c, n); // Alterado aqui
    //LIKWID_MARKER_STOP("Monta_Jacobiana");
	*tempoJ += timestamp() - tempo;
  }
  free(menos_f);
}