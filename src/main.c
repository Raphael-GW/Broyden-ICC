#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <likwid.h>

#include "broyden.h"
#include "utils.h"

int main(int argc, char **argv) {
  int n, max_iter;
  long double x0, tol;

  // Tratamento dos argumentos da linha de comando para o parâmetro -o
  char *output_filename = NULL;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "-o") == 0) {
      if (i + 1 < argc) {
        output_filename = argv[i + 1];
        i++;
      } else {
        fprintf(stderr, "Erro: A opção -o requer o caminho de um arquivo de saída.\n");
        return 1;
      }
    }
  }

  // Cria um arquivo com o filename caso não seja encontrado
  FILE *out = stdout;
  if (output_filename != NULL) {
    out = fopen(output_filename, "w");
    if (out == NULL) {
      perror("Erro ao criar/abrir o arquivo de saída especificado");
      return 1;
    }
  }

  // printf ("Digite a dimensao do sistema (N): ");
  if (scanf ("%d", &n) != 1) {
    fprintf(stderr, "Erro: Entrada inválida para a dimensão (N).\n");
    return 1;
  }
    
  // printf ("Digite o chute inicial (x0): ");
  if (scanf ("%Lf", &x0) != 1) {
      fprintf(stderr, "Erro: Entrada inválida para o chute inicial (x0).\n");
      return 1;
  }
    
  // printf ("Digite a tolerancia (TOL): ");
  if (scanf ("%Lf", &tol) != 1) {
    fprintf(stderr, "Erro: Entrada inválida para a tolerância (TOL).\n");
    return 1;
  }
    
  // printf ("Digite o numero maximo de iteracoes: ");
  if (scanf ("%d", &max_iter) != 1) {
    fprintf(stderr, "Erro: Entrada inválida para o número de iterações.\n");
    return 1;
  }

  long double *x =  malloc(n * sizeof(long double)); // vetor com as variáveis x
  long double *f =  malloc(n * sizeof(long double)); // vetor com cada resultado de f(x)
  long double *s =  malloc(n * sizeof(long double)); // delta usado para atualizar x (i+1)
    
  // Jacobiana otimizada
  long double *a =  malloc(n * sizeof(long double)); // subdiagonal
  long double *b =  malloc(n * sizeof(long double)); // diagonal principal
  long double *c =  malloc(n * sizeof(long double)); // superdiagonal
    

  // Chute inicial para todos
  for (int i = 0; i < n; i++) x[i] = x0;
    
  avaliaF(x, f, n); //monta o vetor F(X)

  LIKWID_MARKER_INIT;

  rtime_t tempo_montaJ = 0.0, tempo_resolucao = 0.0;

  LIKWID_MARKER_START("Newton");
  rtime_t tempo_newton = timestamp();

  // newton(j, f, s, x, tol, max_iter, n, &tempo_montaJ, &tempo_resolucao);

  // Newton otimizado
  newton(a, b, c, f, s, x, tol, max_iter, n, &tempo_montaJ, &tempo_resolucao, out);
  tempo_newton = timestamp() - tempo_newton;
  
  LIKWID_MARKER_STOP("Newton");

    
  fprintf (out, "##########\n");
  fprintf (out, "# Tempo Total: %f\n", tempo_newton);
  fprintf (out, "# Tempo Jacobiana: %f\n", tempo_montaJ);
  fprintf (out, "# Tempo SL: %f\n", tempo_resolucao);
  fprintf (out, "##########\n");

  free(x); 
  free(f); 
  free(s);
  free(a);
  free(b);
  free(c);

  if (out != stdout) {
    fclose(out);
  }

  LIKWID_MARKER_CLOSE;
  return 0;
}