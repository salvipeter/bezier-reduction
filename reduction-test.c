#include <stdio.h>
#include <stdlib.h>

#include "reduction.h"

int main(int argc, char **argv) {
  int n = 7, m = 5, r = 1, s = 1;
  double *Q = (double *)malloc((m + 1) * (n + 1) * sizeof(double));
  reduction_matrix(n, m, r, s, Q);
  for (int i = 0; i <= m; ++i) {
    for (int j = 0; j <= n; ++j) {
      int index = i * (n + 1) + j;
      printf("%.4f\t", Q[index]);
    }
    printf("\n");
  }
  free(Q);
  return 0;
}
