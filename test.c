//test.c
#include "omp.h"
#include <stdio.h>

int main(void) {
  #pragma omp parallel
  printf("thread %d\n", omp_get_thread_num());
}
