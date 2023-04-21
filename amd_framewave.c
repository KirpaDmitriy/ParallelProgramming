#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <float.h>
#include "fwBase.h"
#include "fwSignal.h"

#define A 280


void gnome_sort(double* mas, int size) {
  if(size < 2) return;
  int i = 1;
  int j = 2;
  while(i < size) {
    if(mas[i - 1] < mas[i]) {
      i = j;
      j++;
    }
    else {
      double tmp_copy = mas[i];
      mas[i] = mas[i - 1];
      mas[i - 1] = tmp_copy;
      i--;
      if(i == 0) {
        i = j;
        j++;
      }
    }
  }
}

double aggr_func(double x1, double x2) {
  return fmax(x1, x2);
}


double bounded_random(unsigned int* seed, int lower_bound, int upper_bound) {
  int int_rand = rand_r(seed) % (upper_bound - lower_bound) + lower_bound;
  double not_int_rand = 1.0 / ((double) (rand_r(seed) % 100) + 0.00001);
  return (double) int_rand + not_int_rand;
}


int main(int argc, char* argv[]) {
  unsigned int seed_i, N;
  struct timeval T1, T2;
  long delta_ms;
  N = atoi(argv[1]);
  fwSetNumThreads(atoi(argv[2]));
  printf("%s threads are used\n\n", argv[2]);
  gettimeofday(&T1, NULL);
  for(int i=0; i<100; i++) {
    seed_i = i;
    double array1[N];
    double array2[N / 2];
    double array2_copy[N / 2 + 1];
    array2_copy[0] = 0;

    // filling arrays with random data
    for(int mas_i=0; mas_i<N; mas_i++) {
      array1[mas_i] = bounded_random(&seed_i, 1, A);
    }
    for(int mas_i=0; mas_i<N/2; mas_i++) {
      array2[mas_i] = bounded_random(&seed_i, A, 10 * A);
      array2_copy[mas_i + 1] = array2[mas_i];
    }

    // applying functions to each array
    fwsSqrt_64f_A53(array1, array1, N);
    fwsTan_64f_A53(array1, array1, N);
    fwsInv_64f_A53(array1, array1, N);
    for(int mas_i=0; mas_i<N/2; mas_i++) {
      array2[mas_i] = array2[mas_i] + array2_copy[mas_i];
    }
    fwsTanh_64f_A53(array2, array2, N / 2);
    fwsAbs_64f(array2, array2, N / 2);

    // aggregating
    double min_array2_el = DBL_MAX;
    for(int mas_i=0; mas_i<N/2; mas_i++) {
      array2[mas_i] = aggr_func(array1[mas_i], array2[mas_i]);
      if((min_array2_el > array2[mas_i]) && array2[mas_i]) {
        min_array2_el = array2[mas_i];
      }
    }

    // sorting
    gnome_sort(array2, N / 2);

    // verifying
    double X = 0;
    for(int mas_i=0; mas_i<N/2; mas_i++) {
      double current_el = array2[mas_i];
      if(((int) (current_el / min_array2_el)) % 2 == 0) X += sin(current_el);
    }
    printf("%f ", X);
  }
  gettimeofday(&T2, NULL);
  delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
  printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);
  return 0;
}
