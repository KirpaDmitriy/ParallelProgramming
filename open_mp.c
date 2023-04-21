#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <float.h>
#include <omp.h>

#define A 280


void gnome_sort(double* mas, int size) {
	if(size < 2) return;
	int i = 1;
	int j = 2;
	// can not run in parallel as num iterations is undefined
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


double func1(double x) {
	return 1 / tan(sqrt(x));
}


double func2(double x) {
	return abs(1 / tanh(x));
}


double aggr_func(double x1, double x2) {
	return fmax(x1, x2);
}


double bounded_random(unsigned int* seed, int lower_bound, int upper_bound) {
	// can not be ran in parallel as the seed must be passed sequentially
	unsigned int r1 = rand_r(seed);
	unsigned int r2 = rand_r(seed);

	int int_rand = r1 % (upper_bound - lower_bound) + lower_bound;
	double not_int_rand = 1.0 / ((double) (r2 % 100));
	return (double) int_rand + not_int_rand;
}


int main(int argc, char* argv[]) {
	omp_set_num_threads(4);
	unsigned int seed_i, N;
	struct timeval T1, T2;
	long delta_ms;
	N = atoi(argv[1]);
	gettimeofday(&T1, NULL);
	for(int i=0; i<100; i++) {
		seed_i = i;
		double array1[N];
		double array2[N / 2];
		double array2_copy[N / 2 + 1];
		array2_copy[0] = 0;

		// filling arrays with random data
		// if we run these cycles in parallel, the set of the elements will remain, but their order can break
		for(int mas_i=0; mas_i<N; mas_i++) {
			array1[mas_i] = bounded_random(&seed_i, 1, A);
		}
		// can not run previous and next cycles in parallel as seed will break
		for(int mas_i=0; mas_i<N/2; mas_i++) {
			array2[mas_i] = bounded_random(&seed_i, A, 10 * A);
			array2_copy[mas_i + 1] = array2[mas_i];
		}

		// applying functions to each array
		#pragma omp parallel for default(none) shared(array1, N)
		for(int mas_i=0; mas_i<N; mas_i++) {
			array1[mas_i] = func1(array1[mas_i]);
		}
		#pragma omp parallel for default(none) shared(array2, array2_copy, N)
		for(int mas_i=0; mas_i<N/2; mas_i++) {
			array2[mas_i] = func2(array2[mas_i] + array2_copy[mas_i]);
		}

		// aggregating
		double min_array2_el = DBL_MAX;
		#pragma omp parallel for default(none) shared(array1, array2, min_array2_el, N)
		for(int mas_i=0; mas_i<N/2; mas_i++) {
			array2[mas_i] = aggr_func(array1[mas_i], array2[mas_i]);
			#pragma omp critical
			if((min_array2_el > array2[mas_i]) && array2[mas_i]) {
				min_array2_el = array2[mas_i];
			}
		}

		// sorting
		gnome_sort(array2, N / 2);

		// verifying
		double X = 0;
		#pragma omp parallel for default(none) shared(N, X, array2, min_array2_el)
		for(int mas_i=0; mas_i<N/2; mas_i++) {
			double current_el = array2[mas_i];
			if(((int) (current_el / min_array2_el)) % 2 == 0) {
				#pragma omp atomic
				X += sin(current_el);
			}
		}
		printf("%f ", X);
	}
	gettimeofday(&T2, NULL);
	delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
	printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

	return 0;
}
