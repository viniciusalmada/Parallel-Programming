/* calculatePI.cpp -- Vinicius Almada
 * The purpose of this program is to calculate 
 * PI number with multithread programming.
 * PI can be obtained from integral of 
 * 4 / (1 + x²) from 0 to 1.
 *
 * To obtain this value approximattely we can 
 * use sum of rectangles:
 * PI = sum(F(x_i) * dx) with i equals 0 to N.
 *
 * As larger N the PI will be better approximated.
 *
 * First approach with 5 bilions in single thread
 * 	coasting ~20 seconds
 */
#include <cstdio>
#include <cmath>

/* Calculate function for a given x */
double integralFunction(double x);

/* Calculate PI for a given number of steps */
double calculatePI(double down, double up, long numSteps);

int main() {
	double pi = calculatePI(0.0, 1.0, 5'000'000'000L);

	printf("pi = %.20f\n", pi);
	printf("PI = %.20f\n", M_PI);
}

double integralFunction(double x) {
	return 4.0 / (1.0 + x * x);
}

double calculatePI(double down, double up, long numSteps) {
	double dx = (up - down) / (double) numSteps;

	double sum = 0.0;

	for (long i = 0; i < numSteps; i++) {
		double x = (i + 0.5) * dx;
		double funRes = integralFunction(x);
		sum += funRes * dx;
	}

	return sum;
}
