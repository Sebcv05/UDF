#include <stdio.h>
#include <gsl/gsl_math.h>

int main() {
    double x = gsl_pow_2(3.0);
    printf("Result: %f\n", x);
    return 0;
}
