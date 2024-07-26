//CubicSolver.c

#include "lagrangian/env.h"
#include <user_header.h>

#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <CubicSolver.h>


// FUNCTION TO SOLVE CUBIC EQUATIONS OF THE FORM ax^3 + bx^2 + cx + d = 0
// Outputs roots type structure containg 3 complex doubles with the 3 roots of the cubic. If one or two roots is complex then the real root is outputted as x0. If all roots are complex then a warning is printed
struct roots Cubic_Solver(CONVERGE_precision_t a, CONVERGE_precision_t b, CONVERGE_precision_t c, CONVERGE_precision_t d)
{
   
  //Quadriatc Sovler for a=0 case
   if(a==0)
   {
   double discriminant, root1, root2, realPart, imaginaryPart;

    discriminant = c * c - 4 * b * d;

    // Case for real and different roots
    if (discriminant > 0) {
        root1 = (-c + sqrt(discriminant)) / (2 * b);
        root2 = (-c - sqrt(discriminant)) / (2 * b);

        // Ensure that root1 has the larger real part
        if (root1 < root2) {
            double temp = root1;
            root1 = root2;
            root2 = temp;
        }
    }
    // Case for real and same roots
    else if (discriminant == 0) {
        root1 = root2 = -c / (2 * b);
    }
    // Case for complex roots
    else {
        realPart = -c / (2 * b);
        root1 = root2 = realPart;      
    }
    struct roots r1;
    r1.x0 = root1;
    r1.x1 = root2;
    return(r1);
   }
   
   CONVERGE_precision_t D0, D1;
   float complex x0, x1, x2;
   D0 = pow(b, 2) - 3 * a * c;
   D1 = (2 * pow(b, 3)) - (9 * a * b * c) + (27 * pow(a, 2) * d);
   float complex z1 = csqrt(pow(D1, 2) - 4 * pow(D0, 3));
   float complex aa = 0.5 * ((creal(D1) + creal(z1)) + I * (cimag(z1) + cimag(D1)));
   float complex cc = 0.33333333333333333333333333333333333333333333333333333333333 + I * 0;
   float complex C = cpowf((aa), cc);
   if (creal(C == 0) && cimag(C) == 0)
   {
      if (D1 == 0)
      {
         printf("\n D1 = 0");
         x0 = -b / (3 * a);
         x1 = x0;
         x2 = x0;
      }
      else
      {
         printf("\n C = 0");
         aa = 0.5 * ((creal(D1) - creal(z1)) + I * (cimag(D1) - cimag(z1)));
         C = cpowf(aa, cc);
      }
   }
  
   double complex eps = (-1 + sqrt(3) * I) / 2;
  
   x0 = -(1 / (3 * a)) * (b + C + (D0 / C));
   x1 = -(1 / (3 * a)) * (b + (eps * C) + (D0 / (eps * C)));
   x2 = -(1 / (3 * a)) * (b + (cpow(eps, 2) * C) + (D0 / (cpow(eps, 2) * C)));
   // printf("\nD0 = %f+i%f D1= %f+i%f      C=%f+i%f        x0 = %f+i%f       x1 = %f+i%f         x3 = %f+i%f",creal(D0),cimag(D0),creal(D1),cimag(D1),creal(C),cimag(C),creal(x0),cimag(x0),creal(x1),cimag(x1),creal(x2),cimag(x2));
   //We want largest real root
  /* if(cimag(x0)>0.001){
                float complex xtemp = x0;
                if(cimag(x1)<0.001){
                        x0 = x1;
                        x1=x2;
                        x2=xtemp;
                }
                else if(cimag(x2)<0.001){
                        x0=x2;
                        x2=x1;
                        x1=xtemp;
                }
                else if(cimag(x1)!=0 && cimag(x2)!=0){
                        printf("\nWarning: Cubic has no real roots\n");
                }
        }*/
   struct roots r1;
   r1.x0 = x0;
   r1.x1 = x1;
   r1.x2 = x2;
   

   return (r1);
}