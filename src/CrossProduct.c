//CrossProduct.c

#include "lagrangian/env.h"
#include <user_header.h>
#include <spray_break.h>
#include <CONVERGE/udf.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <complex.h>
#include <CrossProduct.h>

 struct vect cross_product(struct vect a, struct vect b)
   {
      struct vect d;

      d.x = a.y * b.z - a.z * b.y;
      d.y = a.z * b.x - a.x * b.z;
      d.z = a.x * b.y - a.y * b.x;
      return (d);
   }