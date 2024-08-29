//CubicSolver.h

#ifndef CUBIC_SOLVER_H
#define CUBIC_SOLVER_H
struct roots
{
   float complex x0;
   float complex x1;
   float complex x2;
};
struct roots Cubic_Solver(CONVERGE_precision_t a, CONVERGE_precision_t b, CONVERGE_precision_t c, CONVERGE_precision_t d);

#endif // CUBIC_SOLVER_H