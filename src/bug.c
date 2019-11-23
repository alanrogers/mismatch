#include <stdio.h>
#include <math.h>

/****************************************************************
This program demonstrates a gcc compiler bug, which appears when
either the -O2 or -O3 optimization flags are used.  The bug has been
found on the following architectures:

1. Pentium under Red Hat Linux with gcc version 2.7.2.3

To reproduce it, compile this program as follows:

gcc -O2 -o bug bug.c -lm

Executing bug produces:

  top: x=1.5
  bottom: x=5 y=1.5

Note that the value of x at the bottom of the loop should be identical
to that at the top.  The program should make 4 passes through the loop
but instead quits after 1 pass.

****************************************************************/
void main(void)
{
  double x, y;

  for(x=1.5; x <= 3.25; x += 0.5)
  {
    printf("top: x=%g\n", x);
    y = pow(10.0, x);
    printf("bottom: x=%g y=%g\n\n", x, y);
  }
}



