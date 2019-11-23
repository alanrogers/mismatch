#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "eprintf.h"

/* test the eprintf function */
void main(void)
{
  int i=2;
  double x=3.4;
  
  if(fork() == 0) {
    sleep(5);
    eprintf("testing absence of args", "eprintf w/ no args:");
  }

  if(fork() == 0) {
    eprintf(NULL, "eprintf w/ no args:");
  }

  if(fork() == 0) {
    eprintf("testing ordinary args","eprintf w/ args. i=%d x=%lf:", i, x);
  }

  if(fork() == 0) {
    x /= 0.0;
    eprintf(NULL,"eprintf after division by 0. x=%lf:", x);
  }

  if(fork() == 0) {
    x = log(0.0);
    eprintf("testing log(0.0)","eprintf after log(0.0) x=%lf:", x);
  }

  sleep(1);

  printf("Parent exiting with status 0\n");
  exit(0);
}
