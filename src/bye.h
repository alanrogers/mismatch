#ifndef BYE_H
#define BYE_H
#include <stdio.h>
FILE *mustopen(char *name, char *mode);
void error(char *s);
char *mustalloc(unsigned bytes);
FloatVec *newFloatVec(int lo, int hi);
DblVec *newDblVec(int lo, int hi);
IntVec *newIntVec(int lo, int hi);
void initFloatVec(FloatVec *v);
void initDblVec(DblVec *v);
void initIntVec(IntVec *v);
#endif /* BYE_H */
