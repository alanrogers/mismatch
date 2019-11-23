#ifndef MEMSTACK_H
#define MEMSTACK_H

#define USE_MEMSTACK
#define MAX_STACKS 100
typedef struct {
  char *stack;
  int   ndx;
  int   size;
} MEMSTACK ;

MEMSTACK *setmemstack(int size);
char *stackalloc(unsigned size);
void  stackfree(void);
void  deallocstack(void);
int   checkstack(void);
char *stackmustalloc(unsigned size);
char *stackstrdup(char *s);
void  stackstatus(FILE *fp);

#endif /* MEMSTACK_H */
