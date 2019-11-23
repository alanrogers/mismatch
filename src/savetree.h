#ifndef SAVETREE_H
#define SAVETREE_H

/* prototypes of functions defined in savetree.c */

int writetree(NODE *node, FILE *ofp);
NODE *readtree(FILE *ifp);
int treediff(NODE *n1, NODE *n2);
void freetree(NODE *node);

#endif /* SAVETREE_H */
