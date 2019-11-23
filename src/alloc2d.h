#ifndef ALLOC2D_H
#define ALLOC2D_H
/** defined in alloc2d.c **/
void  **alloc2d(unsigned dim1, unsigned dim2, unsigned size);
int     free2d(void **mat);
void ***alloc3d(unsigned dim1, unsigned dim2, unsigned dim3, unsigned size);
int     free3d(void ***mat);
void ****alloc4d(unsigned dim1, unsigned dim2, unsigned dim3, unsigned dim4,
		 unsigned size);
int     free4d(void ****mat);
void  **alloclt(unsigned dim, unsigned size);
void  **allocut(unsigned dim, unsigned size);
#endif  /* ALLOC2D_H */
