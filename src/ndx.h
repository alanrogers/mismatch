#ifndef NDX_H
#define NDX_H
/* Master indices (into ndx_e and ndx_s) */
extern int i_theta0;
extern int i_theta1;
extern int i_tau;
extern int i_mse;
extern int i_mae;
extern int i_roughness;
extern int i_seg;

/****************************************************************
ndx_e[i_tau] is the index of tau in the vector of variables
estimated from real data.
****************************************************************/
extern int   dim_e;                /* number of variables to estimate */
extern char *lbl_e[MAXESTIMATES];  /* labels of the variables         */
extern int   ndx_e[MAXESTIMATES];

/* Indices to variables to be estimated from simulated data */
extern int   dim_s;
extern char *lbl_s[MAXESTIMATES];
extern int   ndx_s[MAXESTIMATES];

void setindices(void);

#endif /* NDX_H */
