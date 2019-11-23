#include <stdio.h>
#include "memstack.h"
#include "mismatch.h"
/****************************************************************
Below, I give define extern to mean "nothing", read in ndx.h, and
then restore the meaning of extern with the #undef command.  This
causes the variables declared in ndx.h to be allocated here.  They
are extern definitions everywhere else.
****************************************************************/
#define extern  /* now extern means "nothing" */
#include "ndx.h"
#undef extern   /* now extern resumes its usual meaning */

void setindices(void)
{
  int i;

  /** initialize index vectors **/
  for(i=0; i<MAXESTIMATES; i++)
    ndx_e[i] = ndx_s[i] = -1;  /* turned off */

  /** master indices **/
  i_theta0    = 0;
  i_theta1    = 1;
  i_tau       = 2;
  i_mse       = 3;
  i_mae       = 4;
  i_roughness = 5;
  i_seg       = 6;

/***** Indices of vector estimated from data ******/
  dim_e = 0;
#if E_THETA0
  ndx_e[i_theta0] = dim_e++;
  lbl_e[ndx_e[i_theta0]] = "theta0";
#endif
#if E_THETA1
  ndx_e[i_theta1] = dim_e++;
  lbl_e[ndx_e[i_theta1]] = "theta1";
#endif
#if E_TAU
  ndx_e[i_tau] = dim_e++;
  lbl_e[ndx_e[i_tau]] = "tau";
#endif
#if E_MSE
  ndx_e[i_mse] = dim_e++;
  lbl_e[ndx_e[i_mse]] = "MSE";
#endif
#if E_MAE
  ndx_e[i_mae] = dim_e++;
  lbl_e[ndx_e[i_mae]] = "MAE";
#endif
#if E_ROUGHNESS
  ndx_e[i_roughness] = dim_e++;
  lbl_e[ndx_e[i_roughness]] = "Rghns";
#endif  
#if E_SEG
  ndx_e[i_seg] = dim_e++;
  lbl_e[ndx_e[i_seg]] = "Seg";
#endif  

/***** Indices of vector estimated from simulated data ******/
  dim_s = 0;
/** Model 2 indices **/
#if S_THETA0 && E_THETA0
  ndx_s[i_theta0] = dim_s++;
  lbl_s[ndx_s[i_theta0]] = "theta0";
#endif
#if S_THETA1 && E_THETA1
  ndx_s[i_theta1] = dim_s++;
  lbl_s[ndx_s[i_theta1]] = "theta1";
#endif
#if S_TAU && E_TAU
  ndx_s[i_tau] = dim_s++;
  lbl_s[ndx_s[i_tau]] = "tau";
#endif
#if S_MSE && E_MSE
  ndx_s[i_mse] = dim_s++;
  lbl_s[ndx_s[i_mse]] = "MSE";
#endif
#if S_MAE && E_MAE
  ndx_s[i_mae] = dim_s++;
  lbl_s[ndx_s[i_mae]] = "MAE";
#endif
#if S_ROUGHNESS && E_ROUGHNESS
  ndx_s[i_roughness] = dim_s++;
  lbl_s[ndx_s[i_roughness]] = "Rghns";
#endif
#if S_SEG && E_SEG
  ndx_s[i_seg] = dim_s++;
  lbl_s[ndx_s[i_seg]] = "Seg";
#endif  

#if 0
  printf("\n%d valid E-indices:", dim_e);
  for(i=0; i<MAXESTIMATES; i++)
  {
    printf("\n%d: ndx=%d", i, ndx_e[i]);
    if(VALIDNDX(ndx_e[i]))
      printf(" lbl=%s", lbl_e[ndx_e[i]]);
  }

  printf("\n%d valid S-indices:", dim_s);
  for(i=0; i<MAXESTIMATES; i++)
  {
    printf("\n%d: ndx=%d", i, ndx_s[i]);
    if(VALIDNDX(ndx_s[i]))
      printf(" lbl=%s", lbl_s[ndx_s[i]]);
  }
#endif
}


