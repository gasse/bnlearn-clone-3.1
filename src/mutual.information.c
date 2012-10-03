
#include "common.h"

/* Number of degrees of freedom adjustment heuristic, to be used in C code. */
int c_adjusted_df(int **xmargins, int **ymargins, int xdims, int ydims, int zdims) {

int dfx = 0, dfy = 0, df = 0;
int i = 0, k = 0;

  /* for each Z dimension */
  for (k = 0 ; k < zdims ; k++) {

    dfx = (xdims - 1);
    dfy = (ydims - 1);

    /* deduce empty cols (X margins) */
    for (i = 0 ; i < xdims; i++)
      if (xmargins[k][i] == 0)
        dfx--;

    /* deduce empty rows (Y margins) */
    for (i = 0 ; i < ydims; i++)
      if (ymargins[k][i] == 0)
        dfy--;

    df += MAX(dfx, 0) * MAX(dfy, 0);

  }/*FOR*/

  return df;

}/*C_ADJUSTED_DF*/

/* unconditional mutual information, to be used in C code. */
double c_mi(int *xx, int *llx, int *yy, int *lly, int *num, int *df) {

int i = 0, j = 0, k = 0;
int  **n = NULL, *ni = NULL, *nj = NULL;
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc2dcont(*llx, *lly);
  ni = alloc1dcont(*llx);
  nj = alloc1dcont(*lly);

  /* compute the joint frequency of x and y. */
  for (k = 0; k < *num; k++) {

    n[xx[k] - 1][yy[k] - 1]++;

  }/*FOR*/

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) {

    ni[i] += n[i][j];
    nj[j] += n[i][j];

  }/*FOR*/

  /* if requested, fill the number of degrees of freedom */
  if (df != NULL) {

    df[0] = (*llx - 1) * (*lly - 1); // standard computation
    df[1] = c_adjusted_df(&ni, &nj, *llx, *lly, 1); // adjusted computation

  }/*THEN*/

  /* compute the mutual information from the joint and marginal frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++) 
      res += MI_PART(n[i][j], ni[i], nj[j], *num);

  return (res)/(*num);

}/*C_MI*/

/* unconditional mutual information, to be used for the asymptotic test. */
SEXP mi(SEXP x, SEXP y, SEXP lx, SEXP ly, SEXP length, SEXP withdf) {

int *llx = INTEGER(lx), *lly = INTEGER(ly), *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y);
int df[2];
double *res = NULL;
SEXP result;

  if (isTRUE(withdf)) {

    PROTECT(result = allocVector(REALSXP, 3));
    res = REAL(result);

    res[0] = c_mi(xx, llx, yy, lly, num, df); // MI measure
    res[1] = df[0]; // degrees of freedom
    res[2] = df[1]; // adjusted degrees of freedom

  }/*THEN*/
  else {

    PROTECT(result = allocVector(REALSXP, 1));
    res = REAL(result);

    res[0] = c_mi(xx, llx, yy, lly, num, NULL); // MI measure

  }/*ELSE*/

  UNPROTECT(1);
  return result;

}/*MI*/

/* conditional mutual information, to be used in C code. */
double c_cmi(int *xx, int *llx, int *yy, int *lly, int *zz, int *llz, int *num, int *df) {

int i = 0, j = 0, k = 0; 
int ***n = NULL, **ni = NULL, **nj = NULL, *nk = NULL;
double res = 0;

  /* initialize the contingency table and the marginal frequencies. */
  n = alloc3dcont(*llz, *llx, *lly);
  ni = alloc2dcont(*llz, *llx);
  nj = alloc2dcont(*llz, *lly);
  nk = alloc1dcont(*llz);

  /* compute the joint frequency of x, y, and z. */
  for (k = 0; k < *num; k++) {

    n[zz[k] - 1][xx[k] - 1][yy[k] - 1]++;

  }/*FOR*/

  /* compute the marginals. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++) {

        ni[k][i] += n[k][i][j];
        nj[k][j] += n[k][i][j];
        nk[k] += n[k][i][j];

      }/*FOR*/

  /* if requested, fill the number of degrees of freedom */
  if (df != NULL) {

    df[0] = (*llx - 1) * (*lly - 1) * (*llz); // standard computation
    df[1] = c_adjusted_df(ni, nj, *llx, *lly, *llz); // adjusted computation

  }/*THEN*/

  /* compute the conditional mutual information from the joint and
     marginal frequencies. */
  for (i = 0; i < *llx; i++)
    for (j = 0; j < *lly; j++)
      for (k = 0; k < *llz; k++) 
        res += MI_PART(n[k][i][j], ni[k][i], nj[k][j], nk[k]);

  res = res/(*num);

  return res;

}/*C_CMI*/

/* conditional mutual information, to be used for the asymptotic test. */
SEXP cmi(SEXP x, SEXP y, SEXP z, SEXP lx, SEXP ly, SEXP lz, SEXP length, SEXP withdf) {

int *llx = INTEGER(lx), *lly = INTEGER(ly), *llz = INTEGER(lz);
int *num = INTEGER(length);
int *xx = INTEGER(x), *yy = INTEGER(y), *zz = INTEGER(z);
int df[2];
double *res = NULL;
SEXP result;

  if (isTRUE(withdf)) {

    /* allocate and initialize result to zero. */
    PROTECT(result = allocVector(REALSXP, 3));
    res = REAL(result);

    res[0] = c_cmi(xx, llx, yy, lly, zz, llz, num, df); // MI measure
    res[1] = df[0]; // degrees of freedom
    res[2] = df[1]; // adjusted degrees of freedom

  }/*THEN*/
  else {

    /* allocate and initialize result to zero. */
    PROTECT(result = allocVector(REALSXP, 1));
    res = REAL(result);

    res[0] = c_cmi(xx, llx, yy, lly, zz, llz, num, NULL); // MI measure

  }/*ELSE*/

  UNPROTECT(1);
  return result;

}/*CMI*/

/* unconditional Gaussian mutual information, to be used in C code. */
double c_mig(double *xx, double *yy, int *num) {

double cor = c_fast_cor(xx, yy, num);

  return - 0.5 * log(1 - cor * cor);

}/*C_MIG*/

/* unconditional Gaussian mutual information, to be used in the asymptotic test. */
SEXP mig(SEXP x, SEXP y, SEXP length) {

double *xx = REAL(x), *yy = REAL(y);
int *num = INTEGER(length);
SEXP result;

  PROTECT(result = allocVector(REALSXP, 1));
  NUM(result) = c_mig(xx, yy, num);
  UNPROTECT(1);

  return result;

}/*MIG*/
