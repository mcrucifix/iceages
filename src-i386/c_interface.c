#include <R.h>
#include <Rdefines.h>  
#include <stdio.h>

typedef void ddt_det ( double [], int *, int *, int *, double [], 
                          double [], double [], double *, double [], 
                          double [], double [], double [] );

typedef void ddt_sto ( double [], double [], double [], double [], double [], double *, 
                       double [], double [], double [], int *, int * );

/* the following mainly to checking at compile time that functions were properly
 * defined. No strict need to get all of them here  */

extern ddt_det vdp_f_ ;
extern ddt_det i80_f_ ;
extern ddt_det i11_f_ ;
extern ddt_det s90_f_ ;
extern ddt_det s91_f_ ;
extern ddt_det t06_f_ ;
extern ddt_det pp04_f_ ;
extern ddt_det pp12_f_ ;
extern ddt_sto vdp_s_ ;


extern void propagate_d_ (int *, double [], double [], double *,  double *, 
                       int *, double [], double [], 
                       int *, int *,                     /* nap, nao */
                       double [], double [], double [],  /* astronomical forcing */
                       double [], double [], double [],  /* astronomical forcing */
                       int *, int *, ddt_det );


/* extern void propagate_s (n, state, par, t0, t1,deltat, ix, nap, nao, ndim, npar) */

extern void propagate_s_ (int *, double [], double [], double *, double *,  double *, 
                          int *,  /* ix */ 
                          int *, int *,                     /* nap, nao */
                          double [], double [], double [],  /* astronomical forcing */
                          double [], double [], double [],  /* astronomical forcing */
                          int *, int *, ddt_sto );
                      
                          


typedef void readfunc (int*, double [], double[], double[], char *, int* ) ;

 
extern readfunc read_obliquity_ ;
extern readfunc read_precession_ ;

 /* read obliquity and precession from file filename */ 
 SEXP c_readinsol (SEXP Snap, SEXP Snao, SEXP Sfpre, SEXP Sfobl)
{
  int nap, nao,ilpre, ilobl;

  ilpre = strlen(CHAR(STRING_ELT(Sfpre, 0)));
  ilobl = strlen(CHAR(STRING_ELT(Sfobl, 0)));

  char * filepre = calloc(ilpre, sizeof(char));
  char * fileobl = calloc(ilpre, sizeof(char));
    
  strcpy(filepre, CHAR(STRING_ELT(Sfpre, 0)));
  strcpy(fileobl, CHAR(STRING_ELT(Sfobl, 0)));

  nap = INTEGER_VALUE(Snap);
  nao = INTEGER_VALUE(Snao);
 
  double *angpre, *omepre, *amppre;
  double *angobl, *omeobl, *ampobl;
  
  SEXP Sangpre, Samppre, Somepre, Sangobl, Sampobl, 
                Someobl, list, list_names;

  PROTECT(Sangpre = NEW_NUMERIC(nap)); angpre = NUMERIC_POINTER(Sangpre);
  PROTECT(Samppre = NEW_NUMERIC(nap)); amppre = NUMERIC_POINTER(Samppre);
  PROTECT(Somepre = NEW_NUMERIC(nap)); omepre = NUMERIC_POINTER(Somepre);
  PROTECT(Sangobl = NEW_NUMERIC(nao)); angobl = NUMERIC_POINTER(Sangobl);
  PROTECT(Sampobl = NEW_NUMERIC(nao)); ampobl = NUMERIC_POINTER(Sampobl);
  PROTECT(Someobl = NEW_NUMERIC(nao)); omeobl = NUMERIC_POINTER(Someobl);

  read_precession_ (&nap, amppre, omepre, angpre, filepre, &ilpre);
  read_obliquity_ (&nao, ampobl, omeobl, angobl, fileobl, &ilobl);

  PROTECT(list = allocVector(VECSXP, 8));
  PROTECT(list_names = allocVector(STRSXP, 8));

  SET_VECTOR_ELT(list, 0, Snap);       SET_STRING_ELT(list_names,0,mkChar("nap"));
  SET_VECTOR_ELT(list, 1, Snao);       SET_STRING_ELT(list_names,1,mkChar("nao"));
  SET_VECTOR_ELT(list, 2, Samppre);    SET_STRING_ELT(list_names,2,mkChar("amppre"));
  SET_VECTOR_ELT(list, 3, Somepre);    SET_STRING_ELT(list_names,3,mkChar("omepre"));
  SET_VECTOR_ELT(list, 4, Sangpre);    SET_STRING_ELT(list_names,4,mkChar("angpre"));
  SET_VECTOR_ELT(list, 5, Sampobl);    SET_STRING_ELT(list_names,5,mkChar("ampobl"));
  SET_VECTOR_ELT(list, 6, Someobl);    SET_STRING_ELT(list_names,6,mkChar("omeobl"));
  SET_VECTOR_ELT(list, 7, Sangobl);    SET_STRING_ELT(list_names,7,mkChar("angobl"));

  setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names

  UNPROTECT(8);
   
  return(list);

};  

/* propagate_deterministic */ 

 SEXP c_propagate_d  (SEXP Sfunc,  SEXP SNs, 
                      SEXP Sstate , SEXP Spar, 
                      SEXP St0, SEXP St1, 
                      SEXP Sicalclyap, SEXP Sds, SEXP Slyap , 
                      SEXP SAstroList)
{
  ddt_det * func = (ddt_det * ) R_ExternalPtrAddr(Sfunc); 
  
  int *n, *nap, *nao, *npar;
  double *amppre, *omepre, *angpre ;
  double *ampobl, *omeobl, *angobl ;
  int  *ndim, *icalclyap;
  double *t0, *t1, *state, *par, *ds, *lyap;

  n    = INTEGER_POINTER (VECTOR_ELT(SNs, 0));
  npar = INTEGER_POINTER (VECTOR_ELT(SNs, 1));
  ndim = INTEGER_POINTER (VECTOR_ELT(SNs, 2));

  state = NUMERIC_POINTER (Sstate);
  par =  NUMERIC_POINTER (Spar);
  t0 =  NUMERIC_POINTER (St0);
  t1 =  NUMERIC_POINTER (St1);
  printf ("to : %f ; t1 : %f", *t0, *t1);
  icalclyap = INTEGER_POINTER (Sicalclyap);
  ds =  NUMERIC_POINTER (Sds);
  lyap =  NUMERIC_POINTER (Slyap);

  nap = INTEGER_POINTER (AS_INTEGER(VECTOR_ELT(SAstroList, 0)));
  nao = INTEGER_POINTER (AS_INTEGER(VECTOR_ELT(SAstroList, 1)));

  amppre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 2));
  omepre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 3));
  angpre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 4));
  ampobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 5));
  omeobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 6));
  angobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 7));

  propagate_d_ (n , state, par , t0, t1, icalclyap, ds, lyap,
            nap, nao, amppre, omepre, angpre, ampobl, omeobl, angobl,
            ndim, npar, func);    

  return(R_NilValue);
} ;

/* stochastic propagator */

 SEXP c_propagate_s  (SEXP Sfunc, SEXP SNs, 
                      SEXP Sstate , SEXP Spar, 
                      SEXP St0, SEXP St1, SEXP Sdeltat, SEXP Six, 
                      SEXP SAstroList)
{

  ddt_sto * func = (ddt_sto * ) R_ExternalPtrAddr(Sfunc); 
  int *n, *nap, *nao, *npar, *ndim, *ix;
  double *amppre, *omepre, *angpre ;
  double *ampobl, *omeobl, *angobl ;
  double *t0, *t1, *deltat, *state, *par;

  n    = INTEGER_POINTER (AS_INTEGER ( VECTOR_ELT(SNs, 0)));
  npar = INTEGER_POINTER (AS_INTEGER (VECTOR_ELT(SNs, 1)));
  ndim = INTEGER_POINTER (AS_INTEGER (VECTOR_ELT(SNs, 2)));

  state = NUMERIC_POINTER (Sstate);
  par =  NUMERIC_POINTER (Spar);

  t0 =  NUMERIC_POINTER (St0);
  t1 =  NUMERIC_POINTER (St1);
  deltat =  NUMERIC_POINTER (Sdeltat);
  ix = INTEGER_POINTER (Six);
 
  nap = INTEGER_POINTER (AS_INTEGER(VECTOR_ELT(SAstroList, 0)));
  nao = INTEGER_POINTER (AS_INTEGER(VECTOR_ELT(SAstroList, 1)));

  amppre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 2));
  omepre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 3));
  angpre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 4));
  ampobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 5));
  omeobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 6));
  angobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 7));

  propagate_s_ (n , state, par , t0, t1, deltat, ix,
            nap, nao, amppre, omepre, angpre, ampobl, omeobl, angobl,
            ndim, npar, func);    

  return(R_NilValue);
} ;


