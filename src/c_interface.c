/* Copyright (c) 2012 Michel Crucifix <michel.crucifix@uclouvain.be>

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject

 the following conditions:

 The above copyright notice and this permission notice shall be
 incluudedin all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND INFRINGEMENT
 IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR

 CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 ------------------------------------------------------------------
 C code  for R
 ------------------------------------------------------------------ */

#include <R.h>
#include <Rdefines.h>  
#include <stdio.h>

typedef void ddt_det ( double [], int *, int *, int *, double [], 
                          double [], double [], double *, double [], 
                          double [], double [], double [] );

/* typedef void ddt_sto ( double [], double [], double [], double [], double [], double *, 
                       double [], double [], double [], int *, int * );  */

 typedef void ddt_sto () ; 

extern void propagate_d_ (int *, double [], double [], double *,  double *, 
                       int *, double [], double [], 
                       int *, int *,                     /* nap, nao */
                       double [], double [], double [],  /* astronomical forcing */
                       double [], double [], double [],  /* astronomical forcing */
                       int *, int *, ddt_det, double *);


/* extern void propagate_s (n, state, par, t0, t1,deltat, ix, nap, nao, ndim, npar) */

extern void propagate_s_ (int *, double [], double [], double[], 
                          double *, double *,  double *, 
                          int *,  /* ix */ 
                          int *, int *,                     /* nap, nao */
                          double [], double [], double [],  /* astronomical forcing */
                          double [], double [], double [],  /* astronomical forcing */
                          int *, int *,  ddt_sto, int* );
                      
                          


typedef void readfunc (int*, double [], double[], double[], char *, int* ) ;

 
extern readfunc read_obliquity_ ;
extern readfunc read_precession_ ;

 /* read obliquity and precession from file filename */ 
 SEXP c_readinsol (SEXP Snap, SEXP Snao, SEXP Sfpre, SEXP Sfobl)
{
  SEXP SnapI, SnaoI;
  int nap, nao,ilpre, ilobl;

  ilpre = strlen(CHAR(STRING_ELT(Sfpre, 0)));
  ilobl = strlen(CHAR(STRING_ELT(Sfobl, 0)));

  char * filepre = calloc(ilpre, sizeof(char));
  char * fileobl = calloc(ilpre, sizeof(char));
    
  strcpy(filepre, CHAR(STRING_ELT(Sfpre, 0)));
  strcpy(fileobl, CHAR(STRING_ELT(Sfobl, 0)));

  /* make sure than nap and nao are returned as 'integers' */ 

  SnapI = AS_INTEGER(Snap); 
  SnaoI = AS_INTEGER(Snao);

  nap = INTEGER_VALUE(SnapI);
  nao = INTEGER_VALUE(SnaoI);
 
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

  SET_VECTOR_ELT(list, 0, SnapI);       SET_STRING_ELT(list_names,0,mkChar("nap"));
  SET_VECTOR_ELT(list, 1, SnaoI);       SET_STRING_ELT(list_names,1,mkChar("nao"));
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
                      SEXP SAstroList, SEXP Sdeltat)
{
  
  ddt_det * func = (ddt_det * ) R_ExternalPtrAddr(Sfunc); 
  
  int *n, *nap, *nao, *npar;
  double *amppre, *omepre, *angpre ;
  double *ampobl, *omeobl, *angobl ;
  int  *ndim, *icalclyap;
  double *t0, *t1, *state, *par, *ds, *lyap, *deltat;

  SEXP list, list_names, OSstate, OSlyap, OSds;

  n    = INTEGER_POINTER (VECTOR_ELT(SNs, 0));
  npar = INTEGER_POINTER (VECTOR_ELT(SNs, 1));
  ndim = INTEGER_POINTER (VECTOR_ELT(SNs, 2));

  PROTECT (OSstate = duplicate(Sstate));
  PROTECT (OSlyap = duplicate(Slyap));
  PROTECT (OSds = duplicate(Sds));

  state = NUMERIC_POINTER (OSstate);
  par =  NUMERIC_POINTER (Spar);
  deltat =  NUMERIC_POINTER (Sdeltat);
  t0 =  NUMERIC_POINTER (St0);
  t1 =  NUMERIC_POINTER (St1);
  icalclyap = INTEGER_POINTER (Sicalclyap);
  ds =  NUMERIC_POINTER (OSds);
  lyap =  NUMERIC_POINTER (OSlyap);

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
            ndim, npar, func, deltat);    

  PROTECT(list = allocVector(VECSXP, 3));
  PROTECT(list_names = allocVector(STRSXP, 3));

  SET_VECTOR_ELT(list, 0, OSstate);    SET_STRING_ELT(list_names,0,mkChar("state"));
  SET_VECTOR_ELT(list, 1, OSds);    SET_STRING_ELT(list_names,1,mkChar("ds"));
  SET_VECTOR_ELT(list, 2, OSlyap);    SET_STRING_ELT(list_names,2,mkChar("lyap"));
 
  setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names
  UNPROTECT(5);
  return(list);
} ;

/* stochastic propagator */

 SEXP c_propagate_s  (SEXP Sfunc, SEXP SNs, 
                      SEXP Sstate , SEXP Spar, SEXP Sscaletime,
                      SEXP St0, SEXP St1, SEXP Sdeltat, SEXP Six, SEXP Sisum,
                      SEXP SAstroList)
{

  ddt_sto * func = (ddt_sto * ) R_ExternalPtrAddr(Sfunc); 
  int *n, *nap, *nao, *npar, *ndim, *ix, *isum;
  double *amppre, *omepre, *angpre ;
  double *ampobl, *omeobl, *angobl ;
  double *t0, *t1, *deltat, *state, *par, *scaletime;

  SEXP list, list_names,  OSstate;

  n    = INTEGER_POINTER (AS_INTEGER ( VECTOR_ELT(SNs, 0)));
  npar = INTEGER_POINTER (AS_INTEGER (VECTOR_ELT(SNs, 1)));
  ndim = INTEGER_POINTER (AS_INTEGER (VECTOR_ELT(SNs, 2)));

  PROTECT (OSstate = duplicate(Sstate));

  state = NUMERIC_POINTER (OSstate);
  par =  NUMERIC_POINTER (Spar);
  scaletime =  NUMERIC_POINTER (Sscaletime);

  t0 =  NUMERIC_POINTER (St0);
  t1 =  NUMERIC_POINTER (St1);
  deltat =  NUMERIC_POINTER (Sdeltat);
  ix   =  INTEGER_POINTER (AS_INTEGER(Six));
  isum =  INTEGER_POINTER (AS_INTEGER(Sisum));
 
  nap = INTEGER_POINTER (AS_INTEGER(VECTOR_ELT(SAstroList, 0)));
  nao = INTEGER_POINTER (AS_INTEGER(VECTOR_ELT(SAstroList, 1)));

  amppre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 2));
  omepre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 3));
  angpre = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 4));
  ampobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 5));
  omeobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 6));
  angobl = NUMERIC_POINTER (VECTOR_ELT(SAstroList, 7));
  propagate_s_ (n , state, par , scaletime, t0, t1, deltat, ix,
            nap, nao, amppre, omepre, angpre, ampobl, omeobl, angobl,
            ndim, npar, func, isum);    

  PROTECT(list = allocVector(VECSXP, 1));
  PROTECT(list_names = allocVector(STRSXP, 1));

  SET_VECTOR_ELT(list, 0, OSstate);    SET_STRING_ELT(list_names,0,mkChar("state"));
 
  setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names
  UNPROTECT (3) ;
  return(list);
} ;


