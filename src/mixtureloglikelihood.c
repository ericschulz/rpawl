#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#define UW exp(p_param[nbparam * icomponent + iparam])
#define PR exp(p_param[nbparam * (2 * nbcomponents + icomponent) + iparam])
#define SQPR exp(0.5 * p_param[nbparam * (2 * nbcomponents + icomponent) + iparam])


SEXP mixtureloglikelihood(SEXP param, SEXP mixturesample){
  int nprotect = 0;
  int * dim;
  dim = INTEGER(getAttrib(param, R_DimSymbol));
  int nbparam = dim[0];
  SEXP results;
  PROTECT(results = allocVector(REALSXP, nbparam)); nprotect ++;
  int dimparam = dim[1];
  int nbcomponents = (dimparam - 1) / 3;
  //printf("n compo: %d\n", nbcomponents);
  double *p_param, *p_sample;
  PROTECT(param = coerceVector(param, REALSXP)); nprotect ++;
  p_param = REAL(param);
  PROTECT(mixturesample = coerceVector(mixturesample, REALSXP)); nprotect ++;
  p_sample = REAL(mixturesample);
  int samplesize = length(mixturesample);
  //printf("sample size: %d\n", samplesize);
  //printf("first element: %f\n", p_sample[0]);

  float SumOfWeights = 0.;
  float centeredobs = 0.;
  float samplesquared = 0.;
  float mixturedensity = 0.;

  //SEXP uw, prec, sqrt;
  /*
  double *unnormalizedweights, *precisions, *sqrtprecisions;
  unnormalizedweights = allocVector(REALSXP, nbcomponents);
  precisions = allocVector(REALSXP, nbcomponents);
  sqrtprecisions = allocVector(REALSXP, nbcomponents);
  */
  for (int iparam = 0; iparam < nbparam; iparam++){
    //printf("param #%d\n", iparam);
    REAL(results)[iparam] = 0.;
    SumOfWeights = 0.;
    for (int icomponent = 0; icomponent < nbcomponents; icomponent++){
      //unnormalizedweights[icomponent] = exp(p_param[nbparam * icomponent + iparam]);
      //printf("%f\n", unnormalizedweights[icomponent]);
      SumOfWeights += UW;
      //precisions[icomponent] = exp(p_param[nbparam * (2 * nbcomponents + icomponent) + iparam]);
      //sqrtprecisions[icomponent] = exp(0.5 * p_param[nbparam * (2 * nbcomponents + icomponent) + iparam]);
    }
    for (int isample = 0; isample < samplesize; isample ++){
      //printf("sample #%d\n", isample);
      mixturedensity = 0.;
      for (int icomponent = 0; icomponent < nbcomponents; icomponent++){
        //printf("component #%d\n", icomponent);
        centeredobs = (p_sample[isample] - p_param[nbparam * (nbcomponents + icomponent) + iparam]);
        samplesquared = centeredobs * centeredobs;
        mixturedensity += 0.3989423 * (UW / SumOfWeights) * SQPR * exp(-0.5 * PR * samplesquared);
        /*
           mixturedensity += 0.3989423 * unnormalizedweights[icomponent] / SumOfWeights *
           sqrtprecisions[icomponent] * exp(-0.5 * precisions[icomponent] * samplesquared);
           */
      }
      REAL(results)[iparam] += log(mixturedensity);
    }
    //printf("result param %d: %f\n", iparam, REAL(results)[iparam]);
  }
  UNPROTECT(nprotect);
  return results;
}


