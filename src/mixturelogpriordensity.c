#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP mixturelogpriordensity(SEXP param, SEXP hyperparam){
  int * dim;
  SEXP results;
  dim = INTEGER(getAttrib(param, R_DimSymbol));
  int nbparam = dim[0];
  PROTECT(results = allocVector(REALSXP, nbparam));
  int dimparam = dim[1];
  int nbcomponents = (dimparam - 1) / 3;
  double *p_hyperparam;
  double *p_param;
  p_hyperparam = REAL(hyperparam);
  p_param= REAL(param);
  float delta = p_hyperparam[0];
  float gammadelta = p_hyperparam[1];
  float xi = p_hyperparam[2];
  float tau = p_hyperparam[3];
  float alpha = p_hyperparam[4];
  float gammaalpha = p_hyperparam[5];
  float g = p_hyperparam[6];
  float gammag = p_hyperparam[7];
  float h = p_hyperparam[8];
  float meancst = - 0.9189385 - 0.5 * log(tau);
  meancst *= nbcomponents;
  float betacst = g * log(h) - log(gammag);
  float weightcst = - log(gammadelta);
  weightcst *= nbcomponents;
  float precisioncst = - log(gammaalpha);
  precisioncst *= nbcomponents;
  float fromMeans;
  float fromBeta;
  float fromWeights;
  float fromPrecs;
  int betaindex, windex, mindex, pindex;
  for(int iparam = 0; iparam < nbparam; iparam ++){
    fromMeans = meancst;
    fromWeights = weightcst;
    fromPrecs = precisioncst;
    betaindex = nbparam * (3 * nbcomponents) + iparam;
    fromBeta = betacst + p_param[betaindex] + (g - 1) * p_param[betaindex] - h * exp(p_param[betaindex]);
    for (int icomponent = 0; icomponent < nbcomponents; icomponent++){
      windex = nbparam * icomponent + iparam;
      mindex = windex + nbparam * nbcomponents;
      pindex = mindex + nbparam * nbcomponents;
      fromWeights += p_param[windex] + (delta - 1) * p_param[windex] - 1 * exp(p_param[windex]);
      fromMeans += - 0.5 / tau * (xi - p_param[mindex]) * (xi - p_param[mindex]);
      fromPrecs += p_param[pindex] + alpha * p_param[betaindex] + (alpha - 1) * p_param[pindex] - exp(p_param[betaindex]) * exp(p_param[pindex]);
    }
    REAL(results)[iparam] = fromBeta + fromMeans + fromWeights + fromPrecs;
  }
  UNPROTECT(1);
  //printf("finished with this prior shit\n");
  return results;
}

