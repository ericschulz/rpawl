#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP compareToData(SEXP chains, SEXP dataimg, SEXP imagesize){
  int imgsize;
  int nprotect = 0;
  int *dim;
  int *p_chains;
  int *p_data;
  dim = INTEGER(getAttrib(chains, R_DimSymbol));
  PROTECT(chains = coerceVector(chains, INTSXP)); nprotect ++;
  PROTECT(dataimg = coerceVector(dataimg, INTSXP)); nprotect ++;
  p_chains = INTEGER(chains);
  p_data = INTEGER(dataimg);
  imgsize = INTEGER_VALUE(imagesize);
  int nchains = dim[0];
  int dimension = dim[1];
  SEXP results;
  PROTECT(results = allocVector(REALSXP, nchains)); nprotect ++;
  int count_similarities = 0; 
  int row, col, indexchain;
  for (indexchain = 0; indexchain < nchains; indexchain ++){
    count_similarities = 0;
    for (row = 0; row < imgsize; row++){
      for (col = 0; col < imgsize; col ++){
           count_similarities += (p_chains[indexchain + nchains * (imgsize * col + row)] 
               == p_data[imgsize * col + row]);
      }
    }
    REAL(results)[indexchain] = count_similarities;
  }
  UNPROTECT(nprotect);
  return results;
}

SEXP countSimilarities(SEXP chains, SEXP imagesize){
  int imgsize;
  int nprotect = 0;
  int *dim;
  int *p_chains;
  dim = INTEGER(getAttrib(chains, R_DimSymbol));
  PROTECT(chains = coerceVector(chains, INTSXP)); nprotect ++;
  p_chains = INTEGER(chains);
  imgsize = INTEGER_VALUE(imagesize);
  int nchains = dim[0];
  int dimension = dim[1];
  SEXP results;
  PROTECT(results = allocVector(REALSXP, nchains)); nprotect ++;
  int count_similarities = 0; 
  int row, col, indexchain;
  int minx, maxx, miny, maxy;
  for (indexchain = 0; indexchain < nchains; indexchain ++){
    count_similarities = 0;
    for (row = 0; row < imgsize; row++){
      for (col = 0; col < imgsize; col ++){
        minx = (int) fmax(col - 1, 0);
        maxx = (int) fmin(col + 1, imgsize - 1);
        miny = (int) fmax(row - 1, 0);
        maxy = (int) fmin(row + 1, imgsize -1 );
        for (int dx = minx; dx <= maxx; dx++){
          for (int dy = miny; dy <= maxy; dy++){
            count_similarities += (
               p_chains[indexchain + nchains * (imgsize * dx + dy)]
               == p_chains[indexchain + nchains * (imgsize * col + row)]);
          }
        }
      }
    }
    count_similarities -= imgsize * imgsize;
    REAL(results)[indexchain] = count_similarities;
  }
  UNPROTECT(nprotect);
  return results;
}

