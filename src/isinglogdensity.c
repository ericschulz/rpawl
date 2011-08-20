#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

SEXP IsingUpdate(SEXP chains, SEXP dataimg, SEXP imagesize, SEXP flippedindex){
  int imgsize;
  int nprotect = 0;
  int *dim;
  int *p_chains;
  int *p_data;
  int *p_flippedindex;

  SEXP list, list_names;
  char *names[2] = {"equalToData", "localSimilarity"};
  PROTECT(list_names = allocVector(STRSXP, 2)); nprotect ++;
  for(int i = 0; i < 2; i++)
      SET_STRING_ELT(list_names, i,  mkChar(names[i]));
  PROTECT(list = allocVector(VECSXP, 2)); nprotect ++;
  setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names
  dim = INTEGER(getAttrib(chains, R_DimSymbol));
  PROTECT(chains = coerceVector(chains, INTSXP)); nprotect ++;
  PROTECT(dataimg = coerceVector(dataimg, INTSXP)); nprotect ++;
  PROTECT(flippedindex = coerceVector(flippedindex, INTSXP)); nprotect ++;
  p_chains = INTEGER(chains);
  p_data = INTEGER(dataimg);
  p_flippedindex = INTEGER(flippedindex);
  imgsize = INTEGER_VALUE(imagesize);
  int nchains = dim[0];
  SEXP equalToData;
  PROTECT(equalToData = allocVector(REALSXP, nchains)); nprotect ++;
  SEXP localSimilarity;
  PROTECT(localSimilarity = allocVector(REALSXP, nchains)); nprotect ++;
  int col, row;
  int indexchain, flip;
  int minx, maxx, miny, maxy;
  for (indexchain = 0; indexchain < nchains; indexchain ++){
      flip = p_flippedindex[indexchain];
      REAL(equalToData)[indexchain] = (p_chains[flip] 
              == p_data[flip]);
      REAL(localSimilarity)[indexchain] = 0;
      col = (int) (flip / imgsize);
      row = flip - col * imgsize;
      minx = (int) fmax(col - 1, 0);
      maxx = (int) fmin(col + 1, imgsize - 1);
      miny = (int) fmax(row - 1, 0);
      maxy = (int) fmin(row + 1, imgsize -1 );
      for (int dx = minx; dx <= maxx; dx++){
          for (int dy = miny; dy <= maxy; dy++){
              REAL(localSimilarity)[indexchain] += (
                      p_chains[indexchain + nchains * (imgsize * dx + dy)]
                      == p_chains[indexchain + nchains * (imgsize * col + row)]);
          }
      }
  }
  UNPROTECT(nprotect);
  SET_VECTOR_ELT(list, 0, equalToData);
  SET_VECTOR_ELT(list, 1, localSimilarity);
  return list;
}

SEXP IsingCompareToData(SEXP chains, SEXP dataimg, SEXP imagesize){
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

SEXP IsingCountSimilarities(SEXP chains, SEXP imagesize){
  int imgsize;
  int nprotect = 0;
  int *dim;
  int *p_chains;
  dim = INTEGER(getAttrib(chains, R_DimSymbol));
  PROTECT(chains = coerceVector(chains, INTSXP)); nprotect ++;
  p_chains = INTEGER(chains);
  imgsize = INTEGER_VALUE(imagesize);
  int nchains = dim[0];
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

