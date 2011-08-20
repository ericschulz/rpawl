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
  char *names[2] = {"equalToData", "SimilarityChange"};
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
  SEXP localSimilarityChange;
  PROTECT(localSimilarityChange = allocVector(REALSXP, nchains)); nprotect ++;
  int col, row;
  int indexchain, flip;
  int minx, maxx, miny, maxy;
  int nneighbors;
  int nsimilar;
  for (indexchain = 0; indexchain < nchains; indexchain ++){
      flip = p_flippedindex[indexchain] - 1;
//      printf("states:%d\n", p_chains[p_chains[indexchain + nchains * flip]]);
//      printf("realdata:%d\n", p_data[flip]);
      REAL(equalToData)[indexchain] = (p_chains[indexchain + nchains * flip] 
              == p_data[flip]);
      nsimilar = -1;
      col = (int) (flip / imgsize);
      row = flip - col * imgsize;
      minx = (int) fmax(col - 1, 0);
      maxx = (int) fmin(col + 1, imgsize - 1);
      miny = (int) fmax(row - 1, 0);
      maxy = (int) fmin(row + 1, imgsize -1 );
      nneighbors = (1 + maxx - minx) * (1 + maxy - miny) - 1;
      for (int dx = minx; dx <= maxx; dx++){
          for (int dy = miny; dy <= maxy; dy++){
              nsimilar += (
                      p_chains[indexchain + nchains * (imgsize * dx + dy)]
                      == p_chains[indexchain + nchains * (imgsize * col + row)]);
          }
      }
      //printf("nsimilar: %d, nneighbors: %d\n", nsimilar, nneighbors);
      REAL(localSimilarityChange)[indexchain] = 2 * (nneighbors - 2 * nsimilar);
  }
  SET_VECTOR_ELT(list, 0, equalToData);
  SET_VECTOR_ELT(list, 1, localSimilarityChange);
  UNPROTECT(nprotect);
  return list;
}

SEXP IsingCounter(SEXP chains, SEXP dataimg, SEXP imagesize){
  int imgsize;
  int nprotect = 0;
  int *dim;
  int *p_chains;
  int *p_data;

  SEXP list, list_names;
  char *names[2] = {"countData", "countSimilarities"};
  PROTECT(list_names = allocVector(STRSXP, 2)); nprotect ++;
  for(int i = 0; i < 2; i++)
      SET_STRING_ELT(list_names, i,  mkChar(names[i]));
  PROTECT(list = allocVector(VECSXP, 2)); nprotect ++;

  setAttrib(list, R_NamesSymbol, list_names); //and attaching the vector names
  dim = INTEGER(getAttrib(chains, R_DimSymbol));
  PROTECT(chains = coerceVector(chains, INTSXP)); nprotect ++;
  PROTECT(dataimg = coerceVector(dataimg, INTSXP)); nprotect ++;
  p_chains = INTEGER(chains);
  p_data = INTEGER(dataimg);
  imgsize = INTEGER_VALUE(imagesize);
  int nchains = dim[0];
  SEXP resultsData;
  SEXP resultsSimilar;
  PROTECT(resultsData = allocVector(REALSXP, nchains)); nprotect ++;
  PROTECT(resultsSimilar = allocVector(REALSXP, nchains)); nprotect ++;
  int countEqualToData = 0; 
  int row, col, indexchain;
  int count_similarities = 0; 
  int minx, maxx, miny, maxy;
  for (indexchain = 0; indexchain < nchains; indexchain ++){
      countEqualToData = 0;
      count_similarities = 0;
      for (row = 0; row < imgsize; row++){
          for (col = 0; col < imgsize; col ++){
              countEqualToData += (p_chains[indexchain + nchains * (imgsize * col + row)] 
                      == p_data[imgsize * col + row]);
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
      REAL(resultsSimilar)[indexchain] = count_similarities;
      REAL(resultsData)[indexchain] = countEqualToData;
  }
  SET_VECTOR_ELT(list, 0, resultsData);
  SET_VECTOR_ELT(list, 1, resultsSimilar);
  UNPROTECT(nprotect);
  return list;
}

//SEXP IsingCompareToData(SEXP chains, SEXP dataimg, SEXP imagesize){
//  int imgsize;
//  int nprotect = 0;
//  int *dim;
//  int *p_chains;
//  int *p_data;
//  dim = INTEGER(getAttrib(chains, R_DimSymbol));
//  PROTECT(chains = coerceVector(chains, INTSXP)); nprotect ++;
//  PROTECT(dataimg = coerceVector(dataimg, INTSXP)); nprotect ++;
//  p_chains = INTEGER(chains);
//  p_data = INTEGER(dataimg);
//  imgsize = INTEGER_VALUE(imagesize);
//  int nchains = dim[0];
//  SEXP results;
//  PROTECT(results = allocVector(REALSXP, nchains)); nprotect ++;
//  int count_similarities = 0; 
//  int row, col, indexchain;
//  for (indexchain = 0; indexchain < nchains; indexchain ++){
//    count_similarities = 0;
//    for (row = 0; row < imgsize; row++){
//      for (col = 0; col < imgsize; col ++){
//           count_similarities += (p_chains[indexchain + nchains * (imgsize * col + row)] 
//               == p_data[imgsize * col + row]);
//      }
//    }
//    REAL(results)[indexchain] = count_similarities;
//  }
//  UNPROTECT(nprotect);
//  return results;
//}
//
//SEXP IsingCountSimilarities(SEXP chains, SEXP imagesize){
//  int imgsize;
//  int nprotect = 0;
//  int *dim;
//  int *p_chains;
//  dim = INTEGER(getAttrib(chains, R_DimSymbol));
//  PROTECT(chains = coerceVector(chains, INTSXP)); nprotect ++;
//  p_chains = INTEGER(chains);
//  imgsize = INTEGER_VALUE(imagesize);
//  int nchains = dim[0];
//  SEXP results;
//  PROTECT(results = allocVector(REALSXP, nchains)); nprotect ++;
//  int count_similarities = 0; 
//  int row, col, indexchain;
//  int minx, maxx, miny, maxy;
//  for (indexchain = 0; indexchain < nchains; indexchain ++){
//      count_similarities = 0;
//      for (row = 0; row < imgsize; row++){
//          for (col = 0; col < imgsize; col ++){
//              minx = (int) fmax(col - 1, 0);
//              maxx = (int) fmin(col + 1, imgsize - 1);
//              miny = (int) fmax(row - 1, 0);
//              maxy = (int) fmin(row + 1, imgsize -1 );
//              for (int dx = minx; dx <= maxx; dx++){
//                  for (int dy = miny; dy <= maxy; dy++){
//                      count_similarities += (
//                              p_chains[indexchain + nchains * (imgsize * dx + dy)]
//                              == p_chains[indexchain + nchains * (imgsize * col + row)]);
//                  }
//              }
//          }
//      }
//      count_similarities -= imgsize * imgsize;
//      REAL(results)[indexchain] = count_similarities;
//  }
//  UNPROTECT(nprotect);
//  return results;
//}
//
