#ifndef RCASPER
#define RCASPER

int lencdf;  //define global variables needed by cumu_fragsta
double *startcdf;
double cumu_fragsta(double x);

extern "C" {

  SEXP calcDenovoSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorprobR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP verboseR);

}

Casper* initCasper(int *exons, int *exonwidth, SEXP transcriptsR, int geneid, int nexons, SEXP pathCountsR, double *fraglen, int *lenvals, int nfraglen, int readLength, SEXP fragstaR, double priorq, int verbose);

Variant* path2Variant(DataFrame *df, Fragment* f, Gene *gene);
#endif
