#ifndef RCASPER
#define RCASPER

int lencdf;  //define global variables needed by cumu_fragsta
double *startcdf;
double cumu_fragsta(double x);

extern "C" {

  SEXP calcKnownSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorqR, SEXP strandR);

  SEXP calcDenovoSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP modelUnifPriorR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP methodR, SEXP niterR, SEXP exactMarginalR, SEXP strandR);

}

#endif
