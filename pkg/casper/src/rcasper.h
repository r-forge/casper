#ifndef RCASPER
#define RCASPER

SEXP fun_fragsta;  //define global variable needed by cumu_fragsta

double cumu_fragsta(double x);

Casper* initCasper(int *exons, int *exonwidth, SEXP transcriptsR, int geneid, int nexons, SEXP pathCountsR, double *fraglen, int *lenvals, int nfraglen, int readLength, SEXP fragstaR);

Variant* path2Variant(DataFrame *df, Fragment* f, Gene *gene);
#endif
