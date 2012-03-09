#ifndef RCASPER
#define RCASPER

int lencdf;  //define global variables needed by cumu_fragsta
double *startcdf;
//SEXP fun_fragsta;  
double cumu_fragsta(double x);

Casper* initCasper(int *exons, int *exonwidth, SEXP transcriptsR, int geneid, int nexons, SEXP pathCountsR, double *fraglen, int *lenvals, int nfraglen, int readLength, SEXP fragstaR, int verbose);

Variant* path2Variant(DataFrame *df, Fragment* f, Gene *gene);
#endif
