#ifndef CASPER
#define CASPER

#include "dataframe.h"
#include "cstat.h"
#include <list>
#include <vector>

using namespace std;

class Casper
{
 public:
  Casper(Model* model, DataFrame* frame);

  // current model (set of variants) that tries to explain the data
  Model* model;
  // data
  DataFrame* frame;

  // gives the mode given the current model and data
  double* calculateMode();

  // gives the integral given the current model and data
  double calculateIntegral();  //uses mode=calculateMode() and n= model->count()
  double calculateIntegral(double* mode, int n);  //do integral with pre-computed mode

  // indicates whether the model can explain all fragments
  bool isValid();

  // random double between 0 and 1
  static double randd();
  // random int between 0 incl and n excl
  static int randi(int n);

  static double priorq;

 private:
  map<Fragment*, map<Variant*, double> > mempprobs;
  map<Variant*, map<Fragment*, double> > memvprobs;

  static const int is_runs;
  static const int em_maxruns;
  static const double mh_gammah;
  
  double priorLn(double* pi);
  double likelihoodLn(double* pi);
  double priorLikelihoodLn(double* pi);

  map<Fragment*, double> fragdist(double* pi);

  double** normapprox(double** G, double*** H, double* mode, double* thmode, int n);
  double* mlogit(double* pi, int n);
  double* milogit(double* theta, int n);
  double** vtGradG(double* th, int n);
  double vtGradLogdet(double** G, int n);
  double*** vtHess(double* th, int n);
  double det(double** a, int n);
};

#endif
