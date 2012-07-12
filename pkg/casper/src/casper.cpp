#include "casper.h"
#include <algorithm>

using namespace std;

const int Casper::is_runs = 100;
const int Casper::em_maxruns = 100;
const double Casper::em_tol = 0.001;
const double Casper::mh_gammah = 2;
double Casper::priorq = 3;

Casper::Casper(Model* model, DataFrame* frame)
{
	this->model = model;
	this->frame = frame;

	vector<Variant* >::const_iterator vi;
	for (vi = model->items.begin(); vi != model->items.end(); vi++)
	{
	  Variant *v = *vi;
		map<Fragment*, double> table = this->frame->probabilities(v);

		memvprobs[v] = table;
        
		map<Fragment*, double>::const_iterator fi;
		for (fi = table.begin(); fi != table.end(); fi++)
		{
			Fragment* f = fi->first;
			mempprobs[f][v] = fi->second;
		}
	}
}

double* Casper::calculateMode() {
  int n = model->count();
  double* pi = new double[n];
  for (int i = 0; i < n; i++) pi[i] = 1.0 / (double)n;
  this->calculateMode(pi);
  return pi;
}

void Casper::calculateMode(double* pi) {

  double normali = (double)memvprobs.size() * (priorq - 1.0);
  map<Fragment*, map<Variant*, double> >::const_iterator ofi;
  for (ofi = mempprobs.begin(); ofi != mempprobs.end(); ofi++) {
    normali += ofi->first->count;
  }

  double err= 1.0, newpi;
  int r;
  for (r = 0; (r < em_maxruns) & (err>em_tol); r++) {
    map<Fragment*, double> mem = fragdist(pi);
		
    map<Variant*, map<Fragment*, double> >::const_iterator vi;
    for (vi = memvprobs.begin(), err=0; vi != memvprobs.end(); vi++) {
      int i = model->indexOf(vi->first);

      double nsum = 0;
      map<Fragment*, double>::const_iterator fi;
      for (fi = vi->second.begin(); fi != vi->second.end(); fi++) {
	nsum += (double)fi->first->count * fi->second / mem[fi->first];
      }

      newpi = (priorq - 1.0 + nsum * pi[i]) / normali;
      err= max_xy(err, fabs(newpi - pi[i]));
      pi[i]= newpi;
    }
  }
}

void Casper::IPMH(double *pi, double *paccept, int niter, int burnin) {
  double* mode = calculateMode();
  IPMH(pi, paccept, niter, burnin, mode);
  delete [] mode;
}

void Casper::IPMH(double *pi, double *paccept, int niter, int burnin, double *mode) {
  double **S;
  int n = model->count();
  S= dmatrix(0,n,0,n);
  normapprox(S, mode, n);
  IPMH(pi, paccept, niter, burnin, mode, S);
  free_dmatrix(S,0,n,0,n);
}

void Casper::IPMH(double *pi, double *paccept, int niter, int burnin, double *mode, double **S) {
  double det, lold, lnew;
  double *thmode, *thold, *thnew, *piold, *pinew, **Gold, **Gnew, **cholS, **cholSinv;
  int n = model->count();

  //Pre-compute useful quantities
  thmode = new double[n - 1];
  mlogit(thmode, mode, n);

  cholS= dmatrix(0,n,0,n);
  cholSinv= dmatrix(0,n,0,n);
  choldc(S,n,cholS);
  choldc_inv(S,n,cholSinv); 
  det= choldc_det(cholSinv,n);

  //MCMC
  thold = new double[n - 1];
  thnew = new double[n - 1];
  piold = new double[n];
  pinew = new double[n];
  Gold = dmatrix(0,n,0,n);
  Gnew = dmatrix(0,n,0,n);

  rmvtC(thold, n, thmode, cholS, 3);
  milogit(piold, thold, n);
  lold= priorLikelihoodLn(piold) - dmvtC(thold, n, thmode, cholSinv, det, 3, 1);;

  vtGradG(Gold,thold, n);
  lold+= vtGradLogdet(Gold, n);

  (*paccept)= 0;
  for (int i=0; i<niter; i++) {
    rmvtC(thnew, n, thmode, cholS, 3);
    milogit(pinew, thnew, n);
    lnew= priorLikelihoodLn(pinew) - dmvtC(thnew, n, thmode, cholSinv, det, 3, 1);
    vtGradG(Gnew,thnew,n);
    lnew+= vtGradLogdet(Gnew, n);
    double p= exp(lnew - lold);
    double u = runif();
    if (u <= p)	{
      (*paccept) += 1;
      double ltemp= lnew;
      lnew = lold; lold = ltemp;
      double *thtemp;
      thtemp= thnew; thnew= thold; thold= thtemp;
      double *pitemp;
      pitemp= pinew; pinew= piold; piold= pitemp;
      double **Gtemp;
      Gtemp= Gnew; Gnew= Gold; Gold= Gtemp;
    } 

    if (i>=burnin) {
      int idx= i-burnin;
      for (int j=0; j<n; j++) pi[idx+j*niter]= piold[j];
    }
  }
  (*paccept) = (*paccept)/(niter+.0);

  delete [] thmode;
  delete [] thold;
  delete [] thnew;
  delete [] piold;
  delete [] pinew;
  free_dmatrix(Gold,0,n,0,n);
  free_dmatrix(Gnew,0,n,0,n);
  free_dmatrix(cholS,0,n,0,n);
  free_dmatrix(cholSinv,0,n,0,n);
}

double Casper::calculateIntegral() {
  int n = model->count();
  double* mode = calculateMode();
  double ans= calculateIntegral(mode, n);
  delete [] mode;
  return ans;
}
double Casper::calculateIntegral(double *mode, int n)
{
  if (n == 1) { return priorLikelihoodLn(mode); }

  double *thmode, ***H, **G, **S;

  thmode = new double[n - 1];
  mlogit(thmode, mode, n);

  H= darray3(n,n,n);
  vtHess(H, thmode, n);

  G = dmatrix(0,n,0,n);
  vtGradG(G,thmode, n);

  S= dmatrix(0,n,0,n);
  normapprox(S, G, H, mode, thmode, n);

  double emlk = priorLikelihoodLn(mode);
  double gdet = vtGradLogdet(G, n);
  double sdet = log(det(S, n - 1));

  double integral = emlk + gdet + (double)(n - 1) / 2.0 * log(2 * M_PI) - 0.5 * sdet;

  delete [] thmode;
  free_darray3(H,n,n,n);
  free_dmatrix(G,0,n,0,n);
  free_dmatrix(S,0,n,0,n);
	
  return integral;
}


bool Casper::isValid()
{
	list<Fragment*>::const_iterator fi;
	for (fi = frame->data.begin(); fi != frame->data.end(); fi++)
	{
		Fragment* f = *fi;
		if (mempprobs.count(f) == 0 || mempprobs[f].size() == 0)
        {
	        return false;
		}
	}

	return true;
}

double Casper::priorLn(double* pi)
{
	int n = model->count();

	double* alpha = new double[n];
	for (int i = 0; i < n; i++) alpha[i] = priorq;

	int log = 1;
	double ans= ddirichlet(pi, alpha, &n, &log);

	delete [] alpha;

	return ans;
}

double Casper::likelihoodLn(double* pi)
{
	double outer = 0;

	map<Fragment*, map<Variant*, double> >::const_iterator fi;
	for (fi = mempprobs.begin(); fi != mempprobs.end(); fi++)
	{
		double inner = 0;

		map<Variant*, double>::const_iterator vi;
		for (vi = fi->second.begin(); vi != fi->second.end(); vi++)
		{
			int i = model->indexOf(vi->first);
			inner += pi[i] * vi->second;
		}
		outer += fi->first->count * log(inner);
	}
	return outer;
}
double Casper::priorLikelihoodLn(double* pi)
{
	return priorLn(pi) + likelihoodLn(pi);
}

map<Fragment*, double> Casper::fragdist(double* pi)
{
	map<Fragment*, double> mem;

	map<Fragment*, map<Variant*, double> >::const_iterator fi;
	for (fi = mempprobs.begin(); fi != mempprobs.end(); fi++)
	{
		double sum = 0;

		map<Variant*, double>::const_iterator vi;
		for (vi = fi->second.begin(); vi != fi->second.end(); vi++)
		{
			int i = model->indexOf(vi->first);
			sum += pi[i] * vi->second;
		}
		mem[fi->first] = sum;
	}

	return mem;
}


void Casper::normapprox(double **S, double *mode, int n) {
  double **G, ***H, *thmode;

  thmode = new double[n - 1];
  mlogit(thmode, mode, n);

  H= darray3(n,n,n);
  vtHess(H, thmode, n);

  G = dmatrix(0,n,0,n);
  vtGradG(G,thmode, n);

  normapprox(S, G, H, mode, thmode, n);

  delete [] thmode;
  free_darray3(H,n,n,n);
  free_dmatrix(G,0,n,0,n);
}


void Casper::normapprox(double **S, double** G, double*** H, double* mode, double* thmode, int n)
{
	map<Fragment*, double> mem = fragdist(mode);

	for (int l = 0; l < n - 1; l++)
	{
		for (int m = l; m < n - 1; m++)
		{
			S[l][m] = 0;
			map<Fragment*, map<Variant*, double> >::const_iterator fi;
			for (fi = mempprobs.begin(); fi != mempprobs.end(); fi++)
			{
				double term1 = 0, term2 = 0, term3 = 0;
				map<Variant*, double>::const_iterator vi;
				for (vi = fi->second.begin(); vi != fi->second.end(); vi++)
				{
					int d = model->indexOf(vi->first);
					double P = vi->second;
					term1 += P * H[d][l][m];
					term2 += P * G[d][l];
					term3 += P * G[d][m];
				}
				S[l][m] -= fi->first->count * (term1 * mem[fi->first] - term2 * term3) / pow(mem[fi->first], 2);
			}
			for (int d = 0; d < n; d++)
			{
				S[l][m] -= (priorq - 1.0) * (H[d][l][m] * mode[d] - G[d][l] * G[d][m]) / pow(mode[d], 2);
			}
			if (l != m)
			{
				S[m][l] = S[l][m];
			}
		}
	}

}
void Casper::mlogit(double *theta, double* pi, int n) {

  for (int i = 0; i < n - 1; i++) theta[i] = log(pi[i + 1] / pi[0]);

}
void Casper::milogit(double *pi, double* theta, int n)
{
	double sum = 1.0;

	for (int i = 0; i < n - 1; i++) sum += exp(theta[i]);

	pi[0] = 1.0 / sum;
	for (int i = 0; i < n - 1; i++) pi[i + 1] = exp(theta[i]) / sum;

}
void Casper::vtGradG(double **G, double* th, int n)
{
	double sum = 1.0;
	for (int i = 0; i < n - 1; i++)
	{
		sum += exp(th[i]);
	}

	for (int l = 0; l < n - 1; l++)
	{
		G[0][l] = -exp(th[l]) / pow(sum, 2);
	}

	for (int d = 1; d < n; d++)
	{
		for (int l = 0; l < n - 1; l++)
		{
			if (l != d - 1)
			{
				G[d][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2);
			}
			else
			{
				G[d][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2) + exp(th[l]) / sum;
			}
		}
	}
	//return G;
}
double Casper::vtGradLogdet(double** G, int n)
{
  double** GS = &G[1];
  double mydet = det(GS, n - 1);
  if (mydet < 0) mydet = -mydet;
  double logdet = log(mydet);
  return logdet;
}
void Casper::vtHess(double ***H, double* th, int n)
{
	double sum = 1.0;
	for (int i = 0; i < n - 1; i++)
	{
		sum += exp(th[i]);
	}

	for (int d = 0; d < n; d++)
	{
		for (int l = 0; l < n - 1; l++)
		{
			for (int m = l; m < n - 1; m++)
			{
				if (d == 0)
				{
					if (m == l)
					{
						H[0][l][l] = -exp(th[l]) / pow(sum, 2) + 2.0 * exp(2.0 * th[l]) / pow(sum, 3);
					}
					else
					{
						H[0][l][m] = H[0][m][l] = 2.0 * exp(th[l] + th[m]) / pow(sum, 3);
					}
				}
				else
				{
					if (m == l)
					{
						if (l != d - 1)
						{
							H[d][l][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2) + 2 * exp(th[d - 1] + 2 * th[l]) / pow(sum, 3);
						}
						else
						{
							H[d][l][l] = -2 * exp(2 * th[d - 1]) / pow(sum, 2) + 2 * exp(3 * th[d - 1]) / pow(sum, 3) + exp(th[d - 1]) / sum - exp(2 * th[d - 1]) / pow(sum, 2);
						}
					}
					else
					{
						if (m == d - 1)
						{
							H[d][l][m] = H[d][m][l] = -exp(th[d - 1] + th[l]) / pow(sum, 2) + 2 * exp(th[d - 1] + th[m] + th[l]) / pow(sum, 3);
						}
						else
						{
							if (l == d - 1)
							{
								H[d][l][m] = H[d][m][l] = 2 * exp(th[d - 1] + th[l] + th[m]) / pow(sum, 3) - exp(th[m] + th[l]) / pow(sum, 2);
							}
							else
							{
								H[d][l][m] = H[d][m][l] = 2 * exp(th[d - 1] + th[l] + th[m]) / pow(sum, 3);
							}
						}
					}
				}
			}
		}
	}

}

double Casper::det(double** a, int n)
{
	double **aout = dmatrix(0, n - 1, 0, n - 1);

	int i,j,k;
	double sum;

	for (i=0;i<n;i++) { for (j=i;j<n;j++) { aout[i][j]= a[i][j]; } }  //copy a into aout
	for (i=0;i<n;i++) {
		for (j=i;j<n;j++) {
			for (sum=aout[i][j],k=i-1;k>=0;k--) sum -= aout[i][k]*aout[j][k];
			if (i == j) {
			  if (sum <= 0.0) {
                            char proc[]="choldc failed", act[]="", what[]="matrix is not positive definite";
                            nrerror(proc,act,what); 
			  }
				aout[i][i]=sqrt(sum);
			} else aout[j][i]=sum/aout[i][i];
		}
	}
	for (i=0;i<n;i++) { for (j=i+1;j<n;j++) { aout[i][j]= 0; } }  //set upper-diagonal elem to 0


	double det;
	for (det=1,i=0; i<n; i++) { det *= aout[i][i]*aout[i][i]; }
	
	free_dmatrix(aout, 0, n - 1, 0, n - 1);
	return(det);
}
