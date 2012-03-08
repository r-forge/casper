#include "casper.h"
#include <algorithm>
#include <time.h>

using namespace std;

const int Casper::is_runs = 100;
const int Casper::em_maxruns = 100;
const double Casper::mh_gammah = 2;
const double Casper::prior_q = 2;

Casper::Casper(Model* model, DataFrame* frame)
{
	this->model = model;
	this->frame = frame;

	vector<Variant* >::const_iterator vi;
	for (vi = model->items.begin(); vi != model->items.end(); vi++)
	{
		Variant* v = *vi;
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
double* Casper::calculateMode()
{
	int n = model->count();

	double* pi = new double[n];
	for (int i = 0; i < n; i++)
	{
		pi[i] = 1.0 / (double)n;
	}

	double normali = (double)memvprobs.size() * (prior_q - 1.0);
	map<Fragment*, map<Variant*, double> >::const_iterator ofi;
	for (ofi = mempprobs.begin(); ofi != mempprobs.end(); ofi++)
	{
		normali += ofi->first->count;
	}

	for (int r = 0; r < em_maxruns; r++)
	{
		map<Fragment*, double> mem = fragdist(pi);
		
		map<Variant*, map<Fragment*, double> >::const_iterator vi;
		for (vi = memvprobs.begin(); vi != memvprobs.end(); vi++)
		{
			int i = model->indexOf(vi->first);

			double nsum = 0;
			map<Fragment*, double>::const_iterator fi;
			for (fi = vi->second.begin(); fi != vi->second.end(); fi++)
			{
				nsum += (double)fi->first->count * fi->second / mem[fi->first];
			}

            pi[i] = (prior_q - 1.0 + nsum * pi[i]) / normali;
		}
	}

	return pi;
}
double Casper::calculateIntegral()
{
	int n = model->count();
	double* mode = calculateMode();
	if (n == 1)
	{
		return priorLikelihoodLn(mode);
	}

	double* thmode = mlogit(mode, n);
    double*** H = vtHess(thmode, n);
    double** G = vtGradG(thmode, n);
    double** S = normapprox(G, H, mode, thmode, n);

	double emlk = priorLikelihoodLn(mode);
	double gdet = vtGradLogdet(G, n);
	double sdet = log(det(S, n - 1));

	double integral = emlk + gdet + (double)(n - 1) / 2.0 * log(2 * M_PI) - 0.5 * sdet;

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

bool Casper::isFragValid(Fragment *f) {
  return(mempprobs.count(f) > 0 && mempprobs[f].size() > 0);
}

double Casper::priorLn(double* pi)
{
	int n = model->count();

	double* alpha = new double[n];
	for (int i = 0; i < n; i++)
    {
		alpha[i] = prior_q;
	}

	int log = 1;
	return ddirichlet(pi, alpha, &n, &log);
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

double** Casper::normapprox(double** G, double*** H, double* mode, double* thmode, int n)
{
	map<Fragment*, double> mem = fragdist(mode);
	double** S = new double*[n - 1];
	for (int i = 0; i < n - 1; i++)
	{
		S[i] = new double[n - 1];
	}
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
				S[l][m] -= (prior_q - 1.0) * (H[d][l][m] * mode[d] - G[d][l] * G[d][m]) / pow(mode[d], 2);
			}
			if (l != m)
			{
				S[m][l] = S[l][m];
			}
		}
	}
	return S;
}
double* Casper::mlogit(double* pi, int n)
{
	double* theta = new double[n - 1];
	for (int i = 0; i < n - 1; i++)
	{
		theta[i] = log(pi[i + 1] / pi[0]);
	}
	return theta;
}
double* Casper::milogit(double* theta, int n)
{
	double sum = 1.0;
	for (int i = 0; i < n - 1; i++)
	{
		sum += exp(theta[i]);
	}
	double* pi = new double[n];
	pi[0] = 1.0 / sum;
	for (int i = 0; i < n - 1; i++)
	{
		pi[i + 1] = exp(theta[i]) / sum;
	}
	return pi;
}
double** Casper::vtGradG(double* th, int n)
{
	double sum = 1.0;
	for (int i = 0; i < n - 1; i++)
	{
		sum += exp(th[i]);
	}

	double** G = new double*[n];
	for (int i = 0; i < n; i++)
	{
		G[i] = new double[n - 1];
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
	return G;
}
double Casper::vtGradLogdet(double** G, int n)
{
	double** GS = &G[1];
    double mydet = det(GS, n - 1);
	if (mydet < 0)
	{
		mydet = -mydet;
	}
	double logdet = log(mydet);
    return logdet;
}
double*** Casper::vtHess(double* th, int n)
{
	double sum = 1.0;
	for (int i = 0; i < n - 1; i++)
	{
		sum += exp(th[i]);
	}

	double*** H = new double**[n];
	for (int i = 0; i < n; i++)
	{
		H[i] = new double*[n - 1];
		for (int j = 0; j < n - 1; j++)
		{
			H[i][j] = new double[n - 1];
		}
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

	return H;
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

double Casper::randd()
{
	return (double)rand() / (double)(RAND_MAX + 1);
}
int Casper::randi(int n)
{
	return (int)((double)n * randd());
}
