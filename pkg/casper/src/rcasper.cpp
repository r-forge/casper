#include <sstream>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include "stdafx.h"
#include "casper.h"
#include "rcasper.h"
#include "seppel.h"

int verbose = 1;

double SIMPLEFRAGSTA(double x)
{
	return x;
}
double cumu_fragsta(double x) 
{
	if (x<=0) return 0;
	if (x>=1) return 1;
	int idx= (int) (x*lencdf);
	double y1= startcdf[idx], x1= (double) idx / (lencdf-1);
	idx++;
	double y2= startcdf[idx], x2= (double) idx / (lencdf-1);
	return y1 + (x-x1) * (y2-y1)/(x2-x1);
}

DataFrame* importDataFrame(SEXP exonsR, SEXP exonwidthR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR) 
{
	int nexons=Rf_length(exonsR), nfraglen=Rf_length(fraglenR), readLength=INTEGER(readLengthR)[0];
	int *exons=INTEGER(exonsR), *exonwidth=INTEGER(exonwidthR), *lenvals=INTEGER(lenvalsR);
	double *fraglen= REAL(fraglenR);

	//Define fragment length/start distributions and create DataFrame
	//DiscreteDF* fraglen_dist = new DiscreteDF(fraglen, lenvals, nfraglen);	
	double* fraglens = new double[1];
        int* lenvals= new int[1];
        fraglens[0]= 1;
        lenvals[0]= 200;
	DiscreteDF* fraglen_dist = new DiscreteDF(fraglens, lenvals, 1);

	lencdf= Rf_length(fragstaR);
	startcdf= REAL(fragstaR);
	//::fun_fragsta = fragstaR;

	//DataFrame* df = new DataFrame(fraglen_dist, cumu_fragsta);
	DataFrame* df = new DataFrame(fraglen_dist, SIMPLEFRAGSTA);
	df->frag_readlen= readLength;

	//Add exons to DataFrame
	for (int i = 0; i < nexons; i++) 
	{
		Exon *ex= new Exon(exons[i], exonwidth[i]);
		df->addExon(ex);
	}

	//Add exon path counts
	int np = LENGTH(pathCountsR);
	SEXP pnames = getAttrib(pathCountsR, R_NamesSymbol);
	int* pathCounts = INTEGER(pathCountsR);  

	for (int i = 0; i < np; i++) 
	{
		//Set counts
		int count = pathCounts[i];

		//Set nb of left & right visited exons
		const char* pname = CHAR(STRING_ELT(pnames, i));
		char* left = new char[strlen(pname) + 1];
		strcpy(left, pname);
		char* mid = strchr(left, '-');
		if (mid == NULL) continue;
		mid[0] = '\0';
		char* right = mid+1;
		int leftc = 0, rightc= 0;
		for (int l = strlen(left)-1; l >= 0; l--) { if (left[l] == '.') leftc++; }  
		for (int r = strlen(right)-1; r >= 0; r--) { if (right[r] == '.') rightc++; } 

		Fragment* f = new Fragment(leftc, rightc, count);

		//Set sequence of visited exons
		left = left+1;
		right[strlen(right)-1] = '\0';
		char* item;
		item = strtok(left, ".");  //split string
		for (int j = 0; item != NULL; j++) {
			int eid = atoi(item);
			f->left[j] = eid;
			item = strtok(NULL, ".");
		}
		item = strtok(right, ".");
		for (int j = 0; item != NULL; j++) {
			int eid = atoi(item);
			f->right[j] = eid;
			item = strtok(NULL, ".");
		}

		df->addData(f);
	}

	return df;
}
set<Variant*, VariantCmp>* importTranscripts(DataFrame* df, SEXP transcriptsR)
{
	int nt = LENGTH(transcriptsR);
	set<Variant*, VariantCmp>* initvars = new set<Variant*, VariantCmp>();
	SEXP tnames = getAttrib(transcriptsR, R_NamesSymbol);
	Variant *v;
	SEXP trow;
	int ntsub, *tvals;
	for (int i = 0; i < nt; i++) {
		trow = VECTOR_ELT(transcriptsR, i);
		ntsub = LENGTH(trow);
		tvals = INTEGER(trow);

		vector<Exon*>* el = new vector<Exon*>();
		for (int s = 0; s < ntsub; s++) {
			int eid = tvals[s];
			Exon* ex = df->id2exon[eid];
			el->push_back(ex);
		}

		v = new Variant(el);
		v->id= i;
		int nbchar= Rf_length(STRING_ELT(tnames,i));
		v->name= string(CHAR(STRING_ELT(tnames, i)), nbchar);

		initvars->insert(v);
	}

	return initvars;
}

extern "C"
{

  //Estimate isoform expression for multiple genes. Calls calcKnownSingle repeadtedly
  //Input args same as for calcDenovoMultiple
  SEXP calcKnownMultiple(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorqR) 
	{
		int i, ngenes=LENGTH(geneidR);
		SEXP ansMultiple, ansSingle;

		PROTECT(ansMultiple= allocVector(VECSXP, ngenes));

		for (i=0; i<ngenes; i++) {
		  ansSingle= calcKnownSingle(VECTOR_ELT(exonsR,i), VECTOR_ELT(exonwidthR,i), VECTOR_ELT(transcriptsR,i), VECTOR_ELT(pathCountsR,i), fragstaR, fraglenR, lenvalsR, readLengthR, VECTOR_ELT(geneidR,i), priorqR);
		  SET_VECTOR_ELT(ansMultiple,i,ansSingle);
		}
		
		UNPROTECT(1);
		return ansMultiple;
	}

  //Estimate isoform expression (known variants case)
  //Input args same as for calcDenovoSingle
  SEXP calcKnownSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP geneidR, SEXP priorqR)
	{
		DataFrame* df = importDataFrame(exonsR, exonwidthR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR);
		set<Variant*, VariantCmp>* initvars = importTranscripts(df, transcriptsR);
		
		double priorq = REAL(priorqR)[0];

		Model* model = new Model(new vector<Variant*>(initvars->begin(), initvars->end()));
		Casper* casp = new Casper(model, df);
		casp->priorq = priorq;

		int vc = model->count();
		double* em = casp->calculateMode();

		SEXP Rc;
		PROTECT(Rc = allocVector(REALSXP, vc));
		double* res = REAL(Rc);
		for (int i = 0; i < vc; i++)
		{
			res[i] = em[i];
		}
		UNPROTECT(1);

		return Rc;
	}

  SEXP calcDenovoMultiple(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP methodR, SEXP verboseR) 
	{
		//De novo isoform discovery and estimate expression for multiple genes. Calls calcDenovoSingle repeadtedly
		int i, ngenes=LENGTH(geneidR);
		SEXP ansMultiple, ansSingle;

		PROTECT(ansMultiple= allocVector(VECSXP, ngenes));

		for (i=0; i<ngenes; i++) {
		  ansSingle= calcDenovoSingle(VECTOR_ELT(exonsR,i), VECTOR_ELT(exonwidthR,i), VECTOR_ELT(transcriptsR,i), VECTOR_ELT(pathCountsR,i), fragstaR, fraglenR, lenvalsR, readLengthR, VECTOR_ELT(geneidR,i), nvarPriorR, nexonPriorR, priorqR, minppR, selectBest, methodR, verboseR);
			SET_VECTOR_ELT(ansMultiple,i,ansSingle);
		}
		
		UNPROTECT(1);
		return ansMultiple;
	}

  SEXP calcDenovoSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP geneidR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP methodR, SEXP verboseR)
	{
		//De novo isoform discovery and estimate expression for a single gene
		//Input
		// - exons: vector with exon ids
		// - exonwidth: vector exon widths
		// - geneid: integer with gene id
		// - transcripts: list of transcripts. Names indicate transcript id. Each element contains list of exon ids. Used to initialize model
		// - pathCounts: vector with path counts. Names indicate series of visited exons e.g. e1.e2-e5
		// - fragsta: function that returns start distrib cdf
		// - fraglen: vector with fragment length distrib, i.e. P(length=lenvals[0]),P(length=lenvals[1]),... up to max length
		// - lenvals: vector with possibles values for length
		// - readLength: read length
		// - nvarPrior: prior prob of nb expressed variants. i^th elem has NegBinom param for genes with i exons
		// - nexonPrior: prior prob of nb exons in a variant. i^th elem has Beta-Binomial param for genes with i exons
		// - priorq: prior on model-specific parameters pi ~ Dirichlet(priorq)
		// - minpp: only models with post prob >= minpp are reported
		// - selectBest: set to !=0 to return only results for model with highest posterior probability
	        // - methodR: set to 1 for exact; 2 for random walk MCMC; 3 for prior proposal MCMC; 0 for auto (exact for genes with <=4 exons, RW-MCMC for >4 exons)
		// - verbose: set verbose != 0 to print warnings & progress info

		int geneid = INTEGER(geneidR)[0];
		int selBest = INTEGER(selectBest)[0];
		double minpp = REAL(minppR)[0];
		double priorq = REAL(priorqR)[0];
                int method= INTEGER(methodR)[0];
		verbose = INTEGER(verboseR)[0];

		SEXP exongenesR;
		PROTECT(exongenesR = allocVector(INTSXP, LENGTH(exonsR)));
		for (int i = 0; i < LENGTH(exonsR); i++) INTEGER(exongenesR)[i] = geneid;

		DataFrame* df = importDataFrame(exonsR, exonwidthR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR);
		set<Variant*, VariantCmp>* initvars = importTranscripts(df, transcriptsR);

		UNPROTECT(1);	

		df->debugprint();
		Model* tmpm = new Model(initvars);
		tmpm->debugprint();

		// END OF INPUT READING

		int discarded = df->fixUnexplFrags(initvars);
		Model* tmp2 = new Model(initvars);
		tmp2->debugprint();
		if (verbose > 0)
		{
			Rprintf("discarded %i fragments\n", discarded);
		}

		Casper::priorq = priorq;

		map<Model*, double, ModelCmp> resProbs;
		map<Model*, double*, ModelCmp> resModes;
		
		Seppel* seppl = new Seppel(df);

		if (method == 1 || method == 0 && df->exons.size() <= 4) 
		{
			seppl->exploreExact();
			resProbs = seppl->resultPPIntegral();
		} 
		else if (method == 3)
		{
			seppl->explorePrior(100000);
			resProbs = seppl->resultPPIntegral();
		}
		else
		{
			Model* startmodel = new Model(initvars);
			seppl->exploreSmart(startmodel, 100000);
			resProbs = seppl->resultPPMCMC();
		}
		resModes = seppl->resultModes();

		vector<Variant*>* allpossvariants = df->allVariants();

		// END OF CALCULATIONS
		
		Model* bestModel;
		double bestModelProb = -1;

		map<Model*, double, ModelCmp>::iterator mvi = resProbs.begin();
		while (mvi != resProbs.end())
		{
			if (mvi->second > bestModelProb)
			{
				bestModelProb = mvi->second;
				bestModel = mvi->first;
			}

			if (mvi->second < minpp) 
			{
				resModes.erase(mvi->first);
				resProbs.erase(mvi++);
			}
			else
			{
				++mvi;
			}
		}

		// END OF FILTERING
		
		SEXP ans;
		PROTECT(ans= allocVector(VECSXP, 6));

		//Report posterior probabilities
		int i;
		int nrowpi=0, nx;
		if (selBest==0) { nx= resProbs.size(); } else { nx=1; }
		SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP,nx,2));
		SET_VECTOR_ELT(ans, 5, allocVector(STRSXP,nx));
		double *resProbsR= REAL(VECTOR_ELT(ans,0));
		
		set<Variant*, VariantCmp>::const_iterator vi;
		map<Model*, double, ModelCmp>::const_iterator mi;
		if (selBest==0) {
			for (mi = resProbs.begin(), i=0; mi != resProbs.end(); mi++, i++) {
				Model* m = mi->first;
				resProbsR[i]= i;
				resProbsR[i+nx]= resProbs[m];
				nrowpi+= m->count();
				SET_STRING_ELT(VECTOR_ELT(ans,5),i,mkChar(m->getCodeStr(allpossvariants)));  //model name
			}
		} else {
			resProbsR[0]= 0;
			resProbsR[nx]= resProbs[bestModel];
			nrowpi+= bestModel->count();
			SET_STRING_ELT(VECTOR_ELT(ans,5),0,mkChar(bestModel->getCodeStr(allpossvariants)));  //model name
		}

		//Report estimated expression
		SET_VECTOR_ELT(ans, 1, allocMatrix(REALSXP,nrowpi,2));  //stores model id and estimated expression
		SET_VECTOR_ELT(ans, 2, allocVector(STRSXP,nrowpi)); //stores variant names
		double *expr= REAL(VECTOR_ELT(ans,1));
		SEXP varnamesR= VECTOR_ELT(ans,2);

		int rowid=0, nrowexons=0;
		set<Variant*, VariantCmp> allvariants;
		if (selBest==0) {
			for (mi = resProbs.begin(), i=0; mi != resProbs.end(); mi++, i++) {
				Model* m = mi->first;
				for (int j=0; j< m->count(); j++) {
					expr[rowid]= i;  //model id
					Variant* v = m->get(j);
					int varidx= m->indexOf(v);
					expr[rowid+nrowpi]= resModes[m][varidx]; //estimated expression
					if (initvars->count(v)>0) v->name= (*initvars->find(v))->name;  //respect initial variant names
					const char *cname= (v->name).c_str();;
					SET_STRING_ELT(varnamesR,rowid,mkChar(cname));  //variant name
					if (allvariants.count(v)==0) {
						allvariants.insert(v);
						nrowexons+= v->exonCount;
					}
					rowid++;
				}
			}
		} else {
			Model* m = bestModel;
			for (int j=0; j< m->count(); j++) {
				expr[rowid]= i;  //model id
				Variant* v = m->get(j);
				int varidx= m->indexOf(v);
				expr[rowid+nrowpi]= resModes[m][varidx]; //estimated expression
				if (initvars->count(v)>0) v->name= (*initvars->find(v))->name;  //respect initial variant names
				const char *cname= (v->name).c_str();
				SET_STRING_ELT(varnamesR,rowid,mkChar(cname));  //variant name
				if (allvariants.count(v)==0) {
					allvariants.insert(v);
					nrowexons+= v->exonCount;
				}
				rowid++;
			}
		}

		//Report exons in each variant
		SET_VECTOR_ELT(ans, 3, allocVector(INTSXP,nrowexons));  //stores exon id
		SET_VECTOR_ELT(ans, 4, allocVector(STRSXP,nrowexons));   //stores variant name
		int *ranges= INTEGER(VECTOR_ELT(ans,3));
		SEXP varnames2R= VECTOR_ELT(ans,4);

		rowid=0;
		for (vi= allvariants.begin(); vi != allvariants.end(); vi++) {
			Variant *v= *vi;
			for (int j=0; j< v->exonCount; j++) {
				ranges[rowid]= v->exons[j]->id ; //exon id
				const char *cname= (v->name).c_str();
				SET_STRING_ELT(varnames2R,rowid,mkChar(cname));  //variant name
				rowid++;
			}
		}

		UNPROTECT(1);
		return(ans);
	}
}
