#include <sstream>
#include <time.h>
#include <R.h>
#include <Rinternals.h>
#include "stdafx.h"
#include "casper.h"
#include "rcasper.h"
#include "seppelsmart.h"

int verbose = 1;

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

int fixUnexplFrags(set<Variant*, VariantCmp>* initvars, DataFrame* df, Gene* gene)
{
	// copy all fragments
	set<Fragment*>* queue = new set<Fragment*>(df->data.begin(), df->data.end());

	// remove the fragments from the queue we can explain with our variants
	set<Variant*, VariantCmp>::iterator vi;
	for (vi = initvars->begin(); vi != initvars->end(); vi++) 
	{
		// remove the fragments that this variant can explain from our queue
		map<Fragment*, double> probs = df->probabilities(*vi);
		map<Fragment*, double>::iterator si;
		for (si = probs.begin(); si != probs.end(); si++) 
		{
			set<Fragment*>::iterator ri = queue->find(si->first);
			if (ri != queue->end())
			{
				queue->erase(ri);
			}
		}
	}

	int discarded = 0;

	// while we still have unexplained fragments
	while (queue->size() > 0)
	{
		// pop the first fragment
		Fragment* frag = *queue->begin();
		queue->erase(queue->begin());

		Variant* nv = path2Variant(df, frag, gene);

		// check if the new variant can explain the fragment
		map<Fragment*, double> probs = df->probabilities(nv);
		if (probs.count(frag) > 0)
		{
			initvars->insert(nv);

			// delete all fragments that this variant can explain
			map<Fragment*, double>::iterator si;
			for (si = probs.begin(); si != probs.end(); si++) 
			{
				set<Fragment*>::iterator ri = queue->find(si->first);
				if (ri != queue->end())
				{
					queue->erase(ri);
				}
			}
		}
		else
		{
			// this fragment cant be explained
			discarded++;
			df->data.remove(frag);
		}
	}

	return discarded;
}
DataFrame* importDataFrame(SEXP exonsR, SEXP exonwidthR, SEXP exongenesR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR) 
{
	int nexons=Rf_length(exonsR), nfraglen=Rf_length(fraglenR), readLength=INTEGER(readLengthR)[0];
	int *exons=INTEGER(exonsR), *exonwidth=INTEGER(exonwidthR), *exongenes=INTEGER(exongenesR), *lenvals=INTEGER(lenvalsR);
	double *fraglen= REAL(fraglenR);

	//Define fragment length/start distributions and create DataFrame
	DiscreteDF* fraglen_dist = new DiscreteDF(fraglen, lenvals, nfraglen);

	lencdf= Rf_length(fragstaR);
	startcdf= REAL(fragstaR);
	//::fun_fragsta = fragstaR;

	DataFrame* df = new DataFrame(fraglen_dist, cumu_fragsta);
	df->frag_readlen= readLength;

	//Add exons & gene to DataFrame
	for (int i = 0; i < nexons; i++) 
	{
		int gid = exongenes[i];

		Gene* g;
		if (df->genes.count(gid) == 0)
		{
			g = new Gene(gid);
			df->addGene(g);
		}
		else
		{
			g = df->genes[gid];
		}

		Exon *ex= new Exon(exons[i], exonwidth[i]);
		df->addExon(ex);
		g->addExon(ex);
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
set<Variant*, VariantCmp>* importTranscripts(DataFrame* df, SEXP transcriptsR, SEXP transgenesR)
{
	int* transgenes = INTEGER(transgenesR);

	int nt = LENGTH(transcriptsR);
	set<Variant*, VariantCmp>* initvars = new set<Variant*, VariantCmp>();
	SEXP tnames = getAttrib(transcriptsR, R_NamesSymbol);
	Gene* gene;
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
			Exon* ex = df->exons[eid];
			el->push_back(ex);
		}

		int geneid = transgenes[i];
		gene = df->genes[geneid];
		v = new Variant(gene, el);
		v->id= i;
		int nbchar= Rf_length(STRING_ELT(tnames,i));
		v->name= string(CHAR(STRING_ELT(tnames, i)), nbchar);

		initvars->insert(v);
	}

	return initvars;
}
Variant* path2Variant(DataFrame *df, Fragment* f, Gene *gene) 
{
	int eid; Exon *ex;
	vector<Exon*>* el = new vector<Exon*>();
	map<int, Exon*>::iterator itexon;
	for (itexon= df->exons.begin(); (*itexon).first != f->left[0]; itexon++) {
		ex= (*itexon).second;
		el->push_back(ex);
	}
	for (int i=0; i< f->leftc; i++) {
		eid = f->left[i];
		ex = df->exons[eid];
		el->push_back(ex);
	}
	if (eid != f->right[0]) {
		eid = f->right[0];
		ex = df->exons[eid];
		el->push_back(ex);
	}
	for (int i=1; i< f->rightc; i++) {
		eid = f->right[i];
		ex = df->exons[eid];
		el->push_back(ex);
	}
	while ((*itexon).first != eid) { itexon++; }
	itexon++;
	while (itexon != df->exons.end()) {
		ex= (*itexon).second;
		el->push_back(ex);
		itexon++;
	}
	Variant* v = new Variant(gene, el);
	return v;
}

void debugdf(DataFrame* df)
{		
	REprintf("Genes:\n");
	map<int, Gene*>::const_iterator gi;
	for (gi = df->genes.begin(); gi != df->genes.end(); gi++)
	{
		Gene* g = gi->second;
		REprintf("%i\t%i\n", g->id, g->exons.size());
		vector<Exon*>::const_iterator ei;
		for (ei = g->exons.begin(); ei != g->exons.end(); ei++)
		{
			REprintf("%i\n", (*ei)->id);
		}
	}
	REprintf("\n");

	// Exons
	REprintf("Exons:\n");
	map<int, Exon*>::const_iterator ei;
	for (ei = df->exons.begin(); ei != df->exons.end(); ei++)
	{
		Exon* e = ei->second;
		REprintf("%i\t%i\n", e->id, e->length);
	}
	REprintf("\n");

	// Fragments
	REprintf("Fragments:\n");
	list<Fragment*>::const_iterator fi;
	for (fi = df->data.begin(); fi != df->data.end(); fi++)
	{
		Fragment* f = *fi;
		REprintf("%i\t%i\t%i\n", f->leftc, f->rightc, f->count);
		for (int l = 0; l < f->leftc; l++)
		{
			REprintf("%i\n", f->left[l]);
		}
		for (int r = 0; r < f->rightc; r++)
		{
			REprintf("%i\n", f->right[r]);
		}
	}
	REprintf("\n");
}
void debugmodel(Model* model)
{
	// Transcripts
	REprintf("Model:\n");
	vector<Variant*>::const_iterator vi;
	for (vi = model->items.begin(); vi != model->items.end(); vi++)
	{
		Variant* v = *vi;
		REprintf("%i\t%i\n", v->id, v->exonCount);
		for (int e = 0; e < v->exonCount; e++)
		{
			REprintf("%i\n", v->exons[e]->id);
		}
	}
	REprintf("\n");
}
char* getmodelcode(vector<Variant*>* allvariants, Model* model)
{
	int n = allvariants->size();
	char* str = new char[n + 1];
	str[n] = '\0';

	for (int i = 0; i < (int)allvariants->size(); i++)
	{
		if (model->contains(allvariants->at(i)))
		{
			str[i] = '1';
		}
		else
		{
			str[i] = '0';
		}
	}

	return str;
}

extern "C"
{
	SEXP calcMode(SEXP exonsR, SEXP exonwidthR, SEXP exongenesR, SEXP transcriptsR, SEXP transgenesR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP priorqR)
	{
		srand((unsigned)time( NULL ));

		DataFrame* df = importDataFrame(exonsR, exonwidthR, exongenesR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR);
		set<Variant*, VariantCmp>* initvars = importTranscripts(df, transcriptsR, transgenesR);
		
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
		srand((unsigned)time( NULL ));

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

		SEXP exongenesR, transgenesR;
		PROTECT(exongenesR = allocVector(INTSXP, LENGTH(exonsR)));
		PROTECT(transgenesR = allocVector(INTSXP, LENGTH(transcriptsR)));
		for (int i = 0; i < LENGTH(exonsR); i++) INTEGER(exongenesR)[i] = geneid;
		for (int i = 0; i < LENGTH(transcriptsR); i++) INTEGER(transgenesR)[i] = geneid;

		DataFrame* df = importDataFrame(exonsR, exonwidthR, exongenesR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR);
		set<Variant*, VariantCmp>* initvars = importTranscripts(df, transcriptsR, transgenesR);

		UNPROTECT(1);		
		UNPROTECT(1);

		// END OF INPUT READING

		Gene* gene = df->genes[geneid];

		int discarded = fixUnexplFrags(initvars, df, gene);
		if (verbose > 0)
		{
			Rprintf("discarded %i fragments\n", discarded);
		}

		Casper::priorq = priorq;

		map<Model*, double, ModelCmp> resProbs;
		map<Model*, double*, ModelCmp> resModes;

		if (method == 1 || method == 0 && df->exons.size() <= 4) 
		{
			SeppelExact* sepex = new SeppelExact(df, gene);
			sepex->calculate();
			resProbs = sepex->resProbs;
			resModes = sepex->resModes;
		} 
		else if (method == 3)
		{
			SeppelPrior* sepsm = new SeppelPrior(df, gene);
			sepsm->calculate();
			resProbs = sepsm->resProbs;
			resModes = sepsm->resModes;
		}
		else
		{
			Model* startmodel = new Model(initvars);
			SeppelSmart* sepsm = new SeppelSmart(df, gene, startmodel);
			sepsm->calculate();
			resProbs = sepsm->resProbs;
			resModes = sepsm->resModes;
		}

		vector<Variant*>* allpossvariants = df->allVariants(gene);
		printf("%s OK DONE", getmodelcode(allpossvariants, new Model(initvars)));

		// END OF CALCULATIONS
		
		Model* bestModel;
		double bestModelProb = -1;

		map<Model*, double, ModelCmp>::iterator mvi;
		for (mvi = resProbs.begin(); mvi != resProbs.end(); mvi++) 
		{
			if (mvi->second < minpp) 
			{
				resModes.erase(mvi->first);
				resProbs.erase(mvi);
			}
			if (mvi->second > bestModelProb)
			{
				bestModelProb = mvi->second;
				bestModel = mvi->first;
			}
		}

		// END OF FILTERING
		
		SEXP ans;
		PROTECT(ans= allocVector(VECSXP, 5));

		//Report posterior probabilities
		int i;
		int nrowpi=0, nx;
		if (selBest==0) { nx= resProbs.size(); } else { nx=1; }
		SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP,nx,2));
		double *resProbsR= REAL(VECTOR_ELT(ans,0));
		
		set<Variant*, VariantCmp>::const_iterator vi;
		map<Model*, double, ModelCmp>::const_iterator mi;
		if (selBest==0) {
			for (mi = resProbs.begin(), i=0; mi != resProbs.end(); mi++, i++) {
				Model* m = mi->first;
				resProbsR[i]= i;
				resProbsR[i+nx]= resProbs[m];
				nrowpi+= m->count();
			}
		} else {
			resProbsR[0]= 0;
			resProbsR[nx]= resProbs[bestModel];
			nrowpi+= bestModel->count();
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
