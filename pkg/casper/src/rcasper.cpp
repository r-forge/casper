#include <sstream>
#include <set>
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
			queue->erase(queue->find(si->first));
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
		if (initvars->count(nv) == 0) 
		{
			// check if the new variant can explain the fragment
			map<Fragment*, double> probs = df->probabilities(nv);
			if (probs.count(frag) > 0)
			{
				initvars->insert(nv);

				// delete all fragments that this variant can explain
				map<Fragment*, double>::iterator si;
				for (si = probs.begin(); si != probs.end(); si++) 
				{
					queue->erase(queue->find(si->first));
				}
			}
			else
			{
				// this fragment cant be explained
				discarded++;
			}
		}
	}

	return discarded;
}
DataFrame* importDataFrame(SEXP exonsR, SEXP exonwidthR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, int geneid) 
{
	int nexons=Rf_length(exonsR), nfraglen=Rf_length(fraglenR), readLength=INTEGER(readLengthR)[0];
	int *exons=INTEGER(exonsR), *exonwidth=INTEGER(exonwidthR), *lenvals=INTEGER(lenvalsR);
	double *fraglen= REAL(fraglenR);

	//Define fragment length/start distributions and create DataFrame
	DiscreteDF* fraglen_dist = new DiscreteDF(fraglen, lenvals, nfraglen);

	lencdf= Rf_length(fragstaR);
	startcdf= REAL(fragstaR);
	//::fun_fragsta = fragstaR;

	DataFrame* df = new DataFrame(fraglen_dist, cumu_fragsta);
	df->frag_readlen= readLength;

	//Add exons & gene to DataFrame
	Gene* g= new Gene(geneid);
	for (int i=0; i<nexons; i++) {
		Exon *ex= new Exon(exons[i], exonwidth[i]);
		df->addExon(ex);
		g->addExon(ex);
	}
	df->addGene(g);

	//Add exon path counts
	int np = LENGTH(pathCountsR);
	SEXP pnames = getAttrib(pathCountsR, R_NamesSymbol);
	int* pathCounts = INTEGER(pathCountsR);  

	for (int i=0; i<np; i++) {
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
set<Variant*, VariantCmp>* importTranscripts(DataFrame* df, SEXP transcriptsR, int geneid)
{
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

extern "C"
{
	SEXP calcDenovoMultiple(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP verboseR) {
		//De novo isoform discovery and estimate expression for multiple genes. Calls calcDenovoSingle repeadtedly
		int i, ngenes=LENGTH(geneidR);
		SEXP ansMultiple, ansSingle;

		PROTECT(ansMultiple= allocVector(VECSXP, ngenes));

		for (i=0; i<ngenes; i++) {
			ansSingle= calcDenovoSingle(VECTOR_ELT(exonsR,i), VECTOR_ELT(exonwidthR,i), VECTOR_ELT(transcriptsR,i), VECTOR_ELT(geneidR,i), VECTOR_ELT(pathCountsR,i), fragstaR, fraglenR, lenvalsR, readLengthR, nvarPriorR, nexonPriorR, priorqR, minppR, selectBest, verboseR);
			SET_VECTOR_ELT(ansMultiple,i,ansSingle);
		}

		UNPROTECT(1);
		return ansMultiple;
	}

	SEXP calcDenovoSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP verboseR)
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
		// - verbose: set verbose != 0 to print warnings & progress info

		int geneid=INTEGER(geneidR)[0];
		int selBest=INTEGER(selectBest)[0];
		double minpp=REAL(minppR)[0], priorq= REAL(priorqR)[0];

		PROTECT(geneidR= allocVector(INTEGERSXP, LENGTH(transcriptsR)));
		for (int i=0; i< LENGTH(transcriptsR); i++) INTEGER(geneidR)[i]= geneid;

		DataFrame* df = importDataFrame(exonsR, exonwidthR, pathCountsR, fragstaR, fraglenR, lenvalsR, readLengthR, geneidR);
		set<Variant*, VariantCmp>* initvars = importTranscripts(df, transcriptsR, geneidR);

		UNPROTECT(1);		

		// EMD OF INPUT READING

		fixUnexplFrags(initvars, df, df->genes[geneid]);
		Model* model = new Model(new vector<Variant*>(initvars->begin(), initvars->end()));
		Casper* casp = new Casper(model, df);
		casp->priorq= priorq;

		SeppelExact* sepex;
		SeppelSmart* sepsm;

		map<Model*, double, ModelCmp> posprob;
		map<Model*, double*, ModelCmp> mode;
		Model* bestModel;
		double* modeBest;

		if (df->exons.size() <= 4) {
			sepex = new SeppelExact(df, df->genes[geneid]);
			sepex->calculate();
			sepex->rmModels(minpp);
			posprob = sepex->posprob;
			mode = sepex->mode;
			bestModel= sepex->bestModel;
			modeBest= sepex->modeBest;
		} else {
			SeppelSmart* sepsm = new SeppelSmart(df, df->genes[geneid]);
			posprob = sepsm->calculate(model);
		}

		// END OF CALCULATIONS
		
		SEXP ans;
		PROTECT(ans= allocVector(VECSXP, 5));

		//Report posterior probabilities
		int i;
		int nrowpi=0, nx;
		if (selBest==0) { nx= posprob.size(); } else { nx=1; }
		SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP,nx,2));
		double *posprobR= REAL(VECTOR_ELT(ans,0));

		map<Model*, double, ModelCmp>::const_iterator mi;
		if (selBest==0) {
			for (mi = posprob.begin(), i=0; mi != posprob.end(); mi++, i++) {
				Model* m = mi->first;
				posprobR[i]= i;
				posprobR[i+nx]= posprob[m];
				nrowpi+= m->count();
			}
		} else {
			posprobR[0]= 0;
			posprobR[nx]= posprob[bestModel];
			nrowpi+= bestModel->count();
		}

		//Report estimated expression
		SET_VECTOR_ELT(ans, 1, allocMatrix(REALSXP,nrowpi,2));  //stores model id and estimated expression
		SET_VECTOR_ELT(ans, 2, allocVector(STRSXP,nrowpi)); //stores variant names
		double *expr= REAL(VECTOR_ELT(ans,1));
		SEXP varnamesR= VECTOR_ELT(ans,2);
		map<int, Variant*> variantHash;
		variantHash= casp->model->getVarHash();

		int rowid=0, nrowexons=0;
		set<Variant*> allvariants;
		if (selBest==0) {
			for (mi = posprob.begin(), i=0; mi != posprob.end(); mi++, i++) {
				Model* m = mi->first;
				for (int j=0; j< m->count(); j++) {
					expr[rowid]= i;  //model id
					Variant* v = m->get(j);
					int varidx= m->indexOf(v);
					expr[rowid+nrowpi]= sepex->mode[m][varidx]; //estimated expression
					if (variantHash.count(v->hashcode)>0) v->name= (variantHash[v->hashcode])->name;  //respect initial variant names
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
				expr[rowid+nrowpi]= sepex->mode[m][varidx]; //estimated expression
				if (variantHash.count(v->hashcode)>0) v->name= (variantHash[v->hashcode])->name;  //respect initial variant names
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

		set<Variant*>::iterator vi;
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
