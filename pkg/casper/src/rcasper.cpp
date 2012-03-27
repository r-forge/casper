#include <sstream>
#include <set>
#include <R.h>
#include <Rinternals.h>
#include "stdafx.h"
#include "casper.h"
#include "rcasper.h"
#include "seppelsmart.h"

double cumu_fragsta(double x) {
  if (x<=0) return 0;
  if (x>=1) return 1;
  int idx= (int) (x*lencdf);
  double y1= startcdf[idx], x1= (double) idx / (lencdf-1);
  idx++;
  double y2= startcdf[idx], x2= (double) idx / (lencdf-1);
  return y1 + (x-x1) * (y2-y1)/(x2-x1);
}

/*double cumu_fragsta(double x) {
  SEXP sval, R_fcall;
  PROTECT(sval = allocVector(REALSXP, 1));
  REAL(sval)[0]= x; //store x in sval

  PROTECT(R_fcall = lang2(fun_fragsta, R_NilValue)); //create executable pairlist
  SETCADR(R_fcall, sval);  //prepare call for value in sval

  SEXP funval = eval(R_fcall, R_NilValue);
  double res = REAL(funval)[0];
  UNPROTECT(2);

  return res;
}
*/


extern "C"
{

  //Compute expression for known variants.
  SEXP calc(SEXP exons, SEXP transcripts, SEXP pathCounts, SEXP fragsta, SEXP fraglen, SEXP lenvals, SEXP readLength) 
	{
	  lencdf= Rf_length(fragsta);
          startcdf= REAL(fragsta);
	  //::fun_fragsta = fragsta;
		
		// DiscreteDF
		SEXP lims;
		PROTECT(lims = getAttrib(fraglen, R_DimSymbol));
		int ln = INTEGER(lims)[0];
		double* ld = REAL(fraglen);
		UNPROTECT(1);
		int* lv= INTEGER(lenvals);

		DiscreteDF* lenfun = new DiscreteDF(ld, lv, ln);

		// Main
		DataFrame* df = new DataFrame(lenfun, cumu_fragsta);
                df->frag_readlen= INTEGER(readLength)[0];
		
		// Exons
		SEXP dims;
		PROTECT(dims = getAttrib(exons, R_DimSymbol));
		int ne = INTEGER(dims)[0];
		UNPROTECT(1);

		int* me = INTEGER(exons);

		REprintf("Reading %i exons..\n", ne);
		for (int i = 0; i < ne; i++)
		{
			int id = me[0*ne+i];
			int length = me[1*ne+i];

			df->addExon(new Exon(id, length));
		}

		// PathCounts
		int np = LENGTH(pathCounts);

		SEXP pnames = getAttrib(pathCounts, R_NamesSymbol);
		int* pvalus = INTEGER(pathCounts);

		REprintf("Reading %i pathCounts...\n", np);
		for (int i = 0; i < np; i++)
		{
			const char* pname = CHAR(STRING_ELT(pnames, i));
			int count = pvalus[i];

			char* left = new char[strlen(pname) + 1];
			strcpy(left, pname);
			char* mid = strchr(left, '-');
			if (mid == NULL)
			{
				continue;
			}
			mid[0] = '\0';
			char* right = mid+1;

			int leftc = 0;
			for (int l = strlen(left)-1; l >= 0; l--)
			{
				if (left[l] == '.')
				{
					leftc++;
				}
			}
			int rightc = 0;
			for (int r = strlen(right)-1; r >= 0; r--)
			{
				if (right[r] == '.')
				{
					rightc++;
				}
			}

			left = left+1;
			right[strlen(right)-1] = '\0';

			Fragment* f = new Fragment(leftc, rightc, count);

			char* item;
			item = strtok(left, ".");
			for (int i = 0; item != NULL; i++)
			{
				int eid = atoi(item);
				f->left[i] = eid;

				item = strtok(NULL, ".");
			}

			item = strtok(right, ".");
			for (int i = 0; item != NULL; i++)
			{
				int eid = atoi(item);
				f->right[i] = eid;

				item = strtok(NULL, ".");
			}

			df->addData(f);
		}

		// Genes
		// TODO: add genes here like this of positive strand
		Gene* g = new Gene(1);
		g->addExon(df->exons[1]);
		g->addExon(df->exons[2]);
		g->addExon(df->exons[3]);
		df->addGene(g);

		// Variants
		int nt = LENGTH(transcripts);
		vector<Variant*>* vm = new vector<Variant*>();
		SEXP tnames = getAttrib(transcripts, R_NamesSymbol);

		REprintf("Reading %i transcripts...\n", nt);
		for (int i = 0; i < nt; i++)
		{
			int tid = atoi(CHAR(STRING_ELT(tnames, i)));

			SEXP trow = VECTOR_ELT(transcripts, i);
			int ntsub = LENGTH(trow);
			int* tvals = INTEGER(trow);

			vector<Exon*>* el = new vector<Exon*>();

			for (int s = 0; s < ntsub; s++)
			{
				int eid = tvals[s];

				Exon* ex = df->exons[eid];
				el->push_back(ex);
			}

			// TODO: read strand, positive or negative?
			bool strand = true;

			// TODO: read gene id of variant here
			int gid = -1;
			Gene* gene = df->genes[gid];

			//Variant* v = new Variant(gene, strand, el);
			Variant* v = new Variant(gene, el);
			v->id = tid;

			vm->push_back(v);
		}
		Model* mdl = new Model(vm);

		
/*
		// Debug
		map<Variant*, map<Fragment*, double> >::const_iterator mi;
		for (mi = c->memvprobs.begin(); mi != c->memvprobs.end(); mi++)
		{
			map<Fragment*, double>::const_iterator ni;
			for (ni = mi->second.begin(); ni != mi->second.end(); ni++) 
			{
				if (ni->first->leftc == 1 && ni->first->rightc == 2)
				{
				//	REprintf("%i\t%i %i %i\t%f\n", mi->first->id, ni->first->left[0], ni->first->right[0], ni->first->right[1], ni->second);
				}
			}
		}
//		REprintf("\n");
		
		map<int, Exon*>::const_iterator ei;
		for (ei = c->exons.begin(); ei != c->exons.end(); ei++)
		{
			Exon* e = ei->second;
//			REprintf("%i\t%i\n", e->id, e->length);
		}
//		REprintf("\n");
		
		list<Fragment*>::const_iterator fi;
		for (fi = c->fragments.begin(); fi != c->fragments.end(); fi++)
		{
			Fragment* f = *fi;
		//	REprintf("%i\t%i\t%i\n", f->leftc, f->rightc, f->count);
			for (int l = 0; l < f->leftc; l++)
			{
		//		REprintf("%i\n", f->left[l]);
			}
			for (int r = 0; r < f->rightc; r++)
			{
		//		REprintf("%i\n", f->right[r]);
			}
		}
	//	REprintf("\n");

		vector<Variant*>::const_iterator vi;
		for (vi = c->variants.begin(); vi != c->variants.end(); vi++)
		{
			Variant* v = *vi;
	//		REprintf("%i\t%i\n", v->id, v->exonCount);
			for (int e = 0; e < v->exonCount; e++)
			{
	//			REprintf("%i\n", v->exons[e]);
			}
		}
	//	REprintf("\n");
*/

		// Output
		Casper* c = new Casper(mdl, df);
		int vc = c->model->count();
		double* em = c->calculateMode();

		SEXP Rc;
		PROTECT(Rc = allocVector(REALSXP, vc));
		double* res = REAL(Rc);
		for (int i = 0; i < vc; i++)
		{
			res[i] = em[i];
		}
		UNPROTECT(1);

		return(Rc);
	}


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

  SEXP calcDenovoSingle(SEXP exonsR, SEXP exonwidthR, SEXP transcriptsR, SEXP geneidR, SEXP pathCountsR, SEXP fragstaR, SEXP fraglenR, SEXP lenvalsR, SEXP readLengthR, SEXP nvarPriorR, SEXP nexonPriorR, SEXP priorqR, SEXP minppR, SEXP selectBest, SEXP verboseR) {
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

    int i, geneid=INTEGER(geneidR)[0], nexons=Rf_length(exonsR), nfraglen=Rf_length(fraglenR), readLength=INTEGER(readLengthR)[0], verbose=INTEGER(verboseR)[0], selBest=INTEGER(selectBest)[0];
    int *exons=INTEGER(exonsR), *exonwidth=INTEGER(exonwidthR), *lenvals=INTEGER(lenvalsR);
    double *fraglen= REAL(fraglenR), minpp=REAL(minppR)[0], priorq= REAL(priorqR)[0]; 
    SEXP ans;

    Casper* casp = initCasper(exons, exonwidth, transcriptsR, geneid, nexons, pathCountsR, fraglen, lenvals, nfraglen, readLength, fragstaR, priorq, verbose);
    SeppelExact* sepex;
    SeppelSmart* sepsm;

    map<Model*, double, ModelCmp> posprob;
    map<Model*, double*, ModelCmp> mode;
    Model* bestModel;
    double* modeBest;

    if (nexons<=4) {
      sepex = new SeppelExact(casp->frame, casp->frame->genes[geneid]);
      sepex->calculate();
      sepex->rmModels(minpp);
      posprob = sepex->posprob;
      mode = sepex->mode;
      bestModel= sepex->bestModel;
      modeBest= sepex->modeBest;
    } else {
      SeppelSmart* sepsm = new SeppelSmart(casp->frame, casp->frame->genes[geneid]);
      posprob = sepsm->calculate(casp->model);
    }

    PROTECT(ans= allocVector(VECSXP, 5));

    //Report posterior probabilities
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


Casper* initCasper(int *exons, int *exonwidth, SEXP transcriptsR, int geneid, int nexons, SEXP pathCountsR, double *fraglen, int *lenvals, int nfraglen, int readLength, SEXP fragstaR, double priorq, int verbose) {

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

  //Set initial variants
  int nt = LENGTH(transcriptsR);
  vector<Variant*>* initvars = new vector<Variant*>();
  map<int, Variant*> hash2var;
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

    hash2var[v->hashcode]= v;
    initvars->push_back(v);
  }

  //Initialize variant * path probabilities (stored in df->cache)
  vector<Variant*>::iterator vi;
  map <Fragment*, double> varprobs;
  for (vi = initvars->begin(); vi != initvars->end(); vi++) {
    Variant *v= *vi;
    varprobs = df->probabilities(v);
  }

  //Add variants to initial set so that all observed exon paths have positive probability
  int deletedPath=0, newct=nt, ndelpaths=0, posprob;
  list<Fragment*>::iterator fi;
  for (fi = df->data.begin(); fi != df->data.end(); fi++) {
    Fragment* f = *fi;
    for (vi = initvars->begin(), posprob=0; vi != initvars->end() && posprob==0; vi++) {
      varprobs = df->probabilities(*vi);
      if (varprobs.count(f)>0 && varprobs[f]>0) { posprob=1; }
    }
    if (posprob==0) {
      Variant* v = path2Variant(df, f, gene);
      v->id= newct;
      if (hash2var.count(v->hashcode)==0) {  //if the new variant not in initvars
	newct++;
        initvars->push_back(v);
        varprobs = df->probabilities(v);
        deletedPath= varprobs.count(f)==0 || !(varprobs[f]>0);
	hash2var[v->hashcode]= v;
      } else {
	deletedPath= true;
	delete v;
      }
      if (deletedPath) {
        fi= df->data.erase(fi);
        ndelpaths++;
	deletedPath= 0;
	fi--;
      }
    }
  }
  if (verbose!=0 && ndelpaths>0) Rprintf("%d exon paths could not be explained by defining new variants. Paths removed\n", ndelpaths); 

  Model* model = new Model(initvars);
  Casper* casp = new Casper(model, df);
  casp->priorq= priorq;

  return casp;
}


Variant* path2Variant(DataFrame *df, Fragment* f, Gene *gene) {
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
