#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "hash.h"
void makeIslands(int **p_exons, int **p_islands, int nex, int tot, int **ex2tx, int **tx2ex, hash_t *ex2txP, hash_t *tx2exP, hash_t *ex2pos);
int connectTxs(int **p_exons, int **p_islands, int i, int allDone, int tot, int **ex2tx, int **tx2ex, hash_t *ex2txP, hash_t *tx2exP, hash_t *ex2posP);
int are_connected(int i, int j, int **p_exons, int **ex2tx, int **tx2ex, hash_t *ex2txP, hash_t *tx2exP);
int connectWithinTx(int **p_exons, int **p_islands, int i, int allDone, int **tx2ex, int **ex2tx, hash_t *ex2txP, hash_t *tx2exP, hash_t *ex2posP);

void makeIslands(int **p_exons, int **p_islands, int nex, int tot, int **ex2tx, int **tx2ex, hash_t *ex2txP, hash_t *tx2exP, hash_t *ex2posP){
  int i=0, islandCount=1, allDone=0;
  char id[30];
  for(i=0; i<nex; i++){
    //printf("%d out/n", i);
    if(p_islands[1][i]==0){
      //printf("%d\n", i);
      p_islands[1][i]=islandCount;
      allDone++;
      allDone = connectWithinTx(p_exons, p_islands, i, allDone, tx2ex, ex2tx, ex2txP, tx2exP, ex2posP);
      //printf("done with i=%d %d %d %d\n", i, allDone, nex, islandCount);
      //      allDone = connectTxs(p_exons, i, allDone, tot, tx2ex, ex2tx, ex2txP, tx2exP, ex2posP);
      if(allDone==tot) break;
      //printf("%d %d %d %d %d %d\n", i, allDone, p_exons[0][i], p_exons[1][i], p_exons[2][i], islandCount);
      islandCount++;
    }
  }
}

int connectWithinTx(int **p_exons, int **p_islands, int i, int allDone, int **tx2ex, int **ex2tx, hash_t *ex2txP, hash_t *tx2exP, hash_t *ex2posP){
  int l, j, m, exi, tx;
  char id[30];

  //get exon id
  //printf("inside %d %d\n", i, allDone);
  sprintf(id, "%d", p_islands[0][i]);
  // Find transcripts for this exon
  exi = hash_lookup(ex2txP, id);
  //  printf("%d %s %d %d %d\n", i, id, exi, ex2tx[exi][0], ex2tx[exi][1]);
  // For each of these transcripts
  for(j=1; j<ex2tx[exi][0]+1; j++){
    // Find all exons
    sprintf(id, "%d", ex2tx[exi][j]);
    tx = hash_lookup(tx2exP, id);
    //printf("before %d %d %s %d %d\n", i, exi, id, tx, tx2ex[tx][0]);
    //    return(0);
    for(m=1; m<tx2ex[tx][0]+1; m++){
      //Find position of exon in islands array
      sprintf(id, "%d", tx2ex[tx][m]);
      l = hash_lookup(ex2posP, id);
      //      printf("inside %d %d %d %s %d\n", l, m, tx2ex[tx][m], id, p_islands[1][l]);
      if(p_islands[1][l] == 0){
	p_islands[1][l] = p_islands[1][i];
	//printf("%d %d %d %s %d\n", i, m, tx, id, p_islands[1][l]);
	allDone++;
	allDone = connectWithinTx(p_exons, p_islands, l, allDone, tx2ex, ex2tx, ex2txP, tx2exP, ex2posP);
      }
    }
  }
  return(allDone);
}

int connectTxs(int **p_exons, int **p_islands, int i, int allDone, int tot, int **ex2tx, int **tx2ex, hash_t *ex2txP, hash_t *tx2exP, hash_t *ex2posP){
  int j;
  for(j=0; j<tot; j++) {
    if(p_exons[2][j]==0) {
      if(are_connected(i, j, p_exons, tx2ex, ex2tx, ex2txP, tx2exP)==1) {
        p_exons[2][j]=p_exons[2][i];
	allDone++;
	allDone = connectWithinTx(p_exons, p_islands, j, allDone, tx2ex, ex2tx, ex2txP, tx2exP, ex2posP);
        allDone = connectTxs(p_exons, p_islands, j, allDone, tot, tx2ex, ex2tx, ex2txP, tx2exP, ex2posP);
      }
    }
  }
  return(allDone);
}

int are_connected(int i, int j, int **p_exons, int **ex2tx, int **tx2ex, hash_t *ex2txP, hash_t *tx2exP){
  int txi, txj, k, m;
  char id[30];
  sprintf(id, "%d", p_exons[0][i]);
  txi = hash_lookup(ex2txP, id);
  sprintf(id, "%d", p_exons[0][j]);
  txj = hash_lookup(ex2txP, id);

  for(k=1; k < ex2tx[txi][0]+1; k++) {
    for(m=1; m < ex2tx[txj][0]+1; m++) {  
      if(ex2tx[txi][k] == ex2tx[txj][m]) {
	return(1);
      }
    }
  }

  return(0);
}

#define getDims(A) INTEGER( coerceVector( getAttrib(A, R_DimSymbol ) , INTSXP) )

SEXP makeGeneIslands(SEXP exons, SEXP isl, SEXP exisl, SEXP txs, SEXP totEx, SEXP nexR){
  //ex2tx is a list of exons with the transcripts they belong to
  //exons is a named vector with exons where names will be the islands

  int **p_exons, totIslands, totExo, l, i, **p_islands, nex;

  p_exons = malloc(3 * sizeof(int *));
  p_islands = malloc(3 * sizeof(int *));
  PROTECT (exons = coerceVector(exons, INTSXP));
  PROTECT (isl = coerceVector(isl, INTSXP)); 
  PROTECT (exisl = coerceVector(exisl, INTSXP));
  PROTECT (txs = coerceVector(txs, INTSXP));
  PROTECT (totEx = coerceVector(totEx, INTSXP));
  p_exons[0] = INTEGER(exons);
  //p_exons[2] = INTEGER(islands);
  p_exons[1] = INTEGER(txs);
  p_islands[1] = INTEGER(isl);
  p_islands[0] = INTEGER(exisl);
  
  totExo = INTEGER(totEx)[0];
  nex = INTEGER(nexR)[0];

  hash_t ex2txH, *ex2txP, tx2exH, *tx2exP, ex2posH, *ex2posP;
  ex2txP = &ex2txH;
  tx2exP = &tx2exH;
  ex2posP = &ex2posH;
  hash_init(ex2txP, totExo*2);
  hash_init(tx2exP, totExo*2);
  hash_init(ex2posP, nex*2);

  int txid=1;
  int exid=1;
  char id[30];
  for(i=0; i<totExo; i++) {
    sprintf(id, "%d", p_exons[1][i]);
    l=hash_lookup(tx2exP,  id);
    if(l==HASH_FAIL) { hash_insert(tx2exP, id, txid); txid++; }
    sprintf(id, "%d", p_exons[0][i]);
    l=hash_lookup(ex2txP, id);
    if(l==HASH_FAIL) { hash_insert(ex2txP, id, exid); exid++; }
   }

  for(i=0; i<nex; i++){
    sprintf(id, "%d", p_islands[0][i]);
    hash_insert(ex2posP, id, i);
    //    printf("%d %d %d %s\n", nex, i, hash_lookup(ex2posP, id), id);
  }

  //printf("Hash ready\n");

  int **ex2tx, **tx2ex;
  ex2tx = malloc((totExo*2) * sizeof(int *));
  tx2ex = malloc((totExo*2) * sizeof(int *));
  for(i=0; i<totExo+1; i++) {
    ex2tx[i] = malloc(500 * sizeof(int));
    ex2tx[i][0]=0;
    tx2ex[i] = malloc(500 * sizeof(int));
    tx2ex[i][0]=0;
  }

  //printf("Arrays ready\n");

   for(i=0; i<totExo; i++){
     sprintf(id, "%d", p_exons[0][i]);
     l=hash_lookup(ex2txP, id);
     ex2tx[l][ex2tx[l][0]+1] = p_exons[1][i];
     ex2tx[l][0]++;
     sprintf(id, "%d", p_exons[1][i]);
     l=hash_lookup(tx2exP, id);
     tx2ex[l][tx2ex[l][0]+1] = p_exons[0][i];
     tx2ex[l][0]++;
 }


   //printf("Making islands\n");
   makeIslands(p_exons, p_islands, nex, totExo, ex2tx, tx2ex, ex2txP, tx2exP, ex2posP);
  
  UNPROTECT(5);
  return(isl);

 }
