#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "hash.h"
#include "header.h"	

void addExon2Frag(int exon, int start, int rid, int pos, path_t *frags, int first);

void addExon2Frag(int exon, int start, int rid, int pos, path_t *frags, int first){
  int *newexons;
  int *newstarts;
  int *newrids;

  if(first==1){
    int exonNum=11;
    frags[pos].nexon=0;
    frags[pos].exons=malloc(exonNum * sizeof(int));
    frags[pos].starts=malloc(exonNum * sizeof(int));
    frags[pos].rids=malloc(exonNum * sizeof(int));	
  }
  else {
  if((frags[pos].nexon % 10)==0) {
    printf("Reallocating %d %d\n", frags[pos].nexon, frags[pos].nexon+10);
        newexons = realloc(frags[pos].exons, (frags[pos].nexon + 10) * sizeof(int));
        if(newexons != NULL) frags[pos].exons=newexons;
	else puts("Error reallocating memory\n");
	newstarts=realloc(frags[pos].starts, (frags[pos].nexon + 10) * sizeof(int));
        if(newstarts != NULL) frags[pos].starts=newstarts;
        else puts("Error reallocating memory\n"); 
	newrids=realloc(frags[pos].rids, (frags[pos].nexon + 10) * sizeof(int));
        if( newrids != NULL) frags[pos].rids = newrids;
	else puts("Error reallocating memory\n");
    }
  }
  frags[pos].exons[frags[pos].nexon]=exon;
  frags[pos].starts[frags[pos].nexon]=start;
  frags[pos].rids[frags[pos].nexon]=rid;
  frags[pos].nexon++;		
}

int buildFrags(hash_t *fragsHashPtr, SEXP reid, int *rid, int *start, int *exon, int nreads, path_t **frags){ 
  int l=0, totF=0, i, totFnow=1000;
  path_t *newF;

  *frags = malloc((totFnow + 1) * sizeof(path_t));

  for(i=0; i<nreads; i++) {
    l=hash_lookup(fragsHashPtr, CHAR(STRING_ELT(reid,i)));
    if(l!=HASH_FAIL) {
      addExon2Frag(exon[i], start[i], rid[i], l, *frags, 2);
    }
    else {
      hash_insert(fragsHashPtr,  CHAR(STRING_ELT(reid,i)), totF);
      addExon2Frag(exon[i], start[i], rid[i], totF, *frags, 1);
      totF++;
      if(totF==(totFnow-1)){
	totFnow=totFnow*2;
	newF= realloc(*frags, totFnow * sizeof(path_t));
	if(newF != NULL) *frags=newF;
	else puts("Error reallocating memory frags\n");
      } 
    }
  }

  return(totF);
}



