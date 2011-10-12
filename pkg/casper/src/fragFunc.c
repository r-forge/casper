#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "hash.h"
#include "header.h"	

void addExon2Frag(int exon, int start, int rid, int pos, path_t *frags, int first);

void addExon2Frag(int exon, int start, int rid, int pos, path_t *frags, int first){
  if(first==1){
    int exonNum=300;
    frags[pos].nexon=0;
    frags[pos].exons=malloc(exonNum * sizeof(int));
    frags[pos].starts=malloc(exonNum * sizeof(int));
    frags[pos].rids=malloc(exonNum * sizeof(int));	
  }
  frags[pos].exons[frags[pos].nexon]=exon;
  frags[pos].starts[frags[pos].nexon]=start;
  frags[pos].rids[frags[pos].nexon]=rid;
  frags[pos].nexon++;		
}

int buildFrags(hash_t *fragsHashPtr, SEXP reid, int *rid, int *start, int *exon, int nreads, path_t *frags){ 
  int l=0, totF=0, i;

  for(i=0; i<nreads; i++) {
    l=hash_lookup(fragsHashPtr, CHAR(STRING_ELT(reid,i)));
    if(l!=HASH_FAIL) {
      addExon2Frag(exon[i], start[i], rid[i], l, frags, 2);
    }
    else {
      hash_insert(fragsHashPtr,  CHAR(STRING_ELT(reid,i)), totF);
      addExon2Frag(exon[i], start[i], rid[i], totF, frags, 1);
      totF++;
    }
  }
 

  return(totF);
}



