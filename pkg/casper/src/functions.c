#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hash.h"
#include "header.h"

//********** Functions ***********

int *procCigar(char *cigar){
	char *pch, *tab="MN";
	int *cigs;
	cigs=malloc(50 * sizeof(char));
	pch=strtok(cigar, tab);
	cigs[0]=0;
	while(pch!=NULL){
		sscanf(pch, "%d", &cigs[cigs[0]+1]);
		pch=strtok(NULL, tab);
		cigs[0]++;
	}
	return(cigs);
}


void addRead2Frag(char *holder, int totF, read_t *frags, int read){
	char *pch, tab[]="\t";
	if(read==1){
		pch=strtok(holder, tab);
        	frags[totF].qname = malloc(100 * sizeof(char));
		strcpy(frags[totF].qname, pch);
		pch=strtok(NULL, tab);
		sscanf(pch, "%d", &frags[totF].flag_1);
		pch=strtok(NULL, tab);
		strcpy(frags[totF].chr_1, pch);
		pch=strtok(NULL, tab);
		sscanf(pch, "%d", &frags[totF].st_1);
		pch=strtok(NULL, tab);
		pch=strtok(NULL, tab);
		strcpy(frags[totF].cigar_1, pch);
		frags[totF].nreads=1;
	} else {    
		pch=strtok(holder, tab);
		pch=strtok(NULL, tab);
		sscanf(pch, "%d", &frags[totF].flag_2);
		pch=strtok(NULL, tab);
		strcpy(frags[totF].chr_2, pch);
		pch=strtok(NULL, tab);
		sscanf(pch, "%d", &frags[totF].st_2);
		pch=strtok(NULL, tab);
		pch=strtok(NULL, tab);
		strcpy(frags[totF].cigar_2, pch);
		frags[totF].nreads=2;
	}   
}
