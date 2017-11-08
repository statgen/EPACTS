#include <R.h>
#include <Rmath.h>
#include <math.h>


/*********************************************************************


************************************************************************/
double SL_runif_double(){

 double val ;
 val = unif_rand();
 return val ;
 
}

/* Generate 0 to max -1 */
int SL_runif_int(int max){

 int val;
 val = (int)(SL_runif_double() * max);
 return val ;
 
}

void  SL_setseed(int * seed){
  
  /*Rprintf("Seed Number %d\n",*seed);*/
  srand ( *seed );

}


/*
	n : number of samples
	k : number of samples to be sampled
	y : output
	x : buffer
*/

void SL_GetSample(int n, int k, int *y, int *x){
    int i, j, temp;
    
    for (i = 0; i < n; i++){
        x[i] = i;
    }
    for (i = 0; i < k; i++) {
        
        temp = rand();
        //j = n * temp / (RAND_MAX+1);
        j = temp % n;
        y[i] = x[j] ;
        x[j] = x[--n];
    }
}

void SL_Binary_Boot1(int n, int ncase, double * pcase, int * buf1, int * buf2, int * z_one, int *err){

	int i, i1,j,k, n1, ncase1;
	double temp;

	SL_GetSample(n, n, buf1, buf2);
	ncase1 = 0;
	n1 = n;
	for(k=0;k<20;k++){
		i1=0;
		for(i=0;i<n1;i++){
			j=  buf1[i];
			temp = SL_runif_double();
			if(temp <= pcase[j]){
				z_one[j]=1;
				ncase1++;
			} else {
				buf2[i1] = j;
				i1++;
			}
			if(ncase1 == ncase)
				break;
		} 
		
		if(ncase1 == ncase){
			break;
		} else if(ncase1 > ncase) {
			*err = -1;
			return;
		} else {
			n1 = n - ncase1;
			memcpy(buf1, buf2, sizeof(int) * n1);		
		}
	}
	
	if(ncase != ncase1){
		*err = -1;
		return;
	}
	
	*err=1;
	return;
	/*Rprintf(" %d %d %d\n", ncase, ncase1, k);*/

}

void SL_Binary_Boot(int * pn, int * pm, int *pncase, double * pcase, int * buf1, int * buf2, int *Z, int *err ){

	int i, n, m, ncase, idx;
	int * z_one;
	n = *pn;
	m = *pm;
	ncase = *pncase;
	
	GetRNGstate();
	
	for(i=0;i<m;i++){
	
		idx = i * n ;
		z_one = &(Z[idx]);
		SL_Binary_Boot1(n, ncase, pcase, buf1, buf2, z_one, err);
		if(*err == -1){
			PutRNGstate();
			return;
		}
	}
	
	PutRNGstate();
	return;
}


void Test1(int * pn, int * pm, int *pncase, double * pcase, int * buf1, int * buf2, int *buf3){

	int i, n, m, ncase, idx;
	int * z_one;
	n = *pn;
	m = *pm;
	ncase = *pncase;
	
	for(i=0;i<m;i++){
	
		idx = i * n ;
		/*SL_Binary_Boot1(n, ncase, pcase, buf1, buf2, buf3);*/
		SL_GetSample(n, n, buf1, buf2);
	}
}

