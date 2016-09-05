#include "fstrack.h"
/*

  sort array by order, returns index

  give input index shiftes (ie.    x-1)
  output will be 0...n-1
  

  based on numerical recipes

  $Id: indexx.c,v 1.3 2004/04/19 18:41:12 becker Exp $
  
*/


#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 5000

void indexx(int n,COMP_PRECISION *arr,int *indx)
{
  int i,indxt,ir=n,itemp,j,k,l=1;
  int jstack=0,*istack;
  COMP_PRECISION a;
  //
  istack=(int *)malloc(sizeof(int)*(NSTACK+1));
  if(!istack)MEMERROR("indexx");
  //
  for (j=1;j<=n;j++) 
    indx[j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	indxt=indx[j];
	a=arr[indxt];
	for (i=j-1;i>=1;i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1]=indx[i];
	}
	indx[i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(indx[k],indx[l+1]);
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	SWAP(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l],indx[ir]);
      }
      if (arr[indx[l+1]] > arr[indx[l]]) {
	SWAP(indx[l+1],indx[l]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l];
      a=arr[indxt];
      for (;;) {
	do i++; while (arr[indx[i]] < a);
	do j--; while (arr[indx[j]] > a);
	if (j < i) break;
	SWAP(indx[i],indx[j]);
      }
      indx[l]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack > NSTACK) {
	fprintf(stderr,"indexx: NSTACK too small");
	exit(-1);
      }
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free(istack);
  // go back to normal numbering
  for(i=1;i<=n;i++)
    indx[i]--;
}
#undef M
#undef NSTACK
#undef SWAP

