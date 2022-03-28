#include "mex.h"
#include <math.h>

/*******************************************************************/
/*                                                                 */
/*          Programme written by Tofik Amara, July 2002            */
/*       Modified by Steeve Zozor, January 2003, July 2013         */
/*                                                                 */
/*      Based on the algorithm proposed by Kaspar and Schuster     */
/*    Physical Review A, vol. 36, no. 2, pp. 842-848, July 1987    */
/*                                                                 */
/*                                                                 */
/*      Bug reports: steeve.zozor@gipsa-lab.grenoble-inp.fr        */
/*                                                                 */
/*******************************************************************/


/* function to interfacage Matlab and C . */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
  void LZComplexity(double*, int ,int, double*); /*  prototype of the function called in C  */

  double *signal;  /*  matrix of the signals   */
  int T;           /*  duration of the signal  */
  int d;           /*  dimension of the signal */
  double *c;       /*  Lempel-Ziv complexity   */


  /*  matlab usage: c = LZComplexity(signal)  */
  /*  check for proper number of arguments    */
  if(nrhs!=1){mexErrMsgTxt(" not proper input arguments ");}
  if(nlhs!=1){mexErrMsgTxt(" not proper output");}


 /*  Create matrix for the return arguments  */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);


  /*  Assign pointers to each input and output  */ 
  signal = mxGetPr(prhs[0]);
  c = (double *) mxGetPr(plhs[0]);


  /*  get the length of the matrix  */
  d = mxGetN(prhs[0]);      /*  column = components     */
  T = mxGetM(prhs[0]);      /*  line = times */

  /*  call of the routine to evaluate the complexity  */
  LZComplexity(signal,d,T,c);
}


/********************************************************************/

 /*  Function to compare test the equality between two vectors  */

int equalvectors(double *vector1, double *vector2, int d)

 /*  Remind that C put the signal in line order contrary to matlab */
 /*  i.e.    in matlab (i,j)  -->  (j-1) * nc + j                  */
 /*          in   C    (i,j)  -->  (i-1) * nc + j-1                */
 /*   !!! it supposed here that the vectors are line vectors !!!   */

{
  int j = 0;
  int eq = 1;/* answer of the equality test : before testng they are assumed equal */

  while ((j < d) && (eq))
  {
    if ( *(vector1+j) == *(vector2+j)) {j++;} // the current components are equal
    else {eq = 0;} // the current component are not equal => the vectors are not equal (unnecessary to compare the next components)
  }
  return(eq);
}

/********************************************************************/

 /*  Function to evaluate the joint Lempel-Ziv complexity  */

void LZComplexity(double *signal, int d, int T, double *c)

     /*  Inputs: signal  =   address of the signal        */
     /*            d     =   dimension of the signal      */
     /*            T     =   duration of the signal       */
     /*  Output:   c     =   complexity of the sequence   */

{
  int m, n, nmax;
  int l = 1;
  int eq;

  /*  initialization of c  */
  *c = 1.0;

 step1 : m = 0; n = 1 ; nmax = 1; // we assume the production is nmax = 1 (minimal possible), m is the pointer

 step2 : eq = equalvectors(signal+m+n-1 , signal+l+n-1 , d); // test of the equality between vectors (-1 because the indices stars from 0 in C)
  if(eq==1) // we have still a reproduction
    {
      n++; // still a reproduction
      if((l+n) > T) // we ended the sequence
	{
	  *c = (*c) + 1.0; // addition of the last word
	  goto output;
	}
      else {goto step2;} // to make the production that is not a reproduction
    }
  else // we face a prodution that is not a reproduction
    {
      if ( n > nmax){nmax = n;}  //the pointer is better since the production is larger -> up to date nmax
      m++; //test of the next point as a possible pointer 
      if(m == l) // all the possible pointer are tested
	{
	  *c = (*c) + 1.0; // the exhaustive production is made : c increases
	  //m = m+1;
	  l = l + nmax; //we skip to the symbol following the exhaustive word
	  if(l+1 > T){goto output;} // the sequence ended
	  else{goto step1;} // we start for seeking the next optimal production step
	}
      else //  we start again testing the next possible pointer
	{
	  n = 1;
	  goto step2;
	}
    }
 output : ;
}
