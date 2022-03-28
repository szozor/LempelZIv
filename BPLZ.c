#define char16_t uint16_T
#include "mex.h"
#include <math.h>

/*******************************************************************/
/*                                                                 */
/*          Programme written Steeve Zozor, July 2013              */
/*                                                                 */
/*       Based on the algorithms proposed by                       */
/*         - Kaspar and Schuster, Phys. Rev. A 36(2): 842, 1987    */
/*         - Bandt and Pompe,                                      */
/*                                                                 */
/*      Bug reports: steeve.zozor@gipsa-lab.grenoble-inp.fr        */
/*                                                                 */
/*******************************************************************/


/* function to interfacage Matlab and C . */


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
  void BPLZ(double*, int ,int, double, int, double*, int*); /*  prototype of the function called in C  */

  double *signal;  /*   matrix of the signals             */
                   /*                                     */
  int d;           /*   embedding dimension               */
  int tau;         /*   delay                             */
  int Tb;          /*   duration of the (bloc of) signal  */
  double conf;     /*   confidence interval               */
                   /*                                     */
  double *c;       /*   LZ complexity on BP vectors       */
  int *Rt;         /*   Sequence of ranks                 */


  /*  matlab usage: [c Rt] = BPLZ(signal,d,tau,conf)  */
  /*  check for proper number of arguments    */
  if((nrhs!=3)*(nrhs!=4)){mexErrMsgTxt(" not proper input arguments ");}
  if((nlhs!=1)*(nlhs!=2)){mexErrMsgTxt(" not proper output");}

 /*  Create matrix for the first return argument  */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);

  /*  Assign pointers to each input and output  */ 
  signal = mxGetPr(prhs[0]);
  d = mxGetScalar(prhs[1]);
  tau = mxGetScalar(prhs[2]);
  c = (double *) mxGetPr(plhs[0]);
  if(nrhs==4){conf = mxGetScalar(prhs[3]);}else{conf = 0.0;}


  /*  get the length of the signal  */
  Tb = mxGetM(prhs[0]);      /*  line = times */
  if(Tb==1){mexErrMsgTxt(" the vector must be one column");}

  if(nlhs==2)
  {
    plhs[1] = mxCreateNumericMatrix(d,Tb-(d-1)*tau,mxINT32_CLASS,mxREAL);
    Rt = (int *) mxGetPr(plhs[1]);
  }
  else{Rt = (int*)malloc(d*(Tb-(d-1)*tau)*sizeof(int));}

  /*  call of the routine to evaluate the complexity  */
  BPLZ(signal,d,tau,conf,Tb,c,Rt);
  if(nlhs==1){free(Rt);}
}


/********************************************************************/

 /*  Function to compare test the equality between two vectors  */

int EqualVectors(int *Vector1, int *Vector2, int d)


     /*  Inputs: Vector1, Vector2  =   address of the vectors to be compared    */
     /*                      d     =   dimension of the vectors                 */
     /*  Output:            eq     =   1 if the vectors are equal, 0 otherwise  */

 /*  Remind that C put the signal in line order contrary to matlab */
 /*  i.e.    in matlab (i,j)  -->  (j-1) * nc + j                  */
 /*          in   C    (i,j)  -->  (i-1) * nc + j-1                */
 /*   !!! it supposed here that the vectors are line vectors !!!   */

{
  int j = 0;
  int eq = 1;/* answer of the equality test : before tested they are assumed equal */

  while ((j < d) && (eq))
  {
    if ( *(Vector1+j) == *(Vector2+j)) {j++;} /* current components equal */
    else {eq = 0;} /* current component not equal */
                   /* => vectors not equal (unnecessary to compare the next components) */
  }
  return(eq);
}

/*****************************************************************/

 /*  Function to built the rank vector (brut force)  */

void Rank(double *Vector, int d, double conf, int *Rk)

     /*  Inputs:   Vector   =    address of the vector to rank                                  */
     /*              d      =    dimension of the vector                                        */
     /*            conf     =    confidence interval (test of the equality)                     */
     /*  Output:    Rk    =    address of the vector of the ranks of the compoents of Vector  */
{
  int i, j;
  double Vi;

  for (i = 0; i < d; i++) /* current component, to compare to that (j) already sorted */
    {
      *(Rk+i) = 0; /* initially the rank is assumed 0 */
      for (j = 0; j < i; j++)
	{
	  if (*(Vector+i) < *(Vector+j)-conf) {(*(Rk+j)) = (*(Rk+j))+1;} /* lower => increase rank j */
	  else {(*(Rk+i)) = (*(Rk+i))+1;}/* higher, increases its own rank */
	}
    }
}

/*****************************************************************/

 /*  Function to up-to-date a phase vector and its rank vector  */

void UtDRank(double *Vector, int d, int tau, double conf, double NewV, int *RkN)

     /*  Inputs:   Vector   =    address of the "old" vector already ranked                           */
     /*              d      =    dimension of the vector                                              */
     /*             tau     =    delay (needed to search for the old rank vector that is up-to-dates) */
     /*             conf    =    confidence interval (test of the equality)                           */
     /*             NewV    =    point that arrive and induces the up-to-date                         */
     /*  Output:    RkN   =    address of the vector of the ranks up-to-dated                       */
{

  int i;
  int nr = d*tau; /* -nr = shift for searching the old rank vector */
  double OldV = *Vector; /* point that will quit the vector of values */
  int pin = 0; /* rank of the new point */

  for (i = 0; i < d-1; i++) /* for component i */
    {
      *(Vector+i) = *(Vector+i+1); /* up to date of the vector of values */
      *(RkN+i) = *(RkN-nr+i+1); /* copy of the old rank */

      if (NewV >= *(Vector+i)-conf) /* if the new point higher: with the margin due to possible numeric noise */
	{
	  pin++; /* rank increases */
	  if (OldV <= *(Vector+i)+conf){*(RkN+i) = *(RkN+i)-1;} /* and if old was lower => rank i decreases */
	}
      else {if (OldV > *(Vector+i)+conf){*(RkN+i) = *(RkN+i)+1;}} /* conversely : new lower, old higher, rank increases */
    }

  *(Vector+d-1) = NewV; /* up-to-date of the last point of the values */
  *(RkN+d-1) = pin; /* pin is teh rank of the new point */
}

/************************************************************************************/

 /*  Function to evaluate the Lempel-Ziv complexity of a sequence of a rank vector  */
 /*                            (Bandt and Pompe construction of the rank vectors)   */

void BPLZ(double *signal, int d, int tau, double conf, int Tb, double *c, int *Rt)

     /*  Inputs: signal  =   address of the signal                      */
     /*            d     =   embedding dimension                        */
     /*           tau    =   delay                                      */
     /*           conf   =   confidence interval (test of the equality) */
     /*            Tb    =   duration of the bloc of signal             */
     /*  Outputs:  c     =   complexity of the sequence                 */
     /*            Rt    =   sequence of rank vectors                   */

{
  int i, j, k, km, l, eq;
  int t, u;
  double *Y;
  int N = Tb-(d-1)*tau;

  Y = (double*)malloc(d*tau*sizeof(double));


  /* Initializations */
  for (j = 0; j < tau; j++)
    {
      for (i = 0; i < d; i++)
	{
	  (*(Y+i+j*d)) = (*(signal+i*tau+j));
	}
      Rank(Y+j*d, d, conf, Rt+j*d);
    }
  t = tau - 1;

 step0 : *c = 1.0; l = 1; if (tau > 1){goto step2;}

 step1 : if ( t >= (Tb-(d-1)*tau) ){goto stop;}
  else
    {
      t = t+1;
      u = t - (t/tau)*tau; /* (euclidean division) */
      UtDRank(Y+u*d, d, tau, conf, *(signal+t+(d-1)*tau), Rt+t*d);
      goto step2;
    }

 step2 : j = 0; k = 0 ; km = 1; /* we assume production is km = 1 (max length km), j = pointer */

 step3 : eq = EqualVectors(Rt+(l+k)*d , Rt+(j+k)*d , d); /* test of the equal between vect. */
                                                         /* (-1 because the indices stars from 0 in C) */
  if (eq == 1) /* we have still a reproduction */
    {
      k++; /* still a reproduction */
      if ( (l+k) > t ) /* we ended the sequence up to the current point */
	{
	  if ( (l+k) >= (Tb-(d-1)*tau) ) {*c = (*c) + 1.0; goto stop;} /* addition of the last word */
	  else
	    {
	      t = t+1;
	      u = t - (t/tau)*tau; /* (euclidean division) */
	      UtDRank(Y+u*d, d, tau, conf, *(signal+t+(d-1)*tau), Rt+t*d);
	      goto step3;
	    } /* to make the production that is not a reproduction */
	}
      else {goto step3;}
    }
  else /* we face a production that is not a reproduction */
    {
      if ( (k+1) > km ) {km = k+1;}  /* the pointer is better since the production is larger -> up to date nmax */
      j++; /* test of the next point as a possible pointer  */
      if(j == l) /* all the possible pointer are tested */
	{
	  *c = (*c) + 1.0; /* the exhaustive production is made : c */
	  l = l + km; /* we skip to the symbol following the exhaustive word */
	  if ( l > t )
	    {
	      if ( l >= (Tb-(d-1)*tau) ){goto stop;} /* the sequence ended */
	      else
		{
		  t = t+1;
		  u =  t - (t/tau)*tau; /* (euclidean division) */
		  UtDRank(Y+u*d, d, tau, conf,  *(signal+t+(d-1)*tau), Rt+t*d);
		  goto step2;
		}
	    }
	  else {goto step2;} /* we start for seeking the next optimal production step */
	}
      else {k = 0; goto step3;} /* we start again testing the next possible pointer */
    }
 stop : free(Y);
}
