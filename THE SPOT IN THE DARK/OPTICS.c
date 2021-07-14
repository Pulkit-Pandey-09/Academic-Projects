#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <string.h>
#ifndef GNU
#include <malloc.h>
#endif

#define MY_MALLOC_THRESH 1
#ifdef MALLOCDEBUG

/*
extern int nmalloc=0, ncalloc=0, nrealloc=0, nfree=0;
*/

#define MALLOC(name,len,type) { \
    /* nmalloc++; */ \
  if ( (len)*sizeof(type) >= MY_MALLOC_THRESH ) \
      printf("MALLOC(  %12s, %12s=%-6d," #type " )\n", \
      #name, #len, len); \
  if ( ( (name) = (type *)malloc( (len) * sizeof(type) ) ) == NULL \
         && (len)>0) { \
      printf( "MALLOC(" #name "," #len "=%d," #type ")" \
      ": cannot allocate space", len); \
      exit(1); \
  } \
}

#define CALLOC(name,len,type) { \
  int i; \
  /* ncalloc++; */ \
  if ( (len)*sizeof(type) >= MY_MALLOC_THRESH ) \
      printf("CALLOC(  %12s, %12s=%-6d," #type " )\n", \
      #name, #len, len); \
  if ( ( (name) = (type *)calloc( (len) , sizeof(type) ) ) == NULL \
         && (len)>0) { \
      printf( "CALLOC(" #name "," #len "=%d," #type ")" \
            ": cannot allocate space", len); \
      exit(1); \
  } \
  for (i=0; i<(len); i++) name[i] = 0; \
}

#define REALLOC(name,len,type) { \
  /* nrealloc++; */ \
  if ( (len)*sizeof(type) >= MY_MALLOC_THRESH ) \
      printf("REALLOC( %12s, %12s=%-6d," #type " )\n", \
      #name, #len, len); \
  if (((name) = (type *)realloc( (name), (len)*sizeof(type) )) == NULL \
         && (len)>0) { \
      printf( "REALLOC(" #name "," #len "=%d," #type ")" \
      ": cannot reallocate space", len); \
      exit(11); \
  } \
}

#define FREE(name) { \
  /* nfree++; */ \
  printf("FREE(    %12s )\n", #name); \
  if ( (name) != NULL ) free( (name) ); \
  (name) = NULL; \
}

#else

#define MALLOC(name,len,type) { \
  if ( ( (name) = (type *)malloc( (len) * sizeof(type) ) ) == NULL \
    && (len)>0) { \
      printf( "MALLOC(" #name "," #len "=%d," #type ")" \
      ": cannot allocate space", len); \
      exit(1); \
  } \
}

#define CALLOC(name,len,type) { \
  int i; \
  if ( ( (name) = (type *)calloc( (len) , sizeof(type) ) ) == NULL \
          && (len)>0) { \
      printf( "CALLOC(" #name "," #len "=%d," #type ")" \
      ": cannot allocate space", len); \
      exit(1); \
  } \
  for (i=0; i<(len); i++) name[i] = 0; \
}

#define REALLOC(name,len,type) { \
  if (((name) = (type *)realloc( (name), (len)*sizeof(type) )) == NULL \
         && (len)>0) { \
      printf( "REALLOC(" #name "," #len "=%d," #type ")" \
      ": cannot reallocate space", len); \
      exit(1); \
  } \
}

#define FREE(name) { \
  if ( (name) != NULL ) free( (name) ); \
  (name) = NULL; \
}

#endif

// #include "myalloc.h"

void spline(double *x, double *y, int n, double yp0, double ypn1, double *y2);
double splint(double *xa, double *ya, double *y2a, int n, double x);
void compute_PSF(double *rt,double drt, double at, double shft,
	    double **A2_real,double **A2_imag,
	    int Npup0, int Nimg,double *rhos, double **P);
void compute_PSF2(double *rt,double drt, double at, double shft,
	    double **A2_real,double **A2_imag,
	    int Npup0, int Nimg,double *rhos, double **P);

int Nocc, Npup;
int Ns = 0, K = 1;
double units = 1;  // 1 = meters, 1e+6 = microns
double pi, pi2;
double lambda;
double f;
double b, c, d;
int n;
double z, R;
double rho1;
double dr;
double j0(double z){
    return sin(z)/z;
}
double b1(double z){
    return sin(z)/(z*z)-cos(z)/z;
}

double jn(double z,int n){
    double out;
    if (n==0){
        out = j0(z);
    }
    else if(n==1){
        out = b1(z);
    }
    /*using recurrence relation */
    else{
        out = (2*n-1)*jn(z,n-1)/z-jn(z,n-2);
    }
    return out;
}

main() 
{
    int i,j,k,lam=5,sgn;

    Nocc=100000;
    Npup=2000;
    z = 21; R = 0.5*0.36*89*1e-3; c = 0.5*89*1e-3, d = 0.5*100*1e-3;
    lambda = 5.32e-7;

    rho1 = d;

    pi = 4*atan(1);
    pi2 = pi/2;

    double *r=NULL, *rho=NULL;
    double *A = NULL;
    char line[80];
    char fname[50];
    FILE *fp;
    
    /**/
    REALLOC(r, Nocc, double);
    REALLOC(A, Nocc, double);
    for (i=0; i<Nocc; i++) {
	r[i] = d*(i+0.5)/Nocc;
	if (r[i] < R || r[i] > c) {
	    A[i] = 0;
	} else {
	    A[i] = 1;
	}
    }
    dr = r[1]-r[0];
    fp = fopen("Aweb","w");
    for (i=0; i<Nocc; i++) { 
	fprintf(fp,"%10.7f %10.7f \n",r[i],A[i]); 
    }
    fclose(fp);
    /**/


    REALLOC(rho, Npup, double);
    for (i=0; i<Npup; i++) {
	rho[i] = (rho1*i)/Npup;
    }

    double **Ereal = NULL;
    REALLOC(Ereal, K, double *);
    for (k=0; k<K; k++) { Ereal[k]=NULL; REALLOC(Ereal[k], Npup, double); }

    double **Eimag = NULL;
    REALLOC(Eimag, K, double *);
    for (k=0; k<K; k++) { Eimag[k]=NULL; REALLOC(Eimag[k], Npup, double); }

      printf("computing Ereal and Eimag for lambda = %e...\n",lambda); fflush(stdout);
      sgn = 1;
      for (k=0; k<K; k++) {
	for (i=0; i<Npup; i++) Ereal[k][i] = 0;
	for (i=0; i<Npup; i++) Eimag[k][i] = 0;
	for (i=0; i<Npup; i++) {
	    for (j=0; j<Nocc; j++) {
		if (k==0) {
		  Ereal[k][i] += sin(pi*(r[j]*r[j]+rho[i]*rho[i])/(lambda*z))*
			      A[j]* j0(-2*pi*r[j]*rho[i]/(lambda*z))*r[j]*dr;
		  Eimag[k][i] += cos(pi*(r[j]*r[j]+rho[i]*rho[i])/(lambda*z))*
			      A[j]* j0(-2*pi*r[j]*rho[i]/(lambda*z))*r[j]*dr;
		} else {
		  Ereal[k][i] += sin(pi*(r[j]*r[j]+rho[i]*rho[i])/(lambda*z))*
			      (sin(pi*k*A[j])/(pi*k))*
			      jn(k*Ns,-2*pi*r[j]*rho[i]/(lambda*z))*r[j]*dr;
		  Eimag[k][i] += cos(pi*(r[j]*r[j]+rho[i]*rho[i])/(lambda*z))*
			      (sin(pi*k*A[j])/(pi*k))*
			      jn(k*Ns,-2*pi*r[j]*rho[i]/(lambda*z))*r[j]*dr;
		}
	    }
	    Ereal[k][i] *=  2*pi*sgn/(lambda*z);
	    Eimag[k][i] *= -2*pi*sgn/(lambda*z);
	    /*
	    if (k == 0) {
		Ereal[k][i] = 1 - Ereal[k][i];
		Eimag[k][i] =   - Eimag[k][i];
	    } else {
		Ereal[k][i] =   - Ereal[k][i];
		Eimag[k][i] =   - Eimag[k][i];
	    }
	    */
	}
	sgn *= -1;
      }

      printf("computing intensity vs. radius for lambda = %e...\n",lambda); fflush(stdout);
      double *Iofr=NULL; REALLOC(Iofr, Npup, double);
      for (i=0; i<Npup; i++) {
	  Iofr[i] = 0;
	  for (k=0; k<K; k++) {
	      Iofr[i] += Ereal[k][i]*Ereal[k][i] + Eimag[k][i]*Eimag[k][i];
	  }
      }

      sprintf(fname,"Iofr%0d",lam);
      fp = fopen(fname,"w");
      for (i=0; i<Npup; i++) { 
	fprintf(fp,"%10.7f %10.7f \n",rho[i], log10(Iofr[i]));
      }
      fclose(fp);
      // FREE(Iofr);

      printf("writing Ereal for lambda = %e...\n",lambda); fflush(stdout);
      sprintf(fname,"Ereal%0d",lam);
      fp = fopen(fname,"w");
      for (i=0; i<Npup; i++) { 
	fprintf(fp,"%10.7f ",rho[i]);
        for (k=0; k<K; k++) {
	    fprintf(fp,"%10.3e ",Ereal[k][i]);
	}
	fprintf(fp,"\n");
	fflush(fp);
      }
      fclose(fp);

      printf("writing Eimag for lambda = %e...\n",lambda); fflush(stdout);
      sprintf(fname,"Eimag%0d",lam);
      fp = fopen(fname,"w");
      for (i=0; i<Npup; i++) { 
	fprintf(fp,"%10.7f ",rho[i]);
        for (k=0; k<K; k++) {
	    fprintf(fp,"%10.3e ",Eimag[k][i]);
	}
	fprintf(fp,"\n");
	fflush(fp);
      }
      fclose(fp);

      printf("computing pupil plane electric field for lambda = %e...\n",lambda); fflush(stdout);
      int px, py, Npup0 = 1000;
      double x, y, rho2, phi, **Er=NULL, **Ei=NULL;
      REALLOC(Er, 2*Npup0+1, double *); Er += Npup0;
      for (px=-Npup0; px<=Npup0; px++) 
          { Er[px]=NULL; REALLOC(Er[px], 2*Npup0+1, double); Er[px] += Npup0;}
      REALLOC(Ei, 2*Npup0+1, double *); Ei += Npup0;
      for (px=-Npup0; px<=Npup0; px++) 
          { Ei[px]=NULL; REALLOC(Ei[px], 2*Npup0+1, double); Ei[px] += Npup0;}
      FILE *fp1, *fp2;
      sprintf(fname,"Er%0d",lam); fp1 = fopen(fname,"w");
      sprintf(fname,"Ei%0d",lam); fp2 = fopen(fname,"w");
      for (px=-Npup0; px<=Npup0; px++) {
	  for (py=-Npup0; py<=Npup0; py++) {
	      x = px*rho1/Npup0;
	      y = py*rho1/Npup0;
	      rho2 = sqrt(x*x+y*y);
	      phi = atan2(y,x);
	      i = (int) (Npup*rho2/rho1);
	      if (i < Npup) {
		  Er[px][py] = Ereal[0][i];
		  Ei[px][py] = Eimag[0][i];
		  for (k=1; k<K; k++) {
		      Er[px][py] += Ereal[k][i]*2*cos(k*Ns*(phi-pi/2));
		      Ei[px][py] += Eimag[k][i]*2*cos(k*Ns*(phi-pi/2));
		  }
	      } else {
		  Er[px][py] = 0;
		  Ei[px][py] = 0;
	      }
	      fprintf(fp1,"%10.3e ", Er[px][py]);
	      fprintf(fp2,"%10.3e ", Ei[px][py]);
	  }
	  fprintf(fp1,"\n"); fprintf(fp2,"\n");
	  fflush(fp1); fflush(fp2);
      }
      fclose(fp1); fclose(fp2);

      int Nimg = 200;
      double **P = NULL; REALLOC(P, 2*Nimg+1, double *); P += Nimg;
      for (px=-Nimg; px<=Nimg; px++) 
          { P[px]=NULL; REALLOC(P[px], 2*Nimg+1, double); P[px] += Nimg;}
      double *rhos = NULL; REALLOC(rhos, 2*Nimg+1, double); rhos += Nimg;

      rho = NULL; REALLOC(rho, 2*Npup0+1, double); rho += Npup0;
      for (px=-Npup0; px<=Npup0; px++) { rho[px] = px*rho1/Npup0; }

      f = 400;
      int s;
      printf("computing image plane intensity for lambda = %e...\n",lambda); fflush(stdout);
      /*
      for (s=0; s<1; s++) {
	  printf("shift = %d \n",s);
	  compute_PSF2(rho,rho[1]-rho[0],2.0, s*R, Er,Ei, Npup0, Nimg,rhos, P);
	  sprintf(fname,"psf%0d_%0d",lam,s);
	  fp = fopen(fname,"w");
	  for (px=-Nimg; px<=Nimg; px++) { 
	     for (py=-Nimg; py<=Nimg; py++) { 
		fprintf(fp,"%15e ", P[px][py]); 
	     }
	     fprintf(fp,"\n");
	     fflush(fp);
	  }
	  fclose(fp);
      }
      */

}

void compute_PSF(double *rt,double dx, double at, double shft,
	    double **Er,double **Ei, int Npup0, int Nimg, double *rhos, double **P)
{
    int i, j, px, py;
    
    // Fourier Transform to image plane electric field
    double drho = 11.2e-6 * units;    // THIS IS 3 TIMES TOO SMALL

    double **E_real=NULL; REALLOC(E_real, 2*Nimg+1, double *); E_real+=Nimg;
    double **E_imag=NULL; REALLOC(E_imag, 2*Nimg+1, double *); E_imag+=Nimg;
    for (px=-Nimg; px<=Nimg; px++) 
        { E_real[px]=NULL; REALLOC(E_real[px], 2*Nimg+1, double); E_real[px]+=Nimg;}
    for (px=-Nimg; px<=Nimg; px++) 
        { E_imag[px]=NULL; REALLOC(E_imag[px], 2*Nimg+1, double); E_imag[px]+=Nimg;}

    int xx = (int) (at/dx);
    int ss = (int) (shft/dx);
    for (px=-Nimg; px<=Nimg; px++) { rhos[px] = px*drho; }
    for (px=-Nimg; px<=Nimg; px++) { 
      for (py=-Nimg; py<=Nimg; py++) { 
	E_real[px][py] = 0;
	E_imag[px][py] = 0;
	for (i=-xx+ss; i<=xx+ss; i++) {
	  for (j=-xx; j<=xx; j++) {
	    if ((i-ss)*(i-ss)+j*j > xx*xx) continue;
	    E_real[px][py] 
		+= cos(-2*pi*(rhos[px]*rt[i]+rhos[py]*rt[j])/(lambda*f))*Er[i][j]
		 - sin(-2*pi*(rhos[px]*rt[i]+rhos[py]*rt[j])/(lambda*f))*Ei[i][j];
	    E_imag[px][py] 
		+= cos(-2*pi*(rhos[px]*rt[i]+rhos[py]*rt[j])/(lambda*f))*Ei[i][j]
		 + sin(-2*pi*(rhos[px]*rt[i]+rhos[py]*rt[j])/(lambda*f))*Er[i][j];
	  }
	}
	E_real[px][py] *= 2*pi*dx*dx/(lambda*f);
	E_imag[px][py] *= 2*pi*dx*dx/(lambda*f);

	P[px][py] = E_real[px][py]*E_real[px][py] + E_imag[px][py]*E_imag[px][py];
      }
    }
}

// Faster n^3 algorithm
void compute_PSF2(double *rt,double dx, double at, double shft,
	    double **Er,double **Ei, int Npup0, int Nimg, double *rhos, double **P)
{
    int i, j, px, py;
    
    // Fourier Transform to image plane electric field
    double drho = 11.2e-6 * units;    // THIS IS 3 TIMES TOO SMALL

    double **E_real=NULL; REALLOC(E_real, 2*Nimg+1, double *); 
    double **E_imag=NULL; REALLOC(E_imag, 2*Nimg+1, double *); 
    E_real+=Nimg;
    E_imag+=Nimg;
    for (px=-Nimg; px<=Nimg; px++) { 
	E_real[px]=NULL;REALLOC(E_real[px], 2*Nimg+1, double); E_real[px]+=Nimg;
    }
    for (px=-Nimg; px<=Nimg; px++) { 
	E_imag[px]=NULL;REALLOC(E_imag[px], 2*Nimg+1, double); E_imag[px]+=Nimg;
    }

    double **Etmp_real=NULL; REALLOC(Etmp_real, 2*Npup0+1, double *); 
    double **Etmp_imag=NULL; REALLOC(Etmp_imag, 2*Npup0+1, double *); 
    Etmp_real+=Npup0;
    Etmp_imag+=Npup0;
    for (i=-Npup0; i<=Npup0; i++) {
	Etmp_real[i]=NULL; 
	REALLOC(Etmp_real[i], 2*Nimg+1, double); Etmp_real[i]+=Nimg;
    }
    for (i=-Npup0; i<=Npup0; i++) {
	Etmp_imag[i]=NULL;
	REALLOC(Etmp_imag[i], 2*Nimg+1, double); Etmp_imag[i]+=Nimg;
    }

    int xx = (int) (at/dx);
    int ss = (int) (shft/dx);
    for (px=-Nimg; px<=Nimg; px++) { rhos[px] = px*drho; }

    for (py=-Nimg; py<=Nimg; py++) { 
	for (i=ss-xx; i<=ss+xx; i++) {
	  if (i<-Npup0 || i>Npup0) {
	      printf("help! ss=%d, xx=%d, Npup0=%d \n", ss,xx,Npup0); 
	      fflush(stdout);
	  }
	  Etmp_real[i][py] = 0;
	  Etmp_imag[i][py] = 0;
	  for (j=-xx; j<=xx; j++) {
	    if ((i-ss)*(i-ss)+j*j > xx*xx) continue;
	    Etmp_real[i][py] 
		+= cos(-2*pi*(rhos[py]*rt[j])/(lambda*f))*Er[i][j]
		 - sin(-2*pi*(rhos[py]*rt[j])/(lambda*f))*Ei[i][j];
	    Etmp_imag[i][py] 
		+= cos(-2*pi*(rhos[py]*rt[j])/(lambda*f))*Ei[i][j]
		 + sin(-2*pi*(rhos[py]*rt[j])/(lambda*f))*Er[i][j];
	  }
	}
    }

    for (px=-Nimg; px<=Nimg; px++) { 
        for (py=-Nimg; py<=Nimg; py++) {
	  E_real[px][py] = 0;
	  E_imag[px][py] = 0;
	  for (i=-Npup0; i<=Npup0; i++) {
	    E_real[px][py] 
		+= cos(-2*pi*(rhos[px]*rt[i])/(lambda*f))*Etmp_real[i][py]
		 - sin(-2*pi*(rhos[px]*rt[i])/(lambda*f))*Etmp_imag[i][py];
	    E_imag[px][py] 
		+= cos(-2*pi*(rhos[px]*rt[i])/(lambda*f))*Etmp_imag[i][py]
		 + sin(-2*pi*(rhos[px]*rt[i])/(lambda*f))*Etmp_real[i][py];
	  }
	  E_real[px][py] *= 2*pi*dx*dx/(lambda*f);
	  E_imag[px][py] *= 2*pi*dx*dx/(lambda*f);
	}
    }

    for (px=-Nimg; px<=Nimg; px++) { 
      for (py=-Nimg; py<=Nimg; py++) { 
	P[px][py] = E_real[px][py]*E_real[px][py] + E_imag[px][py]*E_imag[px][py];
      }
    }
}

void spline(double x[], double y[], int n, double yp0, double ypn1, double y2[])

{
    int i;
    for (i=1; i<n-1; i++) {
        y2[i] = 2*( (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]) )
	        /(x[i+1]-x[i-1]);
    }
    y2[0] = 2*(y[1]-y[0] - yp0*(x[1]-x[0])) / ((x[1]-x[0])*(x[1]-x[0]));
    y2[n-1] = 2*(y[n-2]-y[n-1] - ypn1*(x[n-2]-x[n-1])) 
	      /((x[n-2]-x[n-1])*(x[n-2]-x[n-1]));
}

double splint(double xa[], double ya[], double y2a[], int n, double x)

{
    int klo, khi, k;
    double h, b, a;

    klo = 0;
    khi = n-1;
    while (khi-klo > 1) {
        k = (khi+klo) >> 1;
	if (xa[k] > x) khi = k;
	else klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0) {printf("spline error \n"); return 0;}
    a = (xa[khi]-x)/h;
    b = (x-xa[klo])/h;
    return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6;
}
