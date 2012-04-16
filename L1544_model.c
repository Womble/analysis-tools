/*
 *  model.c
 *  LIME, The versatile 3D line modeling tool 
 *
 *  Created by Christian Brinch on 11/05/07.
 *  Copyright 2006-2011, Christian Brinch, 
 *  <brinch@strw.leidenuniv.nl>
 *  Sterrewacht Leiden, 
 *	Leiden University. 
 *	All rights reserved.
 *
 */

#include "lime.h"
#include "std_mdl.h"
#include "keto_mdl.h"
#include <stdio.h>

#define h2_PO_ratio 0.75

/******************************************************************************/
//L1544 model and accessor functions

double interpol(double lower, double upper, double frac){
  double x;
  if (0<=frac<=1)    x=upper*frac+lower*(1.-frac);
  else x=(upper+lower)/2.;
  return x;
}

int num_picker(double r, int arr_size){
  int i, x, m=arr_size-1;
  if (r>std_array[0][m]) {
    x=m;
  }
  else for (i=0;i<arr_size;i++){
      if (std_array[0][i]>r) {
	x=i-1;
	i=m;
      }
    }
  //  fprintf(stderr,"r:%.3e x:%d\n", r,x);
  return x;
}

int anum_picker(double r, int arr_size){
  int i, x, m=arr_size-1;
  if (r>keto_array[0][m]) {
    x=m;
  }
  else for (i=0;i<arr_size;i++){
      if (keto_array[0][i]>r) {
	x=i-1;
	i=m;
      }
    }
  //  fprintf(stderr,"r:%.3e x:%d\n", r,x);
  return x;
}

double dens(double x, double y, double z){
  double r=sqrt(x*x+y*y+z*z), dr,frac,ans;
  int i=num_picker(r, std_arr_size), m=std_arr_size-1;
  if (i<0) return std_array[1][0];
  else if (i>=m) return std_array[1][m];
  else {
    dr=std_array[0][i+1]-std_array[0][i];
    frac=(r-std_array[0][i])/dr;
    ans=interpol(std_array[1][i],std_array[1][i+1],frac);
    //    fprintf(stderr,"%.3e %.3e %d %.3e\n", sqrt(x*x+y*y+z*z),ans,i,frac);
    return ans;
  }
}

double temp(double x, double y, double z){
  double r=sqrt(x*x+y*y+z*z), dr,frac,ans;
  int i=num_picker(r, std_arr_size),m=std_arr_size-1;
  if (i<0) return std_array[2][0];
  else if (i>=m) return std_array[2][m];
  else {
    dr=std_array[0][i+1]-std_array[0][i];
    frac=(r-std_array[0][i])/dr;
    ans=interpol(std_array[2][i],std_array[2][i+1],frac);
    //    fprintf(stderr,"r:%.3e t:%.3e\n", sqrt(x*x+y*y+z*z),ans);
    return ans;
  }
}


double dtemp(double x, double y, double z){
  double r=sqrt(x*x+y*y+z*z), dr,frac,ans;
  int i=num_picker(r, std_arr_size), m=std_arr_size-1;
  if (i<0) return std_array[3][0];
  else if (i>=m) return std_array[3][m];
  else {
    dr=std_array[0][i+1]-std_array[0][i];
    frac=(r-std_array[0][i])/dr;
    ans=interpol(std_array[3][i],std_array[3][i+1],frac);
    //    fprintf(stderr,"r:%.3e u:%.3e\n", sqrt(x*x+y*y+z*z),ans);
    return ans;
  }
}


double velo(double x, double y, double z){
  double r=sqrt(x*x+y*y+z*z), dr,frac,ans;
  int i=num_picker(r, std_arr_size),m=std_arr_size-1;
  if (i<0) return std_array[4][0];
  else if (i>=m) return std_array[4][m];
  else {
    dr=std_array[0][i+1]-std_array[0][i];
    frac=(r-std_array[0][i])/dr;
    ans=interpol(std_array[4][i],std_array[4][i+1],frac);
    //    fprintf(stderr,"r:%.3e v:%.3e\n", sqrt(x*x+y*y+z*z),ans);
    return ans;
  }
}

double abuns(double x, double y, double z){
  double r=sqrt(x*x+y*y+z*z), dr,frac,ans;
  int i=anum_picker(r, keto_arr_size),m=keto_arr_size-1;
  if (i<0) return 1.2e-9+2.5e-37*(r*r);
  else if (i>=m) return keto_array[5][m];
  else {
    dr=keto_array[0][i+1]-keto_array[0][i];
    frac=(r-keto_array[0][i])/dr;
    ans=interpol(keto_array[5][i],keto_array[5][i+1],frac);
    //    fprintf(stderr,"r:%.3e a:%.3e %d %.3e %.3e\n", sqrt(x*x+y*y+z*z),ans,i,frac, dr);
    return ans;
  }
}



void
input(inputPars *par, image *img){
/*
 * Basic parameters. See cheat sheet for details.
 */
  par->radius		= +1.e16 ;
  par->minScale	   	= 500*AU;
  par->pIntensity    	= 25000;
  par->sinkPoints    	= 5000;
  par->dust		= "jena_thick_e6.tab";
  par->moldatfile[0] 	= "oh2o-h2.dat";
  par->outputfile 	= "populations_std_hugeits.pop"; 
//par->pregrid		= "pregrid.dat";
  par->gridfile		= "grid_std.vtk";

/* 
 * Definitions for image #0. Add blocks for additional images.
 */
  img[0].nchan			= 101;		  // Number of channels
  img[0].velres			= 25.;       // Channel resolution in m/s
  img[0].trans			= 0;          // zero-indexed J quantum number
  //  img[0].freq                   = 97*1e9;           //central freqency
  //  img[0].bandwidth              = 10*1e9;             //full width of the image
  img[0].pxls			= 401;	      // Pixels per dimension
  img[0].imgres			= 1;		  // Resolution in arc seconds
  img[0].theta			= +0.0;		  // 0: face-on, pi/2: edge-on
  img[0].distance		= 140*PC;	  // source distance in m
  img[0].source_vel		= 0.;          // source velocity in m/s
  img[0].unit			= 0;		  // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename		= "h2ox0p75_1-0_std_KappaThick1e6_100TS.fits";	// Output filename
}

/******************************************************************************/

void
density(double x, double y, double z, double *density){	
  double d;
  d=dens(x,y,z);
  density[0]=d;//*h2_PO_ratio;
  //  density[0]=d*(1-h2_PO_ratio);
  //  fprintf(stderr,"%.3e,%.3e,%.3e r:%.3e d:%.3e\n", x,y,z, sqrt(x*x+y*y+z*z),density[0]);
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
  temperature[0]=temp(x,y,z);
  temperature[1]=dtemp(x,y,z);
}

/******************************************************************************/

void
abundance(double x, double y, double z, double *abundance){
  //  double r=sqrt(x*x+y*y+z*z);
  //  if (r<5000*AU) {
  //    abundance[0]=0.95e-11;
  //  }
  //  else {
  //    abundance[0]=0.95e-8;
  //  }
  //  abundance[0]=4.e-9;
  abundance[0]=0.75*abuns(x,y,z);// / h2_PO_ratio;
//    fprintf(stderr,"r:%.3e v:%.3e\n", sqrt(x*x+y*y+z*z),ans);
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  *doppler = 80.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
  double phi,theta,v;

  v=velo(x,y,z);
  theta=atan2(sqrt(x*x+y*y),z);
  phi=atan2(y,x);  

  vel[0]=v*sin(theta)*cos(phi);
  vel[1]=v*sin(theta)*sin(phi);
  vel[2]=v*cos(theta);
}

/******************************************************************************/

