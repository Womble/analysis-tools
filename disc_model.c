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
//#include "disk.h"
//#include "rho_disk_reduced2.h"
#include "vel.h"
#include "CO.h" //create this header file witht he fits2header.py script froma chemistry fits file
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h> //no idea which of these are need from example so just including them all

#define elements 20971520
#define file 3 //ugly ugly hack, for some reason cant open a file and save its int here or manage to extern it from main, but as it's the only file im opening it *SEEMS* to always be 3 when i open it from main.c so I'll run with it, should ask dan how to do this propperly.
//int file=open("disk.bdat",O_RDONLY);
//  double *disk_arr;
//  disk_arr= (double*)mmap(0, elements*sizeof(double), PROT_READ, MAP_SHARED,file,0);
//  if (disk_arr == MAP_FAILED) {
//    close(file);
//    perror("Error mmapping the file");
//    exit(EXIT_FAILURE);
//    } //open and memmap the binary disk array 

/******************************************************************************/

void
input(inputPars *par, image *img){
/*
 * Basic parameters. See cheat sheet for details.
 */
  par->radius		= 95*AU;
  par->minScale	   	= 2*AU;
  par->pIntensity    	= 30000;
  par->sinkPoints    	= 8000;
  par->sampling         = 0;
  par->dust		= "jena_thick_e6.tab";
  par->moldatfile[0] 	= "co.dat"; //replace this with relevent LAMDA file
  par->outputfile 	= "populations.pop";
  //  par->pregrid	= "pregrid.dat";
  par->gridfile		= "grid.vtk";

  //
  //Definitions for image #0. Add blocks for additional images.
  //

  img[0].nchan			= 251;	      // Number of channels
  img[0].velres			= 75.;        // Channel resolution in m/s
  img[0].trans			= 0;          // zero-indexed J quantum number
  //  img[0].freq                   = 75.0e9;  //central freqency
  //  img[0].bandwidth              = 1500e6;
  img[0].pxls			= 201;	      // Pixels per dimension
  img[0].imgres			= 0.05;      // Resolution in arc seconds
  img[0].theta			= 0.0;	      // 0: face-on, pi/2: edge-on
  img[0].distance		= 12.5*PC;      // source distance in m
  img[0].source_vel		= 0;          // source velocity in m/s
  img[0].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename		= "imageCO_1-0F.fits";	// Output filename

  img[1].nchan			= 251;	      // Number of channels
  img[1].velres			= 75.;        // Channel resolution in m/s
  img[1].trans			= 0;          // zero-indexed J quantum number
  //  img[1].freq                   = 75.0e9;  //central freqency
  //  img[1].bandwidth              = 1500e6;
  img[1].pxls			= 201;	      // Pixels per dimension
  img[1].imgres			= 0.05;      // Resolution in arc seconds
  img[1].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[1].distance		= 12.5*PC;      // source distance in m
  img[1].source_vel		= 0;          // source velocity in m/s
  img[1].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[1].filename		= "imageCO_1-0E.fits";	// Output filename

  img[2].nchan			= 251;	      // Number of channels
  img[2].velres			= 75.;        // Channel resolution in m/s
  img[2].trans			= 2;          // zero-indexed J quantum number
  //  img[2].freq                   = 75.0e9;  //central freqency
  //  img[2].bandwidth              = 1500e6;
  img[2].pxls			= 201;	      // Pixels per dimension
  img[2].imgres			= 0.05;      // Resolution in arc seconds
  img[2].theta			= 0.0;	      // 0: face-on, pi/2: edge-on
  img[2].distance		= 12.5*PC;      // source distance in m
  img[2].source_vel		= 0;          // source velocity in m/s
  img[2].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[2].filename		= "imageCO_3-2F.fits";	// Output filename

  img[3].nchan			= 251;	      // Number of channels
  img[3].velres			= 75.;        // Channel resolution in m/s
  img[3].trans			= 2;          // zero-indexed J quantum number
  //  img[3].freq                   = 75.0e9;  //central freqency
  //  img[3].bandwidth              = 1500e6;
  img[3].pxls			= 201;	      // Pixels per dimension
  img[3].imgres			= 0.05;      // Resolution in arc seconds
  img[3].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[3].distance		= 12.5*PC;      // source distance in m
  img[3].source_vel		= 0;          // source velocity in m/s
  img[3].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[3].filename		= "imageCO_3-2E.fits";	// Output filename

  img[4].nchan			= 251;	      // Number of channels
  img[4].velres			= 75.;        // Channel resolution in m/s
  img[4].trans			= 4;          // zero-indexed J quantum number
  //  img[4].freq                   = 75.0e9;  //central freqency
  //  img[4].bandwidth              = 1500e6;
  img[4].pxls			= 201;	      // Pixels per dimension
  img[4].imgres			= 0.05;      // Resolution in arc seconds
  img[4].theta			= 0.0;	      // 0: face-on, pi/2: edge-on
  img[4].distance		= 12.5*PC;      // source distance in m
  img[4].source_vel		= 0;          // source velocity in m/s
  img[4].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[4].filename		= "imageCO_5-4F.fits";	// Output filename

  img[5].nchan			= 251;	      // Number of channels
  img[5].velres			= 75.;        // Channel resolution in m/s
  img[5].trans			= 4;          // zero-indexed J quantum number
  //  img[5].freq                   = 75.0e9;  //central freqency
  //  img[5].bandwidth              = 1500e6;
  img[5].pxls			= 201;	      // Pixels per dimension
  img[5].imgres			= 0.05;      // Resolution in arc seconds
  img[5].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[5].distance		= 12.5*PC;      // source distance in m
  img[5].source_vel		= 0;          // source velocity in m/s
  img[5].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[5].filename		= "imageCO_5-4E.fits";	// Output filename

  img[6].nchan			= 251;	      // Number of channels
  img[6].velres			= 75.;        // Channel resolution in m/s
  img[6].trans			= 8;          // zero-indexed J quantum number
  //  img[6].freq                   = 75.0e9;  //central freqency
  //  img[6].bandwidth              = 1500e6;
  img[6].pxls			= 201;	      // Pixels per dimension
  img[6].imgres			= 0.05;      // Resolution in arc seconds
  img[6].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[6].distance		= 12.5*PC;      // source distance in m
  img[6].source_vel		= 0;          // source velocity in m/s
  img[6].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[6].filename		= "imageCO_9-8E.fits";	// Output filename

  img[7].nchan			= 251;	      // Number of channels
  img[7].velres			= 75.;        // Channel resolution in m/s
  img[7].trans			= 6;          // zero-indexed J quantum number
  //  img[7].freq                   = 75.0e9;  //central freqency
  //  img[7].bandwidth              = 1500e6;
  img[7].pxls			= 201;	      // Pixels per dimension
  img[7].imgres			= 0.05;      // Resolution in arc seconds
  img[7].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[7].distance		= 12.5*PC;      // source distance in m
  img[7].source_vel		= 0;          // source velocity in m/s
  img[7].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[7].filename		= "imageCO_7-6E.fits";	// Output filename

  /*
  img[3].nchan			= 150;	      // Number of channels
  img[3].velres			= 50.;        // Channel resolution in m/s
  img[3].trans			= 8;          // zero-indexed J quantum number
  //  img[1].freq                   = 75.0e9;  //central freqency
  //  img[1].bandwidth              = 1500e6;
  img[3].pxls			= 151;	      // Pixels per dimension
  img[3].imgres			= 0.05;      // Resolution in arc seconds
  img[3].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[3].distance		= 12.5*PC;      // source distance in m
  img[3].source_vel		= 0;          // source velocity in m/s
  img[3].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[3].filename		= "imageCO_8Ex100v.fits";	// Output filename


  img[4].nchan			= 150;	      // Number of channels
  img[4].velres			= 50.;        // Channel resolution in m/s
  img[4].trans			= 9;          // zero-indexed J quantum number
  //  img[0].freq                   = 75.0e9;  //central freqency
  //  img[0].bandwidth              = 1500e6;
  img[4].pxls			= 151;	      // Pixels per dimension
  img[4].imgres			= 0.05;      // Resolution in arc seconds
  img[4].theta			= 0.0;	      // 0: face-on, pi/2: edge-on
  img[4].distance		= 12.5*PC;      // source distance in m
  img[4].source_vel		= 0;          // source velocity in m/s
  img[4].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[4].filename		= "imageCO_9Fx100v.fits";	// Output filename

  img[5].nchan			= 150;	      // Number of channels
  img[5].velres			= 50.;        // Channel resolution in m/s
  img[5].trans			= 9;          // zero-indexed J quantum number
  //  img[1].freq                   = 75.0e9;  //central freqency
  //  img[1].bandwidth              = 1500e6;
  img[5].pxls			= 151;	      // Pixels per dimension
  img[5].imgres			= 0.05;      // Resolution in arc seconds
  img[5].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[5].distance		= 12.5*PC;      // source distance in m
  img[5].source_vel		= 0;          // source velocity in m/s
  img[5].unit			= 0;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[5].filename		= "imageCO_9Ex100v.fits";	// Output filename

  img[6].nchan			= 150;	      // Number of channels
  img[6].velres			= 50.;        // Channel resolution in m/s
  img[6].trans			= 9;          // zero-indexed J quantum number
  //  img[1].freq                   = 75.0e9;  //central freqency
  //  img[1].bandwidth              = 1500e6;
  img[6].pxls			= 151;	      // Pixels per dimension
  img[6].imgres			= 0.05;      // Resolution in arc seconds
  img[6].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[6].distance		= 12.5*PC;      // source distance in m
  img[6].source_vel		= 0;          // source velocity in m/s
  img[6].unit			= 4;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[6].filename		= "imageCO_9E_Taux100v.fits";	// Output filename

  img[7].nchan			= 150;	      // Number of channels
  img[7].velres			= 50.;        // Channel resolution in m/s
  img[7].trans			= 7;          // zero-indexed J quantum number
  //  img[7].freq                   = 75.0e9;  //central freqency
  //  img[7].bandwidth              = 1500e6;
  img[7].pxls			= 151;	      // Pixels per dimension
  img[7].imgres			= 0.05;      // Resolution in arc seconds
  img[7].theta			= PI/2;	      // 0: face-on, pi/2: edge-on
  img[7].distance		= 12.5*PC;      // source distance in m
  img[7].source_vel		= 0;          // source velocity in m/s
  img[7].unit			= 4;	      // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[7].filename		= "imageCO_7E_Taux100v.fits";	// Output filename
  */
}


/******************************************************************************/

void density(double x, double y, double z, double *density){	
  int X=x/AU*2+128.1, Y=y/AU*2+128.1, Z=z/AU*2+32.1; //convert co-ords to AU, then to cell number, offset to 0,0,0 at -x_max,-y_max,-z_max 
  double *disk_arr;
  disk_arr= (double*)mmap(0, elements*sizeof(double), PROT_READ, MAP_SHARED,file,0);
  if (disk_arr == MAP_FAILED) {
    close(file);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
    } //open and memmap the binary disk array 
  if (X>255) X=255;
  else if (X<0) X=0;
  if (Y>255) Y=255;
  else if (Y<0) Y=0;
  if (Z>63)  Z=63;
  else if (Z<0) Z=0;
  int n=X*256*64+Y*64+Z; // grid-ref to arr pos
  //  fprintf(stderr, "D: %.3e %.3e %.3e : %.3e\n",x,y,z,disk_arr[n*5+1]);
  density[0] = disk_arr[n*5+1];
  munmap(disk_arr,elements*sizeof(double));
  //density[0] = 1.0e6;
}

/******************************************************************************/

void temperature(double x, double y, double z, double *temperature){
int X=x/AU*2+127.1, Y=y/AU*2+127.1, Z=z/AU*2+31.1; //convert co-ords to AU, then to cell number, offset to 0,0,0 at -x_max,-y_max,-z_max
  double *disk_arr;
  disk_arr= (double*)mmap(0, elements*sizeof(double), PROT_READ, MAP_SHARED,file,0);
  if (disk_arr == MAP_FAILED) {
    close(file);
    perror("Error mmapping the file");
    exit(EXIT_FAILURE);
    } //open and memmap the binary disk array
  if (X>255) X=255;
  else if (X<0) X=0;
  if (Y>255) Y=255;
  else if (Y<0) Y=0;
  if (Z>63)  Z=63;
  else if (Z<0) Z=0;
  int n=X*256*64+Y*64+Z; // grid-ref to arr pos
  //  fprintf(stderr, "T: %.3e %.3e %.3e : %.3e\n",x,y,z,disk_arr[n*5+0]);
  temperature[0]=disk_arr[n*5+0];
  munmap(disk_arr, elements*sizeof(double));
  //  temperature[0]=50.;
}

/******************************************************************************/

void abundance(double x, double y, double z, double *abundance){
  int X=x/AU/002.1764+25.1, Y=y/AU/002.1764+25.1, Z=z/AU/0.218214+25.1; //convert co-ords to AU, then to cell number, offset to 0,0,0 at -x_max,-y_max,-z_max
  if (X>50) X=50;
  else if (X<0) X=0;
  if (Y>50) Y=50;
  else if (Y<0) Y=0;
  if (Z>50)  Z=50;
  else if (Z<0) Z=0;
  int n=X*51*51+Y*51+Z;
  double ans =chem_arr[n];
  if (ans>=-0.001) ans=-30.0;
  //  fprintf(stderr, "A: %.3e %.3e %.3e : %.3e\n",x,y,z,pow(10,ans));
  abundance[0] = pow(10,ans); 
  //  abundance[0] = 1.e-11;
}

/******************************************************************************/

void doppler(double x, double y, double z, double *doppler){
  *doppler = 100.;
}

/******************************************************************************/

void velocity(double x, double y, double z, double *velocity){
  int X=2*x/AU+128.01, Y=2*y/AU+128.01, Z=2*z/AU+32.01;
  if (X<0) X=0;
  else if (X>255) X=255;
  if (Y<0) Y=0;
  else if (Y>255) Y=255;
  if (Z<0) Z=0;
  else if (Z>63) Z=63;
  //  fprintf(stderr,"%.3d %.3d %.2e %.2e\n",X-128,Y-128,vel_arr[X][Y][Z][0],vel_arr[X][Y][Z][1]);
  velocity[0]=vel_arr[X][Y][Z][0];
  velocity[1]=vel_arr[X][Y][Z][1];
  velocity[2]=vel_arr[X][Y][Z][2];
  /*  double r =sqrt(x*x+y*y), p=10000/sqrt(r/AU), phi=atan2(y,x);
  velocity[0]=p*sin(phi);
  velocity[1]=p*cos(phi);
  velocity[2]=-1000/(z/AU);*/
}

/******************************************************************************/
