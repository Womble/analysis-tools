#undef MAIN
#define USER
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <math.h>
#include "const.h"
#include "struct.h"
#include "model.h"
#include "user.h"
#include "funcs.h"
#include "para.h"
 
#include "starwind_16p2.h"
#include "innerDisc.h"

#include <map>
#include <string>

#include "disc_1e-4_16p2.h"
#include "largeDisc_16p2.h"

#define M_acc     (1e-4*m_sol/(3600*24*365.25))  
#define M_star    (10*m_sol)
#define r_star    (16.2*r_sol)
#define L_star    (8500*L_sol)

typedef std::map<std::string,double> TStrDblMap;
typedef std::pair<std::string,double> TStrDblPair;


TStrDblMap radCacheR;
TStrDblMap radCachez;
   
/* solar values  */
#define rho_sol (1.408e3)   
#define m_sol   (1.9891e30) 
#define r_sol   (6.955e8)   
#define L_sol   (3.839e26)  
#define amu     (1.660538921e-27)
#define MMM     (0.909*1.008*amu+0.088*4.002602*amu+0.003*12*amu)
#define GravConst (6.67384e-11) 
#define C0 (2.998e8/r_star)  
#define pi (3.141592653589793)
#define sigma_e (0.6652458545e-28) 
#define forceScalingFactor (sigma_e/MMM/3e8/r_star)

#define MG (10*m_sol*GravConst/pow(r_star,3))   
#define rho0 (1e-5*pow(r_star*5.5/16.2,3))
#define rho_amb (1e13)
#define kb      (1.38065e-23)
#define cs0 (sqrt(kb*13.8e3/MMM)/r_star)
#define cs_amb (cs0*sqrt(rho0/rho_amb))

#define disc_albedo (1.000e0) 
#define k  (0.79)     
#define alpha (0.523) 

#define Mmax (1.0e3)    
#define continuum (1.0)
#define discP (-2.0/4) 
 
#define rsoft (1e-3)  
#define fixOuterRadius (1.1) 
#define discInnerRadius (1.0) 
#define flatDisc (0)          

#define radflag 1  
#define gravflag 1 
#define radforceupdate 2 
#define caching 1 
#define insDisc 1 
#define insDiscRadius 9.5 
 
#define array_R_size (100)
#define array_R_max (10)
#define array_z_size (100)
#define array_z_max (10)
#define array_log_rd_max (log10(10))

extern bool get_r(float &value);
extern bool get_i(int &value);
extern bool get_n(char key[]);
extern void report();
extern void write_var(char item[], float value);
extern void write_var(char item[], int size, float value[]); 
extern void write_var(char item[], int value);
extern void write_var(char item[], int size, int value[]);
extern void write_var(char item[], bool value);
extern void write_var(char item[], char value[]);
extern void read_var(float *value);
extern void read_var(int *value);
extern void read_var(bool *value);
extern void read_var(char *value);
extern void add_usrvar(char item[], char type[], void *value);
extern void merge_array(float buf[], int nvalues[]);
extern float global_sum(float value);
extern float global_max(float value);
extern void broadcast(bool *value);
extern void broadcast(float buf[], int nvalues);
void usrqis(float xc[], float qa[]);
float usrss(float d, float p);
void usrqp(float pa[], float qa[]);
void usrpq(float qa[], float pa[]);
void usrfxq(cell *cptr);
void usrfxp(float pa[]);
bool usrps(cell *cptr, float &s);
void usrsaf(float dt, cell *cptr, float so[]);
void usrref(cell* cptr, float qref[]);
float usrqsol( cell *cptr);
float usrte(float pa[]);
void usrmu(float pa[], float mu[]);
void usrbd(int id, int bound, join *jptr, float pl[], float pr[]);
bool usrcmd(char key[]);
void usrprt();
void usrwr();
bool usrrd(char vname[]);
void usrsum();
void usrmon();
void usrdef();
void usrkeps();
void usrsetg();
bool usroff();
float usrvel(int id, struct cell *cptr);
float usrbfield(int id, struct cell *cptr);
float usrgfield(int id, struct cell *cptr);
void usrga(int id, struct cell *cptr, struct cell *left,
  struct cell *right);
void usrprompt(char prompt[]);
void usrpath(char path[]);



/*
+ usrqis - classify and initialise cell, given coordinates.
     Description:
        This initialises qa in  the cell with coordinates xc according
        to cas.

pa[iqd] density of fluid
pa[iqu0+1] x velocity of fluid
pa[iqu0+2] y velocity of fluid
pa[iqu0+3] z velocity of fluid
pa[iqe] pressure of fluid
pa[iqie] internal energy per unit volume (only if INTEN command is invoked)
pa[iqal0+1] value of scalar 1
pa[iqal0+2] value of scalar 2 etc
*/
void usrqis(float xc[], float qa[])
{
  int i, id;
  float pa[MAXQ];
  /* Defaults (safety) */
  pa[iqd] = 1.0;
  for (id = 1; id <= nd; id++)
    pa[iqu0+id] = 0.0;
  pa[iqe] = 1.0;
  for (i = iqe + 1; i < nq; i++)
    pa[i] = 0.0;
  /* Sod case */
  if (strcmp(cas, "SOD") == 0){
    pa[iqu0+1] = 0.0;
    if (xc[0] <= 0.0){
      pa[iqd] = 1.0;
      pa[iqe] = 1.0;
    }
    else{
      pa[iqd] = 0.125;
      pa[iqe] = 0.1;
    }
  }
  if (strcmp(cas, "SOD2") == 0){
    pa[iqu0+1] = 0.0;
    if (xc[0] <= 0.0 || xc[0] >= 2){
      pa[iqd] = 1.0;
      pa[iqe] = 1.0;
    }
    else{
      pa[iqd] = 0.125;
      pa[iqe] = 0.1;
    }
  }
  else if (strcmp(cas, "SODXY") == 0){
    if (xc[0] <= 0.0){
      pa[iqd] = 1.0;
      pa[iqe] = 1.0;
    }
    else{
      pa[iqd] = 0.125;
      pa[iqe] = 0.1;
    }
  }
  else if (strcmp(cas, "SODYX") == 0){
    if (xc[1] <= 0.0){
      pa[iqd] = 1.0;
      pa[iqe] = 1.0;
    }
    else{
      pa[iqd] = 0.125;
      pa[iqe] = 0.1;
    }
  }
 else if (strcmp(cas, "EMPTY") == 0){
      pa[iqd] = 1;
      pa[iqe] = 1;
      if (xc[1]*xc[1]+xc[0]*xc[0]<1 || xc[1]/xc[0]<0.25) pa[iqd]*=2;
  }
  else if (strcmp(cas, "CIRCLE") == 0){
    if (xc[0]*xc[0] + xc[1]*xc[1] <= 1.0){
      pa[iqd] = 1.0;
      pa[iqe] = 1.0;
      pa[iqal0+1] = 1.0;
    }
    else{
      pa[iqd] = 0.125;
      pa[iqe] = 0.1;
      pa[iqal0+1] = 0.0;
    }
  }

 else if (strcmp(cas, "RTINST") == 0){
   double r=sqrt(xc[1]*xc[1]+xc[0]*xc[0]);
   if (sqrt(xc[1]*xc[1]+xc[0]*xc[0]) < 0.2){
      pa[iqd] = 0.1;
      pa[iqe] = pa[iqd]*10*10;
      pa[iqal0+1] = 1.0;
    }
    else{
      pa[iqd] = 1.0;
      pa[iqe] = pa[iqd]*1*1;
      pa[iqal0+1] = 0.0;
    }
  }
  else if (strcmp(cas, "GRAVXYZ") == 0){
    float roff= 0.52;
    pa[iqu0+1] = 0.0;
    pa[iqu0+2] = 0.0;
    pa[iqu0+3] = 0.0;
    pa[iqe] = pfloor;
    if (xc[0]*xc[0] + xc[1]*xc[1] + (xc[2] - roff)*(xc[2] -roff)
      < rat*rat)
      
      
      pa[iqd] = 1.0;
    else
      pa[iqd] = dfloor;
  }
  else if (strcmp(cas, "GRAVRZ") == 0){
    pa[iqu0+1] = 0.0;
    pa[iqu0+2] = 0.0;
    pa[iqe] = pfloor;
    float r;
    r = sqrt(xc[0]*xc[0] + xc[1]*xc[1]);
    if (r < rat)
      pa[iqd] = 1.0;
    else
      pa[iqd] = dfloor;
  }
  else if (strcmp(cas, "GRAVX") == 0){
    pa[iqu0+1] = 0.0;
    pa[iqal0+1] = 0.0;
    if (mabs(xc[0]) < rat)
      pa[iqd] = 1.0;
    else
      pa[iqd] = 1.0e-10;
    pa[iqe] = 1.0;
  }
  else if (strcmp(cas, "GRAVXY") == 0){
    pa[iqu0+1] = 0.0;
    pa[iqu0+2] = 0.0;
    if (xc[0]*xc[0] + xc[1]*xc[1] < rat*rat){
      pa[iqd] = 1.0;
      pa[iqe] = 0.1*pa[iqd];
    }
    else{
      pa[iqd] = 1.0e-4;
      pa[iqe] = 0.1*pa[iqd];
    }
  }
  else if (strcmp(cas, "EGRAVRZ") == 0){
    pa[iqu0+1] = 0.0;
    pa[iqu0+2] = 0.0;
    pa[iqe] = 1.0;
    float r;
    r = sqrt(xc[0]*xc[0] + xc[1]*xc[1]);
    if (r < rat){
      pa[iqd] = 1.0;
    if (grav)
      pa[iqpot] = 4*1.0642*gkpi*(r*r/6.0 - 0.5*rat*rat);
    }
    else{
      pa[iqd] = 1.0e-10;
      if (grav)
	pa[iqpot] = - 4*1.0642*gkpi*rat*rat*rat/(3.0*r);
    }
  }
  else if (strcmp(cas, "KHINST") == 0){
    g=1.001;
    if (xc[1]+0.005*sin(3*3.14*xc[0])>=0){
      pa[iqu0+1] = 5.0;
      pa[iqd]=0.1;
      pa[iqal0+1]=1;
    }
    else{
      pa[iqu0+1] = -5.0;
      pa[iqd]=10.0;
    }
  }

  else if (strcmp(cas, "STARWIND") == 0){
    g=1.001;
    if (xc[0]<1){
      pa[iqd]=rho0;
      pa[iqe]=rho0*cs0*cs0;
    } else {
      pa[iqd]=rho_amb;
      pa[iqe]=rho_amb*cs0*cs0;
      pa[iqu0+1]=1000/r_star*sqrt(1-1/xc[0]);
    }
  }
  else if (strcmp(cas, "STARWIND2") == 0){
    g=1.001;
    double R=fabs(xc[0]);
    double z=xc[1];
    double r=sqrt(z*z+R*R+rsoft); 
    double st=R/r; 
    double ct=sqrt(1-st*st);
    int RR = (r/10*2048 + 0.5);
    if (RR>2047) RR=2047;
    if (r<1){
      pa[iqd]=rho0;
      pa[iqe]=rho0*cs0*cs0;
    } else {
      pa[iqd]=star_arr[RR][1];
      pa[iqe]=pa[iqd]*cs0*cs0;
      pa[iqu0+1]=star_arr[RR][2]*st;
      pa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st);
    }
  }

   else if (strcmp(cas, "DISCSTAR") == 0){
    g=1.001;
    double R=fabs(xc[0]);
    double z=xc[1];
    double r=sqrt(z*z+R*R+rsoft); 
    double st=R/r; 
    double ct=sqrt(1-st*st);
    double cs=cs0*pow(R*sqrt(r_star/r_sol/5.5), 0.5 * discP); 
    if (R<1) cs=cs0*pow(1.0/5, 0.5 * discP);
    double H=cs/sqrt(MG/pow(R,3));
    double gaussian=exp(-pow(z,2)/2/H/H);
    if (flatDisc ) {
	if (z<0.02) gaussian=1;
	else gaussian=1e-10;
    }
    int RR = (r/10*2048 + 0.5);
    if (RR>2047) RR=2047;

	pa[iqd]=star_arr[RR][1];
	if (r>10) pa[iqd]*=100/r/r; 
	pa[iqe]=pa[iqd]*cs*cs;
    if (radflag && L_star > L_sol){
      pa[iqu0+1]=star_arr[RR][2]*st; 
      pa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st);
      }else{
      pa[iqu0+1]=0;
      pa[iqu0+2]=0;
    }
    
    if (R >= discInnerRadius && z <(8*H)){
      double rho=rho0*pow(R, 3*discP/2) * gaussian;
      if (rho>pa[iqd]){
	pa[iqu0+1]=star_arr[RR][2]*st*pa[iqd]/rho * gaussian;
	pa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st)*pa[iqd]/rho * gaussian;
	pa[iqd]=rho;
	pa[iqe]=rho*cs*cs;
      }
    }
    if(pa[iqd]<rho_amb) {
      pa[iqd]+=rho_amb;
      pa[iqe]+=rho_amb*cs*cs;
    }
    pa[iqal0+1]=sqrt( R * MG )*st;  
    if (z>(6*H)){ 
      H*=2;
      gaussian=exp(-pow(z,2)/2/H/H);
      pa[iqal0+1]*=gaussian;
    }
  }

  else if (strcmp(cas, "LARGEDISCSTAR") == 0){
    g=1.001;
    double R=fabs(xc[0]);
    double z=xc[1];
    double r=sqrt(z*z+R*R+rsoft); 
    double st=R/r; 
    double ct=sqrt(1-st*st);
    double cs=cs0*pow(R*sqrt(r_star/r_sol/5.5), 0.5 * discP); 
    if (R<1) cs=cs0*pow(1.0/5, 0.5 * discP);
    double H=cs/sqrt(MG/pow(R,3));

    pa[iqd]=rho_amb/r/r;
    pa[iqe]=pa[iqd]*cs*cs;
    if (insDisc){
      pa[iqd]=star_arr[2047][1]*100/r/r;
      pa[iqe]=pa[iqd]*cs*cs;
      pa[iqu0+1]=star_arr[2047][2]*st;
      pa[iqu0+2]=star_arr[2047][2]*ct;
    } 


    int ZZ,RR = (r/10*2048 + 0.5);
    if (RR>2047) RR=2047;
    double gaussian=exp(-pow(z,2)/2/H/H);
    double rho=rho0*pow(R*sqrt(r_star/r_sol/5.5), 3*discP/2) * gaussian;
    if (rho>pa[iqd]){
      pa[iqu0+1]=star_arr[RR][2]*st*pa[iqd]/rho * gaussian;
      pa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st)*pa[iqd]/rho * gaussian;
      pa[iqd]=rho;
      pa[iqe]=rho*cs*cs;
    } 
  
    if(pa[iqd]<(rho_amb/r/r)) {
      pa[iqd]+=rho_amb/r/r;
      pa[iqe]+=rho_amb/r/r*cs*cs;
    }

    pa[iqal0+1]=sqrt( R * MG )*st;  
    if (z>(6*H)){ 
      H*=2;
      gaussian=exp(-pow(z,2)/2/H/H);
      pa[iqal0+1]*=gaussian;
    } 


    RR = (acos(st)/(pi/2)*905);
    /* ZZ = (z/10.0*128+0.5); */

    if (RR>904) RR=904;
    /* if (ZZ>127) ZZ=127; */

    if (insDisc && r<(1.05*insDiscRadius)){ 
      pa[iqd]=largeDisc[0][RR];
      pa[iqe]=largeDisc[1][RR];
      pa[iqu0+1]=largeDisc[2][RR];
      pa[iqu0+2]=largeDisc[3][RR];
      pa[iqal0+1]=largeDisc[4][RR];
    }  else if (insDisc && pa[iqd]<(largeDisc[0][RR]*100/r/r)){
      pa[iqd]=largeDisc[0][RR]*100/r/r;
      pa[iqe]=largeDisc[1][RR]*100/r/r;
      pa[iqu0+1]=largeDisc[2][RR];
      pa[iqu0+2]=largeDisc[3][RR];
    } 


    /*   RR = st*128; */
    /*   ZZ = ct*128; */
    /*   if (RR>126) RR=126; */
    /*   if (ZZ>126) ZZ=126; */
    /*   double Rw = st*128-int(st*128); */
    /*   double zw = ct*128-int(ct*128); */
      
    /*   

    /*   rho=innerDisc_arr[0][RR][ZZ]    *((1-Rw)*(1-zw))+ */
    /* 	  innerDisc_arr[0][RR+1][ZZ]  *((  Rw)*(1-zw))+ */
    /* 	  innerDisc_arr[0][RR][ZZ+1]  *((1-Rw)*(  zw))+ */
    /* 	  innerDisc_arr[0][RR+1][ZZ+1]*((  Rw)*(  zw)); */


    /*   if (rho>pa[iqd]){pa[iqd]=rho;} */
    /*   pa[iqe]=innerDisc_arr[3][RR][ZZ]    *((1-Rw)*(1-zw))+ */
    /* 	      innerDisc_arr[3][RR+1][ZZ]  *((  Rw)*(1-zw))+ */
    /* 	      innerDisc_arr[3][RR][ZZ+1]  *((1-Rw)*(  zw))+ */
    /* 	      innerDisc_arr[3][RR+1][ZZ+1]*((  Rw)*(  zw))*rho/pa[iqd]; */

    /*   pa[iqal0+1]=innerDisc_arr[4][RR][ZZ]    *((1-Rw)*(1-zw))+ */
    /* 	          innerDisc_arr[4][RR+1][ZZ]  *((  Rw)*(1-zw))+ */
    /* 	          innerDisc_arr[4][RR][ZZ+1]  *((1-Rw)*(  zw))+ */
    /* 	          innerDisc_arr[4][RR+1][ZZ+1]*((  Rw)*(  zw)); */


    /*   pa[iqu0+1]=innerDisc_arr[1][RR][ZZ]    *((1-Rw)*(1-zw))+ */
    /*     	 innerDisc_arr[1][RR+1][ZZ]  *((  Rw)*(1-zw))+ */
    /* 	         innerDisc_arr[1][RR][ZZ+1]  *((1-Rw)*(  zw))+ */
    /* 	         innerDisc_arr[1][RR+1][ZZ+1]*((  Rw)*(  zw))*rho/pa[iqd]; */

    /*   pa[iqu0+2]=innerDisc_arr[2][RR][ZZ]    *((1-Rw)*(1-zw))+ */
    /* 	         innerDisc_arr[2][RR+1][ZZ]  *((  Rw)*(1-zw))+ */
    /* 	         innerDisc_arr[2][RR][ZZ+1]  *((1-Rw)*(  zw))+ */
    /* 	         innerDisc_arr[2][RR+1][ZZ+1]*((  Rw)*(  zw))*rho/pa[iqd]; */
    /* } */
    
    
  }
  
  usrqp(pa,qa);
  return;
}






/*
+ usrss - sound speed
*/
float usrss(float d, float p)
{
  return (sqrt(g*p/d));
}



/*
+ usrqp - conserved  variables from primitive
*/
void usrqp(float pa[], float qa[])
{
  float rho, divrho, ke, u;
  int i, id;
  rho = pa[iqd];
  qa[iqd] = rho;
  ke = 0.0;
  for (id = 1; id <= nd; id++){
    u = pa[iqu0+id];
    qa[iqu0+id] = u*rho;
    ke += u*u;
  }
  qa[iqe] = pa[iqe]/(g - 1.0) + 0.5*rho*ke;
  for (i = iqe+1; i <= iqal0 + scad; i++)
    qa[i] = rho*pa[i];
  for (i = iqal0 + scad + 1; i < nq; i++)
    qa[i] = pa[i];
  if (inten)
    qa[iqie] = rho*pa[iqie];
  return;
}



/*
+ usrpq - primitive variables from conserved
qa[iqd] density of fluid
qa[iqu0+1] x momentum of fluid
qa[iqu0+2] y momentum of fluid
qa[iqu0+3] z momentum of fluid
qa[iqe] total energy of fluid per unit volume (internal + kinetic)
qa[iqie] internal energy per unit volume (only if INTEN command is invoked)
qa[iqal0+1] mass of scalar 1 per unit volume (i.e. density times_ pa[iqal0+1]) 
qa[iqal0+2] value of scalar 2 etc
*/
void usrpq(float qa[], float pa[])
{
  float rho, divrho, ke, u;
  int i, id;
  rho = max(qa[iqd], dfloor);
  divrho = 1.0/rho;
  pa[iqd] = rho;
  ke = 0.0;
  for (id = 1; id <= nd; id++){
    u = qa[iqu0+id]*divrho;
    pa[iqu0+id] = u;
    ke += u*u;
  }
  ke *= 0.5*rho;
  pa[iqe] = (g - 1.0)*(qa[iqe] - ke);
  for (i = iqe + 1; i <= iqal0 + scad; i++)
    pa[i] = qa[i]*divrho;
  for (i = iqal0 + scad + 1; i < nq; i++)
    pa[i] = qa[i];
  if (fxp)
    usrfxp(pa);
  else
    pa[iqe] = max(pfloor,pa[iqe]);
  if (inten)
    pa[iqie] = qa[iqie]*divrho;
  return;
}



/*
+ usrfxq - Fix up solution.
     Description :
        This fixes the  solution in a cell to  ensure energy, pressure
        and density are positive. It  is given the cell coordinates in
        xc,  so  that   it  can  also  be  used   for  other  purposes
        e.g. forcing  the solution  to take a  given value  at certain
        points.
*/
void usrfxq(struct cell * cptr)
{
  int i, id;
  float ke, rho, divrho, mu;
  double r =sqrt(cptr->xc[0]*cptr->xc[0]+cptr->xc[1]*cptr->xc[1]);
  double z =cptr->xc[1], R=cptr->xc[0];
  double st = R/r;
  double ct = sqrt(1-st*st);
  double cs=cs0*pow(R*sqrt(r_star/r_sol/5.5), 0.5 * discP); 
  if (R<1) cs=cs0*pow(1.0/5, 0.5 * discP);
  double H=cs/sqrt(MG/pow(R,3));
  int ZZ,RR;


  RR = (R/10.0*128 + 0.5); 
  ZZ = (z/10.0*128 + 0.5); 
  
  if (RR>127) RR=127;
  if (ZZ>127) ZZ=127;

  if (flatDisc || (strcmp(cas, "STARWIND2") == 0) && H>0.05) H=0.05;
  if (R<discInnerRadius) H=0; 
  /* Density */
  
  rho = max(cptr->qa[iqd], dfloor);
  divrho = 1.0/rho;
  cptr->qa[iqd] = rho;
  if (strcmp(cas, "STARWIND2") == 0 && r<1){
    cptr->qa[iqd] = rho0;
    cptr->qa[iqe] = rho0*cs0*cs0;
  }

  if (strcmp(cas, "LARGEDISCSTAR") == 0){
    RR = (acos(st)/(pi/2)*905);
    if (RR>904) RR=904;
    if (insDisc && r<insDiscRadius){
      cptr->pa[iqd]=largeDisc[0][RR];
      cptr->pa[iqe]=largeDisc[1][RR];
      cptr->pa[iqu0+1]=largeDisc[2][RR];
      cptr->pa[iqu0+2]=largeDisc[3][RR];
      cptr->pa[iqal0+1]=largeDisc[4][RR];
      usrqp(cptr->pa, cptr->qa);
    } else if (z<(0.5*H)){
      double gaussian=exp(-pow(z,2)/2/H/H);
      double rho=rho0*pow(R*sqrt(r_star/r_sol/5.5), 3*discP/2) * gaussian;
      cptr->pa[iqd]=rho;
      cptr->pa[iqe]=rho*cs*cs;
      cptr->pa[iqu0+1]=0;
      cptr->pa[iqu0+2]=0;
      cptr->pa[iqal0+1]=sqrt( R * MG )*st;  
      usrqp(cptr->pa, cptr->qa);
    }
  }
  
    
  /* force solution to be spherical CAK solution up to fixOuterRadius where the disc isnt present */
  if (strcmp(cas, "DISCSTAR") == 0){
    if (r<fixOuterRadius || (R>1 &&z<=H &&!flatDisc) || (0)){
      
      double gaussian=exp(-pow(z,2)/2/H/H);
      if (flatDisc ) {
	if (z<0.02) gaussian=1;
	else gaussian=1e-10;
      }
      RR = (r/10*2048 + 0.5);
      if (RR>2047) RR=2047;
      
      cptr->qa[iqd]=star_arr[RR][1];
      cptr->pa[iqe]=cptr->qa[iqd]*cs*cs;
      if (radflag && L_star > L_sol){
	cptr->qa[iqu0+1]=star_arr[RR][2]*st*star_arr[RR][1];
	cptr->qa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st)*star_arr[RR][1];
	cptr->pa[iqu0+1]=star_arr[RR][2]*st;
	cptr->pa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st);
	cptr->qa[iqe] = cptr->pa[iqd]*cs*cs/(g - 1.0) + 0.5*star_arr[RR][1]*pow(star_arr[RR][2],2);
      }else{
	cptr->pa[iqu0+1]=0;
	cptr->pa[iqu0+2]=0;
	cptr->qa[iqu0+1]=0;
	cptr->qa[iqu0+2]=0;
      }
      
      if (R >= discInnerRadius && z <(8*H)){
	rho=rho0*pow(R, 3*discP/2) * gaussian;
	if (rho>cptr->qa[iqd]){
	  cptr->qa[iqd]=rho;
	  cptr->pa[iqd]=rho;
	  cptr->pa[iqe]=rho*cs*cs;
	}
      }
      if(cptr->pa[iqd]<rho_amb) {
	cptr->pa[iqd]+=rho_amb;
	cptr->pa[iqe]+=rho_amb*cs*cs;
      }
      cptr->qa[iqal0+1]=sqrt( R * MG )*st*cptr->qa[iqd];  
      /*  if (z>(6*H)){  */
      /* 	H*=2; */
      /* 	gaussian=exp(-pow(z,2)/2/H/H); */
      /* 	cptr->pa[iqal0+1]*=gaussian; */
      /* } */
      
    }
    if (st<0.1){
      double temp=sqrt(cptr->qa[iqu0+1]*cptr->qa[iqu0+1]+cptr->qa[iqu0+2]*cptr->qa[iqu0+2]);
      cptr->qa[iqu0+1]=temp*st; cptr->qa[iqu0+2]=temp*ct;
    }
      
    /* cptr->qa[iqd]=rho0; */
    /*   cptr->pa[iqe]=cptr->qa[iqd]*cs0*cs0; */
    /*   cptr->qa[iqu0+1]=star_arr[RR][2]*st*cptr->qa[iqd]; */
    /*   cptr->qa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st)*cptr->qa[iqd]; */
    /*   cptr->qa[iqe] = cptr->pa[iqe]/(g - 1.0) + 0.5*cptr->qa[iqd]*pow(star_arr[RR][2],2); */
    /* } */
    /* if(){ */
    /*   cptr->qa[iqd]=star_arr[RR][1]; */
    /*   cptr->pa[iqe]=cptr->qa[iqd]*cs0*cs0; */
    /*   cptr->qa[iqu0+1]=star_arr[RR][2]*st*cptr->qa[iqd]; */
    /*   cptr->qa[iqu0+2]=star_arr[RR][2]*sqrt(1-st*st)*cptr->qa[iqd]; */
    /*   cptr->qa[iqe] = cptr->pa[iqe]/(g - 1.0) + 0.5*cptr->qa[iqd]*pow(star_arr[RR][2],2); */
    /*   cptr->qa[iqal0+1]=st*sqrt( MG )*cptr->qa[iqd]; 
    /*   if (R<0.1) cptr->pa[iqal0+1]=1e-30*cptr->qa[iqd]; */
    /*   cptr->pa[iqal0+1]=cptr->qa[iqd]*sqrt( R * MG )*st;  
    /*   if (z>(6*H)){  */
    /* 	H*=2; */
    /* 	cptr->qa[iqal0+1]*=exp(-pow(z-3*H,2)/2/H/H); */
    /*   } */
    /* } */
    /* if (r<fixOuterRadius && z <=(6*H)){ */
    /*   rho=rho0*pow(R, 3*discP/2) * gaussian; */
    /*   if (rho>cptr->qa[iqd]){ */
    /* 	cptr->qa[iqd]=rho; */
    /* 	cptr->qa[iqu0+1]=0; */
    /* 	cptr->qa[iqu0+2]=0; */
    /* 	cptr->qa[iqal0+1]=st*sqrt(R*MG)*cptr->qa[iqd]; */
    /* 	if (z>(6*H)){  */
    /* 	  cptr->pa[iqal0+1]*=exp(-pow(z-6*H,2)/2/pow(H*2,2)); */
    /* 	}       */
    /*   } */
    /* }else if (){ */
      /* usrqis(cptr->xc, cptr->qa); */
      /* cptr->qa[iqd]=rho0*pow(R, 3*discP/2) * gaussian; */
      /* cptr->qa[iqu0+1]=0; */
      /* cptr->qa[iqu0+2]=0; */
      /* cptr->qa[iqal0+1]=st*sqrt(R*MG )*cptr->qa[iqd]; */
    /* } */
  }

  

  if (isnan(cptr->qa[iqu0+1])) cptr->qa[iqu0+1]=0;
  if (isnan(cptr->qa[iqu0+2])) cptr->qa[iqu0+2]=0;
  if (isnan(cptr->qa[iqd])) cptr->qa[iqd]=dfloor;
  if (isnan(cptr->qa[iqe])) cptr->qa[iqe]=pfloor;

   /* Thermal energy */
  ke = 0.0;
  for (id = 1; id <= nd; id++){
    mu = cptr->qa[iqu0+id];
    ke += mu*mu;
  }
  ke *= 0.5*divrho;
  if (cptr->qa[iqe] - ke <= pfloor/(g - 1.0))
    cptr->qa[iqe] = ke + pfloor/(g - 1.0);
  if (keps){
    cptr->qa[iqk] = max(cptr->qa[iqk], rho*kmin);
    cptr->qa[iqeps] = max(cptr->qa[iqeps], rho*epsmin);
  }
  /* Scalars */
  for (i = iqal0+1; i <= iqal0+scad; i++)
    cptr->qa[i] = min(cptr->qa[iqd], max(cptr->qa[i], 0.0));
  return;
}



/*
+ usrfxp - Fix up primitive solution.
     Description :
        This  fixes the  solution in  pa to  ensure that  pressure and
        density are positive.
*/
void usrfxp(float pa[])
{
  int i;
  
  /* Fix density */
  pa[iqd] = max(pa[iqd], dfloor);
  /* Fix pressure */
  pa[iqe] = max(pa[iqe], pfloor);
  if (isnan(pa[iqu0+1])) pa[iqu0+1]=0;
  if (isnan(pa[iqu0+2])) pa[iqu0+2]=0;
  if (keps){
    pa[iqk] = max(pa[iqk], kmin);
    pa[iqeps] = max(pa[iqeps], epsmin);
  }
  /* Scalars */
  for (i = iqal0+1; i <= iqal0+scad; i++)
    pa[i] = min(1.0, max(pa[i], 0.0));
  return;

}



/* 
+ usrps - construct plot scalar
     Description :
        This  returns  the value  of  a plot  scalar  in  s given  the
        conserved solution  in qa. The  scalar is defined by  ITEM and
        returns TRUE if the plot item is found, FALSE otherwise.
*/
bool usrps(cell *cptr, float &s)
{
  int i, j, id;
  bool found;
  float qref[MAXQ], pa[MAXQ], tv[MAXQ], temp;
  struct cell *left, *right;
  found = false;
  if (strncmp(item, "QA", 2) ==0){
    sscanf(&item[2], "%d", &i);
    if ((i < 1) || (i > nq))
      found = false;
    else{
      s = cptr->qa[i-1];
      found = true;
    }
  }
  else if (strncmp(item, "QR", 2) ==0){
    sscanf(&item[2], "%d", &i);
    if ((i < 1) || (i > nq))
      found = false;
    else{
      s = cptr->dqdt[i-1]/cptr->vol;
      found = true;
    }
  }
  else if (strcmp(item, "RHO") == 0){
    s = cptr->qa[iqd];
    found = true;
  }
  else if (strcmp(item, "U1") == 0){
    s = cptr->qa[iqu0+1]/cptr->qa[iqd];
    found = true;
  }
  else if (strcmp(item, "U2") == 0){
    if (nd >= 2){
      s = cptr->qa[iqu0+2]/cptr->qa[iqd];
      found = true;
    }
    else
      found = false;
  }
  else if (strcmp(item, "U3") == 0){
    if (nd == 3){
      s = cptr->qa[iqu0+3]/cptr->qa[iqd];
      found = true;
    }
    else
      found = false;
  }
  else if (strcmp(item, "DU01") == 0){
    s = cptr->dudx[0][1];
    found = true;
  }
  else if (strcmp(item, "DU10") == 0){
    s = cptr->dudx[1][0];
    found = true;
  }
  else if (strcmp(item, "DU00") == 0){
    s = (cptr->pa[iqu0+1]-cptr->rjoin[0]->right->pa[iqu0+1])/(cptr->xc[1]-cptr->rjoin[0]->right->xc[0]);
    found = true;
  }
  else if (strcmp(item, "DU11") == 0){
    s = cptr->dudx[1][1];
    found = true;
  }
  else if (strcmp(item, "UT") == 0){
    s = 0.0;
    for (id = 1; id <= nd; id++){
      temp = cptr->qa[iqu0+id]/cptr->qa[iqd];
      s += temp*temp;
    }
    s = sqrt(s);
    found = true;
  }
  else if (strcmp(item, "MOMT") == 0){
    s = 0.0;
    for (id = 1; id <= nd; id++){
      temp = cptr->qa[iqu0+id];
      s += temp*temp;
    }
    s = sqrt(s);
    found = true;
  }
  else if (strcmp(item, "KE") == 0){
    s = 0.0;
    for (id = 1; id <= nd; id++){
      temp = cptr->qa[iqu0+id];
      s += temp*temp;
    }
    s *= 0.5/cptr->qa[iqd];
    found = true;
  }

  else if (strcmp(item, "MASSFLOW") == 0){
    s = 0.0;
    for (id = 1; id <= nd; id++){
      temp = cptr->qa[iqu0+id];
      s += temp*temp;
    }
    s = sqrt(s)*(cptr->xc[0]*cptr->xc[0]+cptr->xc[1]-cptr->xc[1])/m_sol*3600*24*365.25;
    found = true;
  }
  else if (strcmp(item, "POT") == 0){
    if (grav){
      s = gconst*cptr->qa[iqpot];
      found = true;
    }
    else
      found = false;
  }
  else if (strcmp(item, "G") == 0){
    double R=cptr->xc[0], z=cptr->xc[1];
    double r=sqrt(R*R+z*z); 
    s=MG/r/r;
      found = true;
  }
  else if (strncmp(item, "GA", 2) ==0){
    sscanf(&item[2], "%d", &id);
    if ((!grav) || (id < 1) || (id > nd))
      found = false;
    else{
      s = gconst*cptr->pa[iqu0+id];
      found = true;
    }
  }
  else if (strcmp(item, "QSOL") == 0){
    s = cptr->qsol;
    found = true;
  }
  else if (strcmp(item, "ICELL") == 0){
    s = cptr->icell;
    found = true;
  }
  else if (strcmp(item, "PROC") == 0){
    s = (float)mype;
    found = true;
  }
  else if (strcmp(item, "WORK") == 0){
    i = (int)((cptr->xc[nd-1] - xl[nd-1])/lev[nlevs-1]->dxbar[nd-1]);
    s = (float)wka[i];
    found = true;
  }
  else{
    usrpq(cptr->qa, pa);
    if (strncmp(item, "PA", 2) ==0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nq))
      found = false;
      else{
	s = pa[i-1];
	found = true;
      }
    }
    else if (strncmp(item, "AL", 2) ==0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nal))
      found = false;
      else{
	s = pa[iqal0+i];
	found = true;
      }
    }
    else if (strncmp(item, "AD", 2) ==0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nal))
      found = false;
      else{
	s = pa[iqal0+i]/pa[iqd];
	found = true;
      }
    }
    else if (strcmp(item, "PG") == 0){
      s = pa[iqe];
      found = true;
    }
    else if (strcmp(item, "PGINT") == 0){
      if (inten){
	s = pa[iqie]*powf(pa[iqd], g);
	found = true;
      }
      else
	found = false;
    }
    else if (strcmp(item, "ENT") == 0){
      s = pa[iqe]/powf(pa[iqd], g);
      found = true;
    }
    else if (strcmp(item, "ENTIE") == 0){
      if (inten){
	s = pa[iqie];
	found = true;
      }
    }
    else if (strcmp(item, "C") == 0){
      s = usrss(pa[iqd], pa[iqe]);
      found = true;
    }
    else if (strcmp(item, "MACH") == 0){
      s = 0.0;
      for (id = 1; id <= nd; id++)
	s += pa[iqu0+id]*pa[iqu0+id];
      s = sqrt(s)/usrss(pa[iqd], pa[iqe]);
      found = true;
    }
    else if (strcmp(item, "TE") == 0){
      s = usrte(pa);
      found = true;
    }
    else if (strcmp(item, "K") == 0){
      if (!keps)
	found = false;
      else{
	s = pa[iqk];
	found = true;
      }
    }
    else if (strcmp(item, "EPS") == 0){
      if (!keps)
	found = false;
      else{
	s = pa[iqeps];
	found = true;
      }
    }
    else if (strcmp(item, "TL") == 0){
      if (!keps)
	found = false;
      else{
	s = ceps*pa[iqk]*sqrt(pa[iqk])/max(pa[iqeps],epsmin);
	found = true;
      }
    }
    else if (strncmp(item, "MU", 2) ==0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nq))
      found = false;
      else{
	usrmu(pa, tv);
	s = tv[i];
	found = true;
      }
    }
    else if (strncmp(item, "SA", 2) ==0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nq))
      found = false;
      else{
	for (j = 0; j < nq; j++)
	  cptr->pa[j] = pa[j];
	usrsaf(0.0, cptr, tv);
	s = tv[i-1];
	found = true;
      }
    }
    else if (strncmp(item, "SD", 2) ==0){
      sscanf(&item[2], "%d", &i);
      if ((i < 1) || (i > nq))
      found = false;
      else{
	for (j = 0; j < nq; j++)
	  cptr->pa[j] = pa[j];
	usrsaf(0.0, cptr, tv);
	s = tv[i-1]/cptr->pa[iqd];
	found = true;
      }
    }
    else if (strcmp(item, "lvl") ==0){
      s = 1.0*cptr->icell;
      
    }
  }
  return(found);
}



/*
+ usrsaf - source terms.

so[iqd]____ mass source for fluid
so[iqu0+1]_ x momentum source for fluid
so[iqu0+2]_ y momentum source for fluid
so[iqu0+3]_ z momentum source for fluid
so[iqe]____ energy source for fluid
so[iqal0+1] source for scalar 1
so[iqal0+2] source for scalar 2 etc

*/
void usrsaf(float dt, cell *cptr, float so[])
{
  int i, j;
  float rc, dudx[MAXD][MAXD], pa[MAXQ], mu[MAXQ], vcof, tv, divv, pp;
  bool rzp, zrp, rtp;
  /* Zero so (safety) and copy primitives */
  for (i = 0; i < nq; i++){
    so[i] = 0.0;
    pa[i] = cptr->pa[i];
  }
  /* Geometry */
  rzp = ((strcmp(geom, "RZP") == 0) && (strcmp(coord[0], "R") == 0));
  zrp = ((strcmp(geom, "RZP") == 0) && (strcmp(coord[1], "R") == 0));
  rtp = (strcmp(geom, "RTP") == 0);
  /* Geometric source term (if any) */
  if (rzp){
    /* Cylindrical polars xc[0] = R */
    rc = cptr->xc[0];
    so[iqu0+1] += pa[iqe]/(rc+1/32);
  }
  else if (zrp){
    /* Cylindrical polars xc[1] = R */
    rc = cptr->xc[1];
    so[iqu0+2] += pa[iqe]/(rc+1/32);
  }
  else if (rtp){
    /* Spherical polars */
    rc = cptr->vol/(cptr->rjoin[0]->area - cptr->ljoin[0]->area);
    so[iqu0+1] += pa[iqe]/(rc+1/32);
  }
  if (viscous){
    /* Diffusion coefficients etc */
    usrmu(pa, mu);
    vcof = mu[iqu0+1];
    divv = 0.0;
    /*Cartesian divergence */
    for (i = 0; i < nd; i++){
      divv += cptr->dudx[i][i];
    }
    /* Geometric source term (if any) */
    if (rzp){
      /* Cylindrical polars xc[0] = R */
      divv += pa[iqu0+1]/(rc+1/32);
      so[iqu0+1] += 2.0*vcof*(divv/3.0 - pa[iqu0+1]/(rc+1/32))/(rc+1/32);
    }
    else if (zrp){
      /* Cylindrical polars xc[1] = R */
      divv += pa[iqu0+2]/(rc+1/32);
      so[iqu0+2] += 2.0*vcof*(divv/3.0 - pa[iqu0+2]/(rc+1/32))/(rc+1/32);
    }
    else if (rtp){
      /* Spherical polars */
      divv += 2.0*pa[iqu0+1]/(rc+1/32);
      so[iqu0+1] += 4.0*vcof*(divv/3.0 - pa[iqu0+1]/(rc+1/32))/(rc+1/32);
    }
    if (keps){
      /* Velocity gradients */
      for (i = 0; i < nd ; i++)
	for (j = 0; j < nd; j++)
	  dudx[i][j] = cptr->dudx[i][j];
      /* Production term */
      if (rzp){
	/* Turbulent pressure r-momentum source */
	so[iqu0+1] += 2.0*pa[iqd]*pa[iqk]/(3.0*rc);
	/* Geometric production term */
	tv = pa[iqu0+1]/(rc+1/32);
	pp = 2.0*tv*tv;
      }
      else if (zrp){
	/* Turbulent pressure r-momentum source */
	so[iqu0+2] += 2.0*pa[iqd]*pa[iqk]/(3.0*rc);
	/* Geometric production term */
	tv = pa[iqu0+2]/(rc+1/32);
	pp = 2.0*tv*tv;
      }
      else if (rtp){
	/* Turbulent pressure r-momentum source */
	so[iqu0+1] += 4.0*pa[iqd]*pa[iqk]/(3.0*rc);
	/* Geometric production term */
	tv = pa[iqu0+1]/(rc+1/32);
	pp = 4.0*tv*tv;
      }
      else
	pp = 0.0;
      for (i = 0; i < nd; i++){
	for (j = 0; j < nd; j++)
	  pp += dudx[i][j]*(dudx[i][j] + dudx[j][i]);
      }
      pp = vcof*pp - 2.0*divv*(pa[iqd]*pa[iqk] + divv*vcof)/3.0;
      /* k source */
      so[iqk] += pp - pa[iqd]*pa[iqeps];
      /* Energy source */
      so[iqe] += (pa[iqd]*pa[iqeps] - pp);
      /* Epsilon source */
      so[iqeps] += pa[iqeps]*(ce1*pp - pa[iqd]*ce2*pa[iqeps])/pa[iqk]; 
    }
  }

  double R=cptr->xc[0], z=cptr->xc[1];
  double modR=fabs(R)+0.01;     
  double r=sqrt(R*R+z*z); 
  double st=R/r; 
  double ct=sqrt(1-st*st);
  double vR=cptr->pa[iqu0+1];
  double vz=cptr->pa[iqu0+2];
  double H=cs0/sqrt(MG/pow(R,3));
  double gg=0;

  
  
/* gravity source term */

 if(strcmp(cas, "RTINST") == 0){
   so[iqu0+2] = -pa[iqd];

 }
 else if (strcmp(cas, "STARWIND") == 0 || strcmp(cas, "STARWIND2") == 0 ){
   if (strcmp(cas, "STARWIND2") == 0){ R=r; }
   if (R>1) {
     if (gravflag){
       gg =MG/R/R;
     } else {
     gg=0;
     }
     double a,b,c,d ,  theta,phi , IR, area, Multiplier, gradV, vtherm;
     int RR=R*(array_R_size-1)/array_R_max;
     if (RR > (array_R_size-2)) RR=(array_R_size-2);
     
     if (cptr->rjoin[0]->right != NULL && cptr->ljoin[0]->left == NULL){
       gradV=(cptr->rjoin[0]->right->pa[iqu0+1]-cptr->pa[iqu0+1])/(cptr->rjoin[0]->right->xc[0]-cptr->xc[0]);
     } else if (cptr->ljoin[0]->left != NULL && cptr->rjoin[0]->right ==NULL){
       gradV=(cptr->ljoin[0]->left->pa[iqu0+1]-cptr->pa[iqu0+1])/(cptr->ljoin[0]->left->xc[0]-cptr->xc[0]);
     } else{
       gradV=((cptr->ljoin[0]->left->pa[iqu0+1]-cptr->pa[iqu0+1])/(cptr->ljoin[0]->left->xc[0]-cptr->xc[0])+
	      (cptr->rjoin[0]->right->pa[iqu0+1]-cptr->pa[iqu0+1])/(cptr->rjoin[0]->right->xc[0]-cptr->xc[0]))/2;
     }
     
     
     IR=L_star; 
     
     area=0;
     double MU,vr,muMin=cos(atan(1/R));

       if(strcmp(cas, "STARWIND2") == 0) {
	 vr=fabs(pa[iqu0+1]*st+pa[iqu0+2]*ct);
       } else {
	 vr=fabs(pa[iqu0+1]);
       }
     d=4*pi*pow(R*r_star,2);  
     

     if (gradV <1e-10) gradV=1e-10;
     double multi=sigma_e/MMM * cs0*0.3 *pa[iqd] /r_star/r_star; 
     Multiplier =1+k * pow(multi/fabs(gradV), -alpha); 
     
     for (i=0;i<20;i++){
       MU=1-(i*muMin/19.0);
       area+=pow(((1-MU*MU)*vr/R+MU*MU*gradV)/gradV, alpha)*MU*(1-muMin)/20.0;
     }
     Multiplier*=area*(2/(1-muMin*muMin));
     if (Multiplier>Mmax) Multiplier=Mmax;
     
     if (strcmp(cas, "STARWIND2") == 0){
       so[iqu0+1]+=pa[iqd] * L_star/d*forceScalingFactor*Multiplier*st;
       so[iqu0+2]+=pa[iqd] * L_star/d*forceScalingFactor*Multiplier*ct;
       so[iqe]   +=pa[iqd] * forceScalingFactor*Multiplier*(st*pa[iqu0+1]+ct*pa[iqu0+2])*L_star/d;
     } else {
       so[iqu0+1]+=pa[iqd] * (L_star/d*forceScalingFactor*Multiplier-gg);
       so[iqe]   +=pa[iqd] * (L_star/d*forceScalingFactor*Multiplier-gg) * pa[iqu0+1];
     }

   /*     
   /*     Multiplier=1+k *pow(sigma_e/MMM*cs0/3.4*pa[iqd]/fabs(gradV)/r_star/r_star, -alpha); */
   /*     
   /*     if (Multiplier >Mmax) { */
   /* 	 Multiplier=Mmax; */
   /*     } */
   /*     cptr->radR=IR/(4*pi*R*R*r_star*r_star)*sigma_e/3e8/MMM/r_star * Multiplier; */
   /*     if (strcmp(cas, "STARWIND2") == 0) { */
   /* 	 so[iqu0+1]+=pa[iqd]*st*(cptr->radR-MG/r/r); */
   /* 	 so[iqu0+2]+=pa[iqd]*ct*(cptr->radR-MG/r/r); */
   /* 	 so[iqe]   +=pa[iqd]*(pa[iqu0+1] * st*(cptr->radR-MG/r/r) + pa[iqu0+2] * ct*(cptr->radR-MG/r/r)); */
   /*     }else{ */
   /* 	 so[iqu0+1]+=pa[iqd]*(cptr->radR-MG/R/R); */
   /* 	 so[iqe]   +=pa[iqd]*(st*pa[iqu0+1]*(cptr->radR-MG/R/R)); */
   /*     } */
   /*   } else { */
   /*     if (strcmp(cas, "STARWIND2") == 0) { */
   /* 	 so[iqu0+1]+=st*(pa[iqd]*cptr->radR-pa[iqd]*MG/r/r); */
   /* 	 so[iqu0+2]+=ct*(pa[iqd]*cptr->radR-pa[iqd]*MG/r/r); */
   /* 	 so[iqe]+=(so[iqu0+1] * st*(cptr->radR-MG/r/r) + so[iqu0+2] * ct*(cptr->radR-MG/r/r)); */
   /*     }else{ */
   /* 	 so[iqu0+1]+=pa[iqd]*cptr->radR-pa[iqd]*MG/R/R; */
   /* 	 so[iqe]+=st*pa[iqu0+1]*(pa[iqd]*cptr->radR-pa[iqd]*MG/R/R); */
   /*     } */
   /*   } */
   /* } */
   /* cptr->numtimesteps+=1; */
   }
 }
 else if (strcmp(cas, "DISCSTAR2") == 0){
   gg =MG/r/r;
   if (gravflag){
     if (!(R>0 && z>0 && r>=R && r>=z && r>0 && ct>=0 && st>=0 && st<=1 && ct<=1)){
       fprintf(stderr,"%.2f %.2f   %.3e %.3e   %.2f %.2f \n",R, z, pow(pa[iqal0+1]/modR,2)/modR,gg,ct,st);
     }
     
     so[iqu0+1]+= -pa[iqd] *  (gg*st-pow(pa[iqal0+1]/modR,2)/modR);  
     so[iqu0+2]+= -pa[iqd] *  gg*ct;                                 
     so[iqe]   += -gg*pa[iqd]*(pa[iqu0+1]*st+pa[iqu0+2]*ct);         
   }
 }

 else if (strcmp(cas, "LARGEDISCSTAR") == 0){
   gg =MG/r/r;

   so[iqu0+1]+= -pa[iqd] * (gg*st-pow(pa[iqal0+1]/modR,2)/modR);
   so[iqu0+2]+= -pa[iqd] *  gg *ct;                                 
   so[iqe]   += -gg*pa[iqd]*(pa[iqu0+1]*st+pa[iqu0+2]*ct);         
   /* if (!(R>0 && z>0 && r>=R && r>=z && r>0 && ct>=0 && st>=0 && st<=1 && ct<=1)){ */
   /*   fprintf(stderr,"%.2f %.2f   %.3e %.3e   %.2f %.2f \n",R, z, pow(pa[iqal0+1]/modR,2)/modR,gg,ct,st); */
   /* } */

   if (R<9 && z<9) {
     so[iqu0+1]=0;
     so[iqu0+2]=0;
     so[iqe]   =0;
   }
   else if (gravflag){

     double gradvzz=cptr->dudx[1][1],
       gradvRR=cptr->dudx[0][0],
       gradvzR=cptr->dudx[1][0],
       gradvRz=cptr->dudx[0][1];
     

     
     double d=4*pi*pow(r*r_star,2);  
     
     double gradvl=st*(st*fabs(gradvRR)+ct*fabs(gradvRz))+ct*(st*fabs(gradvzR)+ct*fabs(gradvzz));
     if (gradvl <1e-10) gradvl=1e-10;
     double Multiplier =k * pow(sigma_e/MMM * cs0*pow(R/5, 0.5 * -3/4)*0.3 *pa[iqd] /r_star/r_star /fabs(gradvl), -alpha); 
     if (Multiplier>Mmax) Multiplier=Mmax;     

     gg=((1+ct)*0.5*L_star + ct*6.67e-11*M_star*M_acc/2/r_star)/d*forceScalingFactor*Multiplier; 

     so[iqu0+1]+= pa[iqd] * gg * st;
     so[iqu0+2]+= pa[iqd] * gg * ct;                                 
     so[iqe]   += gg*pa[iqd]*(pa[iqu0+1]*st+pa[iqu0+2]*ct);         
   }
 }


 else if (strcmp(cas, "DISCSTAR") == 0 && (r>fixOuterRadius || (r>1 && z<(6*H)) )){
   if (!(r<1 || pa[iqd]<2e15)){
     gg =MG/r/r;
     if (gravflag){
       /* if (!(R>0 && z>0 && r>=R && r>=z && r>0 && ct>=0 && st>=0 && st<=1 && ct<=1)){ */
       /*   fprintf(stderr,"%.2f %.2f   %.3e %.3e   %.2f %.2f \n",R, z, pow(pa[iqal0+1]/modR,2)/modR,gg,ct,st); */
       /* } */
       
       so[iqu0+1]+= -pa[iqd] * (gg*st-pow(pa[iqal0+1]/modR,2)/modR);   
       so[iqu0+2]+= -pa[iqd] *  gg*ct;                                 
       so[iqe]   += -gg*pa[iqd]*(pa[iqu0+1]*st+pa[iqu0+2]*ct);         
     }
     
     /*    /\* /\* end gravity source tem  *\/ *\/ */
     
     
     /*    /\* /\* rad drving force term  *\/ *\/ */
     if (radflag){
       double gradvzz=cptr->dudx[1][1],
       	 gradvRR=cptr->dudx[0][0],
       	 gradvzR=cptr->dudx[1][0],
       	 gradvRz=cptr->dudx[0][1];
	 /* double gradvzz=cptr->dudx[1][1], */
       	 /* gradvRR=cptr->dudx[0][0], */
       	 /* gradvzR=cptr->dudx[1][0], */
       	 /* gradvRz=cptr->dudx[0][1]; */ 

       if ((cptr->rjoin[0]->right)!=NULL){
       	 gradvRR=(cptr->rjoin[0]->right->pa[iqu0+1]-cptr->pa[iqu0+1])/(cptr->rjoin[0]->right->xc[0]-cptr->xc[0]);
       	 gradvzR=(cptr->rjoin[0]->right->pa[iqu0+2]-cptr->pa[iqu0+2])/(cptr->rjoin[0]->right->xc[0]-cptr->xc[0]);
       }else {
       	 gradvRR=(cptr->pa[iqu0+1]-cptr->ljoin[0]->left->pa[iqu0+1])/(cptr->xc[0]-cptr->ljoin[0]->left->xc[0]);
       	 gradvzR=(cptr->pa[iqu0+2]-cptr->ljoin[0]->left->pa[iqu0+2])/(cptr->xc[0]-cptr->ljoin[0]->left->xc[0]);
       }
       if ((cptr->rjoin[1]->right)!=NULL){
       	 gradvRz=(cptr->rjoin[1]->right->pa[iqu0+1]-cptr->pa[iqu0+1])/(cptr->rjoin[1]->right->xc[1]-cptr->xc[1]);
       	 gradvzz=(cptr->rjoin[1]->right->pa[iqu0+2]-cptr->pa[iqu0+2])/(cptr->rjoin[1]->right->xc[1]-cptr->xc[1]);
       }else {
       	 gradvRz=(cptr->pa[iqu0+1]-cptr->ljoin[1]->left->pa[iqu0+1])/(cptr->xc[1]-cptr->ljoin[1]->left->xc[1]);
       	 gradvzz=(cptr->pa[iqu0+2]-cptr->ljoin[1]->left->pa[iqu0+2])/(cptr->xc[1]-cptr->ljoin[1]->left->xc[1]);
       }

       double multi=sigma_e/MMM * cs0*0.3 *pa[iqd] /r_star/r_star; 

       if (r>10 && (z>10 || R>10)){
	 double d=4*pi*pow(r*r_star,2);  
	 double gradvl=st*(st*fabs(gradvRR)+ct*fabs(gradvRz))+ct*(st*fabs(gradvzR)+ct*fabs(gradvzz));	 
	 if (fabs(gradvl) <1e-10) gradvl=1e-10;

	 double Multiplier =k * pow(multi/fabs(gradvl), -alpha); 
	 
	 const double x= L_star/(6.67e-11*M_star*M_acc/(2*r_star));
	 double L= (1+ct)/2*L_star+ (3*6.67e-11*M_star*M_acc / (8*pi*pi*pow(r_star,3))  * (0.1+(4-pi)*x/6/pi))*ct/2; 

	 so[iqu0+1]+=pa[iqd] * L/d*forceScalingFactor*Multiplier*st;
	 so[iqu0+2]+=pa[iqd] * L/d*forceScalingFactor*Multiplier*ct;
	 so[iqe]   +=pa[iqd] * forceScalingFactor*Multiplier*(st*pa[iqu0+1]+ct*pa[iqu0+2])*L/d;
       }else {


       
       
       
	 const static double cosp[15]={1.00000000e+00,  9.74927912e-01,  9.00968868e-01,  7.81831482e-01,  6.23489802e-01,  
				       4.33883739e-01,  2.22520934e-01,  0.00000000e-0,  -2.22520934e-01, -4.33883739e-01, 
				       -6.23489802e-01, -7.81831482e-01, -9.00968868e-01, -9.74927912e-01, -1};
	 const static double sinp[15]={0.        ,  0.22252093,  0.43388374,  0.6234898 ,  0.78183148,
				       0.90096887,  0.97492791,  1.        ,  0.97492791,  0.90096887,
				       0.78183148,  0.6234898 ,  0.43388374,  0.22252093, 0}; 
	 const static double rd[15]={1.000001    ,   1.17876863,   1.38949549,   1.63789371,
				     1.93069773,   2.27584593,   2.6826958 ,   3.16227766,
				     3.72759372,   4.39397056,   5.17947468,   6.1054023 ,
				     7.19685673,   8.48342898,  10.       };
	 
	 static char buff[60];
	 snprintf(buff,52,"%.2e%.2e%+.2e%+.2e%+.2e%+.2e",R,z,gradvRR,gradvRz,gradvzR,gradvzz);
	 std::string s=buff;
	 if (caching && st<0.7 && radCacheR.find(s)!=radCacheR.end()){
	   multi = pow(pa[iqd], -alpha);
	   so[iqu0+1]=pa[iqd]*radCacheR[s]*multi;
	   so[iqu0+2]=pa[iqd]*radCachez[s]*multi;
	   so[iqe]=(radCacheR[s]*pa[iqu0+1]+radCachez[s]*pa[iqu0+2])*pa[iqd]*multi;
	 }
	 else{ 
	   int i,j,di, flag;
	   double IR,Iz, rr, a,b,c,d, area, A,B,C;
	   double gradvl,Multiplier = 0;
	   
	   int RR=R*(array_R_size-1)/array_R_max;
	   int ZZ=z*(array_z_size-1)/array_z_max; 
	   
	   
	   if (RR > (array_R_size-2)) RR=(array_R_size-2);
	   if (ZZ > (array_z_size-2)) ZZ=(array_z_size-2);
	   
	   
	   for (i=0;i<15;i++){
	     
	   
	     for (j=0;j<15;j++){
	       
	       /* disc force  */
	       a=R-rd[i]*cosp[j];    
	       c=rd[i]*sinp[j];      
	       d=sqrt(a*a+z*z+c*c);
	       a/=d;
	       b=z/d;              
	       
	       /* X= (rd[i]*cos(phi), rd[i]*sin(phi), 0)  */
	       /* r= (R, 0 ,Z) */
	       /* dX =(R-rd[i]*cos(phi), -rd[i]*sin(phi), Z)  */
	       /* (X_1+k*dX_1)^2+(X_2+k*dX_2)^2+(k*dX_3)^2 = 1 */
	       /* rearange -> |dX|^2 k^2 + 2(X.dX) k + X^2 - 1 = 0 */
	       /* star intersects the disc-wind vector if 0<k<1 */
	       
	       A=a*a+b*b+c*c;
	       B=2*(a*rd[i]*cosp[j]+c*rd[i]*sinp[j]);
	       C=rd[i]*rd[i]-1;
	       if (4*A*C > B*B) flag=1;
	       else C=sqrt(B*B-4*A*C); 
	       
	       if ((flag || !(-1e-3<(-B-C)/2/A<1) || !(-1e-3<(C-B)/2/A<1))){  
		 gradvl=fabs(a*(a*gradvRR+b*gradvRz)+b*(a*gradvzR+b*gradvzz));
		 if (gradvl <1e-10) gradvl=1e-10;
		 Multiplier = k * pow(multi/gradvl, -alpha); 
		 if (Multiplier>Mmax)  Multiplier=Mmax;
		 
		 IR=disc[RR][ZZ][i][j][0]*disc_albedo;
		 Iz=disc[RR][ZZ][i][j][1]*disc_albedo;
		 
		 so[iqu0+1]+=pa[iqd] * IR*forceScalingFactor*Multiplier;
		 so[iqu0+2]+=pa[iqd] * Iz*forceScalingFactor*Multiplier;
		 so[iqe]   +=pa[iqd] *    forceScalingFactor*Multiplier * 
		   (IR*pa[iqu0+1]+Iz*pa[iqu0+2]);
	       }
	       
	     }
	   }
	   /* end disc force */
	   
	   
	   /* star star force */
	   
	   
	   
	   double vr=fabs(pa[iqu0+1]*st+pa[iqu0+2]*ct);
	   d=4*pi*pow(r*r_star,2);  
	   
	   gradvl=st*(st*fabs(gradvRR)+ct*fabs(gradvRz))+ct*(st*fabs(gradvzR)+ct*fabs(gradvzz));
	   if (gradvl <1e-10) gradvl=1e-10;
	   Multiplier =k * pow(multi/fabs(gradvl), -alpha); 
	   
	   area=0;
	   double MU,muMin=cos(atan(1/r));
	   
	   for (i=0;i<20;i++){
	     MU=1-(i*muMin/19.0);
	     area+=pow(((1-MU*MU)*vr/r+MU*MU*gradvl)/gradvl, alpha)*MU*(1-muMin)/20.0;
	   }
	   Multiplier*=area*(2/(1-muMin*muMin));
	   if (Multiplier>Mmax) Multiplier=Mmax;
	 
	   so[iqu0+1]+=pa[iqd] * (1+ct)*0.5*L_star/d*forceScalingFactor*Multiplier*st;
	   so[iqu0+2]+=pa[iqd] * (1+ct)*0.5*L_star/d*forceScalingFactor*Multiplier*ct;
	   so[iqe]   +=pa[iqd] * forceScalingFactor*Multiplier*(st*pa[iqu0+1]+ct*pa[iqu0+2])*L_star/d;
	   
	   multi=pow(pa[iqd], 1-alpha);	 
	   if (caching && st<0.75){
	     radCacheR[s]= so[iqu0+1]/multi;
	     radCachez[s]= so[iqu0+2]/multi;
	   }
	   
	   
	 }
	 /* else{ */
	 /* 	 so[iqu0+1]=pa[iqd]*cptr->radR; */
	 /* 	 so[iqu0+2]=pa[iqd]*cptr->radz; */
	 /* 	 so[iqe]=(cptr->radR*pa[iqu0+1]+cptr->radz*pa[iqu0+2])*pa[iqd]; 
	 /* 	 
	 /* } */
       }
     }
   }
 }
 else if (strcmp(cas, "KHINST") == 0){
   
 }
 double esrf=1.0;
 
 
 /*
   Physical  sources  (per unit  volume).   Note  primitives are  in
   pa. If the entropy is used, you also need to add in its sources
 */
 return;
 /*new*/}



/*
  + usrref - reference solution.
  Description:
  Provides a reference solution  for cell cptr. 
*/
void usrref(cell *cptr, float qref[])
{
  int i, id;
  float pa[MAXQ], rho, ss;
  for (i = 0; i < nq; i++)
    qref[i] = mabs(cptr->qa[i]);
  ss = sqrt(cptr->qa[iqe]*cptr->qa[iqd]);
  for (id = 1; id <= nd; id++)
    qref[iqu0+id] = mabs(cptr->qa[iqu0+id]) + ss;
  rho = cptr->qa[iqd];
  if (keps){
    qref[iqk] = 0.01*qref[iqe]/rho;
    qref[iqeps] = qref[iqk];
  }
  for (i = iqal0 + 1; i <= iqal0 + scad; i++)
    qref[i] = rho;
  for (i = iqal0 + scad + 1; i <= iqal0 + nal; i++)
    qref[i] = 1.0;
  if (grav)
    qref[iqpot] = big;
  return;
}



/*
+ usrqsol - Solution quality.
     Description:
        This returns  qsol = badcell if  the solution is  poor, qsol =
        good if  not. This is  determined by comparing  cptr->qa (fine
        cell   value),  with   cptr->pa  (mapped   down   coarse  cell
        value). usrref provides a reference solution.
*/
float usrqsol(int ilev, struct cell *cptr)
{
  int i;
  float qsol, qref[MAXQ];
  usrref(cptr, qref);
  double X=1, R,z;
  R=cptr->xc[0];
  z=cptr->xc[1];
  double r=sqrt(R*R+z*z);
  double st=R/r;
  double cs=cs0*pow(R, 0.5 * discP); 
  if (R<1) cs=cs0*pow(1.0/5, 0.5 * discP);
  double H=cs/sqrt(MG/pow(R,3)) * 1-flatDisc;
  if (strcmp(cas, "DISCSTAR") == 0){
    if (r>.975) X=1+st*st*5/r;
    else       X=pow(r,6);
  }else if (strcmp(cas, "LARGEDISCSTAR") == 0){
    if (r>.975) X=1+5 * sqrt(50/r)*st*st;
  }
  qsol = good;
  if (keps){
    for (i = 0; i < iqk; i++)
      if (mabs(cptr->qa[i] - cptr->pa[i]) > rtol*qref[i])
	qsol = badcell;
    for (i = iqal0+1; i < nq; i++)
      if (mabs(cptr->qa[i] - cptr->pa[i]) > rtol*qref[i])
	qsol = badcell;
  }
  else{
    /* if (strcmp(cas, "DISCSTAR") == 0){ */
    /*   if (  (mabs(cptr->qa[iqd] - cptr->pa[iqd]) > rtol*qref[iqd]/X) */
    /* 	  ||(mabs(cptr->qa[iqu0+1] - cptr->pa[iqu0+1]) > 5*rtol*qref[iqu0+1]/X) */
    /* 	  ||(mabs(cptr->qa[iqu0+2] - cptr->pa[iqu0+2]) > 5*rtol*qref[iqu0+2]/X) */
    /* 	  ||(mabs(cptr->qa[iqe] - cptr->pa[iqe]) > 3*rtol*qref[iqe]/X)) */
    /* 	qsol = badcell; */
    /* }else { */
    
    if (mabs(cptr->qa[iqd] - cptr->pa[iqd]) > rtol*qref[iqd]/X ||
	mabs(cptr->qa[iqu0+1] - cptr->pa[iqu0+1]) > rtol*qref[iqu0+1]/X ||
	mabs(cptr->qa[iqu0+2] - cptr->pa[iqu0+2]) > rtol*qref[iqu0+2]/X){
      qsol = badcell;
    }
    if ((strcmp(cas, "DISCSTAR") == 0) &&
	((ilev>=3 && r>10.0)||
	 (ilev>=4 && r>5.0 )||
	 (ilev>=5 && r>2.5 )||
	 (ilev>=6 && r>1.25 ))){

      qsol=good;
    }    if ((strcmp(cas, "LARGEDISCSTAR") == 0) &&
	((ilev>=5  && r>1600)||
	 (ilev>=6  && r>800)||
	 (ilev>=7  && r>400)||
	 (ilev>=8  && r>200)||
	 (ilev>=9  && r>100)||
	 (ilev>=10 && r>50 )||
	 (ilev>=11 && r>25 )||
	 (ilev>=12 && r>12.5 ))){

      qsol=good;
    }
	
  
  }
  if (strcmp(cas, "LARGEDISCSTAR") == 0 && (r < (0.9*insDiscRadius)) || (st<0.1 && r>20) ){
    qsol=good;
  }
  return(qsol);
}



/*
+ usrte - Temperature.
*/
float usrte(float pa[])
{
  return(pa[iqe]/pa[iqd]);
}



/*
+ usrmu - Diffusion coefficients.
*/
void usrmu(float pa[], float mu[])
{
  float rho, tv, tv1, tl, mut, c2, tm2;
  int i;
  /* Zero for safety */
  for (i = 0; i < nq; i++)
    mu[i] = 0.0;
  rho = pa[iqd];
  if (viscous){
    /* Laminar diffusion coefficients */
    for (i = iqu0 + 1; i <= iqu0 + nd; i++)
      mu[i] = rho*vis;
    mu[iqe] = rho*tcond;
    tv = rho*sdiff;
    for (i = iqal0 + 1; i <= iqal0 + scad; i++)
      mu[i] = tv;
    for (i = iqal0 + scad + 1; i <= iqal0 + nal; i++)
      mu[i] = sdiff;
  }
  if (keps){
    /* k-epsilon */
    tv = max(pa[iqk], kmin);
    tv1 = sqrt(tv);
    tl = tv*tv1/max(pa[iqeps],epsmin);
    tl = min(ceps*tl, mtl);
    mut = rho*tv1*tl;
    /* Turbulent viscosity */
    tv = cs*mut;
    for (i = iqu0 + 1; i <= iqu0 + nd; i++)
      mu[i] += tv;
    /* Turbulent thermal conduction */
    mu[iqe] += ct*mut;
    /* k diffusion */
    mu[iqk] = tv;
    /* epsilon diffusion */
    mu[iqeps] = ce*mut;
    /* Turbulent scalar diffusion */
    tv = cd*mut;
    for (i = iqal0 + 1; i <= iqal0 + scad; i++)
      mu[i] += tv;
  }
  return;
}



/*
+ usrbd - Boundary conditions.
     Description:
        Sets  the primitive arrays,  pl(left), pr(right)  according to
        the boundary condition specified by the boundary flag, bound.
*/
void usrbd(int id, int bound, join *jptr, float pl[], float pr[])
{
  int i;
  float qa[MAXQ], pa[MAXQ], xc[MAXQ];
  if (bound == lebound){
    /* Left boundary -- default is pl = pr */
    for (i = 0; i < nq; i++)
      pl[i] = pr[i];
    if (bcl[id-1] == bcref)
      /* Symmetry -- reverse normal velocity */
      pl[iqu0+id] = - pr[iqu0+id];
    else if (bcl[id-1] == bcwall)
      /* Wall==outflow velocities forced to be felt at left right at right */
      pl[iqu0+id] = - fabs(pr[iqu0+id]);
    else if (bcl[id-1] == bcfix){
      /* Fix -- fix solution */
      for (i = 0; i < nd; i++)
	xc[i] = jptr->xj[i];
      xc[id-1] = 2.0*xc[id-1] - jptr->right->xc[id-1];
      usrqis(xc, qa);
      usrpq(qa, pa);
      for (i = 0; i < nq; i++)
	pl[i] = pa[i];
    }
    if ((grav) &&(bcl[id-1] != bcref))
    /* Gravitational potential is zero unless symmetry */
      pl[iqpot] = -pr[iqpot];
  }
  else{
    /* Right boundary -- default is pr = pl */
    for (i = 0; i < nq; i++)
      pr[i] = pl[i];
    if (bcr[id-1] == bcref)
      /* Symmetry -- reverse normal velocity */
      pr[iqu0+id] = - pl[iqu0+id];
    else if (bcr[id-1] == bcwall)
      /* Wall==outflow velocities forced to be felt at left right at right */
      pr[iqu0+i] = fabs(pl[iqu0+i]);
    else if (bcr[id-1] == bcfix){
      /* Fix -- fix solution */
      for (i = 0; i < nd; i++)
	xc[i] = jptr->xj[i];
      xc[id-1] = 2.0*xc[id-1] - jptr->left->xc[id-1];
      usrqis(xc, qa);
      usrpq(qa, pa);
      for (i = 0; i < nq; i++)
	pr[i] = pa[i];
    }
    if ((grav) &&(bcr[id-1] != bcref))
    /* Gravitational potential is zero unless symmetry */
      pr[iqpot] = -pl[iqpot];
  }
  return;
}



/*
+ usrcmd - User commands.
    Description :
     This checks to see if it knows the command string in key, and and
     if necessary  performs the associated actions,  returning true to
     stop other command menus from seeing it.
*/
bool usrcmd(char key[])
{
  float tv;
  char filename[SZNAM];
  bool cmd, found, fail;

  cmd = true;
  if (strcmp(key, "PFLOOR") == 0){
    if (!get_r(pfloor)){
      cout << "Error : reading PFLOOR" << '\n';
      report();
    }
  }
  else if (strcmp(key, "DFLOOR") == 0){
    if (!get_r(dfloor)){
      cout << "Error : reading DFLOOR" << '\n';
      report();
    }
  }
  else if (strcmp(key, "AVISP") == 0){
    if (!get_r(avisp)){
      cout << "Error : reading AVISP" << '\n';
      report();
    }
  }
  else if (strcmp(key, "AVISE") == 0){
    if (!get_r(avise)){
      cout << "Error : reading AVISE" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GAMMA") == 0){
    if (!get_r(g)){
      cout << "Error : reading G" << '\n';
      report();
    }
    else
      usrsetg();
  }
  else if (strcmp(key, "KEPS") == 0){
    keps = true;
    usrkeps();
  }
  else if (strcmp(key, "KMIN") == 0){
    if (get_r(tv) == 0){
      cout << "Error : reading KMIN" << '\n';
      report();
    }
    else{
      kmin = tv;
      if (!keps)
      cout << "Warning : not turbulent" << '\n';
    }
  }
  else if (strcmp(key, "EPSMIN") == 0){
    if (get_r(tv) == 0){
      cout << "Error : reading EPSMIN" << '\n';
      report();
    }
    else{
      epsmin = tv;
      if (!keps)
      cout << "Warning : not turbulent" << '\n';
    }
  }
  else if (strcmp(key, "MTL") == 0){
    if (get_r(tv) == 0){
      cout << "Error : reading MTL" << '\n';
      report();
    }
    else{
      mtl = tv;
      if (!keps)
      cout << "Warning : not turbulent" << '\n';
    }
  }
  else if (strcmp(key, "MON") == 0){
    if (get_n(filename) == 0){
      if (mype == 0)
	cout << "Error : reading monfile name" << '\n';
      report();
    }
    else{
      if (mype == 0){
	strcat(filename,".xq");
	fail = (monfile=fopen(filename, "r"));
	if (fail)
	  fclose(monfile);
      }
      broadcast(&fail);
      if (fail){
	if (mype == 0)
	  cout << "Error: file exists \n";
	report();
      }
      else{
	if (mype == 0)
	  fail = (!(monfile = fopen(filename, "w")));
	broadcast(&fail);
	if (fail){
	  if (mype == 0)
	    cout << "Error: unable to open monitor file \n";
	  report();
	}
	else
	  mon = true;
      }
    }
  }
  else if (strcmp(key, "NOMON") == 0){
    mon = false;
    if ((monfile != NULL) && (mype == 0))
      fclose(monfile);
  }
  else if (strcmp(key, "FXP") == 0)
    fxp = true;
  else if (strcmp(key, "NOFXP") == 0)
    fxp = false;
  else if (strcmp(key, "VISCOUS") == 0)
    viscous = true;
  else if (strcmp(key, "NOVISCOUS") == 0)
    viscous = false;
  else if (strcmp(key, "INTEN") == 0)
    inten = true;
  else if (strcmp(key, "NOINTEN") == 0)
    inten = false;
  else if (strcmp(key, "TCOND") == 0){
    if (get_r(tcond) == 0){
      cout << "Error : reading TCOND" << '\n';
      report();
    }
  }
  else if (strcmp(key, "SDIFF") == 0){
    if (get_r(sdiff) == 0){
      cout << "Error : reading SDIFF" << '\n';
      report();
    }
  } 
  else if (strcmp(key, "VIS") == 0){
    if (get_r(vis) == 0){
      cout << "Error : reading VIS" << '\n';
      report();
    }
  }
  else if (strcmp(key, "GCONST") == 0){
    if (get_r(gconst) == 0){
      cout << "Error : reading GCONST" << '\n';
      report();
    }
  }
  else if (strcmp(key, "RAT") == 0){
    if (get_r(rat) == 0){
      cout << "Error : reading RAT" << '\n';
      report();
    }
  } 
  else if (strcmp(key, "HEAT") == 0){
    if (get_r(heat) == 0){
      cout << "Error : reading HEAT" << '\n';
      report();
    }
  }
  else
    /* Command not known here, pass it back to caller */
    cmd = false;

  return(cmd);
}



/*
+ usrprt - Print user information.
    Description :
     Prints user information. It is called by the USERS command.
*/
void usrprt()
{
  cout.setf(ios::scientific);
  cout << setprecision(4);
  cout << "no of scalars " << nal << " no of of advecting scalars "
       << scad << '\n';
  cout << "g " << g << '\n';
  cout << "dfloor " << dfloor << '\n';
  cout << "pfloor " << pfloor << '\n';
  cout << "avisp " << avisp << '\n';
  cout << "avise " << avise << '\n';
  cout << "rat " << rat << '\n';
  cout << "heat " << heat << '\n';
  if (viscous)
    cout << "viscous " << " vis " << vis << " tcond " << tcond << 
      " sdiff " << sdiff <<'\n';
  else 
    cout << "inviscid " << '\n';
  if (keps)
    cout << "k-epsilon: kmin " << kmin << "  epsmin " << epsmin <<
      " max turb length scale mtl " << mtl << '\n';
  if (inten)
    cout << "calculating entropy " << '\n';
  if (grav){
    cout << "gravity nrel " << nrel << " normrat " << 
      normrat << " gconst " << gconst;
    if (boundpot)
      cout << " boundpot levb " <<  levb << '\n';
    else
      cout << " noboundpot " << '\n';
  }
  else
    cout << "no gravity " << '\n';
  if (mon)
    cout << "monitoring " << '\n';
  else 
    cout << "no monitoring " << '\n';
  return;
}



/*
+ usrwr - Write user information.
    Description :
     Writes  user information to  a model  file. It  is called  by the
     WRITE command.
*/
void usrwr()
{
  write_var("iqd", iqd);
  write_var("iqu0", iqu0);
  write_var("iqe", iqe);
  write_var("iqal0", iqal0);
  write_var("iqie", iqie);
  write_var("iqpot", iqpot);
  write_var("keps", keps);
  if (keps){
    write_var("iqk", iqk);
    write_var("iqeps", iqeps);
    write_var("ce1", ce1);
    write_var("ce2", ce2);
    write_var("cs", cs);
    write_var("ct", ct);
    write_var("cd", cd);
    write_var("ce", ce);
    write_var("ceps", ceps);
    write_var("kmin", kmin);
    write_var("epsmin", epsmin);
    write_var("mtl", mtl);
  }
  write_var("g", g);
  write_var("pfloor", pfloor);
  write_var("dfloor", dfloor);
  write_var("avisp", avisp);
  write_var("avise", avise);
  write_var("mon", mon);
  write_var("fxp", fxp);
  write_var("vis", vis);
  write_var("inten", inten);
  write_var("tcond", tcond);
  write_var("sdiff", sdiff);
  write_var("viscous", viscous);
  write_var("gconst", gconst);
  write_var("rat", rat);
  write_var("heat", heat);
  return;
}



/*
+ usrrd - Read user information.
    Description :
     Reads user information to a model  file. It is called by the READ
     command.
*/
bool usrrd(char vname[])
{
  bool found;
  found = true;
  if (strcmp(vname, "iqd") == 0)
    read_var(&iqd);
  else if (strcmp(vname, "iqu0") == 0)
    read_var(&iqu0);
  else if (strcmp(vname, "iqe") == 0)
    read_var(&iqe);
  else if (strcmp(vname, "keps") == 0)
    read_var(&keps);
  else if (strcmp(vname, "iqk") == 0)
    read_var(&iqk);
  else if (strcmp(vname, "iqeps") == 0)
    read_var(&iqeps);
  else if (strcmp(vname, "iqal0") == 0)
    read_var(&iqal0);
  else if (strcmp(vname, "iqie") == 0)
    read_var(&iqie);
  else if (strcmp(vname, "iqpot") == 0)
    read_var(&iqpot);
  else if (strcmp(vname, "ce1") == 0)
    read_var(&ce1);
  else if (strcmp(vname, "ce2") == 0)
    read_var(&ce2);
  else if (strcmp(vname, "cs") == 0)
    read_var(&cs);
  else if (strcmp(vname, "ct") == 0)
    read_var(&ct);
  else if (strcmp(vname, "cd") == 0)
    read_var(&cd);
  else if (strcmp(vname, "ce") == 0)
    read_var(&ce);
  else if (strcmp(vname, "ceps") == 0)
    read_var(&ceps);
  else if (strcmp(vname, "kmin") == 0)
    read_var(&kmin);
  else if (strcmp(vname, "epsmin") == 0)
    read_var(&epsmin);
  else if (strcmp(vname, "mtl") == 0)
    read_var(&mtl);
  else if (strcmp(vname, "g") == 0)
      read_var(&g);
  else if (strcmp(vname, "pfloor") == 0)
    read_var(&pfloor);
  else if (strcmp(vname, "dfloor") == 0)
    read_var(&dfloor);
  else if (strcmp(vname, "avisp") == 0)
    read_var(&avisp);
  else if (strcmp(vname, "avise") == 0)
    read_var(&avise);
  else if (strcmp(vname, "mon") == 0)
    read_var(&mon);
  else if (strcmp(vname, "fxp") == 0)
    read_var(&fxp);
  else if (strcmp(vname, "vis") == 0)
    read_var(&vis);
  else if (strcmp(vname, "tcond") == 0)
    read_var(&tcond);
  else if (strcmp(vname, "sdiff") == 0)
    read_var(&sdiff);
  else if (strcmp(vname, "gconst") == 0)
    read_var(&gconst);
  else if (strcmp(vname, "viscous") == 0)
    read_var(&viscous);
  else if (strcmp(vname, "inten") == 0)
    read_var(&inten);
  else if (strcmp(vname, "rat") == 0)
    read_var(&rat);
  else if (strcmp(vname, "heat") == 0)
    read_var(&heat);
  else
    /* Unknown here  */
    found = false;
  /* Check for turbulence */
  keps = iqk > 0;
  /* Set gamma constants for nonlinear Riemann solver */
  usrsetg();
  return(found);
}



/*
+ usrsum - sum variables over the grid.
     Description :
        Calculates the sum of variables over the grid. The sum is over
        levels plmin to plmax.
*/
void usrsum()
{
  int i, ilev;
  float cvol, tot[MAXQ], pa[MAXQ]; 
  struct cell *cptr;
  /* Return if no model */
  if (nom)
    return;
  /* Zero sum, */
  for (i = 0; i < nq; i++)
    tot[i] = 0.0;
  /* Sum of all variables */
  for (ilev = plmin; ilev <= plmax; ilev++){
    cptr = lev[ilev]->fcell;
    while (cptr != NULL){
      if (!cptr->refined){
	cvol = cptr->vol;
	for (i = 0; i < nq; i++)
	  tot[i] += cvol*cptr->qa[i];
      }
      cptr = cptr->nxt;
    }
  }
  /* Get global sum */
  for (i = 0; i < nq; i++)
    tot[i] = global_sum(tot[i]);
  /* Output sum of all variables */
  if (mype == 0){
    cout.setf(ios::scientific);
    for (i = 0; i < nq; i++)
      cout << "variable " << i << " " << setprecision(6) << tot[i]
	   << '\n';
  }
  return;
}



/*
+ usrmon - Monitors solution.
     Description :
        This is called  every timestep if mon=true. It  can be used to
        output quantities to stdout or a file at selected intervals.
*/
void usrmon(){
  static int moni=0;
  moni+=1;
  if (!mon || (moni%10)==0){
    /* No monitor - do nothing */
      return;
  } else{
    /* Sum of all variables */
    int i,ilev;
    double cvol, tot[MAXQ];
    struct cell *cptr;
    for (i = 0; i < nq; i++)
      tot[i] = 0.0;
    for (ilev = plmin; ilev <= plmax; ilev++){
      cptr = lev[ilev]->fcell;
      while (cptr != NULL){
	if (!cptr->refined){
	  cvol = cptr->vol;
	  for (i = 0; i < nq; i++)
	    tot[i] += cvol*cptr->qa[i];
	}
	cptr = cptr->nxt;
      }
    }
    fprintf(monfile, "%.3e  ",t);
    for (i = 0; i < nq; i++){
      fprintf(monfile, "%.4e  ",tot[i]);
    }
    fprintf(monfile, "\n",tot[i]);
    return;
  }
}



/*
+ usrdef - set defaults.
     Description :
        Sets the user defaults at start up.
*/
void usrdef()
{
  g = 5.0/3.0;
  usrsetg();
  mon = false;
  fxp = true;
  inten = false;
  dfloor = 0.001;
  pfloor = 0.001;
  avisp = 0.2;
  avise = 0.2;
  viscous = false;
  grav = false;
  gconst = 1.0;
}



/*
+ usrkeps - set k-epsilon constants.
     Description :
        Sets the k-epsilon constants.
*/
void usrkeps()
{
  ce1 = 1.4;
  ce2 = 1.94;
  cs = 0.09;
  ct = 0.225;
  cd = 0.09;
  ce = 0.0692;
  ceps = 1.0;
  kmin = 0.0;
  epsmin = 0.0;
  mtl = 1.0;
  viscous = true;
  return;
}



/*
+ usrsetg - set gamma constants.
     Description :
        Sets   the  gamma   constants  for   the   non-linear  Riemann
        solver. Must be called whenever gamma is changed.
*/
void usrsetg()
{
  g1 = g - 1.0;
  g2 = g/g1;
  g3 = 0.5*(g + 1.0)/g;
  g4 = 0.5*g1/g;
  g5 = 1.0/g4;
  g6 = 1.0/g;
  g7 = 2.0/g1;
  g8 = g1/(g + 1.0);
  return;
}



/*
+ usroff - set ofsets.
     Description :
        Sets the offsets for the variables at start up.
*/
bool usroff()
{
  iqd = 0;
  iqu0 = 0;
  iqe = 1 + nd;
  /* k- epsilon */
  if (keps){
    iqk = iqe + 1;
    iqeps = iqk + 1;
    iqal0 = iqeps;
  }
  else{
    iqk = -1;
    iqal0 = iqe;
  }
  nq = iqal0 + nal + 1;
  if (inten){
    iqie = nq;
    nq = iqie + 1;
  }
  if (grav){
    /* Gravity */
    iqpot = nq;
    nq++;
  }
  else
    iqpot = 0;
  if (nq > MAXQ){
    cout << "Error : nq = " << nq << " > MAXQ = " << MAXQ <<
      " Change MAXQ in struct.h" << '\n';
    return(false);
  }
  return(true);
}



/*
+ usrvel - Returns velocity.
     Description :
        This returns the id direction  velocity in the cell pointed to
        by cptr. Used for plotting velocity vectors and streamlines.
*/
float usrvel(int id, struct cell *cptr)
{
  return(cptr->qa[iqu0+id]/cptr->qa[iqd]);
}



/*
+ usrbfield - Returns velocity.
     Description :
        This returns the id direction  velocity in the cell pointed to
        by   cptr.   Used    for   plotting   velocity   vectors   and
        streamlines. This is just a copy of usrvel.
*/
float usrbfield(int id, struct cell *cptr)
{
  return(cptr->qa[iqu0+id]);
}



/*
+ usrgfield - Returns gravitational acceleration.
     Description :

        This  returns the id  direction gravitational  acceleration in
        the cell  pointed to by cptr.  Used  for plotting acceleration
        vectors and streamlines.
*/
float usrgfield(int id, struct cell *cptr)
{
  return(cptr->pa[iqu0+id]);
}



/*
+ usrga - Returns gravitational acceleration.
     Description :
        This calculates the id direction gravitational acceleration in
        the cell  pointed to by cptr,  with cells left,  right on left
        and right.  The result is  put in cptr->pa[iqu0+id].
*/
void usrga(int id, struct cell *cptr, struct cell *left,
  struct cell *right)
{
  int count;
  float ga;
  /* Initialise */
  count = 0;
  ga = 0.0;
  if (left != NULL){
    ga += cptr->qa[iqpot] - left->qa[iqpot];
    count++;
  }
  if (right != NULL){
    ga +=  right->qa[nq-1] - cptr->qa[nq-1];
    count++;
  }
  if (count == 0)
    cptr->pa[iqu0+id] = 0.0;
  else{
    ga *= - gconst/((float)(count)*
      (cptr->rjoin[id-1]->xj[id-1] - cptr->ljoin[id-1]->xj[id-1]));
    cptr->pa[iqu0+id] = ga;
  }
}



/*
+ usrprompt - set prompt.
     Description:
        This sets the prompt at start up.
*/
void usrprompt(char prompt[])
{
  strcpy(prompt, "MG_G>");
}



/*
+ usrpath - set path.
     Description:
        This sets the path at start up.
*/
void usrpath(char path[])
{
  strcpy(path, "mg_g_data");
}



/*
+ usrquit - clean up user dynamic arrays.
     Description:
        This  deletes any dynamic arrays defined by the user.
*/
void usrquit()
{
for (TStrDblMap::iterator it=radCacheR.begin(); it!=radCacheR.end(); ++it){
            cerr << it->first << "  => " << it->second << '\n';
    }
  return;
}

 

