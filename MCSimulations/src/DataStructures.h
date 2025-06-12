//
//  DataStructures.h
//  TriangulatedMonteCarloMembrane
//

#ifndef DataStructures_h
#define DataStructures_h

#include "PreprocessorDeclarations.h"
#include <stdio.h>
#include "uthash.h"

typedef struct{

  double gamma; // surface tension
  double kappa; // bending modulus

  double sidex,sidey,sidez;
  double Isidex,Isidey,Isidez;
  double side2x,side2y,side2z;

  double vol,Ivol,volC;
  int Lcelx,Lcely,Lcelz,Ncel;
  double celsidex,celsidey,Icelsidex,Icelsidey;
  double celsidez,Icelsidez;

  double etaC,etaR,etaR2;
  double sig2,sig1_sq,sum_sig,um_sig_sq,q;
  double rv,rc,rb,rv2,rb2,Vgap;
  double rb2_2D;
  double C0,C1,C2,CC;
  double eps_mc,eps_mv;

  double rad;
  double Svol;
  double SurfArea;
  double eqBondDist;

  double EbendChange;
  double EbindChange_colloid_bead;
  double EbindChange_bead_colloid;
  double Evolchange;
  double Esurfchange;


} SYSTEM;

typedef struct{
  double x,y,z;     
  double Vx,Vy,Vz;
  double vx,vy,vz;  
  int verl[2700],bond[2700],clus[2700];
  int list[9100+1];
  int nlist;
  int type;
  int Nverl,Nbond;
  int cellID;
  int before,after;
  int Ng,ng[15];
  int ngR[15];
  int ngL[15];
} PARTICLE;

typedef struct{
  int Nv,v[3+1];
  int Nt,t[3+1];
} TRIANGLE;

typedef struct{
  int ngb[30];
  int begin;
} CELL;

typedef struct{
  double x,y;
} vec2D;

typedef struct{
  double x,y,z;
} vec3D;

typedef struct{
  double theta,phi;
} ANGOLO;

typedef struct{
  double x,y,z; 
} MEMBRANE;

typedef struct{
 double sig;
 double sig2;
 double rad;
 double dis_mc,rv,rv2,rc, rc2;
 double update_bcoll_list;
} COLL;

typedef struct {
    int key;            
    UT_hash_handle hh; 
} HashSet;

// ***************************************************
// GLOBAL STRUCTURES
// ***************************************************
SYSTEM S;
PARTICLE *Elem;
TRIANGLE *Tr;
MEMBRANE *MEM;
CELL *Cella;
COLL Coll;
int VERTATEDGE[5000];

// for checking the edges
HashSet *set = NULL;

double V0,A0;
int LX, LY, LZ;

// ***************************************************
// GLOBAL VARIABLES
// ***************************************************
int Ntri, N, Ncoll, Nedges, Nbound, counter, seed, Nbonds;
int Tmax=10; // [ON]: never allow more then Tmax-coordinated points.
int is_mem_fluid;
int start_from_config;
int L_system;
int ID_force;
int inorout;

double D0, phi1, phi2, theta1, theta2;
double lattice,lattice_y;
double cutoff,cutoff2;
double cm_phi,cm_theta;

// factors for harmonic system
double k_harmonic;
double r_eq;
double speed_pulling;
double position_bead;
double equilibrium;
double eq_old;
double factor_LJ;

char readname[500];   
char errorfile[500];  
char Roname[500];    
char Prob_name[500];

FILE *read,*wr1,*wr2,*prob;
vec3D Centre;

int accepted_moves;

#endif /* DataStructures_h */
