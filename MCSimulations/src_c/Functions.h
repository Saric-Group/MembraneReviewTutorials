// *****************************************************************************
//  DYNAMICALLY TRIANGULATED MEMBRANE -- SARIC GROUP (ISTA)
//  Read main.c and the README.md file to get further details about the code
// *****************************************************************************

#ifndef Functions_h
#define Functions_h

#include <stdio.h>
#include "DataStructures.h"
#include "PreprocessorDeclarations.h"
#include "uthash.h"

//..
////////////////////////////////////////////////////////////////////////////////
// *****************************************************************************
// FUNCTION DECLARATION
// *****************************************************************************
////////////////////////////////////////////////////////////////////////////////
// ..

// *****************************
// Functions for INITIALIZATION
// ****************************
void set_all(void);
void set_Elements(void);
void set_nb(void);
void set_triangles(void);
void set_Cells(void);
void set_neighbors(void);
void set_Verlet(void);
void set_vertices_edge(void);
void update_vertices_edge(void);
void b_coll_vlist(void);
void read_Elements(void);

// *******************************
// Functions INVOLVED IN MC moves
// *******************************
int Interaction_1p(int);
int Interaction_ALL(void);
int Interaction_1coll(int);
double bending_1p(int,int);
double bending_edge(int,int);
double bind_1p_colloid_bead(int);
double bind_1p_bead_colloid(int);
double interaction_coll(int);
void MC_xyz(int);
void switch_bond(int);
int check_if_edge(int k);

// ***********************************
// PAINTERS - functions to print data
// ***********************************
void painter(int SWEEP);
void painter_en(int, int, int, int, int);
void painter_nanoparticles(void);
void write_config(void);
void write_config_error(void);

// **************************
// BASIC FUNCTION/OPERATIONS
// **************************
vec3D normale(int,int,int);
vec3D rotate(double , double ,double  ,double, double);
vec3D rotx(vec3D,double);
vec3D rotz(vec3D,double);
ANGOLO angle(double , double ,double);

// ***************************************
// FUNCTIONS ADDED TO CONSTRAIN AREA
// ***************************************
double area_change(int,int);
double area_change_2tr(int,int,int,int,int,int);

// ************************
// Miscellaneous functions
// ************************
double compute_total_surface(void);
double compute_pulling_force(void);


// ..
////////////////////////////////////////////////////////////////////////////////
// *****************************************************************************
// FUNCTION DEFINITION
// *****************************************************************************
////////////////////////////////////////////////////////////////////////////////
//..

void add(HashSet **set, int key) {
    HashSet *entry;
    HASH_FIND_INT(*set, &key, entry);  // Check if it exists
    if (entry == NULL) {
        entry = malloc(sizeof(HashSet));
        entry->key = key;
        HASH_ADD_INT(*set, key, entry);  // Add it to the hash table
    }
}

int contains(HashSet *set, int key) {
    HashSet *entry;
    HASH_FIND_INT(set, &key, entry);
    return entry != NULL;
}

void free_all(HashSet *set) {
    HashSet *curr, *tmp;
    HASH_ITER(hh, set, curr, tmp) {
        HASH_DEL(set, curr);
        free(curr);
    }
}

// -----------------------------
// FUNCTIONS FOR INITIALIZATION
// -----------------------------
void set_all(void){

    // initialize random seed
    srand48(seed);

    // maximum bond length for membrane beads
    cutoff=sqrt(3);
    cutoff2=cutoff*cutoff;

    // KEEP the '+2'
    S.sidex=(double)(LX+2); // length X simulation box
    S.sidey=(double)(LY+2); // length Y simulation box
    S.sidez=(double)(LZ+2); // length Z simulation box

    S.side2x=S.sidex/2.0;
    S.side2y=S.sidey/2.0;
    S.side2z=S.sidez/2.0;
    S.Isidex=1.0/S.sidex;
    S.Isidey=1.0/S.sidey;
    S.Isidez=1.0/S.sidez;
    S.vol=S.sidex*S.sidey*S.sidez;  // volume of the simulation box
    S.Ivol=1./S.vol;

    S.rv=2.2;           // verlet minimum
    S.rb=S.rv;          // bond order minimum
    S.rc=1.78;          // cluster minimum
    S.rb2_2D=1.5*1.5;   // Bond Order Cut-OFF for 2D order parameter
    S.rv2=S.rv*S.rv;
    S.rb2=S.rb*S.rb;
    S.Vgap=0.5*(S.rv-1);  // put the largest instead of 1.

    S.Lcelx=(int)(S.sidex/S.rb);    // num cells along x
    S.Lcely=(int)(S.sidey/S.rb);    // num cells along y
    S.Lcelz=(int)(S.sidez/S.rb);    // num cells along z
    S.Ncel=S.Lcelx*S.Lcely*S.Lcelz; // total number of cells

    printf("Lcelx: %d Lcely: %d Lcelz: %d\n", S.Lcelx, S.Lcely, S.Lcelz);

    S.celsidex=S.sidex/S.Lcelx;
    S.celsidey=S.sidey/S.Lcely;
    S.celsidez=S.sidez/S.Lcelz;

    printf("Celsidex: %lf Celsidey: %lf Celsidez: %lf\n", S.celsidex, S.celsidey, S.celsidez);

    S.Icelsidex=1.0/S.celsidex;
    S.Icelsidey=1.0/S.celsidey;
    S.Icelsidez=1.0/S.celsidez;

    // ---- set of warnings on the dimensions of the box
    if (S.Lcelx>LX){
        // number of cells along x more than length of direction x
        printf("Lcel=%d> SIDEX OF BOX=%d System too sparse increase rv and rb\n",S.Lcelx,LX);
        exit(-1);
    }
    if (S.Lcely>LY){
        // number of cells along y more than length of direction y
        printf("Lcel=%d> SIDEY OF BOX=%d System too sparse increase rv and rb\n",S.Lcely,LY);
        exit(-1);
    }
    if (S.Lcelz>LZ){
        // number of cells along z more than length of direction z
        printf("Lcelz=%d> SIDEZ OF BOX=%d System too sparse increase rv and rb\n",S.Lcelz,LZ);
        exit(-1);
    }
    if (S.celsidex<=S.rb){
        // length of the cell along x less than bond order
        printf("Lcel=%lf<rb decrease Rb\n",S.celsidex);
        exit(-1);
    }
    if (S.celsidey<=S.rb){
        // length of the cell along y less than bond order
        printf("Lcel=%lf<rb decrease Rb\n",S.celsidey);
        exit(-1);
    }
    if (S.celsidez<=S.rb){
        // length of the cell along z less than bond order
        printf("Lcelz=%lf<rb decrease Rb\n",S.celsidez);
        exit(-1);
    }

    // --- MC parameters
    S.eps_mc=0.10;
    S.eps_mv=0.0008;

    // --- allocate space for data structures (the +1 comes bc the zero-th index is not used in this code)
    Elem  =(PARTICLE *) malloc( (N+Ncoll+1)*sizeof(PARTICLE) );
    Tr    =(TRIANGLE *) malloc( (8*N+1)*sizeof(TRIANGLE) ); 
    Cella =(CELL *) malloc( (LX*LY*LZ+1)*sizeof(CELL) ); 
    MEM   =(MEMBRANE *) malloc( (N+1)*sizeof(MEMBRANE) );

    sprintf(readname , "conf_N_%d_Sig_%f_kappa_%lf_gamma_%lf_seed_%d_.dump",N, 2*Coll.rad, S.kappa, S.gamma, seed);
    sprintf(errorfile, "errorconf_N_%d_Sig_%f_kappa_%lf_gamma_%lf_seed_%d_.dump",N, 2*Coll.rad, S.kappa, S.gamma, seed);

    set_Elements();
    
    // --- prepare verlet lists, linked lists etc...
    // set cells in the system
    set_Cells();
    // verlet list for interaction with membranes
    set_Verlet();
    // verlet list for interaction with colloidal particles
    b_coll_vlist();
    // introduce the edge vertices
    set_vertices_edge();

}

void set_vertices_edge(void){

    int vertedge_counter = 1;

    // for the particles at the edges
    for (int k = 1; k<=Nedges; k++){

        // for their neighbours
        for (int w=1; w<=Elem[k].Ng; w++){
            int j=Elem[k].ng[w];
            
            //printf("Adding as edge... %d \n", j);
            // include them in the hash table so that
            // it is more efficient to find them
            add(&set, j);
        }
    }
}

void update_vertices_edge(void){

    // for the particles at the edges
    for (int k = 1; k<=Nedges; k++){
        // for their neighbours
        for (int w=1; w<=Elem[k].Ng; w++){
            int j=Elem[k].ng[w];
            // for their neighbours
            for (int s = 1; s<Elem[j].Ng; s++){
                int m=Elem[j].ng[s];
                add(&set, m);
            }            
        }
    }
}

void set_Elements(void){

    // Function gets called when initial conditions file does not exist
    // This function produces said file by reading from your init file
    // (which might have different names depending on the type of sim)
    // There is a python file that produces these initial configurations
    // (called ProduceInitialConfig.py)
    //
    // Steps in the function:
    // 1. Read in particle coordinates from init file
    // 2. Compute theta and phi angles
    // 3. Compute center of mass of the system
    // 4. Launch nb
    // 5. Set triangles
    // 6. Call painter2

    int i,j,k,p,t=0;
    double ssside,disp,vvv;
    double dx,dy,dz,d2;
    double Ex,Ey,Ez;
    double radius_membrane;
    double dummy_CX, dummy_CY, dummy_CZ;
    char to[500];
    FILE *in;

    // --- read file with initial configuration
    sprintf(to,"flatpatch_N_%d_Sigma_%d_.dat",N,(int)(Coll.sig));
    
    printf("Reading initial conditions from %s\n",to);
    in=fopen(to,"r");

    // --- read in line by line and save the coordinates of each particle (Nbeads and colloids)
    for (t=1;t<=N+Ncoll;t++){
        fscanf(in,"%d %lf %lf%lf\n",&Elem[t].type,&Elem[t].x,&Elem[t].y,&Elem[t].z);

        // slightly randomize initial positions? why?
        Elem[t].x+=.001*(.5-drand48());
        Elem[t].y+=.001*(.5-drand48());
        Elem[t].z+=.001*(.5-drand48());
    }

    // --- compute center of mass of the system
    Centre.x=0.0;Centre.y=0.0;Centre.z=0.0;
    for (t=1;t<=N;t++){
        Centre.x+=Elem[t].x;
        Centre.y+=Elem[t].y;
        Centre.z+=Elem[t].z;
    }

    Centre.x/=(double)(N);
    Centre.y/=(double)(N);
    Centre.z/=(double)(N);

    printf("Center of mass of the system (%2.3f %2.3f %2.3f)\n",Centre.x,Centre.y,Centre.z);

    // --------------------- we must do this so that we can compute the triangles properly (it is important for non-trivial sims)
    // displace first membrane to center
    dummy_CX = 0.0; dummy_CY = 0.0; dummy_CZ = 0.0;
    for(t=1; t<=N; t++){
        dummy_CX+= Elem[t].x;
        dummy_CY+= Elem[t].y;
        dummy_CZ+= Elem[t].z;
    }

    dummy_CX /= (double)(N);
    dummy_CY /= (double)(N);
    dummy_CZ /= (double)(N);

    double max_extension =0;
    for(t=1; t<=N; t++){
        MEM[t].x = Elem[t].x-dummy_CX;
        MEM[t].y = Elem[t].y-dummy_CY;
        MEM[t].z = Elem[t].z-dummy_CZ;

        if(max_extension<MEM[t].x){
            max_extension = MEM[t].x;
        }
    }


    S.rad = max_extension;
    if(max_extension>LX*0.5){
        printf("Problem: The box is not big enough for the membrane\n");
        exit(1);
    }

    set_nb();
    set_triangles();
    printf("All set.");

    S.SurfArea = compute_total_surface();

}

void set_Cells(void){

    int i, x, y, z;
    double facx,facy,facz;
    double xtemp, ytemp, ztemp;

    S.Lcelx=(int)(S.sidex/S.rb);
    S.Lcely=(int)(S.sidey/S.rb);
    S.Lcelz=(int)(S.sidez/S.rb);

    S.Ncel=S.Lcelx*S.Lcely*S.Lcelz;
    S.celsidex=S.sidex/S.Lcelx;
    S.celsidey=S.sidey/S.Lcely;
    S.celsidez=S.sidez/S.Lcelz;
    S.Icelsidex=1./S.celsidex;
    S.Icelsidey=1./S.celsidey;
    S.Icelsidez=1./S.celsidez;
    facx=(double)(S.Lcelx)/S.sidex;
    facy=(double)(S.Lcely)/S.sidey;
    facz=(double)(S.Lcelz)/S.sidez;

    // initialize cella
    for (i=0;i<S.Ncel;i++){
        Cella[i].begin=0;
    }

    // determine what particles are in what cells
    for (i=1;i<=N;i++){

        xtemp = (Elem[i].x+S.side2x)*facx;
        ytemp = (Elem[i].y+S.side2y)*facy;
        ztemp = (Elem[i].z+S.side2z)*facz;

        x = (int)(xtemp);
        y = (int)(ytemp);
        z = (int)(ztemp);

        Elem[i].cellID=x+y*S.Lcelx+z*S.Lcelx*S.Lcely;
        Elem[i].after=Cella[Elem[i].cellID].begin;
        Elem[ Elem[i].after ].before=i;
        Cella[ Elem[i].cellID ].begin=i;
    }

    set_neighbors();
}
void set_neighbors(void){

    // Function that creates cell neighbours

    int k,s,i,j;
    int x,y,z;
    int x_P,x_M;
    int y_P,y_M;
    int z_P,z_M;

    for (x=0;x<S.Lcelx;x++){
        for (y=0;y<S.Lcely;y++){
            for (z=0;z<S.Lcelz;z++){

                k = x + y*S.Lcelx + z*S.Lcelx*S.Lcely;

                if (x==S.Lcelx-1){
                  x_P=0;
                }
                else{
                  x_P=x+1;
                }
                if (x==0){
                  x_M=S.Lcelx-1;
                }
                else{
                  x_M=x-1;
                }
                if (y==S.Lcely-1){
                  y_P=0;
                }
                else{
                  y_P=y+1;
                }
                if (y==0){
                  y_M=S.Lcely-1;
                }
                else{
                  y_M=y-1;
                }
                if (z==S.Lcelz-1){
                  z_P=0;
                }
                else{
                  z_P=z+1;
                }
                if (z==0){
                  z_M=S.Lcelz-1;
                }
                else{
                  z_M=z-1;
                }

                Cella[k].ngb[0]=x+(y)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[1]=x_P+(y)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[2]=x_M+(y)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[3]=x+(y_P)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[4]=x+(y_M)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[5]=x_P+(y_P)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[6]=x_M+(y_P)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[7]=x_M+(y_M)*S.Lcelx+(z)*S.Lcelx*S.Lcely;
                Cella[k].ngb[8]=x_P+(y_M)*S.Lcelx+(z)*S.Lcelx*S.Lcely;

                Cella[k].ngb[9 ]=x_P+(y)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[10]=x_M+(y)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[11]=x+(y_P)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[12]=x+(y_M)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[13]=x_P+(y_P)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[14]=x_M+(y_P)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[15]=x_M+(y_M)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[16]=x_P+(y_M)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;

                Cella[k].ngb[17]=x_P+(y)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
                Cella[k].ngb[18]=x_M+(y)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
                Cella[k].ngb[19]=x+(y_P)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
                Cella[k].ngb[20]=x+(y_M)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
                Cella[k].ngb[21]=x_P+(y_P)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
                Cella[k].ngb[22]=x_M+(y_P)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
                Cella[k].ngb[23]=x_M+(y_M)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;
                Cella[k].ngb[24]=x_P+(y_M)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;

                Cella[k].ngb[25]=x+(y)*S.Lcelx+(z_P)*S.Lcelx*S.Lcely;
                Cella[k].ngb[26]=x+(y)*S.Lcelx+(z_M)*S.Lcelx*S.Lcely;

            }
        }
    }
}

void set_Verlet(void){

    // Function that builds the Verlet list for the interactions
    // between membrane beads

    int i,j,neigh,k;
    double r2,dx,dy,dz;

    // coordinates of the verlet list
    for (i=1;i<=N;i++){
        Elem[i].Nverl=0;
        Elem[i].Nbond=0;

        // center of the verlet list for membrane beads
        Elem[i].Vx=Elem[i].x;
        Elem[i].Vy=Elem[i].y;
        Elem[i].Vz=Elem[i].z;
    }

    // fill in the verlet list (linked list)
    for (i=1;i<=N;i++){

        for (j=0;j<=26;j++){  // 26 is the number of neighboring cells to 1 cell

            neigh=Cella[ Elem[i].cellID ].ngb[j];
            k=Cella[neigh].begin;

            while (k!=0){

                if (k!=i && k>=1){
                    dx=Elem[i].x-Elem[k].x;
                    dy=Elem[i].y-Elem[k].y;
                    dz=Elem[i].z-Elem[k].z;

                    r2=(dx*dx+ dy*dy+ dz*dz);

                    // they are in the same verlet list (they are neighbours)
                    if (r2<S.rv2){
                        Elem[i].Nverl++;
                        Elem[i].verl[Elem[i].Nverl]=k;
                    }

                    // they are connected by a bond
                    if (r2<S.rb2){
                        Elem[i].Nbond++;
                        Elem[i].bond[Elem[i].Nbond]=k;
                    }
                }

                k=Elem[k].after;
            }
        }
    }
}

void set_nb(void){

    // In this function we get the neighbours of beads in the membrane
    // and we set the bonds for the triangles therein

    int i,j,k,w,s,keepL,keepR, keepIDL, keepIDR;
    double dx,dy,dz,dr, dummy;
    double COS,SIN,tm_L,tm_R;
    double low=100000;
    vec3D rr[N+1];
    vec3D n,n_jk;
    ANGOLO cm;

    char iiiii[500];
    FILE *o;
    sprintf(iiiii,"setnb_init.dat");
    o=fopen(iiiii,"w");

    Nbonds = 0;

    // --- iterate over membrane beads to get neighbours (triangular lattice first neighbours are 6)
    for (i=1;i<=N;i++){

        Elem[i].Ng=0;
        rr[i].x=0;
        rr[i].y=0;
        rr[i].z=0;

        for (j=1;j<=N;j++){

            if (i!=j){

                dx=Elem[i].x-Elem[j].x;
                dy=Elem[i].y-Elem[j].y;
                dz=Elem[i].z-Elem[j].z;

                dr=sqrt(dx*dx + dy*dy + dz*dz);

                // if these two membrane beads are close to each other
                if (dr<lattice){

                    // control to find out smallest distance
                    if (dr<low){
                        low=dr;
                    }

                    // update the total number of bonds
                    Nbonds+=1;

                    Elem[i].Ng++;               // neighbours of particle i -- neighbours in the lattice
                    Elem[i].ng[ Elem[i].Ng ]=j; // saving the specific neighbour of i, which is j

                    // write down the neighbours
                    fprintf(o, "%d\t%d\n", i, j);

                }
            }
        }

        // if i has more than 7 --> error (edges in flat membrane can have less than five neighbours)
        if (Elem[i].Ng>7){
            printf("Problem: Triangular lattice seems to be poorly initialized\n");
            exit(-1);
        }
    }

    fclose(o);

    // --- CONTROL: if two particles are overlapping, flag error
    if (low<1.0){
        printf("Problem: Min. dist bw two particles <1 =%f (overlapping)\n",low);
        exit(-1);
    }

    // PLEASE KEEP IN MIND THAT THIS LOOP ONLY WORKS WELL IF THE CM OF THE VESICLE IS AT 0
    // This is why one has to use the MEM object
    // --- set bonds triagles -- k IS particle index
    for (k=1;k<=N;k++){

        // get position of membrane bead and transform into angles (used centered coords)
        cm=angle(MEM[k].x,MEM[k].y,MEM[k].z);
        // rotate vector so that it aligns with z-axis (sort of 'new base')
        rr[k]=rotate(MEM[k].x,MEM[k].y,MEM[k].z,cm.theta,cm.phi);
        rr[k].x=MEM[k].x;
        rr[k].y=MEM[k].y;
        rr[k].z=MEM[k].z;

        // iterate around the neighbours of k
        for (w=1; w<=Elem[k].Ng; w++){
            j=Elem[k].ng[w];
            rr[j].x=MEM[j].x;
            rr[j].y=MEM[j].y;
            rr[j].z=MEM[j].z;
        }

        // for neighbours of k (called j)
        for (w=1;w<=Elem[k].Ng;w++){
            j=Elem[k].ng[w];

            dx=rr[j].x-rr[k].x; 
            dy=rr[j].y-rr[k].y; 

            dr=sqrt(dx*dx+dy*dy);
            n_jk.x=dx/dr;
            n_jk.y=dy/dr;

            tm_L=-1;
            tm_R=-1;

            keepL=-1;
            keepR=-1;
            keepIDL=-1;
            keepIDR=-1;

            // for neighbours of k that are not j (called s)
            for (i=1;i<=Elem[k].Ng;i++){

                if (i!=w){

                    s=Elem[k].ng[i];

                    dx=rr[s].x - rr[k].x; 
                    dy=rr[s].y - rr[k].y; 

                    dr=sqrt(dx*dx+dy*dy);
                    n.x=dx/dr;
                    n.y=dy/dr;

                    // cosine taken from scalar product
                    COS=n.x*n_jk.x + n.y*n_jk.y;

                    // sine taken from vector product
                    SIN=n_jk.x*n.y - n_jk.y*n.x;

                    if (SIN>0.0){ // s is "up" --> counterclockwise w.r.t. jk (right hand rule)
                        if (COS>tm_L){
                            keepL=i; // keep in mind we save the neighbor index ID, not the ID itself
                            keepIDL=s;
                            tm_L=COS;
                        }
                    }

                    if (SIN<0.0){ // s is "down" --> clockwise w.r.t. jk
                        if (COS>tm_R){
                            keepR=i; // keep in mind we save the neighbor index ID, not the ID itself
                            keepIDR=s;
                            tm_R=COS;
                        }
                    }
                }
            }

            // THIS WORKS FOR THE FLAT PATCH BUT NOTE THAT IT IS A VERY SPECIFIC INITIALIZATION
            if(tm_L>0){
                Elem[k].ngL[w]=keepL;
            }
            else{
                Elem[k].ngL[w]=-1;
            }
            if(tm_R>0){
                Elem[k].ngR[w]=keepR;
            }
            else{
                Elem[k].ngR[w]=-1;
            }
            
            printf("k: %d j: %d keepIDR: %d keepIDL: %d\n", k, j, keepIDR, keepIDL);
        }
    }
}

void set_triangles(void){

    int i,j,k,flag,t;
    int tp1,tp2,tp3;

    t=0;

    char iiiii[500];
    FILE *o;
    sprintf(iiiii,"settriangles_init.dat");
    o=fopen(iiiii,"w");

    // for particles in membrane
    for (i=1;i<=N;i++){

        // for neighbours of particle i
        for (j=1;j<=Elem[i].Ng;j++){

            tp1=i;
            tp2=Elem[i].ng[j];  // getting a particle ID
            tp3=Elem[i].ngL[j]; // getting an index

            if (tp3>0){

                // tp3 tells me where is the Left neighbor in the list of site "i"
                tp3=Elem[i].ng[tp3]; // getting a particle ID
                fprintf(o, "%d\t%d\t%d\n", tp1, tp2, tp3);

                //---------------
                // for triangles already added
                for (k=1;k<=t;k++){
                    // [ON] verify that the triangle is not already taken
                    flag=0;
                    if (tp1==Tr[k].v[0] ||tp1==Tr[k].v[1] ||tp1==Tr[k].v[2]){
                        flag++;
                    }
                    if (tp2==Tr[k].v[0] ||tp2==Tr[k].v[1] ||tp2==Tr[k].v[2]){
                        flag++;
                    }
                    if (tp3==Tr[k].v[0] ||tp3==Tr[k].v[1] ||tp3==Tr[k].v[2]){
                        flag++;
                    }
                    if (flag==3){
                        // we have found a triangle 'k' with all vertices agreeing with our system
                        break;
                    }
                }

                // saving the triangles (checks that has not been saved yet)
                if (flag!=3){
                    t++;
                    // add new triangle t with vertices v
                    Tr[t].v[0]=tp1;
                    Tr[t].v[1]=tp2;
                    Tr[t].v[2]=tp3;
                    Tr[t].Nv=3; // needed?
                }
            }
            //---------------
        }
    }

    fclose(o);

    // counts and prints total number of triangles
    Ntri=t;
    printf("T=%d Ntri=%d\n",t,Ntri);

    // [ON] Now define the neighs of each triangle
    //      the neigh of v[1] must be the triangle t[1]
    //      that is opposite to it; (they must share 2 points = 1 line)

    for (i=1;i<=Ntri;i++){

        // [ON] Ngb 0 opposite to v[0]
        tp2=Tr[i].v[1];
        tp3=Tr[i].v[2];

        for (j=1;j<=Ntri;j++){
            if (i!=j){
                flag=0;
                if (tp2==Tr[j].v[0] ||tp2==Tr[j].v[1] ||tp2==Tr[j].v[2]){
                    flag++;
                }
                if (tp3==Tr[j].v[0] ||tp3==Tr[j].v[1] ||tp3==Tr[j].v[2]){
                    flag++;
                }
                if (flag==2){// accept
                    Tr[i].t[0]=j;
                    break;
                }
            }
        }

        // [ON] Ngb 1 opposite to v[1]
        tp1=Tr[i].v[0];
        tp3=Tr[i].v[2];

        for (j=1;j<=Ntri;j++){
            if (i!=j){
                flag=0;
                if (tp1==Tr[j].v[0] ||tp1==Tr[j].v[1] ||tp1==Tr[j].v[2]){
                    flag++;
                }
                if (tp3==Tr[j].v[0] ||tp3==Tr[j].v[1] ||tp3==Tr[j].v[2]){
                    flag++;
                }
                if (flag==2) {// accept
                    Tr[i].t[1]=j;
                    break;
                }
            }
        }

        // Ngb 2 opposite to v[2]
        tp1=Tr[i].v[0];
        tp2=Tr[i].v[1];

        for (j=1;j<=Ntri;j++){
            if (i!=j){
                flag=0;
                if (tp1==Tr[j].v[0] ||tp1==Tr[j].v[1] ||tp1==Tr[j].v[2]){
                    flag++;
                }
                if (tp2==Tr[j].v[0] ||tp2==Tr[j].v[1] ||tp2==Tr[j].v[2]){
                    flag++;
                }
                if (flag==2) {// accept
                    Tr[i].t[2]=j;
                    break;
                }
            }
        }

        Tr[i].Nt=3;

        // print neighbouring triangles/triangle information
        printf("i=%d) %d %d %d\n",i,Tr[i].v[0],Tr[i].v[1],Tr[i].v[2]);
        printf("t1=%d) %d %d %d\n",Tr[i].t[0],Tr[ Tr[i].t[0] ].v[0],Tr[ Tr[i].t[0] ].v[1],Tr[ Tr[i].t[0] ].v[2]);
        printf("t2=%d) %d %d %d\n",Tr[i].t[1],Tr[ Tr[i].t[1] ].v[0],Tr[ Tr[i].t[1] ].v[1],Tr[ Tr[i].t[1] ].v[2]);
        printf("t3=%d) %d %d %d\n",Tr[i].t[2],Tr[ Tr[i].t[2] ].v[0],Tr[ Tr[i].t[2] ].v[1],Tr[ Tr[i].t[2] ].v[2]);
        printf("\n\n");
         
  }
}


void b_coll_vlist(void){

    // Function created neighbour list for colloidal particles
    int j,i;
    double dx,dy,dz,r2;

    // verlet list of the colloidal particle
    for (i=1;i<=N+Ncoll;i++) {
        Elem[i].nlist=0;

        // center of verlet colloidal lists
        Elem[i].vx=Elem[i].x;
        Elem[i].vy=Elem[i].y;
        Elem[i].vz=Elem[i].z;
    }

    for (i=1;i<=N;i++) {
        for (j=N+1;j<=N+Ncoll;j++) {

            dx=Elem[i].x-Elem[j].x;
            dy=Elem[i].y-Elem[j].y;
            dz=Elem[i].z-Elem[j].z;

            r2=(dx*dx+dy*dy+dz*dz);

            if (r2<Coll.rv2) {
                Elem[i].nlist+=1;
                Elem[j].nlist+=1;

                Elem[i].list[Elem[i].nlist]=j;
                Elem[j].list[Elem[j].nlist]=i;
            }
        }
    }
}

// -----------------------------
// PAINTERS
// -----------------------------
void painter(int SWEEP){

    // Function that writes particle types plus coordinates
    // This is the file that ovito reads (the 'time evolution' of the system)

    int p,j,s,l,r, bondnum = 0;
    char iiiii[500], bbbbb[500];
    FILE *o, *bo;
    
    sprintf(iiiii,"out.dump");
    o=fopen(iiiii,"a+");
    sprintf(bbbbb,"bonds.dump");
    bo=fopen(bbbbb,"a+");

    // getting the file ready for ovito
    fprintf(o,"ITEM: TIMESTEP\n");
    fprintf(o,"%d\n",(SWEEP));
    fprintf(o,"ITEM: NUMBER OF ATOMS\n");
    fprintf(o,"%d\n",(N+Ncoll));
    fprintf(o,"ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(o,"%f %f\n",-S.side2x, S.side2x);
    fprintf(o,"%f %f\n",-S.side2y, S.side2y);
    fprintf(o,"%f %f\n",-S.side2z, S.side2z);
    fprintf(o,"ITEM: ATOMS id type x y z\n");

    for (p=1;p<=N;p++){
        // assumes that membrane beads have diameter 1.0
        fprintf(o,"%d %d %f %f %f\n",p, 1,Elem[p].x,Elem[p].y,Elem[p].z);
    }

    for(p=N+1; p<=N+Ncoll; p++){
        fprintf(o,"%d %d %f %f %f\n",p, 2,Elem[p].x,Elem[p].y,Elem[p].z);
    }
    fclose(o);

    // getting the file ready for ovito
    fprintf(bo,"ITEM: TIMESTEP\n");
    fprintf(bo,"%d\n",(SWEEP));
    fprintf(bo,"ITEM: NUMBER OF ENTRIES\n");
    fprintf(bo,"%d\n",(Nbonds));
    fprintf(bo,"ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(bo,"%f %f\n",-S.side2x, S.side2x);
    fprintf(bo,"%f %f\n",-S.side2y, S.side2y);
    fprintf(bo,"%f %f\n",-S.side2z, S.side2z);
    fprintf(bo,"ITEM: ENTRIES index particle1 particle2\n");

    for (int k = 1; k<=N; k++){
        for (int w=1; w<=Elem[k].Ng; w++){
            j=Elem[k].ng[w];
            fprintf(bo,"%d %d %d\n", bondnum, k, j);
            bondnum+=1;
        }
    }

    fclose(bo);

}

void painter_en(int i, int counter, int N, int Ncoll, int Ntri){

    // OLD VERSION OF THE ENERGY PAINTER
    // Function that writes energy of the configuration

    int p, pe,pc;
    double EnCol, EnBend;

    //EnCol=0;
    EnBend=0;

    // can there be redundancies here
    // (given that we might go over the same edge several times?)
    for (p=1;p<=Ntri;p++){
        for (pe=0;pe<=2;pe++){
            if (bending_edge(p,pe) == bending_edge(p,pe)){ //getting rid of nan values
                EnBend += S.kappa*bending_edge(p,pe);
            }
        }
    }

    // compute harmonic bond energy
    double dx=Elem[Nbound].x-Elem[N+Ncoll].x;
    double dy=Elem[Nbound].y-Elem[N+Ncoll].y;
    double dz=Elem[Nbound].z-Elem[N+Ncoll].z;
    double r2=(dx*dx+dy*dy+dz*dz);
    double sqrtr = sqrt(r2);
    double harmenergy= 0.5*k_harmonic *(sqrtr - r_eq)*(sqrtr - r_eq);

    char iiiii[500];
    FILE *o;

    sprintf(iiiii,"energy.dump");
    o=fopen(iiiii,"a+");
    fprintf(o,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",i,S.kappa,S.gamma,EnBend,S.SurfArea,harmenergy, S.EbendChange,S.EbindChange_bead_colloid,S.Esurfchange);

    fclose(o);

}

void write_config(void){

    // Function writes last configuration

    int i,j,npart,trash;
    double tmp1,tmp2,tmp3,tmp4,tmp5;
    FILE *wr1;

    wr1=fopen(readname,"w");
    for (i=1;i<=N;i++){
        fprintf(wr1,"%4.15lf %4.15lf %4.15lf %d\n",Elem[i].x,Elem[i].y,Elem[i].z,Elem[i].Ng);
        for (j=1;j<=Elem[i].Ng;j++){
            fprintf(wr1,"%d %d %d\n",Elem[i].ng[j],Elem[i].ngR[j],Elem[i].ngL[j]);
        }
    }

    // including the colloidal particles
    for(i=N+1; i<=N+Ncoll; i++){
      fprintf(wr1,"%4.15lf %4.15lf %4.15lf %d\n",Elem[i].x,Elem[i].y,Elem[i].z,-1);
    }

    fclose(wr1);
}

void write_config_error(void){

    // Function writes last configuration

    int i,j,npart,trash;
    double tmp1,tmp2,tmp3,tmp4,tmp5;
    FILE *wr1;

    wr1=fopen(errorfile,"w");
    for (i=1;i<=N;i++){
        fprintf(wr1,"%4.15lf %4.15lf %4.15lf %d\n",Elem[i].x,Elem[i].y,Elem[i].z,Elem[i].Ng);
        for (j=1;j<=Elem[i].Ng;j++){
            fprintf(wr1,"%d %d %d\n",Elem[i].ng[j],Elem[i].ngR[j],Elem[i].ngL[j]);
        }
    }

    // including the colloidal particles
    for(i=N+1; i<=N+Ncoll; i++){
      fprintf(wr1,"%4.15lf %4.15lf %4.15lf %d\n",Elem[i].x,Elem[i].y,Elem[i].z,-1);
    }

    fclose(wr1);
}

void painter_nanoparticles(void){
    int i;
    char iiii[500];
    FILE *o;

    sprintf(iiii,"position_pulling_bead.dat");
    o=fopen(iiii,"a+");

    for(i=N+1; i<=N+Ncoll; i++){
        fprintf(o, "%d %d %lf %lf %lf %lf %lf %lf \n",counter, i, Elem[i].x, Elem[i].y, Elem[i].z, Elem[Nbound].x, Elem[Nbound].y, Elem[Nbound].z);
    }
 
    fclose(o);
    return;
}

// -----------------------------
// MC SWEEPS/INTERACTIONS/ENERGY
// -----------------------------

void MC_xyz(int oo){

    double dx,dy,dz,dx1,dy1,dz1,facx,facy,facz;
    double oldX,oldY,oldZ;
    double r2=0,r21=0,rs,rs1,DA1,DA2;
    int touch=0,ee,k;
    int i,s,x,y,z,oCel,nCel;
    double en,eo,bind_o,bind_n;
    double xtemp,ytemp,ztemp;
    double AOLD;
    double pull_old, extension_old, pull_new, extension_new;
    double sus_X, sus_Y, sus_Z;
    double correct_centre_X, correct_centre_Y, correct_centre_Z;

    // [ON] pick a random vertex "ee" from the random triangle "oo"
    ee=(int)(drand48()*3.0);
    k=Tr[oo].v[ee];

    // do not integrate vertices at the edges
    int is_edge = check_if_edge(k);
    if (is_edge==1){
        return;
    }
    
    // see if new verlet list has to be computed
    dx=Elem[k].x-Elem[k].Vx;
    dy=Elem[k].y-Elem[k].Vy;
    dz=Elem[k].z-Elem[k].Vz;

    r2=(dx*dx+ dy*dy+ dz*dz);
    rs=sqrt(r2);

    if (rs>S.Vgap){
        set_Verlet();
    }

    // ---------------------------------
    // OLD CONFIGURATION: positions + energy
    // ---------------------------------
    oldX=Elem[k].x;
    oldY=Elem[k].y;
    oldZ=Elem[k].z;
    oCel=Elem[k].cellID;

    // current bending energy
    eo=bending_1p(oo,ee); 

    // current local area
    DA1=area_change(oo,ee);

    // if it is the particle that is bound to the pulling bead
    // we also need to compute the interaction energy
    if(k==Nbound){
        bind_o=bind_1p_bead_colloid(k);
    }
    else{
        bind_o=0.0;
    }

    // ---------------------------------
    // ATTEMPT SMALL UPDATE: can we accept the move?
    // ---------------------------------
    Elem[k].x+=S.eps_mc*(0.5-drand48());
    Elem[k].y+=S.eps_mc*(0.5-drand48());
    Elem[k].z+=S.eps_mc*(0.5-drand48());

    // check no bonds are broken by this move
    for (i=1;i<=Elem[k].Ng;i++){
        s=Elem[k].ng[i];

        dx=Elem[s].x-Elem[k].x;
        dy=Elem[s].y-Elem[k].y;
        dz=Elem[s].z-Elem[k].z;

        r2=(dx*dx+ dy*dy+ dz*dz);

        // reject here (bond broken, only if not edge)
        if (r2 > cutoff2 ){ 
            is_edge = check_if_edge(s);
            if(is_edge==0){
                // non-edge bonds cannot stretch, but the others yes
                Elem[k].x=oldX;
                Elem[k].y=oldY;
                Elem[k].z=oldZ;
                return;
            }
        }
    }

    // compute new verlet list if needed
    dx=Elem[k].x-Elem[k].Vx;
    dy=Elem[k].y-Elem[k].Vy;
    dz=Elem[k].z-Elem[k].Vz;

    r2=(dx*dx+ dy*dy+ dz*dz);

    rs=sqrt(r2);
    if (rs>S.Vgap){
        set_Verlet();
    }

    // are there overlaps?
    touch=Interaction_1p(k);

    // reject if there are overlaps
    if (touch==1){  // reject
        Elem[k].x=oldX;
        Elem[k].y=oldY;
        Elem[k].z=oldZ;
        return;
    }
    else{

        if(k==Nbound){
            bind_n=bind_1p_bead_colloid(k);
        }
        else{
            bind_n=0.0;
        }

        // reject
        if (bind_n==999){
            Elem[k].x=oldX;
            Elem[k].y=oldY;
            Elem[k].z=oldZ;
            return ;
        }

        // maybe accept
        else{

            en=bending_1p(oo,ee);
            DA2=area_change(oo,ee);

            // save old surface area
            AOLD = S.SurfArea;
            // update surface area after change
            S.SurfArea+= (DA2-DA1);

            // ---------------------------------
            // METROPOLIS HASTINGS criterion to reject move
            // ---------------------------------

            if (drand48()>exp(-S.kappa*(en-eo)-S.gamma*(DA2-DA1)-(bind_n-bind_o))){
                // reject move
                Elem[k].x=oldX;
                Elem[k].y=oldY;
                Elem[k].z=oldZ;
                // return surface area to previous value
                S.SurfArea = AOLD;
                return;
            }
            else{
                // accept move
                accepted_moves+=1;

                S.EbendChange               += S.kappa*(en-eo);
                S.EbindChange_bead_colloid  += (bind_n-bind_o);
                S.Esurfchange               += S.gamma*(DA2-DA1);

                // update position of the center of mass of the system
                Centre.x=Centre.x-oldX/(double)(N)+Elem[k].x/(double)(N);
                Centre.y=Centre.y-oldY/(double)(N)+Elem[k].y/(double)(N);
                Centre.z=Centre.z-oldZ/(double)(N)+Elem[k].z/(double)(N);

                facx=(double)(S.Lcelx)/S.sidex;
                facy=(double)(S.Lcely)/S.sidey;
                facz=(double)(S.Lcelz)/S.sidez;

                x = (int)((Elem[k].x+S.side2x)*facx);
                y = (int)((Elem[k].y+S.side2y)*facy);
                z = (int)((Elem[k].z+S.side2z)*facz);

                nCel=x + y*S.Lcelx + z*S.Lcelx*S.Lcely;

                if (nCel!=oCel){
                    if (Cella[oCel].begin==k){
                        Cella[oCel].begin = Elem[k].after;
                    }
                    else{
                        Elem[Elem[k].before].after=Elem[k].after;
                        Elem[Elem[k].after].before=Elem[k].before;
                    }
                    Elem[k].cellID=nCel;
                    Elem[k].after=Cella[nCel].begin;
                    Cella[nCel].begin=k;
                    Elem[Elem[k].after].before=k;
                }
            }
        }
    }

}

double bind_1p_bead_colloid(int k){

    // Function for LJ-type of interaction

    int j,m;
    double dx,dy,dz,invr2,invr6,invr12;
    double sqrtr, invr;
    double r2=0, ee=0.0;

    // interaction with the pulling bead
    dx=Elem[k].x-Elem[N+Ncoll].x;
    dy=Elem[k].y-Elem[N+Ncoll].y;
    dz=Elem[k].z-Elem[N+Ncoll].z;
    r2=(dx*dx+dy*dy+dz*dz);
    sqrtr = sqrt(r2);
    ee+= 0.5*k_harmonic *(sqrtr - r_eq)*(sqrtr - r_eq);
    return ee;

}

int Interaction_1coll(int k){

    // This function checks whether after the move of colloidal particle k
    // there are overlaps with other colloidal particles

    int j;
    double dx,dy,dz,rs,dd;
    double r2=0;

    for (j=N+1;j<=N+Ncoll;j++){
        if (k!=j) {
            dx=Elem[k].x-Elem[j].x;
            dy=Elem[k].y-Elem[j].y;
            dz=Elem[k].z-Elem[j].z;

            r2=(dx*dx+ dy*dy+ dz*dz);

            if (r2<Coll.sig2){
                // contact-reject
                return 1;
            }
        }
    }

    // no overlap with other colloids
    return 0;

}

int check_if_edge(int k){

    // if the bead if part of the edges, do not move
    if(k<=Nedges){
        return 1;
    }

    if(contains(set, k)){
        return 1;
    }

    return 0;
}

int Interaction_1p(int k){

    int j, partid, is_edge;
    double dx,dy,dz,r2,rs,dd;

    for (j=1;j<=Elem[k].Nverl;j++){

        partid = Elem[k].verl[j];
        dx=Elem[k].x-Elem[ partid ].x;
        dy=Elem[k].y-Elem[ partid ].y;
        dz=Elem[k].z-Elem[ partid ].z;

        r2=(dx*dx+ dy*dy+ dz*dz);

        if (r2<1){

            // never allow compression because then it might block the system
            return 1;
        }
    }
    return 0;
}

int Interaction_ALL(void){
  int k,j;
  double dx,dy,dz,r2,rs,dd;

  for (k=1;k<=N;k++){
    for (j=1;j<=Elem[k].Nverl;j++){

      dx=Elem[k].x-Elem[ Elem[k].verl[j] ].x;
      dy=Elem[k].y-Elem[ Elem[k].verl[j] ].y;
      dz=Elem[k].z-Elem[ Elem[k].verl[j] ].z;

      r2=(dx*dx+ dy*dy+ dz*dz);

      if (r2<1){
        return 1;  // contact-reject
      }
    }
  }
  return 0;
}
void switch_bond(int k){

    // Function responsible for the fluidity of the membrane

    int i,j,Ei,Es,p,w1,w2;
    int r1,r2,r3,r4;
    int tmpK[3];
    int tmpP[3];
    int FixA,FixB,Em,En;
    int G1,G2;
    double d2,scalar,s1,s2,s3,s4,s5;
    double Energy_O,Energy_N;
    double AOLD,DA1,DA2;
    vec3D n2A,n4A,n2B,n4B;
    vec3D n1,n3,n5,n6;
    int r1indr3, r1indr4, r2indr3, r2indr4;

    Ei=(int)(drand48()*3);  // [ON] pick a rand vertex "s" for triangle "K"
    p=Tr[k].t[Ei];          // [ON] this define the second triangle "P" given "K"

    // if there is no triangle defined, we do not attempt to switch bonds
    if (p == 0){
        return;
    }
    
    if (Tr[p].t[0]==k){
        Es=0;
    }
    if (Tr[p].t[1]==k){
        Es=1;
    }
    if (Tr[p].t[2]==k){
        Es=2;}

    // [ON] detect what is the opposite point in P to Ei
    // [ON] now I know Ei & Es  (see figure "triangles.eps")
    // Make sure we are not trying to change the connectivity
    // of the edges at the vertices
    
    r1=Tr[k].v[Ei];
    // if vertex belong to edges, do not update connectivity
    if (r1<=Nedges){
        return;
    }

    r2=Tr[p].v[Es];
    // if vertex belong to edges, do not update connectivity
    if (r2<=Nedges){
        return;
    }

    r3=Tr[k].v[(Ei+1)%3];
    // if vertex belong to edges, do not update connectivity
    if (r3<=Nedges){
        return;
    }

    r4=Tr[k].v[(Ei+2)%3];
    // if vertex belong to edges, do not update connectivity
    if (r4<=Nedges){
        return;
    }

    if ((Elem[r1].Ng>=8)||(Elem[r2].Ng>=8)){
        return;
    }
    if ((Elem[r3].Ng<=4)||(Elem[r4].Ng<=4)){
        return ;
    }

    d2=((Elem[r1].x-Elem[r2].x)*(Elem[r1].x-Elem[r2].x)+
       (Elem[r1].y-Elem[r2].y)*(Elem[r1].y-Elem[r2].y)+
       (Elem[r1].z-Elem[r2].z)*(Elem[r1].z-Elem[r2].z));
    
    // if the candidate link is too long, check if one of the particles is an edge bond
    if (d2>cutoff2){
        return;
    }

    // Check for the convexity condition of the triangle's pair to be changed
    n2A=normale(Tr[k].v[Ei],Tr[k].v[(Ei+1)%3],Tr[k].v[(Ei+2)%3]);
    n4A=normale(Tr[p].v[Es],Tr[p].v[(Es+1)%3],Tr[p].v[(Es+2)%3]);
    scalar=n2A.x*n4A.x + n2A.y*n4A.y + n2A.z*n4A.z;

    // reject
    if (scalar<0){
        return;
    }
    // accept
    else{

        n2B=normale(Tr[k].v[Ei],Tr[k].v[(Ei+1)%3],Tr[p].v[Es]);
        n4B=normale(Tr[p].v[Es],Tr[p].v[(Es+1)%3],Tr[k].v[Ei]);
        scalar=n2B.x*n4B.x + n2B.y*n4B.y + n2B.z*n4B.z;

        // reject
        if (scalar<0){
            return;
        }
        // accept
        else{

            // [ON] Check for metropolis -- energy (Check Fig. triangles2.eps)
            FixA=Tr[k].t[(Ei+1)%3];
            G1=Tr[k].t[(Ei+2)%3];

            FixB=Tr[p].t[(Es+1)%3];
            G2=Tr[p].t[(Es+2)%3];

            n1=normale(Tr[FixA].v[0],Tr[FixA].v[1],Tr[FixA].v[2]);
            n3=normale(Tr[G1].v[0],Tr[G1].v[1],Tr[G1].v[2]);
            n5=normale(Tr[G2].v[0],Tr[G2].v[1],Tr[G2].v[2]);
            n6=normale(Tr[FixB].v[0],Tr[FixB].v[1],Tr[FixB].v[2]);

            s1=n1.x*n2A.x+ n1.y*n2A.y+ n1.z*n2A.z;
            s2=n3.x*n2A.x+ n3.y*n2A.y+ n3.z*n2A.z;
            s3=n4A.x*n2A.x+ n4A.y*n2A.y+ n4A.z*n2A.z;
            s4=n5.x*n4A.x+ n5.y*n4A.y+ n5.z*n4A.z;
            s5=n6.x*n4A.x+ n6.y*n4A.y+ n6.z*n4A.z;

            Energy_O=(5.0-(s1+s2+s3+s4+s5));

            s1=n1.x*n4B.x+ n1.y*n4B.y+ n1.z*n4B.z;
            s2=n4B.x*n2B.x+ n4B.y*n2B.y+ n4B.z*n2B.z;
            s3=n3.x*n2B.x+ n3.y*n2B.y+ n3.z*n2B.z;
            s4=n6.x*n2B.x+ n6.y*n2B.y+ n6.z*n2B.z;
            s5=n5.x*n4B.x+ n5.y*n4B.y+ n5.z*n4B.z;

            Energy_N=(5.0-(s1+s2+s3+s4+s5));
            AOLD=S.SurfArea;
            DA1=area_change_2tr(r1,r3,r4, r2,r4,r3);
            DA2=area_change_2tr(r1,r3,r2, r2,r4,r1);

            // change in surface area
            S.SurfArea+=(DA2-DA1);

            // Metropolis-Hastings criterion for rejecting a move
            if (drand48()>exp(-S.kappa*(Energy_N-Energy_O)-S.gamma*(DA2-DA1))){
                S.SurfArea=AOLD;
                return; //reject
            }

            // accepting the move
            else{
                S.EbendChange    += S.kappa*(Energy_N-Energy_O);
                S.Esurfchange    += S.gamma*(DA2-DA1);

                // [ON] these definitions will be useful later
                if (Tr[FixA].t[0]==k){
                    Em=0;
                }
                if (Tr[FixA].t[1]==k){
                    Em=1;
                }
                if (Tr[FixA].t[2]==k){
                    Em=2;
                }

                if (Tr[FixB].t[0]==p){
                    En=0;
                }
                if (Tr[FixB].t[1]==p){
                    En=1;
                }
                if (Tr[FixB].t[2]==p){
                    En=2;
                }

                // [ON] SWITCH TRIANGLES //

                // [on] first switch vertices
                Tr[k].v[(Ei+2)%3]=r2;
                Tr[p].v[(Es+2)%3]=r1;

                // [on] then triangle neighbors
                tmpK[Ei]=Tr[k].t[Ei];
                tmpK[(Ei+1)%3]=Tr[k].t[(Ei+1)%3];
                tmpK[(Ei+2)%3]=Tr[k].t[(Ei+2)%3];

                tmpP[Es]=Tr[p].t[Es];
                tmpP[(Es+1)%3]=Tr[p].t[(Es+1)%3];
                tmpP[(Es+2)%3]=Tr[p].t[(Es+2)%3];

                Tr[k].t[Ei]=tmpP[(Es+1)%3];
                Tr[k].t[(Ei+1)%3]=p;

                Tr[p].t[Es]=tmpK[(Ei+1)%3];
                Tr[p].t[(Es+1)%3]=k;

                Tr[FixA].t[Em]=p;
                Tr[FixB].t[En]=k;

                // finally add/subtruct bonds to the vertices
                // AND ADJUST NEIGHBORS
                int r1indr3, r1indr4, r2indr3, r2indr4;

                for(i=1;i<=Elem[r1].Ng;i++) {
                    if(Elem[r1].ng[i] == r3)
                    r1indr3 = i;
                    if(Elem[r1].ng[i] == r4)
                    r1indr4 = i;
                }
                for(i=1;i<=Elem[r2].Ng;i++) {
                    if(Elem[r2].ng[i] == r3)
                    r2indr3 = i;
                    if(Elem[r2].ng[i] == r4)
                    r2indr4 = i;
                }

                // add
                Elem[r1].Ng++;
                Elem[r1].ng[Elem[r1].Ng]=r2;
                //Adjust neighbors
                if(Elem[r1].ngL[r1indr4] == r1indr3) {
                    Elem[r1].ngL[r1indr4] = Elem[r1].Ng;
                    Elem[r1].ngR[r1indr3] = Elem[r1].Ng;
                    Elem[r1].ngL[Elem[r1].Ng] = r1indr3;
                    Elem[r1].ngR[Elem[r1].Ng] = r1indr4;
                }
                else if(Elem[r1].ngR[r1indr4] == r1indr3) {
                    Elem[r1].ngR[r1indr4] = Elem[r1].Ng;
                    Elem[r1].ngL[r1indr3] = Elem[r1].Ng;
                    Elem[r1].ngR[Elem[r1].Ng] = r1indr3;
                    Elem[r1].ngL[Elem[r1].Ng] = r1indr4;
                }
                else{
                    printf("\nMAJOR PROBLEM! (adjusting neighbors)\n");
                    exit(-1);
                }

                Elem[r2].Ng++;
                Elem[r2].ng[Elem[r2].Ng]=r1;
                if(Elem[r2].ngL[r2indr4] == r2indr3) {
                    Elem[r2].ngL[r2indr4] = Elem[r2].Ng;
                    Elem[r2].ngR[r2indr3] = Elem[r2].Ng;
                    Elem[r2].ngL[Elem[r2].Ng] = r2indr3;
                    Elem[r2].ngR[Elem[r2].Ng] = r2indr4;
                }
                else if(Elem[r2].ngR[r2indr4] == r2indr3) {
                    Elem[r2].ngR[r2indr4] = Elem[r2].Ng;
                    Elem[r2].ngL[r2indr3] = Elem[r2].Ng;
                    Elem[r2].ngR[Elem[r2].Ng] = r2indr3;
                    Elem[r2].ngL[Elem[r2].Ng] = r2indr4;
                }
                else{
                    printf("\nMAJOR PROBLEM! (adjusting neighbors)\n");
                    exit(-1);
                }

                // subtract
                int w1L, w1R, w2L, w2R, r3NgL, r3NgR, r4NgL, r4NgR;

                w1=0;
                for (i=1;i<=Elem[r3].Ng;i++){
                    if (Elem[r3].ng[i]==r4){
                        w1=i;break
                        ;
                    }
                }
                if (w1==0){
                    printf("Flag1. Problem.\n"); exit(-1);
                }

                w1L = Elem[r3].ngL[w1];
                w1R = Elem[r3].ngR[w1];
                Elem[r3].ngL[w1R] = w1L;
                Elem[r3].ngR[w1L] = w1R;

                if (w1==Elem[r3].Ng){
                    Elem[r3].Ng--;
                    goto via1;
                }

                Elem[r3].ng[w1]=Elem[r3].ng[ Elem[r3].Ng ];
                r3NgL = Elem[r3].ngL[Elem[r3].Ng];
                r3NgR = Elem[r3].ngR[Elem[r3].Ng];
                Elem[r3].ngL[w1] = r3NgL;
                Elem[r3].ngR[w1] = r3NgR;
                Elem[r3].ngR[r3NgL] = w1;
                Elem[r3].ngL[r3NgR] = w1;
                Elem[r3].Ng--;

                via1:
                    w2=0;
                    for (i=1;i<=Elem[r4].Ng;i++){
                        if (Elem[r4].ng[i]==r3){
                            w2=i;
                            break;
                            }
                        }
                if (w2==0){
                    printf("Flag 2. Problem.\n"); 
                    exit(-1);
                }

                w2L = Elem[r4].ngL[w2];
                w2R = Elem[r4].ngR[w2];
                Elem[r4].ngL[w2R] = w2L;
                Elem[r4].ngR[w2L] = w2R;
                
                if (w2==Elem[r4].Ng){
                    Elem[r4].Ng--;
                    goto via2;
                }

                Elem[r4].ng[w2]=Elem[r4].ng[ Elem[r4].Ng ];
                r4NgL = Elem[r4].ngL[Elem[r4].Ng];
                r4NgR = Elem[r4].ngR[Elem[r4].Ng];
                Elem[r4].ngL[w2] = r4NgL;
                Elem[r4].ngR[w2] = r4NgR;
                Elem[r4].ngR[r4NgL] = w2;
                Elem[r4].ngL[r4NgR] = w2;
                Elem[r4].Ng--;

                via2:{}
            }
        }
    }
}

// -----------------------------
// MEMBRANE MODEL PROPERTIES
// -----------------------------

double bending_1p(int oo,int ee){

    int i,j,k,vo,op,nx, here;
    double en=0.0,tmp;
    vec3D nn,nA,nB;

    k=Tr[oo].v[ee];

    here=oo;
    vo=ee;

    for (i=1;i<=Tmax;i++){
        
        // normal current triangle
        nn=normale(Tr[here].v[vo], Tr[here].v[(vo+1)%3],Tr[here].v[(vo+2)%3]);
        
        // triangle opposite to point "k"
        op=Tr[here].t[vo]; 
        nA=normale(Tr[op].v[0],Tr[op].v[1],Tr[op].v[2]);
        // triangle next (C-Clwise) to point "k"
        nx=Tr[here].t[(vo+1)%3]; 
        nB=normale(Tr[nx].v[0],Tr[nx].v[1],Tr[nx].v[2]);

        // discretized form of the bending energy
        tmp=2.-( (nn.x*nA.x+nn.y*nA.y+nn.z*nA.z)+(nn.x*nB.x+nn.y*nB.y+nn.z*nB.z) );
        en+=tmp;

        here=nx;

        vo=4;

        if (Tr[here].v[0]==k){
            vo=0;
        }
        if (Tr[here].v[1]==k){
            vo=1;
        }
        if (Tr[here].v[2]==k){
            vo=2;
        }
        
        if (vo==4){
            printf("BIG PROBLEM, CHECK TRIANGLES\n Triangle %d) p1=%d p2=%d p3=%d looking for(%d)\n\n", here,Tr[here].v[0],Tr[here].v[1],Tr[here].v[2],k );
            exit(-1);
        }

        if (here==oo){
            break;
        }
    }

    return en;
}

double bending_edge(int oo,int ee){

    int i,j,k,vo,op,nx;
    int here;
    double en=0.,tmp;
    vec3D nn,nA;

    k=Tr[oo].v[ee];
    here=oo;
    vo=ee;

    nn=normale(Tr[here].v[vo], Tr[here].v[(vo+1)%3],Tr[here].v[(vo+2)%3]);
    op=Tr[here].t[vo]; 
    nA=normale(Tr[op].v[0],Tr[op].v[1],Tr[op].v[2]);

    tmp=1.-(nn.x*nA.x+nn.y*nA.y+nn.z*nA.z);
    en+=tmp/2;
    return en;
}

// --------------------------------
// BASIC FUNCTIONS
// --------------------------------

vec3D normale(int o, int i,int j){

    // Function that computes the normal vector to the
    // plane defined by o, i and j

    double lg;
    vec3D a,b,n;

    // distance i, o
    a.x=Elem[i].x-Elem[o].x;
    a.y=Elem[i].y-Elem[o].y;
    a.z=Elem[i].z-Elem[o].z;

    if (a.x>S.side2x)    a.x-=S.sidex;
    if (a.x<-S.side2x)   a.x+=S.sidex;
    if (a.y>S.side2y)    a.y-=S.sidey;
    if (a.y<-S.side2y)   a.y+=S.sidey;
    if (a.z>S.side2z)    a.z-=S.sidez;
    if (a.z<-S.side2z)   a.z+=S.sidez;

    // distance j, o
    b.x=Elem[j].x-Elem[o].x;
    b.y=Elem[j].y-Elem[o].y;
    b.z=Elem[j].z-Elem[o].z;

    if (b.x>S.side2x)    b.x-=S.sidex;
    if (b.x<-S.side2x)   b.x+=S.sidex;
    if (b.y>S.side2y)    b.y-=S.sidey;
    if (b.y<-S.side2y)   b.y+=S.sidey;
    if (b.z>S.side2z)    b.z-=S.sidez;
    if (b.z<-S.side2z)   b.z+=S.sidez;

    // vector product?
    n.x=a.y*b.z-a.z*b.y;
    n.y=a.z*b.x-a.x*b.z;
    n.z=a.x*b.y-a.y*b.x;

    // the length of the vector
    lg=sqrt(n.x*n.x+n.y*n.y+n.z*n.z);

    // normalize it
    n.x/=lg;
    n.y/=lg;
    n.z/=lg;

    return n;
}

vec3D rotate(double x,double y,double z,double theta, double phi){

    // Function that rotates x, y, z coordinates along angles theta and phi

    vec3D first,th;
    vec3D second;

    th.x=x;
    th.y=y;
    th.z=z;

    // rotation along z axis
    first=rotz(th,phi);

    // rotation along x axis
    second=rotx(first,theta);

    return second;
}

vec3D rotx(vec3D a,double Aangle){

    // Function that performs a rotation around z axis of angle Aangle

    vec3D rx;

    rx.x=a.x;
    rx.y=a.y*cos(Aangle)-a.z*sin(Aangle);
    rx.z=a.y*sin(Aangle)+a.z*cos(Aangle);

    return rx;
}

vec3D rotz(vec3D a,double Aangle){

    // Function that performs a rotation around z axis of angle Aangle

    vec3D rz;

    rz.x=a.x*cos(Aangle)-a.y*sin(Aangle);
    rz.y=a.x*sin(Aangle)+a.y*cos(Aangle);
    rz.z=a.z;

    return rz;
}

ANGOLO angle(double x,double y, double z){

    // Function gets coordinates of a membrane bead
    // and return the theta and phi angles of this
    // membrane bead on the surface of a sphere

    int i;
    double length,d0,sign;
    ANGOLO pp;
    vec3D nor;

    // --- get normalized position vector
    nor.x=x;
    nor.y=y;
    nor.z=z;

    length=sqrt(nor.x*nor.x+nor.y*nor.y+nor.z*nor.z);

    nor.x/=length;
    nor.y/=length;
    nor.z/=length;

    // --- get the theta angle
    // the acos function in C takes an argumen -1 <= x <= 1 and return the arc cosine
    // the if functions below make sure we are within the region where the function is defined
    if( nor.z>1.0-EPSILON ){
        nor.z=1.-EPSILON;
    }
    if ( nor.z<-1.+EPSILON ){
        nor.z=-1.+EPSILON;
    }

    pp.theta=acos(nor.z);

    // --- get the phi angle
    if (fabs(nor.y)<EPSILON && fabs(nor.x)<EPSILON ){
        pp.phi=0.;
        goto here;
    }
    if (fabs(nor.y)<EPSILON && fabs(nor.x)>EPSILON ){
        pp.phi=M_PI/2.;
    }

    // get the phi angle -- but it needs correction to see in what quadrant
    pp.phi=atan(nor.x/nor.y);

    // not sure this 'here' makes a lot of sense
    here:{
        if (pp.phi!=0.){
            if (nor.x<0 && nor.y<0){
                pp.phi+=M_PI;
            }
            if (nor.x>0 && nor.y<0){
                pp.phi+=M_PI;
            }
        }
    }

    // returns angle with theta and phi
    return pp;
}

// --------------------------------
// PATCHES AND ADDITIONAL FUNCTIONS
// --------------------------------
double compute_total_surface(void){

  double total_surface = 0;
  int k, p, q;
  double h,area,ss;
  double l1,l2,l3,ln;
  vec3D r,u1,u2,u3,n;
  vec3D dir;

  for(int i = 1; i<=Ntri; i++){
    
    //total_surface+= area_change(i, 1);

    // MAYBE NOT AS FAST BUT WE GO TRIANGLE BY TRIANGLE
    k=Tr[i].v[0];
    p=Tr[i].v[1];
    q=Tr[i].v[2];
        
    u1.x=Elem[p].x-Elem[k].x;
    u1.y=Elem[p].y-Elem[k].y;
    u1.z=Elem[p].z-Elem[k].z;

    u2.x=Elem[q].x-Elem[k].x;
    u2.y=Elem[q].y-Elem[k].y;
    u2.z=Elem[q].z-Elem[k].z;

    u3.x=Elem[q].x-Elem[p].x;
    u3.y=Elem[q].y-Elem[p].y;
    u3.z=Elem[q].z-Elem[p].z;

    l1=sqrt(u1.x*u1.x +u1.y*u1.y +u1.z*u1.z);
    l2=sqrt(u2.x*u2.x +u2.y*u2.y +u2.z*u2.z);
    l3=sqrt(u3.x*u3.x +u3.y*u3.y +u3.z*u3.z);

    ss=.5*(l1+l2+l3);
    area=(sqrt(ss*(ss-l1)*(ss-l2)*(ss-l3)));
    total_surface+=area;
  }

  return total_surface;

}

double area_change(int oo,int ee){

    int p,q,j,k,i;
    int here,vo;
    double h,area,ss,areaOut=0.;
    double l1,l2,l3,ln;
    vec3D r,u1,u2,u3,n;
    vec3D dir;

    // CALCULATION OF TRIANGLE AREA USING HERON'S FORMULA
    // (that is, knowing the side lenghts of the triangles)
    // [ON] "oo" is the triangle and "ee" is the
    // vertex of that triangle that has been moved

    int is_p_edge, is_q_edge;

    k=Tr[oo].v[ee];

    here=oo;
    vo=ee;
    for (i=1;i<=Tmax;i++){

            p=Tr[here].v[(vo+1)%3];
            q=Tr[here].v[(vo+2)%3];
            
            u1.x=Elem[p].x-Elem[k].x;
            u1.y=Elem[p].y-Elem[k].y;
            u1.z=Elem[p].z-Elem[k].z;

            u2.x=Elem[q].x-Elem[k].x;
            u2.y=Elem[q].y-Elem[k].y;
            u2.z=Elem[q].z-Elem[k].z;

            u3.x=Elem[q].x-Elem[p].x;
            u3.y=Elem[q].y-Elem[p].y;
            u3.z=Elem[q].z-Elem[p].z;

            l1=sqrt(u1.x*u1.x +u1.y*u1.y +u1.z*u1.z);
            l2=sqrt(u2.x*u2.x +u2.y*u2.y +u2.z*u2.z);
            l3=sqrt(u3.x*u3.x +u3.y*u3.y +u3.z*u3.z);

            // semiperimeter of the triangle
            ss=.5*(l1+l2+l3);
            // application of heron's formula
            area=(sqrt(ss*(ss-l1)*(ss-l2)*(ss-l3)));
            areaOut+=area;

        // iterate over neighbouring triangles
        here=Tr[here].t[(vo+1)%3];;
        vo=4;
        if (Tr[here].v[0]==k) {vo=0;}
        if (Tr[here].v[1]==k) {vo=1;}
        if (Tr[here].v[2]==k) {vo=2;}

        if (vo==4){
        printf("ERROR: Problem with triangles in area_change function\n");
        exit(-1);
        }

        // if we go back to the initial triangle
        if (here==oo) {break;}
    }

    return areaOut;

}

double area_change_2tr(int a1,int b1,int c1,int a2,int b2,int c2){

  // Compute the total area limited by two triangles using heron's formula

  double h,area,ss,sArea=0.;
  double l1,l2,l3,ln;
  vec3D r,u1,u2,u3,n;
  vec3D dir;

    u1.x=Elem[b1].x-Elem[a1].x;
    u1.y=Elem[b1].y-Elem[a1].y;
    u1.z=Elem[b1].z-Elem[a1].z;

    u2.x=Elem[c1].x-Elem[a1].x;
    u2.y=Elem[c1].y-Elem[a1].y;
    u2.z=Elem[c1].z-Elem[a1].z;

    u3.x=Elem[c1].x-Elem[b1].x;
    u3.y=Elem[c1].y-Elem[b1].y;
    u3.z=Elem[c1].z-Elem[b1].z;

    l1=sqrt(u1.x*u1.x +u1.y*u1.y +u1.z*u1.z);
    l2=sqrt(u2.x*u2.x +u2.y*u2.y +u2.z*u2.z);
    l3=sqrt(u3.x*u3.x +u3.y*u3.y +u3.z*u3.z);

    ss=.5*(l1+l2+l3);
    area=(sqrt(ss*(ss-l1)*(ss-l2)*(ss-l3)));

    sArea+=area;

    /// second part

    u1.x=Elem[b2].x-Elem[a2].x;
    u1.y=Elem[b2].y-Elem[a2].y;
    u1.z=Elem[b2].z-Elem[a2].z;

    u2.x=Elem[c2].x-Elem[a2].x;
    u2.y=Elem[c2].y-Elem[a2].y;
    u2.z=Elem[c2].z-Elem[a2].z;

    u3.x=Elem[c2].x-Elem[b2].x;
    u3.y=Elem[c2].y-Elem[b2].y;
    u3.z=Elem[c2].z-Elem[b2].z;

    l1=sqrt(u1.x*u1.x +u1.y*u1.y +u1.z*u1.z);
    l2=sqrt(u2.x*u2.x +u2.y*u2.y +u2.z*u2.z);
    l3=sqrt(u3.x*u3.x +u3.y*u3.y +u3.z*u3.z);

    ss=.5*(l1+l2+l3);
    area=(sqrt(ss*(ss-l1)*(ss-l2)*(ss-l3)));

    sArea+=area;

    return sArea;
}

#endif /* Functions_h */
