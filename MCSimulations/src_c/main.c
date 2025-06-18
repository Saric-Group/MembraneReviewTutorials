// *****************************************************************************
//
//  DYNAMICALLY TRIANGULATED MEMBRANE -- SARIC GROUP (ISTA)
//
//  This code simulates a triangulated membrane
//  using a Monte Carlo algorithm (w/ Metropolis-Hastings criterion for the moves).
//  There are two update moves implemented for the membrane:
//  1. Move the vertex
//  2. Switch bond connectivity (to allow for fluidity)
//      - One can switch off this move and simulate an elastic membrane
//        (use is_mem_fluid = 0 flag in in.ves)
//
//  The code below is designed to simulate the pulling of a tether from a flat 
//  membrane patch by moving a bead that is connected to a single membrane vertex. 
//  The vertices at the edges of the membrane patch do not move, 
//  and the bonds that connect them to the rest of the lattice are not updated.
//  See the README.md associated to this code to learn further details about the boundary
//  conditions of the system.
//
//  ** HOW TO COMPILE **
//
//  To compile this code:
//  - Go to directory where *.c and *.h files are
//  - Run gcc -o EXENAME *.c
//  (Additional flags, like -O3, can be used -- RECOMMENDED FOR SPEED)
//  To compile at a HPC cluster, you might have to run:
//  - gcc *.c -lm -o EXENAME
//  (You can also include the -O3 flag for speed)
//
// *****************************************************************************
//  - The current clean-up and development of the code is maintained by 
//    Maitane Munoz-Basagoiti (maitane.munoz-basagoiti@ista.ac.at)
// *****************************************************************************

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#include "PreprocessorDeclarations.h"
#include "DataStructures.h"
#include "Functions.h"

// *****************************************************************************
// MAIN
// *****************************************************************************
int main(int argc, const char * argv[]) {

    printf("** MC code for pulling a tether from a membrane patch **\n");
    printf("|| Parameters read from in.ves file || \n");

    int i =0;
    int j,k,o,coll;
    double stv_ratio;

    // set all the energy changes to zero
    S.EbendChange = 0;
    S.EbindChange_colloid_bead = 0;
    S.EbindChange_bead_colloid = 0;
    S.Evolchange = 0;
    S.Esurfchange = 0;

    // reading the input data
    FILE *in0;
    in0=fopen("in.ves","r");    
    fscanf(in0,"%lf %lf %d %d %d %d %lf %d %d %d %lf %lf %lf %lf %d %d %d\n",
        &S.kappa, &S.gamma, &N, &Ncoll, &Nedges, &Nbound, &Coll.sig,
        &seed, &is_mem_fluid, &start_from_config, &speed_pulling, &position_bead,
        &k_harmonic, &r_eq, &LX, &LY, &LZ);

    printf("** Total number of particles: %d **\n", N+Ncoll);
    printf("Parameters:\n - Bending modulus: %lf\n - Surface tension eq.: %lf\n - Nmembrane: %d\n - Nnanops: %d\n - Nedges: %d\n - Nbound: %d\n - Sigma nanop: %lf\n - Int. cutoff: %lf\n - Seed: %d\n - Is membrane fluid? 1 = YES, 0 = NO --> %d\n - k_harmonic: %lf\n - r_eq: %lf\n - speed_pulling: %lf\n - LX: %d - LY: %d - LZ: %d\n", S.kappa,S.gamma,N,Ncoll, Nedges, Nbound, Coll.sig, Coll.rc, seed, is_mem_fluid, k_harmonic, r_eq, speed_pulling, LX, LY, LZ);

    // ----------- WHEN COLLOIDS PRESENT: define MC move for colloids, interaction cutoff (assumes sig_BEAD = 1), verlet cutoff...
    // (not entirely applicable for the pulling experiment)
    Coll.dis_mc =0.1;
    Coll.rad    =0.5*Coll.sig;
    Coll.rc     =Coll.rc*(0.5+0.5*Coll.sig);
    Coll.rc2    =Coll.rc*Coll.rc;
    Coll.rv     =1.6*Coll.rc;
    Coll.rv2    =Coll.rv*Coll.rv;
    Coll.sig2   =Coll.sig*Coll.sig;
    Coll.update_bcoll_list = (Coll.rv-Coll.rc)*0.5;
    factor_LJ = 0.5*(1+Coll.sig)*0.5*(1+Coll.sig);

    // ----------- define system lattice
    lattice=1.69;                     // It is used to decide the neighbours at the beginning
    lattice_y=lattice*sqrt(3.0)/2.0;  // Important that you make sure your bonds are shorter than this!

    FILE *acceptfile;

    // ----------- IF DESIRED clean output files to avoid appending to them
    if(start_from_config==0){

        // trajectory of the system
        char iiiii[500];
        sprintf(iiiii,"out.dump");
        FILE *clean;
        clean = fopen(iiiii, "w");
        fclose(clean);

        // bonds of the system
        sprintf(iiiii,"bonds.dump");
        FILE *bonds;
        bonds = fopen(iiiii, "w");
        fclose(bonds);

        // energy changes etc in the system
        sprintf(iiiii,"energy.dump");
        FILE *energyfile;
        energyfile = fopen(iiiii, "w");
        fclose(energyfile);

        // position pulling bead
        sprintf(iiiii,"position_pulling_bead.dat");
        FILE *pullingbeadfile;
        pullingbeadfile = fopen(iiiii, "w");
        fclose(pullingbeadfile);

    }

    // ----------- set configuration of the system (read initial conditions, neighbour lists...)
    printf("- Setting initial configuration...\n");
    set_all();
    printf("! Done.\n");

    // ----------- target area and volume (computed inside set_all; not fully applicable to current simulation)
    A0 = S.SurfArea;
    V0 = S.Svol;

    // ----------- printing initial details
    counter = 0;
    printf("C = (%2.3f, %2.3f, %2.3f) \n",Centre.x,Centre.y,Centre.z);
    printf("Initial membrane surface : %lf\n", S.SurfArea);
    printf("Initial vesicle volume  : %lf\n", S.Svol);
    printf("Initial vesicle radius  : %lf\n", S.rad);
    
    // get the first snapshot
    painter(1);
    write_config();
    painter_en(1, 0, N, Ncoll,Ntri);

    //counting the accepted moves
    accepted_moves = 0;

    // ----------- MC sweeps
    printf("- Starting MC sweeps...\n");
    printf("Step: %d || Surface: %2.3f || Ebend: %2.3f || Esurface: %2.3f || Ebindpulling: %2.3f || C=(%2.3f, %2.3f, %2.3f)\n", i, S.SurfArea, S.EbendChange, S.Esurfchange, S.EbindChange_bead_colloid, Centre.x,Centre.y,Centre.z);
    double max_position = Elem[N+Ncoll].z+position_bead;
    
    for (i=1;i<=MCSTEPS;i++){

        if(i>MCEQUILIBRATE){
            // only move pulling bead after membrane equilibration
            // in pulling experiments, bead moves in a direction at a specific speed
            if (Elem[N+Ncoll].z < (max_position)){
                Elem[N+Ncoll].z += speed_pulling;
            }
        }

        // update position membrane beads
        for (j=1;j<=N;j++){
              MC_xyz( 1+(int)( drand48()*Ntri ) );
        }

        // swap bonds membrane beads to achieve fluidity
        for (j=1;j<=N;j++){
            switch_bond(1+(int)(drand48()*Ntri));
        }

        counter+=1;

        // write down configuration/evolution of the system
        if (i%WRITE_CONF==0 && i>1){
            printf("Step: %d || Surface: %2.3f || Ebend: %2.3f || Esurface: %2.3f || Ebindpulling: %2.3f || C=(%2.3f, %2.3f, %2.3f)\n", i, S.SurfArea, S.EbendChange, S.Esurfchange, S.EbindChange_bead_colloid, Centre.x,Centre.y,Centre.z);
            // writes configuration so that system can be restarted
            write_config();
            // writes "time"-evolution of system (along MC sweeps)
            painter(i);
            // writes the energy changes of the system
            painter_en(i, counter, N, Ncoll,Ntri);
            // print trajectories of the pulling bead
            painter_nanoparticles();
        }
    }

    printf("! Done.\n");

    // writes down the last configuration
    write_config();
    return 0;
}
