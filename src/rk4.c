/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    rk4.c                                                       ********/
/******** @brief   Manages different aspects of the numerical integration      ********/
/******** @author  Jérémy COUTURIER <jeremycouturier.com>                      ********/
/********                                                                      ********/
/******** @section 	LICENSE                                                ********/
/******** Copyright (c) 2023 Jérémy COUTURIER                                  ********/
/********                                                                      ********/
/******** This file is part of NcorpiON                                        ********/
/********                                                                      ********/
/******** NcorpiON is free software. You can redistribute it and/or modify     ********/
/******** it under the terms of the GNU General Public License as published by ********/
/******** the Free Software Foundation, either version 3 of the License, or    ********/
/******** (at your option) any later version.                                  ********/
/********                                                                      ********/
/******** NcorpiON is distributed in the hope that it will be useful,          ********/
/******** but WITHOUT ANY WARRANTY; without even the implied warranty of       ********/
/******** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ********/
/******** GNU General Public License for more details.                         ********/
/********                                                                      ********/
/******** You should have received a copy of the GNU General Public License    ********/
/******** along with NcorpiON.  If not, see <http://www.gnu.org/licenses/>.    ********/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rk4.h"
#include "parameters.h"
#include "structure.h"
#include "physics.h"
#include "ffm.h"
#include "collision.h"
#include "display.h"
#include "spring.h"
#include <errno.h>
#include <math.h>
#include <string.h>

/******** Declaring external global files ********/
FILE ** files;

struct node * FlatTree;


void kick(struct moonlet * X, struct moonlet * C, void (*F)(struct moonlet *)){

      /******** Performs the kick of the leapfrog method   ********/
      /******** F is the vector field such that dX/dt=F(X) ********/
      
      int i;
      
      
      for(i = 0; i <= largest_id; i ++){            
            if(*(exists + i)){ //Checking whether or not there is a body in the i^th cell of the body array
            
                  /******** xx = X ********/
                  *(xx + i) = *(X + i);
            }
      }
      
      /******** xx = F(X) = dX/dt ********/
      (*F)(xx);
      
      /******** Taking care of tides and other dissipation ********/
      if (central_tides_bool){
            tides(X);
      }
      if (viscoelastic_bool){
            KelvinVoigtDamping(X);
      }
      
      /******** Applying the kick ********/
      for(i = 0; i <= largest_id; i ++){            
            if(*(exists + i)){ //Checking whether or not there is a body in the ith cell of the body array      
                  (X + i) -> vx += timestep * (xx + i) -> vx;
                  (X + i) -> vy += timestep * (xx + i) -> vy;
                  (X + i) -> vz += timestep * (xx + i) -> vz;
            }
      }
      if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){ //Applying the kick to the central body or the perturbing body
            C -> vx += timestep * CM_acc[0];  C -> vy += timestep * CM_acc[1];  C -> vz += timestep * CM_acc[2];
      }
}


void drift(struct moonlet * X, struct moonlet * C){

      /******** Performs the drift of the leapfrog method                          ********/
      /******** When this function is called, collision have already been resolved ********/
      
      int i;
      
      
      for(i = 0; i <= largest_id; i ++){            
            if(*(exists + i)){ //Checking whether or not there is a body in the i^th cell of the body array        
                  (X + i) -> x += timestep * (X + i) -> vx;
                  (X + i) -> y += timestep * (X + i) -> vy;
                  (X + i) -> z += timestep * (X + i) -> vz;
            }
      }
      if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){ //Applying the drift to the central body or to the perturbing body
            C -> x += timestep * C -> vx;  C -> y += timestep * C -> vy;  C -> z += timestep * C -> vz;
      }
}


void end_of_timestep(struct moonlet * moonlets, int progressed){

      /******** Takes care of the end of the timestep by reinitializing the hash table, trees and everything   ********/
      /******** that needs to be. Updates the mesh-size gam to match the desired expected number of neighbours ********/
      /******** Updates the timestep to match the mesh-size                                                    ********/
      /******** Prints useful informations in the terminal for the user to read                                ********/
      
      
      int index;
      int j;
      typ X, Y, Z, vX, vY, vZ, m, M;
      typ XX, YY, ZZ, DX, DY, DZ;
      typ com[6];
      typ alkhqp[6];
      
      first_passage = 1;

      /******** Removing from the simulation the bodies that need to be removed ********/
      how_many_moonlets   = 0;
      typ total_mass      = 0.0;
      typ M1, M2, R1, R2;                                 //Masses and radii of the two most massive bodies
      typ a1, a2, e1, e2, I1, I2, l1, l2, o1, o2, O1, O2; //Orbital elements of the two most massive bodies
      int i1, i2;                                         //Indexes          of the two most massive bodies
      M1 = 0.; M2 = 0.; i1 = 0; i2 = 0;
      int losing_that_one = 0;
      XX                  = (central_mass_bool ? CM.x : 0.);
      YY                  = (central_mass_bool ? CM.y : 0.);
      ZZ                  = (central_mass_bool ? CM.z : 0.);
      for (j = 0; j <= largest_id; j ++){
            if (*(exists + j)){
                  X  = (moonlets + j) -> x;
                  Y  = (moonlets + j) -> y;
                  Z  = (moonlets + j) -> z;
                  m  = (moonlets + j) -> mass;
                  DX = X - XX;  DY = Y - YY;  DZ = Z - ZZ;
                  if (DX*DX + DY*DY + DZ*DZ < disruption_threshold*disruption_threshold && central_mass_bool){ //Merging the body with the central body or the inner fluid disk
                        vX       = (moonlets + j) -> vx;
                        vY       = (moonlets + j) -> vy;
                        vZ       = (moonlets + j) -> vz;
                        CM.mass += m;
                        M        = (inner_fluid_disk_bool ? CM.mass + fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) : CM.mass);
                        CM.x    *= (M - m)/M;  CM.y *= (M - m)/M;  CM.z *= (M - m)/M;  CM.vx *= (M - m)/M;  CM.vy *= (M - m)/M;  CM.vz *= (M - m)/M;
                        CM.x    += m*X/M;      CM.y += m*Y/M;      CM.z += m*Z/M;      CM.vx += m*vX/M;     CM.vy += m*vY/M;     CM.vz += m*vZ/M;
                        if (inner_fluid_disk_bool && willMergeWithDisk(moonlets, j)){ //The mass of the dumped body is added to the inner fluid disk only if its periapsis
                              CM.mass -= m;                                           //is above the surface or if it will cross the xy plane before hitting the surface
                        }                                                             //Mass and outer radius are updated during call to function willMergeWithDisk
                        lose_moonlet(j);
                        losing_that_one = 1;
                  }
                  else if (DX*DX + DY*DY + DZ*DZ > high_dumping_threshold*high_dumping_threshold){ //Too far away
                        lose_moonlet(j);
                        losing_that_one = 1;
                        need_to_reduce_COM_bool = 1;
                  }
                  if (m > M1 && !losing_that_one){
                        M2 = M1;
                        i2 = i1;
                        M1 = m;
                        i1 = j;
                  }
                  else if (m > M2 && !losing_that_one){
                        M2 = m;
                        i2 = j;
                  }
                  if(!losing_that_one){
                        how_many_moonlets ++;
                        total_mass += m;
                  }
                  losing_that_one = 0;
            }
      }
      
      /******** Reducing the center of mass if needed (a super-catastrophic disruption happened or a body crossed the high_dumping_threshold) ********/
      if ((need_to_reduce_COM_bool && reduce_to_COM_bool) || progressed){
            total_momentum(moonlets, com); //Getting the speed and position of the center of mass
            if (need_to_reduce_COM_bool && reduce_to_COM_bool){
                  if (central_mass_bool){
                        CM.x -= com[0];  CM.y -= com[1];  CM.z -= com[2];  CM.vx -= com[3];  CM.vy -= com[4];  CM.vz -= com[5];
                  }
                  for (j = 0; j <= largest_id; j ++){
                        if (*(exists + j)){
                              (moonlets + j) -> x  -= com[0];
                              (moonlets + j) -> y  -= com[1];
                              (moonlets + j) -> z  -= com[2];
                              (moonlets + j) -> vx -= com[3];
                              (moonlets + j) -> vy -= com[4];
                              (moonlets + j) -> vz -= com[5];
                        }
                  }
                  need_to_reduce_COM_bool = 0;
                  com[0] = 0.;  com[1] = 0.;  com[2] = 0.;  com[3] = 0.;  com[4] = 0.;  com[5] = 0.;
            }
      }
      
      /******** If the simulation progressed by at least 0.1%, I display useful informations ********/
      if (progressed){
            
            /******** Drifting backward half a timestep before displaying. To be adapted to mesh at some point ********/
            if (central_mass_bool){
                  CM_buffer = CM;
                  CM_buffer.x -= 0.5*timestep*CM_buffer.vx;  CM_buffer.y -= 0.5*timestep*CM_buffer.vy;  CM_buffer.z -= 0.5*timestep*CM_buffer.vz;
            }
            
            /******** Displaying ********/
            if (collision_bool){
                  printf("                  N = %d  |  Collisions = %d  |  t = %.13lf   (largest id in the body array = %d)\n",
                  how_many_moonlets, collision_count, t_init + time_elapsed, largest_id);
            }
            else{
                  printf("                  N = %d  |  t = %.13lf   (largest id in the body array = %d)\n", how_many_moonlets, t_init + time_elapsed, largest_id);
            }
            if (!viscoelastic_bool){
                  if (collision_bool && fragmentation_bool && collision_count != 0){
                        typ cll = (typ) collision_count;
                        typ mrg = 100.0*((typ) merger_count)/cll;
                        typ spc = 100.0*((typ) super_catastrophic_count)/cll;
                        typ hfr = 100.0*((typ) half_fragmentation_count)/cll;
                        typ ffr = 100.0*((typ) full_fragmentation_count)/cll;
                        printf("                  Merger = %.2lf %% | Super-catastrophic = %.2lf %% | Partially fragmented = %.2lf %% | Fully fragmented = %.2lf %%\n",
                        mrg, spc, hfr, ffr);
                  }
                  if (inner_fluid_disk_bool){
                        typ disk_mass = M_PI*(Rout*Rout - R_unit*R_unit)*fluid_disk_Sigma;
                        printf("                  Central mass = %.8lf,  Bodies mass = %.8lf,  Fluid disk mass = %.8lf,  Bodies + disk = %.8lf\n",
                        CM.mass, total_mass, disk_mass, disk_mass + total_mass);
                  }
                  else{
                        if (central_mass_bool){
                              printf("                  Central mass = %.8lf,  Bodies mass = %.8lf\n", CM.mass, total_mass);
                        }
                        else{
                              printf("                  Bodies total mass = %.8lf\n", total_mass);
                        }
                  }
                  typ mu = central_mass_bool ? CM.mass : M_unit;  mu += inner_fluid_disk_bool ? fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) : 0.;  mu *= G;
                  if (how_many_moonlets > 1){
                        typ vp1, vp2;
                        struct moonlet mlt1 = {(moonlets + i1)->x, (moonlets + i1)->y, (moonlets + i1)->z, (moonlets + i1)->vx, (moonlets + i1)->vy, (moonlets + i1)->vz, 0., 0.};
                        struct moonlet mlt2 = {(moonlets + i2)->x, (moonlets + i2)->y, (moonlets + i2)->z, (moonlets + i2)->vx, (moonlets + i2)->vy, (moonlets + i2)->vz, 0., 0.};
                        mlt1.x -= 0.5*timestep*mlt1.vx;  mlt1.y -= 0.5*timestep*mlt1.vy;  mlt1.z -= 0.5*timestep*mlt1.vz;
                        mlt2.x -= 0.5*timestep*mlt2.vx;  mlt2.y -= 0.5*timestep*mlt2.vy;  mlt2.z -= 0.5*timestep*mlt2.vz;
                        R1  = (moonlets + i1) -> radius;  R2 = (moonlets + i2) -> radius;
                        cart2ell(&mlt1, 0, alkhqp, mu + G*M1);
                        a1  = *alkhqp; e1 = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]); I1 = 360.*asin(sqrt(alkhqp[4]*alkhqp[4] + alkhqp[5]*alkhqp[5]))/M_PI;
                        vp1 = 180.*atan2(alkhqp[3], alkhqp[2])/M_PI;
                        l1  = 180.*alkhqp[1]/M_PI - vp1; O1 = 180.*atan2(alkhqp[5], alkhqp[4])/M_PI; o1 = vp1 - O1;
                        cart2ell(&mlt2, 0, alkhqp, mu + G*M2);
                        a2  = *alkhqp; e2 = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]); I2 = 360.*asin(sqrt(alkhqp[4]*alkhqp[4] + alkhqp[5]*alkhqp[5]))/M_PI;
                        vp2 = 180.*atan2(alkhqp[3], alkhqp[2])/M_PI;
                        l2  = 180.*alkhqp[1]/M_PI - vp2; O2 = 180.*atan2(alkhqp[5], alkhqp[4])/M_PI; o2 = vp2 - O2;
                        printf("                  Most massive : (m, R; a, e, i; M, o, O) = (%.7lf, %.4lf; %.3lf, %.3lf, %.3lf°; %.2lf°, %.2lf°, %.2lf°)\n",
                        M1, R1, a1, e1, I1, l1, o1, O1);
                        printf("                  Second most  : (m, R; a, e, i; M, o, O) = (%.7lf, %.4lf; %.3lf, %.3lf, %.3lf°; %.2lf°, %.2lf°, %.2lf°)\n",
                        M2, R2, a2, e2, I2, l2, o2, O2);
                  }
                  else if (how_many_moonlets == 1){
                        typ vp1;
                        struct moonlet mlt1 = {(moonlets + i1)->x, (moonlets + i1)->y, (moonlets + i1)->z, (moonlets + i1)->vx, (moonlets + i1)->vy, (moonlets + i1)->vz, 0., 0.};
                        mlt1.x -= 0.5*timestep*mlt1.vx;  mlt1.y -= 0.5*timestep*mlt1.vy;  mlt1.z -= 0.5*timestep*mlt1.vz;
                        R1  = (moonlets + i1) -> radius;
                        cart2ell(&mlt1, 0, alkhqp, mu + G*M1);
                        a1  = *alkhqp; e1 = sqrt(alkhqp[2]*alkhqp[2] + alkhqp[3]*alkhqp[3]); I1 = 360.*asin(sqrt(alkhqp[4]*alkhqp[4] + alkhqp[5]*alkhqp[5]))/M_PI;
                        vp1 = 180.*atan2(alkhqp[3], alkhqp[2])/M_PI;
                        l1  = 180.*alkhqp[1]/M_PI - vp1; O1 = 180.*atan2(alkhqp[5], alkhqp[4])/M_PI; o1 = vp1 - O1;
                        printf("                  Most massive : (m, R; a, e, i; M, o, O) = (%.7lf, %.4lf; %.3lf, %.3lf, %.3lf°; %.2lf°, %.2lf°, %.2lf°)\n",
                        M1, R1, a1, e1, I1, l1, o1, O1);
                  }
                  if (central_tides_bool && inner_fluid_disk_bool){
                        printf("                  Sideral rotation = %.13lf,  Inner fluid disk outer edge = %.13lf\n", SideralOmega, Rout);
                  }
                  else if (central_tides_bool){
                        printf("                  Sideral rotation = %.13lf\n", SideralOmega);
                  }
                  else if (inner_fluid_disk_bool){
                        printf("                  Inner fluid disk outer edge = %.13lf\n", Rout);
                  }
                  if (pert_mass != 0. && central_mass_bool && pert_ecc <= 0.2){ //Getting the position of the evection resonance
                        typ Xin, Xout, MM, K;
                        int convergedIn, convergedOut;
                        MM   = (inner_fluid_disk_bool ? CM.mass + fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) : CM.mass);
                        K    = pert_sma*pow(1./(1. + pert_mass/MM), 1./3.);
                        Xout = pow(9./2., 1./3.);
                        if (J2_bool){
                              convergedIn  = evectionResonance(&Xin,  1);
                              convergedOut = evectionResonance(&Xout, 0);
                              if (convergedIn && convergedOut){
                                    printf("                  Inner evection resonance at a = %.3lf  |  Outer evection resonance at a = %.1lf\n", K/(Xin*Xin), K/(Xout*Xout));
                              }
                              else if (convergedIn && !convergedOut){
                                    printf("                  Inner evection resonance at a = %.3lf  |  Cannot find outer evection resonance\n", K/(Xin*Xin));
                              }
                              else if (!convergedIn && convergedOut){
                                    printf("                  Cannot find inner evection resonance   |  Outer evection resonance at a = %.1lf\n", K/(Xout*Xout));
                              }
                              else{
                                    printf("                  Cannot find inner evection resonance   |  Cannot find outer evection resonance\n");
                              }
                        }
                        else{
                              printf("                  Outer evection resonance at a = %.2lf\n", K/(Xout*Xout));
                        }
                  }
            }
            printf("                  COM = (%.9lf, %.9lf,  %.9lf,  %.9lf,  %.9lf,  %.9lf)\n", com[0], com[1], com[2], com[3], com[4], com[5]);
      }
      
      /******** Resetting global variables relative to the boxdot tree ********/
      if ((mutual_bool || collision_bool) && !brute_force_bool && !force_naive_bool && (falcON_bool || standard_tree_bool)){
            for (j = 0; j < cell_id; j ++){
                  free((FlatTree + j) -> dots);
                  (FlatTree + j) -> dots = NULL;
            }
            free(FlatTree);
            FlatTree       = NULL;
            how_many_cells = 0;
            cell_id        = 0;
      }
      
      /******** Reinitializing the array did_collide ********/
      if (collision_bool && (one_collision_only_bool || fragmentation_bool)){
            for (j = 0; j <= largest_id; j ++){
                  *(did_collide + j) = 0;
            }
      }
      
      /******** Reinitializing data relative to the mesh algorithm ********/
      if (mesh_bool && (collision_bool || mutual_bool) && !force_naive_bool){
      
            how_many_pairs = 0; //It is unnecessary to reinitialize the array pairs
            /******** Reinitializing the hash table ********/
            int p;
            for (p = 0; p < how_many_modified; p ++){
                  index = *(modified_cells + p);
                  clear_chain(hash + index); //Reinitializing to NULL the index^th cell of the hash table
            }
      
            /******** I update the size of the mesh ********/
            average_neighbours = 2.0*((typ) total_neighbours)/((typ) how_many_moonlets);
            gam *= pow(how_many_neighbours/average_neighbours, 1.0/3.0);
            if (gam < gam_min){
                  gam = gam_min;
            }
            collision_cube = gam*((typ) collision_cube_cells);
      
            /******** If the simulation progressed by at least 0.1% since last display, I display useful informations ********/
            if (progressed){
                  printf("                  average number of neighbours = %.3lf\n", average_neighbours);
                  printf("                  gamma = %.6lf\n", gam);
            }
      
            /******** Reinitializing some global variables ********/
            how_many_modified = 0; //It is unnecessary to reinitialize the array modified_cells
            how_many_big      = 0;
            how_many_small    = 0;
            total_neighbours  = 0;
      }
      
      /******** If there are very few bodies, I switch to the brute-force algorithm ********/
      if      (!brute_force_bool && !force_naive_bool && how_many_moonlets < switch_to_brute_force){ //Switching to the brute-force algorithm for mutual interactions
            force_naive_bool = 1;
      }
      else if (!brute_force_bool &&  force_naive_bool && how_many_moonlets > switch_to_brute_force){ //Switching back to a more efficient algorithm
            force_naive_bool = 0;
      }
      
      /******** Reordering the array moonlets ********/
      if (how_many_moonlets < 9*largest_id/10){ //I push the bodies to the left of the array if there are too many holes
            tidy_up(moonlets);
            if ((mutual_bool || collision_bool) && (falcON_bool || standard_tree_bool) && !brute_force_bool && !force_naive_bool){
                  IndexPeanoHilbertOrder = largest_id + 1;
                  for (j = 0; j < IndexPeanoHilbertOrder; j ++){
                        PeanoHilbertOrder[j] = j;
                  }
            }
      }
      
      /******** Spawning bodies from the inner fluid disk ********/
      if (inner_fluid_disk_bool){
            M             = CM.mass + fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit);
            typ x         = Rout/R_unit;
            typ x52       = fast_pow(sqrt(x), 5) - 1.;
            typ x32       = fast_pow(sqrt(x), 3);
            typ x2        = x*x - 1.;
            typ gx        = (4.*x*x52 - 5.*x2*x32)/(4.*x52 - 5.*x2);
            typ Omega_in  = sqrt(G*M/(R_unit*R_unit*R_unit));
            typ dotMinner = M_PI*M_PI*M_PI*G*G/(Omega_in*Omega_in*Omega_in*(x - 1.)*(1. - gx)); //From Salmon & Canup 2012
            typ dotMouter = -fast_pow(sqrt(x), 11)*dotMinner/gx;      
            /******** Removing from the disk the mass that flows out from the inner edge ********/
            typ flowed_mass   = timestep*dotMinner*fluid_disk_Sigma*fluid_disk_Sigma*fluid_disk_Sigma;
            fluid_disk_Sigma -= flowed_mass/(M_PI*(Rout*Rout - R_unit*R_unit));
            CM.mass          += flowed_mass;
            /******** If the mass that flowed out of the inner fluid disk by the outer edge since the last time a body spawned exceeds the spawn mass, I spawn a body ********/
            flowed_since_last_spawn += timestep*dotMouter*fluid_disk_Sigma*fluid_disk_Sigma*fluid_disk_Sigma;
            typ spawned_mass         = 16.*fast_pow(M_PI, 4)*f_tilde*f_tilde*fast_pow(fluid_disk_Sigma*Rout*Rout, 3)/(M_unit*M_unit);
            spawned_mass             = spawned_mass > 1.0e-7*M_unit ? spawned_mass : 1.0e-7*M_unit; //Threshold needed for very long simulations
            if (flowed_since_last_spawn > spawned_mass){
                  typ cart[6];
                  flowed_since_last_spawn  -= spawned_mass;
                  typ rad                   = pow(3.0*spawned_mass/(4.0*M_PI*spawned_density), 1.0/3.0);
                  int index                 = get_free_index(0); //Retrieving an index in the bodies array for the new body
                  typ nu                    = rdm(0.0, 2.0*M_PI);
                  typ omega                 = rdm(0.0, 2.0*M_PI);
                  typ Omega                 = rdm(0.0, 2.0*M_PI);
                  ell2cart(Rout, 0., 0., nu, omega, Omega, G*(M + spawned_mass), cart);   //Cartesian coordinate of the spawned body in the geocentric reference frame
                  typ disk_mass             = fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit);
                  typ disk_angular_momentum = 0.8*fluid_disk_Sigma*M_PI*(Rout*Rout*sqrt(Rout) - R_unit*R_unit*sqrt(R_unit));
                  innerFluidDiskAngularMomentum(-spawned_mass, Rout, 0., 1., disk_angular_momentum, disk_mass);
                  X  = cart[0] + CM.x;
                  Y  = cart[1] + CM.y;
                  Z  = cart[2] + CM.z;
                  vX = cart[3] + CM.vx;
                  vY = cart[4] + CM.vy;
                  vZ = cart[5] + CM.vz;
                  /******* Spawning a body on an equatorial and circular orbit with respect to the central body ********/
                  (moonlets + index) -> x      = X;
                  (moonlets + index) -> y      = Y;
                  (moonlets + index) -> z      = Z;
                  (moonlets + index) -> vx     = vX;
                  (moonlets + index) -> vy     = vY;
                  (moonlets + index) -> vz     = vZ;
                  (moonlets + index) -> mass   = spawned_mass;
                  (moonlets + index) -> radius = rad;
                  /********  Taking the body's momentum away from Earth for conservation ********/
                  typ mf = spawned_mass;
                  CM.x  *= (M + mf)/M;  CM.y *= (M + mf)/M;  CM.z *= (M + mf)/M;  CM.vx *= (M + mf)/M;  CM.vy *= (M + mf)/M;  CM.vz *= (M + mf)/M;
                  CM.x  -= mf*X/M;      CM.y -= mf*Y/M;      CM.z -= mf*Z/M;      CM.vx -= mf*vX/M;     CM.vy -= mf*vY/M;     CM.vz -= mf*vZ/M;
            }
      }
      
      /******** Updating the J2 and the position of the evection resonance ********/
      if (J2_bool){
            if (J2_value == 0.){
                  J2 = 0.5*SideralOmega*SideralOmega*R_unit*R_unit*R_unit/(G*CM.mass);
            }
            if (pert_mass > 0.){
                  evection_resonance = pow(1.5*sqrt(CM.mass/pert_mass)*J2, 2./7.)*pow(pert_sma/R_unit, 3./7.)*R_unit;
            }
      }
      
      /******** Communicating with REBOUND for 3D visualization ********/
      if (openGL_bool){
            rebound(moonlets);
      }
}


void integration_tree(typ t){


      /******** Performs the numerical integration using a tree algorithm (falcON or standard tree) ********/
      /******** for mutual interactions treatment. t is the total time.                             ********/


      /******** Initializing the array of bodies ********/
      struct moonlet * moonlets = (viscoelastic_bool ? generate_visco_elastic_body() : populate());
      struct moonlet * moonlet_buffer = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //Buffer for output steps
      if (moonlet_buffer == NULL){
            fprintf(stderr, "Error : Can't allocate buffer for array of bodies in function integration_tree.\n");
            abort();
      }
      
      /******** Numerical integration ********/
      int error             = 0;
      bigint iter           = 0;
      int progressed        = 0;
      typ progress          = 0.0;
      typ previous_progress = 0.0;
      int j;
      struct boxdot * root = NULL;
      
      printf("progress = %.1lf %%\n", 0.0);
      if (write_to_files_bool){
            readme();
            if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
                  CM_buffer = CM;
            }
            display(moonlets);
      }
      
      /******** Performing half a drift ********/
      if (t > 0.0){
            timestep /= 2.0;
            if (collision_bool){ //Finding and resolving collisions and going backward
                  root     = root_cell(moonlets);
                  FlatTree = flattree_init(root);
                  clear_boxdot(&root);
                  center_and_maxR_flattree(FlatTree, moonlets);
                  rmax_and_rcrit_flattree (FlatTree, moonlets);
                  if (falcON_bool){
                        collision_flattree(FlatTree, moonlets);
                  }
                  else if (standard_tree_bool){
                        for (j = 0; j <= largest_id; j++){
                              if (exists[j] && !(did_collide[j])){
                                    standard_tree_collision(FlatTree, moonlets, j);
                              }
                        }
                  }
            }
            drift(moonlets, &CM); //Drifting
            timestep *= 2.0;

            /******** Resetting data relative to trees ********/
            if (collision_bool){
                  for (j = 0; j < cell_id; j ++){
                        free((FlatTree + j) -> dots);
                        (FlatTree + j) -> dots = NULL;
                  }
                  free(FlatTree);
                  FlatTree       = NULL;
                  how_many_cells = 0;
                  cell_id        = 0;
                  /******** Reinitializing the array did_collide ********/
                  if (one_collision_only_bool){
                        for (j = 0; j <= largest_id; j ++){
                              *(did_collide + j) = 0;
                        }
                  }
            }
      }

      while(time_elapsed < t && error == 0){
      
            progressed = 0;

            /******** Writing the results of the numerical integration to the output files ********/
            if (write_to_files_bool && iter % output_step == 0 && iter > 0){

                  /******** I drifted too far at the previous timestep, I have to undrift by half a timestep ********/
                  for (j = 0; j <= largest_id; j++){
                        if (*(exists + j)){
                              *(moonlet_buffer + j) = *(moonlets + j);
                        }
                  }
                  if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
                        CM_buffer = CM;
                  }
                  timestep /= -2.0;
                  drift(moonlet_buffer, &CM_buffer); //Drifting backwards half a timestep
                  timestep *= -2.0;
                  display(moonlet_buffer);
            }

            /******** Performing a full kick. Integrator is SABA1 ********/
            if (!force_naive_bool && mutual_bool){
                  root     = root_cell(moonlets);
                  FlatTree = flattree_init(root);
                  clear_boxdot(&root);
                  com_flattree(FlatTree, moonlets);
                  Mtot     = FlatTree -> M0;
                  rmax_flattree(FlatTree, moonlets);
                  rcrit_flattree(FlatTree);
                  tensor_initialization(FlatTree);
                  #if expansion_order > 2
                        multipole_flattree(FlatTree, moonlets);
                  #endif
            }        
            kick(moonlets, &CM, vector_field);
            time_elapsed += 0.5*timestep;

            /******** Performing a full drift ********/
            if (collision_bool){
                  if (force_naive_bool){ //Finding and resolving collisions and going backward brute-forcely
                        brute_force(moonlets);
                  }
                  else{
                        if (!mutual_bool){ //The tree was not previously made and must be made now
                              root     = root_cell(moonlets);
                              FlatTree = flattree_init(root);
                              clear_boxdot(&root);
                        }
                        center_and_maxR_flattree(FlatTree, moonlets);
                        rmax_and_rcrit_flattree (FlatTree, moonlets);
                        if (falcON_bool){ //Finding and resolving collisions and going backward with falcON                        
                              collision_flattree(FlatTree, moonlets);
                        }
                        else if (standard_tree_bool){ //Finding and resolving collisions and going backward with the standard tree code
                              for (j = 0; j <= largest_id; j ++){
                                    if (exists[j] && !(did_collide[j])){
                                          standard_tree_collision(FlatTree, moonlets, j);
                                    }
                              }
                        }
                  }
            }            
            drift(moonlets, &CM);

            iter ++;
            time_elapsed += 0.5*timestep;
            progress = time_elapsed/t;
            
            /******** Displaying informations in the terminal, and reinitializing what needs to be reinitialized after every timestep ********/
            if (progress - previous_progress > 0.001){
                  previous_progress = progress;
                  printf("progress = %.1lf %%\n", 100.0*progress);
                  progressed = 1;
            }
            
            /******** Taking care of the end of the timestep ********/
            end_of_timestep(moonlets, progressed);
            
      }
      
      /******** Last file saving ********/
      if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
            CM_buffer = CM;
      }
      if (t > 0.0){
            timestep /= -2.0;
            drift(moonlets, &CM_buffer); //Performing half a backward drift
            timestep *= -2.0;
      }
      if (write_to_files_bool){
            display(moonlets);
      }
      if (resume_simulation_bool){
            resume(moonlets);
            resumeFile();
      }
      
      /******** Displaying the number of timesteps performed ********/
      printf("progress = %.1lf %%\n", 100.0);
      printf("total number of timestep performed : %ld\n", iter);
      
      /******** Deallocating the array of bodies ********/
      free(moonlets);
      moonlets = NULL;
      free(moonlet_buffer);
      moonlet_buffer = NULL;
}


void integration_mesh(typ t){


      /******** Performs the numerical integration using a mesh algorithm ********/
      /******** for mutual interactions treatment. t is the total time.   ********/


      /******** Initializing the array of bodies ********/
      struct moonlet * moonlets = (viscoelastic_bool ? generate_visco_elastic_body() : populate());
      struct moonlet * moonlet_buffer = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //Buffer for output steps
      if (moonlet_buffer == NULL){
            fprintf(stderr, "Error : Can't allocate buffer for array of bodies in function integration_brute_force.\n");
            abort();
      }
      
      
      /******** Numerical integration ********/
      int error = 0;
      bigint iter = 0;
      int progressed = 0;
      typ progress = 0.0;
      typ previous_progress = 0.0;
      int index, j;
      
      printf("progress = %.1lf %%\n", 0.0);
      if (write_to_files_bool){
            readme();
            if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
                  CM_buffer = CM;
            }
            display(moonlets);
      }
      
      /******** Performing half a kick ********/
      if (t > 0.0){
            timestep /= 2.0;
            three_largest_moonlets(moonlets);
            three_largest_three_first(moonlets);
            get_neighbours_mesh(moonlets);
            kick(moonlets, &CM, vector_field);
            timestep *= 2.0;
      
            /******** Reinitializing data relative to the hash table ********/
            for (j = 0; j < how_many_modified; j++){
                  index = *(modified_cells + j);
                  clear_chain(hash + index); //Reinitializing to NULL the index^th cell of the hash table
            }
            average_neighbours = 2.0*((typ) total_neighbours)/((typ) how_many_moonlets);
            gam *= pow(how_many_neighbours/average_neighbours, 1.0/3.0); //Updating the mesh-size gamma
            if (gam < gam_min){
                  gam = gam_min;
            }
            collision_cube    = gam*((typ) collision_cube_cells);
            how_many_pairs    = 0; //It is unnecessary to reinitialize the array pairs
            how_many_modified = 0; //It is unnecessary to reinitialize the array modified_cells
            total_neighbours  = 0;
      }

      while(time_elapsed < t && error == 0){

            progressed = 0;
                   
            
            /******** Writing the results of the numerical integration to the output files ********/
            if (write_to_files_bool && iter % output_step == 0 && iter > 0){

                  /******** I kicked too far at the previous timestep, I have to unkick by half a timestep ********/
                  for (j = 0; j <= largest_id; j ++){
                        if (*(exists + j)){
                              *(moonlet_buffer + j) = *(moonlets + j);
                        }
                  }
                  if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
                        CM_buffer = CM;
                  }
                  timestep /= -2.0;
                  if (!force_naive_bool){
                        three_largest_moonlets(moonlet_buffer);
                        three_largest_three_first(moonlet_buffer);
                        get_neighbours_mesh(moonlet_buffer);
                  }
                  kick(moonlet_buffer, &CM_buffer, vector_field); //Performing half a backward kick
                  timestep *= -2.0;
                  display(moonlet_buffer);
                  
                  /******** Reinitializing data relative to the hash table ********/
                  if (!force_naive_bool){
                        for (j = 0; j < how_many_modified; j ++){
                              index = *(modified_cells + j);
                              clear_chain(hash + index); //Reinitializing to NULL the index^th cell of the hash table
                        }
                        how_many_pairs    = 0; //It is unnecessary to reinitialize the array pairs
                        how_many_modified = 0; //It is unnecessary to reinitialize the array modified_cells
                        total_neighbours  = 0;
                  }
            }

            /******** Integrator is SBAB1 ********/
            if (collision_bool){
                  if (force_naive_bool){
                        brute_force(moonlets); //Resolving collisions and going backward
                  }
                  else{
                        three_largest_moonlets(moonlets);
                        three_largest_three_first(moonlets);
                        mesh(moonlets);
                  }
            }
            drift(moonlets, &CM); //Performing a full drift
            time_elapsed += 0.5*timestep;
            
            if (!collision_bool && !force_naive_bool){
                  three_largest_moonlets(moonlets);
                  three_largest_three_first(moonlets);
                  get_neighbours_mesh(moonlets);
            }
            kick(moonlets, &CM, vector_field); //Performing a full kick

            iter ++;
            time_elapsed += 0.5*timestep;
            progress = time_elapsed/t;
            
            /******** Displaying informations in the terminal, and reinitializing what needs to be reinitialized after every timestep ********/
            if (progress - previous_progress > 0.001){
                  previous_progress = progress;
                  printf("progress = %.1lf %%\n", 100.0*progress);
                  progressed = 1;
            }
            
            /******** Taking care of the end of the timestep ********/
            end_of_timestep(moonlets, progressed);
            
            
      }
      
      /******** Last file saving ********/
      if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
            CM_buffer = CM;
      }
      if (t > 0.0){
            timestep /= -2.0;
            if (!force_naive_bool){
                  three_largest_moonlets(moonlets);
                  three_largest_three_first(moonlets);
                  get_neighbours_mesh(moonlets);
            }
            kick(moonlets, &CM_buffer, vector_field); //Performing half a backward kick
            timestep *= -2.0;
      }
      if (write_to_files_bool){
            display(moonlets);
      }
      if (resume_simulation_bool){
            resume(moonlets);
            resumeFile();
      }           
      if (!force_naive_bool && t > 0.0){
            for (j = 0; j < how_many_modified; j ++){
                  index = *(modified_cells + j);
                  clear_chain(hash + index); //Reinitializing to NULL the index^th cell of the hash table
            }
            how_many_pairs    = 0; //It is unnecessary to reinitialize the array pairs
            how_many_modified = 0; //It is unnecessary to reinitialize the array modified_cells
            total_neighbours  = 0;
      }
      
      /******** Displaying the number of timesteps performed ********/
      printf("progress = %.1lf %%\n", 100.0);
      printf("total number of timestep performed : %ld\n", iter);

      /******** Deallocating the array of bodies ********/
      free(moonlets);
      moonlets = NULL;
      free(moonlet_buffer);
      moonlet_buffer = NULL;
}


void integration_brute_force_SABA1(typ t){


      /******** Performs the numerical integration with a brute force method. t is the total time ********/


      /******** Initializing the array of bodies ********/
      struct moonlet * moonlets = (viscoelastic_bool ? generate_visco_elastic_body() : populate());
      struct moonlet * moonlet_buffer = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //Buffer for output steps
      if (moonlet_buffer == NULL){
            fprintf(stderr, "Error : Can't allocate buffer for array of bodies in function integration_brute_force.\n");
            abort();
      }
      
      
      /******** Numerical integration ********/
      int error = 0;
      bigint iter = 0;
      int progressed = 0;
      typ progress = 0.0;
      typ previous_progress = 0.0;
      int j;
      
      printf("progress = %.1lf %%\n", 0.0);
      if (write_to_files_bool){
            readme();
            if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
                  CM_buffer = CM;
            }
            display(moonlets);
      }
      
      /******** Performing half a drift ********/
      if (t > 0.0){
            timestep /= 2.0;
            if (collision_bool){
                  brute_force(moonlets);  //Resolving collisions and going backward
            }
            drift(moonlets, &CM);
            timestep *= 2.0;
      }

      while(time_elapsed < t && error == 0){

            progressed = 0;
                   
            
            /******** Writing the results of the numerical integration to the output files ********/
            if (write_to_files_bool && iter % output_step == 0 && iter > 0){
                  
                  /******** I drifted too far at the previous timestep, I have to undrift by half a timestep ********/
                  for (j = 0; j <= largest_id; j ++){
                        if (*(exists + j)){
                              *(moonlet_buffer + j) = *(moonlets + j);
                        }
                  }
                  if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
                        CM_buffer = CM;
                  }
                  timestep /= -2.0;
                  drift(moonlet_buffer, &CM_buffer); //Performing half a backward drift
                  timestep *= -2.0;
                  display(moonlet_buffer);
            }

            /******** Integrator is SBAB1 ********/
            kick(moonlets, &CM, vector_field); //Performing a full kick            
            time_elapsed += 0.5*timestep;            
            if (collision_bool){
                  brute_force(moonlets);  //Resolving collisions and going backward
            }           
            drift(moonlets, &CM);         //Performing a full drift

            iter ++;
            time_elapsed += 0.5*timestep;
            progress = time_elapsed/t;
            
            /******** Displaying informations in the terminal, and reinitializing what needs to be reinitialized after every timestep ********/
            if (progress - previous_progress > 0.001){
                  previous_progress = progress;
                  printf("progress = %.1lf %%\n", 100.0*progress);
                  progressed = 1;
            }
            
            /******** Taking care of the end of the timestep ********/
            end_of_timestep(moonlets, progressed);      
      }
      
      /******** Last file saving ********/
      if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
            CM_buffer = CM;
      }
      if (t > 0.0){
            timestep /= -2.0;
            drift(moonlets, &CM_buffer); //Performing half a backward drift
            timestep *= -2.0;
      }
      if (write_to_files_bool){
            display(moonlets);
      }
      if (resume_simulation_bool){
            resume(moonlets);
            resumeFile();
      }
      
      /******** Displaying the number of timesteps performed ********/
      printf("progress = %.1lf %%\n", 100.0);
      printf("total number of timestep performed : %ld\n", iter);
      
      /******** Deallocating the array of bodies ********/
      free(moonlets);
      moonlets = NULL;
      free(moonlet_buffer);
      moonlet_buffer = NULL;
}
