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
      typ X, Y, Z, vX, vY, vZ, m, R, M;
      typ XX, YY, ZZ, DX, DY, DZ;
      typ com[6];
      
      first_passage = 1;

      /******** Removing from the simulation the bodies that need to be removed ********/
      how_many_moonlets = 0;
      typ total_mass    = 0.0;
      typ maxR          = 0.0;
      XX                = (central_mass_bool ? CM.x : 0.);
      YY                = (central_mass_bool ? CM.y : 0.);
      ZZ                = (central_mass_bool ? CM.z : 0.);
      for (j = 0; j <= largest_id; j++){
            if (*(exists + j)){
                  X  = (moonlets + j) -> x;
                  Y  = (moonlets + j) -> y;
                  Z  = (moonlets + j) -> z;
                  R  = (moonlets + j) -> radius;
                  m  = (moonlets + j) -> mass;
                  DX = X - XX;  DY = Y - YY;  DZ = Z - ZZ;
                  if (DX*DX + DY*DY + DZ*DZ < low_dumping_threshold*low_dumping_threshold && central_mass_bool){ //Merging the body with the central body or the inner fluid disk
                        vX       = (moonlets + j) -> vx;
                        vY       = (moonlets + j) -> vy;
                        vZ       = (moonlets + j) -> vz;
                        CM.mass += m;
                        M        = (inner_fluid_disk_bool ? CM.mass + fluid_disk_Sigma*M_PI*(Rroche*Rroche - R_unit*R_unit) : CM.mass);
                        CM.x    *= (M - m)/M;  CM.y *= (M - m)/M;  CM.z *= (M - m)/M;  CM.vx *= (M - m)/M;  CM.vy *= (M - m)/M;  CM.vz *= (M - m)/M;
                        CM.x    += m*X/M;      CM.y += m*Y/M;      CM.z += m*Z/M;      CM.vx += m*vX/M;     CM.vy += m*vY/M;     CM.vz += m*vZ/M;
                        if (inner_fluid_disk_bool){
                              fluid_disk_Sigma += m/(M_PI*(Rroche*Rroche - R_unit*R_unit)); //The mass of the dumped body is added to the inner fluid disk
                              CM.mass -= m;
                        }
                        lose_moonlet(j);                      
                  }
                  else if (DX*DX + DY*DY + DZ*DZ > high_dumping_threshold*high_dumping_threshold){
                        lose_moonlet(j);
                        need_to_reduce_COM_bool = 1;
                  }
                  if (R > maxR){
                        maxR = R;
                  }
                  how_many_moonlets ++;
                  total_mass += m;
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
            printf("                  N = %d   (largest id in the body array = %d)\n", how_many_moonlets, largest_id);
            if (inner_fluid_disk_bool && central_mass_bool){
                  typ disk_mass = M_PI*(Rroche*Rroche - R_unit*R_unit)*fluid_disk_Sigma;
                  printf("                  Central mass = %.8lf,  Bodies mass = %.8lf,  Fluid disk mass = %.8lf,  Bodies + disk = %.8lf\n",
                  CM.mass, total_mass, disk_mass, disk_mass + total_mass);
            }
            else{
                  if (central_mass_bool){
                        printf("                  Central mass = %.8lf,  Bodies mass = %.8lf\n", CM.mass, total_mass);
                  }
                  else{
                        printf("                  Bodies mass = %.8lf\n", total_mass);
                  }
            }
            printf("                  Largest body radius = %.8lf\n", maxR);
            if (central_tides_bool && central_mass_bool){
                  printf("                  Length of day = %.13lf\n", 2.0*M_PI/SideralOmega);
                  if (J2_bool && Sun_bool){
                        printf("                  Evection resonance at a = %.13lf\n", evection_resonance);
                  }
            }
            if (collision_bool){
                  printf("                  Total collisions = %d\n", collision_count);
                  if (fragmentation_bool && collision_count != 0){
                        typ cll = (typ) collision_count;
                        typ mrg = 100.0*((typ) merger_count)/cll;
                        typ spc = 100.0*((typ) super_catastrophic_count)/cll;
                        typ hfr = 100.0*((typ) half_fragmentation_count)/cll;
                        typ ffr = 100.0*((typ) full_fragmentation_count)/cll;
                        printf("                  Merger = %.2lf %% | Super-catastrophic = %.2lf %% | Partially fragmented = %.2lf %% | Fully fragmented = %.2lf %%\n", 
                        mrg, spc, hfr, ffr);
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
            FlatTree = NULL;
            how_many_cells = 0;
            cell_id = 0;
            tensor_free();
      }
      
      /******** Reinitializing the array did_collide ********/
      if (collision_bool && one_collision_only_bool){
            for (j = 0; j <= largest_id; j++){
                  *(did_collide + j) = 0;
            }
      }
      
      /******** Reinitializing data relative to the mesh algorithm ********/
      if (mesh_bool && (collision_bool || mutual_bool) && !force_naive_bool){
      
            how_many_pairs = 0; //It is unnecessary to reinitialize the array pairs
            /******** Reinitializing the hash table ********/
            int p;
            for (p = 0; p < how_many_modified; p++){
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
                  printf("                  gamma = %.6lf = %.2lf km\n", gam, gam*6371.0);
            }
      
            /******** Reinitializing some global variables ********/
            how_many_modified = 0; //It is unnecessary to reinitialize the array modified_cells
            how_many_big = 0;
            how_many_small = 0;
            total_neighbours = 0;
      }
      
      /******** If there are very few bodies, I switch to the brute-force algorithm ********/
      if (!brute_force_bool && !force_naive_bool && how_many_moonlets < switch_to_brute_force){
            force_naive_bool = 1;
      }
      
      /******** Reordering the array moonlets ********/
      if (how_many_moonlets < 9*largest_id/10 && largest_id > 50){
            tidy_up(moonlets);
            if ((mutual_bool || collision_bool) && (falcON_bool || standard_tree_bool) && !brute_force_bool && !force_naive_bool){
                  IndexPeanoHilbertOrder = largest_id + 1;
                  for (j = 0; j < IndexPeanoHilbertOrder; j++){
                        PeanoHilbertOrder[j] = j;
                  }
            }
      }
      
      /******** Spawning bodies from the inner fluid disk ********/
      if (inner_fluid_disk_bool){
            time_since_last_spawn     += timestep;
            typ quotient               = floor(time_since_last_spawn/time_between_spawn);
            time_since_last_spawn     -= quotient*time_between_spawn;
            int how_many_to_be_spawned = (int) quotient;
            while (how_many_to_be_spawned){ //If the time elapsed since the last spawn if above the characteristic timescale of body spawning, I spawn a body
                  typ mf            = 16.0*fast_pow(M_PI,4)*f_tilde*f_tilde*fast_pow(fluid_disk_Sigma*Rroche*Rroche,3)/(M_unit*M_unit); //Mass of the spawned body
                  typ rad           = pow(3.0*mf/(4.0*M_PI*spawned_density), 1.0/3.0);
                  fluid_disk_Sigma -= mf/(M_PI*(Rroche*Rroche - R_unit*R_unit)); //New surface density of the inner fluid disk
                  int index         = get_free_index(0); //Retrieving an index in the bodies array for the new body
                  typ nu, omega, Omega;
                  nu              = rdm(0.0, 2.0*M_PI);
                  omega           = rdm(0.0, 2.0*M_PI);
                  Omega           = rdm(0.0, 2.0*M_PI);
                  *(moonlets + index) = init(Rroche, 0.0, 0.0, nu, omega, Omega, spawned_density, rad); //Spawning a body at the Roche radius on a circular and equatorial orbit
                  how_many_to_be_spawned --;
                  X  = (moonlets + index) -> x;  //Taking the body's momentum away from Earth for conservation
                  Y  = (moonlets + index) -> y;
                  Z  = (moonlets + index) -> z;
                  vX = (moonlets + index) -> vx;
                  vY = (moonlets + index) -> vy;
                  vZ = (moonlets + index) -> vz;
                  M  = CM.mass + fluid_disk_Sigma*M_PI*(Rroche*Rroche - R_unit*R_unit);
                  CM.x  *= (M + mf)/M;  CM.y *= (M + mf)/M;  CM.z *= (M + mf)/M;  CM.vx *= (M + mf)/M;  CM.vy *= (M + mf)/M;  CM.vz *= (M + mf)/M;
                  CM.x  -= mf*X/M;      CM.y -= mf*Y/M;      CM.z -= mf*Z/M;      CM.vx -= mf*vX/M;     CM.vy -= mf*vY/M;     CM.vz -= mf*vZ/M;
            }
      }

      /******** Updating the position of the star or companion star ********/
      if (Sun_bool){
            * sun_vector      = star_semi_major*cos(time_elapsed*star_mean_motion);
            *(sun_vector + 1) = star_semi_major*sin(time_elapsed*star_mean_motion)*cos(obliquity);
            *(sun_vector + 2) = star_semi_major*sin(time_elapsed*star_mean_motion)*sin(obliquity);
      }
      
      /******** Updating the J2 and the position of the evection resonance ********/
      if (central_tides_bool && J2_bool){
            if (J2_value == 0.0){
                  J2 = 0.5*SideralOmega*SideralOmega*R_unit*R_unit*R_unit/(G*M_unit);
            }
            else{
                  J2 = J2_value*SideralOmega*SideralOmega*Tearth*Tearth/(4.0*M_PI*M_PI);
            }
            if (Sun_bool){
                  evection_resonance = pow(1.5*sqrt(M_unit/star_mass)*J2, 2.0/7.0)*pow(star_semi_major/R_unit, 3.0/7.0)*R_unit;
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
      int error = 0;
      bigint iter = 0;
      int progressed = 0;
      typ progress = 0.0;
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
                  FlatTree = NULL;
                  how_many_cells = 0;
                  cell_id = 0;
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
                  tensor_initialization();
                  if (expansion_order >= 3){
                        multipole_flattree(FlatTree, moonlets);
                  }
            }        
            kick(moonlets, &CM, vector_field);

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
            time_elapsed += timestep;
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
            if (!collision_bool && !force_naive_bool){
                  three_largest_moonlets(moonlets);
                  three_largest_three_first(moonlets);
                  get_neighbours_mesh(moonlets);
            }
            kick(moonlets, &CM, vector_field); //Performing a full kick

            iter ++;
            time_elapsed += timestep;
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
            if (collision_bool){
                  brute_force(moonlets);  //Resolving collisions and going backward
            }           
            drift(moonlets, &CM);         //Performing a full drift

            iter ++;
            time_elapsed += timestep;
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



