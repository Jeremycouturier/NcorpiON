/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    rk4.c                                                       ********/
/******** @brief   Unlike name-suggested, implements a Leap-frog integrator    ********/
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
/******** along with rebound.  If not, see <http://www.gnu.org/licenses/>.     ********/
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
#include <errno.h>
#include <math.h>
#include <string.h>

/******** Declaring external global files ********/
FILE ** files;

struct node * FlatTree;


void kick(struct moonlet * X, void (*F)(struct moonlet *)){

      /******** Performs the kick of the leapfrog method   ********/
      /******** F is the vector field such that dX/dt=F(X) ********/
      
      int i;
      
      
      for(i=0; i<=largest_id; i++){            
            if(*(exists+i)){ //Checking whether or not there is a moonlet in the i^th cell of the moonlet array
            
                  /******** xx=X ********/
                  *(xx+i) = *(X+i);
            }
      }     
      
      /******** xx=F(X)=dX/dt ********/
      (*F)(xx);      
      
      /******** Applying the kick ********/
      for(i=0; i<=largest_id; i++){            
            if(*(exists+i)){ //Checking whether or not there is a moonlet in the ith cell of the moonlet array
                  (X+i) -> vx += timestep * (xx+i) -> vx;
                  (X+i) -> vy += timestep * (xx+i) -> vy;
                  (X+i) -> vz += timestep * (xx+i) -> vz;
            }
      }
}


void drift(struct moonlet * X){

      /******** Performs the drift of the leapfrog method                          ********/
      /******** When this function is called, collision have already been resolved ********/
      
      int i;
      
      
      for(i=0; i<=largest_id; i++){            
            if(*(exists+i)){ //Checking whether or not there is a moonlet in the i^th cell of the moonlet array
                  (X+i) -> x += timestep * (X+i) -> vx;
                  (X+i) -> y += timestep * (X+i) -> vy;
                  (X+i) -> z += timestep * (X+i) -> vz;
            }
      }
}


FILE ** file_opening(){

      /******** Opens the files given as arguments to "display" ********/      
      
      
      /******** Initializing the paths towards the files ********/
      char filex_path[100]; char filey_path[100]; char filez_path[100]; char filevx_path[100]; char filevy_path[100]; char filevz_path[100]; char filerad_path[100];
      char filea_path[100]; char filee_path[100]; char filei_path[100]; char filestat_path[100];
      
      strcpy(filex_path,   path); //The string "path" is defined in structure.h
      strcpy(filey_path,   path);
      strcpy(filez_path,   path);
      strcpy(filevx_path,  path);
      strcpy(filevy_path,  path);
      strcpy(filevz_path,  path);
      strcpy(filerad_path, path);
      strcpy(filea_path,   path);
      strcpy(filee_path,   path);
      strcpy(filei_path,   path);
      strcpy(filestat_path,path);
      
      strcat(filex_path,       "x.txt");
      strcat(filey_path,       "y.txt");
      strcat(filez_path,       "z.txt");
      strcat(filevx_path,     "vx.txt");
      strcat(filevy_path,     "vy.txt");
      strcat(filevz_path,     "vz.txt");
      strcat(filerad_path,"radius.txt");
      strcat(filea_path,       "a.txt");
      strcat(filee_path,       "e.txt");
      strcat(filei_path,       "i.txt");
      strcat(filestat_path, "stat.txt");
      
      
      /******** Creating and opening the files ********/
      FILE * filex;  FILE * filey;  FILE * filez;  FILE * filevx;  FILE * filevy;  FILE * filevz;  FILE * filerad;  FILE * filea;  FILE * filee;  FILE * filei;  FILE * filestat;
      filex   = fopen(filex_path,   "w");
      filey   = fopen(filey_path,   "w");
      filez   = fopen(filez_path,   "w");
      filevx  = fopen(filevx_path,  "w");
      filevy  = fopen(filevy_path,  "w");
      filevz  = fopen(filevz_path,  "w");
      filerad = fopen(filerad_path, "w");
      filea   = fopen(filea_path,   "w");
      filee   = fopen(filee_path,   "w");
      filei   = fopen(filei_path,   "w");
      filestat= fopen(filestat_path,"w");
      if (filex==NULL || filey==NULL || filez==NULL || filevx==NULL || filevy==NULL || filevz==NULL || filerad==NULL || filea==NULL || filee==NULL || filei==NULL || filestat==NULL){
            fprintf(stderr, "Error : Can't create or open output files.\n");
            abort();
      }
      
      /******** Defining and initializing the array of 6 files ********/
      files = (FILE **)malloc(11*sizeof(FILE *));
      if (files==NULL){
            fprintf(stderr, "Error : Can't allocate file array.\n");
            abort();
      }
      *files=filex;
      *(files+1)=filey;
      *(files+2)=filez;
      *(files+3)=filevx;
      *(files+4)=filevy;
      *(files+5)=filevz;
      *(files+6)=filerad;
      *(files+7)=filea;
      *(files+8)=filee;
      *(files+9)=filei;
      *(files+10)=filestat;
      

      return files;

}


void file_closing(){

      /******** Closes the files passed as argument to "display" ********/


      /******** Closing the files ********/
      fclose(*files);
      fclose(*(files+1));
      fclose(*(files+2));
      fclose(*(files+3));
      fclose(*(files+4));
      fclose(*(files+5));
      fclose(*(files+6));
      fclose(*(files+7));
      fclose(*(files+8));
      fclose(*(files+9));
      fclose(*(files+10));
      free(files);
      files=NULL;

}


void display(struct moonlet * moonlets, typ * aei){

      /******** Outputs the moonlets to the output files                                                                      ********/
      /******** The output occurs on 11 files named x.txt, y.txt, ... vz.txt, radius.txt, a.txt, e.txt, i.txt and stat.txt    ********/
      /******** x.txt to vz.txt contain the cartesian coordinates, while a.txt, e.txt and i.txt contains the orbital elements ********/
      /******** In each file, a time step is printed on a single line. For example, the third line of e.txt contains the      ********/
      /******** eccentricities of the N_max moonlets. If the j^th moonlet did not exist at the n^th output, then the j^th     ********/
      /******** column of the n^th line contains 0. The columns of stat.txt are time, total number of moonlets, cumulative    ********/
      /******** number of collisions, radius of the largest moonlet and total mass of the moonlets                            ********/
      
      
      FILE * filex;
      FILE * filey;
      FILE * filez;
      FILE * filevx;
      FILE * filevy;
      FILE * filevz;
      FILE * filerad;
      FILE * filea;
      FILE * filee;
      FILE * filei;
      FILE * filestat;
  
      
      /******** Defining the 7 files ********/
      filex    = * files;
      filey    = *(files+1);
      filez    = *(files+2);
      filevx   = *(files+3);
      filevy   = *(files+4);
      filevz   = *(files+5);
      filerad  = *(files+6);
      filea    = *(files+7);
      filee    = *(files+8);
      filei    = *(files+9);
      filestat = *(files+10);

      
      /******** Writing to the files ********/
      int p;
      typ X,Y,Z,vX,vY,vZ,m,R;
      typ total_mass = 0.0;
      typ maxR = 0.0;
      for (p = 0; p <= largest_id; p++){
            if(*(exists+p)){
                  cart2aei(moonlets, p, aei);
                  X  = (moonlets+p) ->      x;
                  Y  = (moonlets+p) ->      y;
                  Z  = (moonlets+p) ->      z;
                  vX = (moonlets+p) ->     vx;
                  vY = (moonlets+p) ->     vy;
                  vZ = (moonlets+p) ->     vz;
                  m  = (moonlets+p) ->   mass;
                  R  = (moonlets+p) -> radius;
                  fprintf(filex,   "%.13lf ",        X);
                  fprintf(filey,   "%.13lf ",        Y);
                  fprintf(filez,   "%.13lf ",        Z);
                  fprintf(filevx,  "%.13lf ",       vX);
                  fprintf(filevy,  "%.13lf ",       vY);
                  fprintf(filevz,  "%.13lf ",       vZ);
                  fprintf(filerad, "%.13lf ",        R);
                  fprintf(filea,   "%.13lf ",     *aei);
                  fprintf(filee,   "%.13lf ", *(aei+1));
                  fprintf(filei,   "%.13lf ", *(aei+2));
                  total_mass += m;
                  if (R > maxR){
                        maxR = R;
                  }
            }
      }
      fprintf(filestat, "%.13lf %d %d %.13lf %.13lf", time_elapsed, how_many_moonlets, collision_count, maxR, total_mass);
      fprintf(filex,   "\n");
      fprintf(filey,   "\n");
      fprintf(filez,   "\n");
      fprintf(filevx,  "\n");
      fprintf(filevy,  "\n");
      fprintf(filevz,  "\n");
      fprintf(filerad, "\n");
      fprintf(filea,   "\n");
      fprintf(filee,   "\n");
      fprintf(filei,   "\n");
      fprintf(filestat,"\n");  

}


void end_of_timestep(struct moonlet * moonlets, int progressed){

      /******** Takes care of the end of the timestep by reinitializing the hash table and everything that needs to be. ********/
      /******** Updates the mesh-size gam to match the desired expected number of neighbours                            ********/
      /******** Updates the timestep to match the mesh-size                                                             ********/
      /******** Prints useful informations in the terminal for the user to read                                         ********/
      
      
      int index;
      int j;
      typ X,Y,Z,vX,vY,vZ,m,R;
      
      first_passage = 1;

      /******** Removing from the simulation the moonlets that need to be removed ********/
      how_many_moonlets = 0;
      typ total_mass = 0.0;
      typ maxR       = 0.0;
      for (j = 0; j <= largest_id; j++){
            if (*(exists+j)){
                  X = (moonlets+j) -> x;
                  Y = (moonlets+j) -> y;
                  Z = (moonlets+j) -> z;
                  R = (moonlets+j) -> radius;
                  if (X*X+Y*Y+Z*Z < low_dumping_threshold*low_dumping_threshold || X*X+Y*Y+Z*Z > high_dumping_threshold*high_dumping_threshold){
                        lose_moonlet(j);
                  }
                  if (R > maxR){
                        maxR = R;
                  }
                  how_many_moonlets ++;
                  total_mass += (moonlets+j) -> mass;
            }
      }
      
      /******** If the simulation progressed by at least 0.1%, I display useful informations ********/
      if (progressed){
            printf("                  N = %d\n",how_many_moonlets);
            printf("                  largest id = %d\n",largest_id);
            printf("                  Total moonlet mass = %.8lf Mearth\n", total_mass);
            printf("                  Largest moonlet radius = %.8lf Rearth\n", maxR);
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
      }
      
      /******** Resetting global variables relative to the boxdot tree ********/
      if ((mutual_bool || collision_bool) && !brute_force_bool && !force_naive_bool && (falcON_bool || standard_tree_bool)){
            for (j = 0; j < cell_id; j++){
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
      if (collision_bool){
            for (j = 0; j <= largest_id; j++){
                  *(did_collide+j) = 0;
            }
      }
      

      if (mesh_bool && (collision_bool || mutual_bool) && !force_naive_bool){
      
            how_many_pairs = 0; //It is unnecessary to reinitialize the array pairs
            /******** Reinitializing the hash table ********/
            int p;
            for (p = 0; p < how_many_modified; p++){
                  index = *(modified_cells+p);
                  clear_chain(hash+index); //Reinitializing to NULL the index^th cell of the hash table
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
      
      /******** If there are very few moonlets, I switch to the brute-force algorithm ********/
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
}


int integration_brute_force(typ t){


      /******** Performs the numerical integration with a brute force method. t is the final time ********/


      /******** Initializing the array of moonlets ********/
      struct moonlet * moonlets = populate(M_0, radius_stddev);
      struct moonlet * moonlet_buffer = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //Buffer for output steps
      if (moonlet_buffer == NULL){
            fprintf(stderr, "Error : Can't allocate buffer for array of moonlets in function integration_brute_force.\n");
            abort();
      }
      
      
      /******** Numerical integration ********/
      int error = 0;
      int iter = 0;
      int progressed = 0;
      typ progress = 0.0;
      typ previous_progress = 0.0;
      typ * aei = (typ *)malloc(3*sizeof(typ));
      int index, j;
      
      printf("progress = %.1lf %%\n", 0.0);
      if (write_to_files_bool){
            display(moonlets, aei);
      }
      
      timestep /= 2.0;
      kick(moonlets, vector_field); //Performing half a kick
      timestep *= 2.0;
      
      /******** To be removed ********/
      int sample_size = (int) (t_end/timestep);
      typ * time_per_timestep = (typ *)malloc(sample_size * sizeof(typ));
      clock_t t0, t1;
      typ timestep_time, average, stdd;

      while(time_elapsed < t && error == 0){

            progressed = 0;
                   
            
            /******** Writing the results of the numerical integration to the output files ********/
            if (write_to_files_bool && iter % output_step == 0 && iter > 0){
                  
                  /******** I kicked too far at the previous timestep, I have to unkick by half a timestep ********/
                  for (j = 0; j <= largest_id; j++){
                        if (*(exists + j)){
                              *(moonlet_buffer + j) = *(moonlets + j);
                        }
                  }
                  timestep /= -2.0;
                  kick(moonlet_buffer, vector_field); //Performing half a backward kick
                  timestep *= -2.0;
                  display(moonlet_buffer, aei);
            }

            /******** Integrator is SBAB1 ********/
            t0 = clock(); //To be removed
            if (collision_bool){
                  brute_force(moonlets);        //Resolving collisions and going backward
            }
            drift(moonlets);              //Performing a full drift
            kick(moonlets, vector_field); //Performing a full kick
            t1 = clock(); //To be removed
            timestep_time = ((typ) (t1 - t0)) / CLOCKS_PER_SEC; //To be removed
            time_per_timestep[iter] = timestep_time; //To be removed

            iter ++;
            time_elapsed += timestep;
            progress = time_elapsed/t;
            
            /******** Displaying informations in the terminal, and reinitializing what needs to be reinitialized after every timestep ********/
            if (progress-previous_progress > 0.001){
                  previous_progress = progress;
                  printf("progress = %.1lf %%\n", 100.0*progress);
                  progressed = 1;
            }
            
            /******** Taking care of the end of the timestep ********/
            end_of_timestep(moonlets, progressed);
            
            
      }
      printf("progress = %.1lf %%\n", 100.0);
      printf("total number of timestep performed : %d\n",iter);
      
      /******** To be removed ********/
      average = 0.0;
      for (j = 0; j < sample_size; j++){
            average += time_per_timestep[j];
      }
      average /= (typ) sample_size;
      stdd = 0.0;
      for (j = 0; j < sample_size; j++){
            stdd += (time_per_timestep[j] - average) * (time_per_timestep[j] - average);
      }
      stdd /= (typ) sample_size;
      stdd  = sqrt(stdd);
      printf("average            = %.13lf\n", average);
      printf("standard deviation = %.13lf\n", stdd);
      free(time_per_timestep);
      time_per_timestep = NULL;
      
      
      /******** Deallocating the array of moonlets ********/
      free(moonlets);
      moonlets = NULL;
      free(moonlet_buffer);
      moonlet_buffer = NULL;
      free(aei);
      aei = NULL;
       
      return 0;
}


int integration_tree(typ t){


      /******** Performs the numerical integration using a tree algorithm (falcON or standard tree) ********/
      /******** for mutual interactions treatment. t is the final time.                             ********/


      /******** Initializing the array of moonlets ********/
      struct moonlet * moonlets = populate(M_0, radius_stddev);
      struct moonlet * moonlet_buffer = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //Buffer for output steps
      if (moonlet_buffer == NULL){
            fprintf(stderr, "Error : Can't allocate buffer for array of moonlets in function integration_tree.\n");
            abort();
      }
      
      /******** Numerical integration ********/
      int error = 0;
      int iter = 0;
      int progressed = 0;
      typ progress = 0.0;
      typ previous_progress = 0.0;
      typ * aei = (typ *)malloc(3*sizeof(typ));
      int index, j;
      struct boxdot * root = NULL;
      
      printf("progress = %.1lf %%\n", 0.0);
      if (write_to_files_bool){
            display(moonlets, aei);
      }
      
      /******** Performing half a drift ********/
      timestep /= 2.0;
      if (collision_bool){ //Finding and resolving collisions and going backward
            root = root_cell(moonlets);
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
      drift(moonlets); //Drifting
      timestep *= 2.0;

      /******** Resetting data relative to trees ********/
      if (collision_bool){
            for (j = 0; j < cell_id; j++){
                  free((FlatTree + j) -> dots);
                  (FlatTree + j) -> dots = NULL;
            }
            free(FlatTree);
            FlatTree = NULL;
            how_many_cells = 0;
            cell_id = 0;
            /******** Reinitializing the array did_collide ********/
            for (j = 0; j <= largest_id; j++){
                  *(did_collide+j) = 0;
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
                  timestep /= -2.0;
                  drift(moonlet_buffer); //Drifting backwards half a timestep
                  timestep *= -2.0;
                  display(moonlet_buffer, aei);
            }

            /******** Performing a full kick. Integrator is SABA1 ********/
            if (!force_naive_bool && mutual_bool){
                  root = root_cell(moonlets);
                  FlatTree = flattree_init(root);
                  clear_boxdot(&root);
                  com_flattree(FlatTree, moonlets);
                  Mtot = FlatTree -> M0;
                  rmax_flattree(FlatTree, moonlets);
                  rcrit_flattree(FlatTree, moonlets);
                  tensor_initialization();
                  if (expansion_order >= 3){
                        multipole_flattree(FlatTree, moonlets);
                  }
            }
            kick(moonlets, vector_field);
            
            /******** Performing a full drift ********/
            if (collision_bool){
                  if (force_naive_bool){ //Finding and resolving collisions and going backward brute-forcely
                        brute_force(moonlets);
                  }
                  else{
                        if (!mutual_bool){
                              root = root_cell(moonlets);
                              FlatTree = flattree_init(root);
                              clear_boxdot(&root);
                        }
                        center_and_maxR_flattree(FlatTree, moonlets);
                        rmax_and_rcrit_flattree (FlatTree, moonlets);
                        if (falcON_bool){ //Finding and resolving collisions and going backward with falcON
                              collision_flattree(FlatTree, moonlets);
                        }
                        else if (standard_tree_bool){ //Finding and resolving collisions and going backward with the standard tree code
                              for (j = 0; j <= largest_id; j++){
                                    if (exists[j] && !(did_collide[j])){
                                          standard_tree_collision(FlatTree, moonlets, j);
                                    }
                              }
                        }
                  }
            }
            drift(moonlets);

            iter ++;
            time_elapsed += timestep;
            progress = time_elapsed/t;
            
            /******** Displaying informations in the terminal, and reinitializing what needs to be reinitialized after every timestep ********/
            if (progress-previous_progress > 0.001){
                  previous_progress = progress;
                  printf("progress = %.1lf %%\n", 100.0*progress);
                  progressed = 1;
            }
            
            /******** Taking care of the end of the timestep ********/
            end_of_timestep(moonlets, progressed);
            
            
      }
      printf("progress = %.1lf %%\n", 100.0);
      printf("total number of timestep performed : %d\n",iter);
      
      
      /******** Deallocating the array of moonlets ********/
      free(moonlets);
      moonlets = NULL;
      free(moonlet_buffer);
      moonlet_buffer = NULL;
      free(aei);
      aei = NULL;
     
      return 0;
}


int integration_mesh(typ t){


      /******** Performs the numerical integration using a mesh algorithm ********/
      /******** for mutual interactions treatment. t is the final time.   ********/


      /******** Initializing the array of moonlets ********/
      struct moonlet * moonlets = populate(M_0, radius_stddev);
      struct moonlet * moonlet_buffer = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //Buffer for output steps
      if (moonlet_buffer == NULL){
            fprintf(stderr, "Error : Can't allocate buffer for array of moonlets in function integration_brute_force.\n");
            abort();
      }
      
      
      /******** Numerical integration ********/
      int error = 0;
      int iter = 0;
      int progressed = 0;
      typ progress = 0.0;
      typ previous_progress = 0.0;
      typ * aei = (typ *)malloc(3*sizeof(typ));
      int index, j;
      
      printf("progress = %.1lf %%\n", 0.0);
      if (write_to_files_bool){
            display(moonlets, aei);
      }
      
      timestep /= 2.0;
      three_largest_moonlets(moonlets);
      three_largest_three_first(moonlets);
      get_neighbours_mesh(moonlets);
      kick(moonlets, vector_field); //Performing half a kick
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
      collision_cube = gam*((typ) collision_cube_cells);
      how_many_pairs    = 0; //It is unnecessary to reinitialize the array pairs
      how_many_modified = 0; //It is unnecessary to reinitialize the array modified_cells
      total_neighbours  = 0;
      
      /******** To be removed ********/
      int sample_size = (int) (t_end/timestep);
      typ * time_per_timestep = (typ *)malloc(sample_size * sizeof(typ));
      clock_t t0, t1;
      typ timestep_time, average, stdd;

      while(time_elapsed < t && error == 0){

            progressed = 0;
                   
            
            /******** Writing the results of the numerical integration to the output files ********/
            if (write_to_files_bool && iter % output_step == 0 && iter > 0){

                  /******** I kicked too far at the previous timestep, I have to unkick by half a timestep ********/
                  for (j = 0; j <= largest_id; j++){
                        if (*(exists + j)){
                              *(moonlet_buffer + j) = *(moonlets + j);
                        }
                  }
                  timestep /= -2.0;
                  if (!force_naive_bool){
                        three_largest_moonlets(moonlet_buffer);
                        three_largest_three_first(moonlet_buffer);
                        get_neighbours_mesh(moonlet_buffer);
                  }
                  kick(moonlet_buffer, vector_field); //Performing half a backward kick
                  timestep *= -2.0;
                  display(moonlet_buffer, aei);
                  
                  /******** Reinitializing data relative to the hash table ********/
                  if (!force_naive_bool){
                        for (j = 0; j < how_many_modified; j++){
                              index = *(modified_cells + j);
                              clear_chain(hash + index); //Reinitializing to NULL the index^th cell of the hash table
                        }
                        how_many_pairs    = 0; //It is unnecessary to reinitialize the array pairs
                        how_many_modified = 0; //It is unnecessary to reinitialize the array modified_cells
                        total_neighbours  = 0;
                  }
            }

            /******** Integrator is SBAB1 ********/
            t0 = clock(); //To be removed
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
            drift(moonlets); //Performing a full drift
            if (!collision_bool && !force_naive_bool){
                  three_largest_moonlets(moonlets);
                  three_largest_three_first(moonlets);
                  get_neighbours_mesh(moonlets);
            }
            kick(moonlets, vector_field); //Performing a full kick
            t1 = clock(); //To be removed
            timestep_time = ((typ) (t1 - t0)) / CLOCKS_PER_SEC; //To be removed
            time_per_timestep[iter] = timestep_time; //To be removed

            iter ++;
            time_elapsed += timestep;
            progress = time_elapsed/t;
            
            /******** Displaying informations in the terminal, and reinitializing what needs to be reinitialized after every timestep ********/
            if (progress-previous_progress > 0.001){
                  previous_progress = progress;
                  printf("progress = %.1lf %%\n", 100.0*progress);
                  progressed = 1;
            }
            
            /******** Taking care of the end of the timestep ********/
            end_of_timestep(moonlets, progressed);
            
            
      }
      printf("progress = %.1lf %%\n", 100.0);
      printf("total number of timestep performed : %d\n",iter);
      
      /******** To be removed ********/
      average = 0.0;
      for (j = 0; j < sample_size; j++){
            average += time_per_timestep[j];
      }
      average /= (typ) sample_size;
      stdd = 0.0;
      for (j = 0; j < sample_size; j++){
            stdd += (time_per_timestep[j] - average) * (time_per_timestep[j] - average);
      }
      stdd /= (typ) sample_size;
      stdd  = sqrt(stdd);
      printf("average            = %.13lf\n", average);
      printf("standard deviation = %.13lf\n", stdd);
      free(time_per_timestep);
      time_per_timestep = NULL;
      
      
      /******** Deallocating the array of moonlets ********/
      free(moonlets);
      moonlets = NULL;
      free(moonlet_buffer);
      moonlet_buffer = NULL;
      free(aei);
      aei = NULL;

      return 0;
}


int integration_brute_force_SABA1(typ t){


      /******** Performs the numerical integration with a brute force method. t is the final time ********/


      /******** Initializing the array of moonlets ********/
      struct moonlet * moonlets = populate(M_0, radius_stddev);
      struct moonlet * moonlet_buffer = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //Buffer for output steps
      if (moonlet_buffer == NULL){
            fprintf(stderr, "Error : Can't allocate buffer for array of moonlets in function integration_brute_force.\n");
            abort();
      }
      
      
      /******** Numerical integration ********/
      int error = 0;
      int iter = 0;
      int progressed = 0;
      typ progress = 0.0;
      typ previous_progress = 0.0;
      typ * aei = (typ *)malloc(3*sizeof(typ));
      int index, j;
      
      printf("progress = %.1lf %%\n", 0.0);
      if (write_to_files_bool){
            display(moonlets, aei);
      }
      
      timestep /= 2.0;
      drift(moonlets); //Performing half a drift
      timestep *= 2.0;

      /******** To be removed ********/
      int sample_size = (int) (t_end/timestep);
      typ * time_per_timestep = (typ *)malloc(sample_size * sizeof(typ));
      clock_t t0, t1;
      typ timestep_time, average, stdd;

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
                  timestep /= -2.0;
                  drift(moonlet_buffer); //Performing half a backward drift
                  timestep *= -2.0;
                  display(moonlet_buffer, aei);
            }

            /******** Integrator is SBAB1 ********/
            t0 = clock(); //To be removed
            kick(moonlets, vector_field); //Performing a full kick
            if (collision_bool){
                  brute_force(moonlets);  //Resolving collisions and going backward
            }
            drift(moonlets);              //Performing a full drift
            t1 = clock(); //To be removed
            timestep_time = ((typ) (t1 - t0)) / CLOCKS_PER_SEC; //To be removed
            time_per_timestep[iter] = timestep_time; //To be removed
            

            iter ++;
            time_elapsed += timestep;
            progress = time_elapsed/t;
            
            /******** Displaying informations in the terminal, and reinitializing what needs to be reinitialized after every timestep ********/
            if (progress-previous_progress > 0.001){
                  previous_progress = progress;
                  printf("progress = %.1lf %%\n", 100.0*progress);
                  progressed = 1;
            }
            
            /******** Taking care of the end of the timestep ********/
            end_of_timestep(moonlets, progressed);
            
            
      }
      printf("progress = %.1lf %%\n", 100.0);
      printf("total number of timestep performed : %d\n",iter);
      
      /******** To be removed ********/
      average = 0.0;
      for (j = 0; j < sample_size; j++){
            average += time_per_timestep[j];
      }
      average /= (typ) sample_size;
      stdd = 0.0;
      for (j = 0; j < sample_size; j++){
            stdd += (time_per_timestep[j] - average) * (time_per_timestep[j] - average);
      }
      stdd /= (typ) sample_size;
      stdd  = sqrt(stdd);
      printf("average            = %.13lf\n", average);
      printf("standard deviation = %.13lf\n", stdd);
      free(time_per_timestep);
      time_per_timestep = NULL;
      
      
      /******** Deallocating the array of moonlets ********/
      free(moonlets);
      moonlets = NULL;
      free(moonlet_buffer);
      moonlet_buffer = NULL;
      free(aei);
      aei = NULL;
       
      return 0;
}



