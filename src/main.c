/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    main.c                                                      ********/
/******** @brief   The main file of NcorpiON                                   ********/
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
#include "collision.h"
#include "ffm.h"
#include "display.h"
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#if openGL_bool
      #include <mpi.h>
#endif


int main(int argc, char ** argv){

      (void) argc;
      (void) argv;

      /*********************************************************************/
      /******** Verifying that the parameter file is properly setup ********/
      /*********************************************************************/
      printf("\nChecking if the parameter file is properly setup\n");
      verify();


      /***************************************************************/
      /******** Defining the seed for the random function rdm ********/
      /***************************************************************/
      printf("Defining the seed for random number generation\n");
      time_t t;
      if (seed_bool){
            t = seed;
            srand((unsigned) t);
      }
      else{
            time(&t);
            srand((unsigned) t);
      }
      printf("Seed = %d\n", (int) t);
      
      
      /************************************************/
      /******** Initializing the global memory ********/
      /************************************************/
      printf("Allocating and initializing global memory\n");
      variable_initialization();
      array_initialization();
      k_from_s2s3_init();
      s1s2s3_from_kn_init();
      s2s3_from_k_init();
      k_from_ijklmn_init();
      ijklmn_from_kn_init();
      perm_from_kn_init();
      q1fromq2q3_init();
      
      
      /***************************************************/
      /******** Creating and opening output files ********/
      /***************************************************/
      if (write_to_files_bool){
            printf("Opening output files\n");
            file_opening();
      }
      
      
      /************************************************************/
      /******** Initializing the MPI execution environment ********/
      /************************************************************/
      #if openGL_bool
            printf("Initializing MPI environment\n");
            MPI_Init(&argc, &argv);
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            typ * dtt = (typ *)malloc(3*sizeof(typ));
            int * att = (int *)malloc(2*sizeof(int));
            dtt[0] = timestep;  dtt[1] = t_end;  dtt[2] = radius_blow_up_factor;
            att[0] = central_mass_bool;  att[1] = browser_port;
            MPI_Send(dtt, 3, MPI_TYP, !rank, 0, MPI_COMM_WORLD);
            MPI_Send(att, 2, MPI_INT, !rank, 0, MPI_COMM_WORLD);
            free(dtt);  dtt = NULL;
            free(att);  att = NULL;
      #endif
      
      /*****************************************************/
      /******** Launching the numerical integration ********/
      /*****************************************************/
      printf("Launching numerical integration\n");
      typ duration;
      clock_t t0, t1;
      t0 = clock();
      if (brute_force_bool){
            integration_brute_force_SABA1(t_end);
      }
      else if (falcON_bool || standard_tree_bool){
            integration_tree(t_end);
      }
      else if (mesh_bool){
            integration_mesh(t_end);
      }
      t1       = clock();
      duration = ((typ) (t1 - t0))/(typ) CLOCKS_PER_SEC;
      printf("The integration took %.3lf seconds to complete\n", duration);
      
      
      /***********************************************************/
      /******** Terminating the MPI execution environment ********/
      /***********************************************************/  
      #if openGL_bool
            printf("Terminating MPI environment\n");
            MPI_Finalize();
      #endif
      
      
      /******************************************/
      /******** Closing the output files ********/
      /******************************************/
      if (write_to_files_bool){
            printf("Closing output files\n");
            file_closing();
      }
      

      /****************************************************************/
      /******** Deallocating the arrays that are used globally ********/
      /****************************************************************/
      printf("Deallocating global memory\n");
      deallocation();
      printf("Simulation over\n");
      
      
      /*********************************************************/
      /******** Making animations out of the simulation ********/
      /*********************************************************/
      if (make_animation_bool && write_to_files_bool){
            /******** Producing the images      ********/
            /******** A python script is called ********/
            printf("Producing the images for the animations\n");
            char frag_bl[10]; 
            char inner_bl[10];
            char J2_bl[10];
            char Sun_bl[10];
            char tides_bl[10];
            if (J2_bool){ //The second zonal harmonic is passed to the python script
                  sprintf(J2_bl, "%d", 1);
            }
            else{
                  sprintf(J2_bl, "%d", 0);
            }
            if (Sun_bool){ //The user's decision on taking into account perturbations from the sun is passed to the python script
                  sprintf(Sun_bl, "%d", 1);
            }
            else{
                  sprintf(Sun_bl, "%d", 0);
            }
            if (fragmentation_bool && collision_bool){ //The user's decision to use NcorpiON's fragmentation model is passed to the python script
                  sprintf(frag_bl, "%d", 1);
            }
            else{
                  sprintf(frag_bl, "%d", 0);
            }
            if (inner_fluid_disk_bool){ //The user's decision to feature an inner fluid disk is passed to the python script
                  sprintf(inner_bl, "%d", 1);
            }
            else{
                  sprintf(inner_bl, "%d", 0);
            }
            if (central_tides_bool){ //The user's decision to consider tides on the central body is passed to the python script
                  sprintf(tides_bl, "%d", 1);
            }
            else{
                  sprintf(tides_bl, "%d", 0);
            }
            
            int n_images = (int) (t_end/(timestep*(typ) output_step)); //The number of images to produce
            char argument[20];
            sprintf(argument, "%d", n_images); //Transforming it to a string
            char to_system[500];
            strcpy(to_system, "./image_creation.sh ");
            strcat(to_system, argument);
            strcat(to_system, " ");
            strcat(to_system, pth);
            strcat(to_system, " ");
            strcat(to_system, inner_bl);
            strcat(to_system, " ");
            strcat(to_system, frag_bl);
            strcat(to_system, " ");
            strcat(to_system, J2_bl);
            strcat(to_system, " ");
            strcat(to_system, Sun_bl);
            strcat(to_system, " ");
            strcat(to_system, tides_bl);
            (void) (system(to_system) + 1); //Just a weid manoeuver to avoid the annoying "set but not used" compiler warning
            
            /******** Assembling the images into a gif and a mp4 ********/
            /******** Ffmpeg is called                           ********/
            printf("\nAssembling the images\n\n");
            strcpy(to_system, "./ncorpion_animation.sh ");
            strcat(to_system, pth);
            (void) (system(to_system) + 1); //Just a weid manoeuver to avoid the annoying "set but not used" compiler warning
      }

      return 0;
}



