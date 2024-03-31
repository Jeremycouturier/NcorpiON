/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    main.c                                                      ********/
/******** @brief   The main file of the software NcorpiON                      ********/
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
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>


int main(){


      /******** Defining the seed for the random function rdm ********/
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
      
      /******** Initializing the variables that are used globally ********/
      printf("Allocating and initializing global memory\n");
      variable_initialization();
      
      
      /******** Initializing the arrays that are used globally ********/
      array_initialization();
      k_from_s2s3_init();
      s1s2s3_from_kn_init();
      s2s3_from_k_init();
      k_from_ijklmn_init();
      ijklmn_from_kn_init();
      perm_from_kn_init();
      q1fromq2q3_init();
      
      
      /******** Creating and opening the output files ********/
      if (write_to_files_bool){
            printf("Opening output files\n");
            file_opening();
      }
      
      
      /******** Launching the numerical integration ********/
      printf("Launching numerical integration\n");
      int integ;
      if (brute_force_bool){
            integ = integration_brute_force_SABA1(t_end);
      }
      else if (falcON_bool || standard_tree_bool){
            integ = integration_tree(t_end);
      }
      else if (mesh_bool){
            integ = integration_mesh(t_end);
      }
      else{
            fprintf(stderr, "Error : Exactly one of the booleans (brute_force_bool, falcON_bool, standard_tree_bool, mesh_bool) must be 1.\n");
            abort();
      }

      
      /******** Closing the output files ********/
      if (write_to_files_bool){
            printf("Closing output files\n");
            file_closing();
      }
      

      /******** Deallocating the arrays that are used globally ********/
      printf("Deallocating global memory\n");
      deallocation();
      printf("Simulation over\n");
      
      
      /******** Making animations out of the simulation ********/
      if (make_animation_bool && write_to_files_bool){
            /******** Producing the images      ********/
            /******** A python script is called ********/
            printf("Producing the images for the animations\n");
            char frag_bl[10]; 
            char inner_bl[10];
            char sideralPeriod[30];
            char evectionResonance[30];
            if (J2_bool){ //The second zonal harmonic is passed to the python script
                  sprintf(sideralPeriod, "%.13lf", Tearth);
                  if (Sun_bool){ //The position of the evection resonance is passed to the python script 
                        sprintf(evectionResonance, "%.13lf", evection_resonance);
                  }
                  else{
                        sprintf(evectionResonance, "%.13lf", 0.0);
                  }
            }
            else{
                  sprintf(sideralPeriod, "%.13lf", 999999999.0);
            }
            if (fragmentation_bool && collision_bool){ //The user decision to use NcorpiON's fragmentation model is passed to the python script
                  sprintf(frag_bl, "%d", 1);
            }
            else{
                  sprintf(frag_bl, "%d", 0);
            }
            if (inner_fluid_disk_bool){ //The user decision to feature an inner fluid disk is passed to the python script
                  sprintf(inner_bl, "%d", 1);
            }
            else{
                  sprintf(inner_bl, "%d", 0);
            }
            
            int n_images = (int) (t_end/(timestep*(typ) output_step)); //The number of images to produce
            char argument[20];
            sprintf(argument, "%d", n_images); //Transforming it to a string
            char to_system[500];
            strcpy(to_system, "./image_creation.sh ");
            strcat(to_system, argument);
            strcat(to_system, " ");
            strcat(to_system, path);
            strcat(to_system, " ");
            strcat(to_system, inner_bl);
            strcat(to_system, " ");
            strcat(to_system, frag_bl);
            strcat(to_system, " ");
            strcat(to_system, sideralPeriod);
            strcat(to_system, " ");
            strcat(to_system, evectionResonance);      
            int status = system(to_system);
            
            /******** Assembling the images into a gif and a mp4 ********/
            /******** Ffmpeg is called                           ********/
            printf("\nAssembling the images\n\n");
            strcpy(to_system, "./ncorpion_animation.sh ");
            strcat(to_system, path);
            status = system(to_system);
      }

      return integ;
}








