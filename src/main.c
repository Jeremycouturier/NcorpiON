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

      
      /******** Closing the output files ********/
      if (write_to_files_bool){
            printf("Closing output files\n");
            file_closing();
      }
      

      /******** Deallocating the arrays that are used globally ********/
      printf("Deallocating global memory\n");
      deallocation();
      printf("Over\n");

      return integ;
}








