/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    3D.c                                                        ********/
/******** @brief   3D visualization of NcorpiON simulations using REBOUND      ********/
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


/***********************************************************************************************************/
/***********************************************************************************************************/
/******** In this file, I use REBOUND to visualize in 3D with webGL a simulation run by NcorpiON.   ********/
/******** The timesteps of the simulation are sent by NcorpiON to REBOUND using MPI in the          ********/
/******** heartbeat function. I use the integrator REB_INTEGRATOR_NONE of REBOUND that does nothing.********/
/******** Since REB_INTEGRATOR_NONE does nothing, REBOUND does not modify the state of the system   ********/
/******** during the integration but receives the new timesteps from NcorpiON. Both NcorpiON and    ********/
/******** REBOUND are run at the same time. After each timestep, NcorpiON sends to REBOUND the new  ********/
/******** coordinates, masses and radii of the bodies and REBOUND updates the state of the system.  ********/
/******** From REBOUND's point of view, it is as if it was actually integrating the system.         ********/
/******** This unconventional approach allows me to enjoy the nice visualization tool of REBOUND in ********/
/******** NcorpiON without having to go through the hassle of compute shaders.                      ********/
/***********************************************************************************************************/
/***********************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include "rebound.h"

#define typ double                   //Renaming double as typ
#define MPI_TYP MPI_DOUBLE           //Renaming double as typ for mpi

int rank  = 0;
int central_mass_bool;
typ radius_blow_up_factor;

void heartbeat(struct reb_simulation* r);


int main(int argc, char* argv[]){

      /******** Initializing ********/
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Status status;

      /******** Preparing the simulation ********/
      struct reb_simulation * r = reb_simulation_create();
      r -> integrator = REB_INTEGRATOR_NONE; //No integrator. REBOUND is only used for visualization
      r -> heartbeat = heartbeat;
      r -> exact_finish_time = 0;
      
      /******** Receiving some informations from NcorpiON ********/
      typ * dtt = (typ *)malloc(4*sizeof(typ));
      int * att = (int *)malloc(2*sizeof(int));
      MPI_Recv(dtt, 4, MPI_TYP, !rank, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(att, 2, MPI_INT, !rank, 0, MPI_COMM_WORLD, &status);
      r -> dt               = dtt[0];
      typ until             = dtt[1];
      radius_blow_up_factor = dtt[2];
      r -> t                = dtt[3];
      central_mass_bool     = att[0];
      int browser_port      = att[1];
      free(dtt);  dtt = NULL;
      free(att);  att = NULL;

      /******** Starting the fake integration, so the visualization can start ********/
      reb_simulation_start_server(r, browser_port);
      printf("The 3D visualization will start in 5 seconds.\n");
      usleep(5000000.);                     //Waiting 5 seconds to give the user enough time to open a web browser tab
      reb_simulation_integrate(r, until);   //Pseudo - integrating. Only what is done in the heartbeat function serves a purpose

      /******** Finishing ********/
      MPI_Finalize();
      reb_simulation_stop_server(r);
      reb_simulation_free(r);
      return 0;
}


void heartbeat(struct reb_simulation * r){

      /******** Receives the next timestep from NcorpiON and updates the particle array ********/
      
      int N, j;

      /******** Communicating with NcorpiON ********/
      MPI_Status status;
      MPI_Recv(&N, 1, MPI_INT, !rank, 0, MPI_COMM_WORLD, &status);       //Receiving the number of bodies
      typ * buffer = (typ *)malloc(8*N*sizeof(typ));                     //Allocating receiving buffer
      MPI_Recv(buffer, 8*N, MPI_TYP, !rank, 0, MPI_COMM_WORLD, &status); //Receiving the bodies' coordinates

      /******** Updating the simulation ********/
      reb_simulation_remove_all_particles(r);                            //I first remove all the bodies from the simulation
      for (j = 0; j < N; j ++){                                          //I now recreate the bodies with their new coordinates and add them back to the simulation
            struct reb_particle body;
            body.x  = buffer[8*j];
            body.y  = buffer[8*j + 1];
            body.z  = buffer[8*j + 2];
            body.vx = buffer[8*j + 3];
            body.vy = buffer[8*j + 4];
            body.vz = buffer[8*j + 5];
            body.m  = buffer[8*j + 6];
            body.r  = buffer[8*j + 7];
            body.r *= radius_blow_up_factor; //Making bodies other than the central mass bigger than their actual size for better visualization
            body.r /= (!j && central_mass_bool ? radius_blow_up_factor : 1.);
            reb_simulation_add(r, body);
      }
      free(buffer);  buffer = NULL;
      if (N < 4000){
            usleep(1000000.0/80.0); //Sleeping a few milliseconds to slow down the simulation
      }
}
