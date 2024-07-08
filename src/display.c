/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    display.c                                                   ********/
/******** @brief   Exports the simulation (to files or to REBOUND's 3D tool)   ********/
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
#include "parameters.h"
#include "structure.h"
#include "collision.h"
#include "physics.h"
#include "ffm.h"
#include "rk4.h"
#include "display.h"
#include "spring.h"
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#if openGL_bool
      #include <mpi.h>
#endif

struct collisionData * collisionDatas;
int indexCollision;


FILE ** file_opening(){

      /******** Opens the output files used by the function display ********/      
      
      
      /******** Initializing the paths towards the files ********/
      char filex_path[800]; char filey_path[800]; char filez_path   [800]; char filevx_path [800]; char filevy_path[800]; char filevz_path[800];
      char filea_path[800]; char filel_path[800]; char filek_path   [800]; char fileh_path  [800]; char fileq_path [800]; char filep_path [800];
      char filer_path[800]; char filem_path[800]; char filestat_path[800]; char filecol_path[800];
      
      strcpy(filex_path,   pth); //The string "pth" is defined in structure.h
      strcpy(filey_path,   pth);
      strcpy(filez_path,   pth);
      strcpy(filevx_path,  pth);
      strcpy(filevy_path,  pth);
      strcpy(filevz_path,  pth);
      strcpy(filer_path,   pth);
      strcpy(filem_path,   pth);
      strcpy(filea_path,   pth);
      strcpy(filel_path,   pth);
      strcpy(filek_path,   pth);
      strcpy(fileh_path,   pth);
      strcpy(fileq_path,   pth);
      strcpy(filep_path,   pth);
      strcpy(filestat_path,pth);
      strcpy(filecol_path, pth);
      
      strcat(filex_path,           "x.txt");
      strcat(filey_path,           "y.txt");
      strcat(filez_path,           "z.txt");
      strcat(filevx_path,         "vx.txt");
      strcat(filevy_path,         "vy.txt");
      strcat(filevz_path,         "vz.txt");
      strcat(filer_path,      "radius.txt");
      strcat(filem_path,        "mass.txt");
      strcat(filea_path,           "a.txt");
      strcat(filel_path,      "lambda.txt");
      strcat(filek_path,           "k.txt");
      strcat(fileh_path,           "h.txt");
      strcat(fileq_path,           "q.txt");
      strcat(filep_path,           "p.txt");
      strcat(filestat_path,     "stat.txt");
      strcat(filecol_path, "collision.txt");
      
      /******** Creating and opening the files ********/
      FILE * filex;  FILE * filey;  FILE * filez;  FILE * filevx;  FILE * filevy;  FILE * filevz;  FILE * filer;     FILE * filem;
      FILE * filea;  FILE * filel;  FILE * filek;  FILE * fileh;   FILE * fileq;   FILE * filep;   FILE * filestat;  FILE * filecol;
      
      filestat = fopen(filestat_path, "w");
      if (write_cartesian_bool || write_elliptic_bool){
            filer = fopen(filer_path, "w");
            filem = fopen(filem_path, "w");
      }
      if (collision_bool && write_collisions_bool){
            filecol = fopen(filecol_path, "w");
      }
      if (write_cartesian_bool){
            filex  = fopen(filex_path,  "w");
            filey  = fopen(filey_path,  "w");
            filez  = fopen(filez_path,  "w");
            filevx = fopen(filevx_path, "w");
            filevy = fopen(filevy_path, "w");
            filevz = fopen(filevz_path, "w");
      }
      if (write_elliptic_bool){
            filea = fopen(filea_path, "w");
            filel = fopen(filel_path, "w");
            filek = fopen(filek_path, "w");
            fileh = fopen(fileh_path, "w");
            fileq = fopen(fileq_path, "w");
            filep = fopen(filep_path, "w");
      }
      
      if (filestat == NULL){
            fprintf(stderr, "Error : Cannot create or open output files. Did you specify the path 'pth' in the parameter file ?\n");
            abort();
      }
      if (write_cartesian_bool && (filex == NULL || filey == NULL || filez == NULL || filevx == NULL || filevy == NULL || filevz == NULL || filer == NULL || filem == NULL)){
            fprintf(stderr, "Error : Cannot create or open output files. Did you specify the path 'pth' in the parameter file ?\n");
            abort();
      }
      if (write_elliptic_bool && (filea == NULL || filel == NULL || filek == NULL || fileh == NULL || fileq == NULL || filep == NULL || filer == NULL || filem == NULL)){
            fprintf(stderr, "Error : Cannot create or open output files. Did you specify the path 'pth' in the parameter file ?\n");
            abort();
      }
      if (collision_bool && write_collisions_bool && filecol == NULL){
            fprintf(stderr, "Error : Cannot create or open output files. Did you specify the path 'pth' in the parameter file ?\n");
            abort();
      }
      
      /******** Defining and initializing the array of files ********/
      int collision = collision_bool && write_collisions_bool ? 1 : 0;
      int coordinates = write_cartesian_bool || write_elliptic_bool ? 1 : 0;
      int nfiles = 1 + collision + 6*(write_cartesian_bool + write_elliptic_bool) + 2*coordinates;
      files = (FILE **)malloc(nfiles*sizeof(FILE *));
      if (files == NULL){
            fprintf(stderr, "Error : Cannot allocate file array in function file_opening.\n");
            abort();
      }
      * files = filestat;
      if (collision){
            *(files + 1) = filecol;
      }
      if (write_cartesian_bool){
            *(files + 1 + collision) = filem;
            *(files + 2 + collision) = filer;
            *(files + 3 + collision) = filex;
            *(files + 4 + collision) = filey;
            *(files + 5 + collision) = filez;
            *(files + 6 + collision) = filevx;
            *(files + 7 + collision) = filevy;
            *(files + 8 + collision) = filevz;
            if (write_elliptic_bool){
                  *(files + 9  + collision) = filea;
                  *(files + 10 + collision) = filel;
                  *(files + 11 + collision) = filek;
                  *(files + 12 + collision) = fileh;
                  *(files + 13 + collision) = fileq;
                  *(files + 14 + collision) = filep;
            }
      }
      else if (write_elliptic_bool){
            *(files + 1 + collision) = filem;
            *(files + 2 + collision) = filer;
            *(files + 3 + collision) = filea;
            *(files + 4 + collision) = filel;
            *(files + 5 + collision) = filek;
            *(files + 6 + collision) = fileh;
            *(files + 7 + collision) = fileq;
            *(files + 8 + collision) = filep;
      }
      return files;
}


void file_closing(){

      /******** Closes the output files ********/

      int collision = collision_bool && write_collisions_bool ? 1 : 0;
      /******** Closing the files ********/
      fclose(* files);
      if (collision){
            fclose(*(files + 1));
      }
      if (write_cartesian_bool || write_elliptic_bool){
            fclose(*(files + 1 + collision));
            fclose(*(files + 2 + collision));
            fclose(*(files + 3 + collision));
            fclose(*(files + 4 + collision));
            fclose(*(files + 5 + collision));
            fclose(*(files + 6 + collision));
            fclose(*(files + 7 + collision));
            fclose(*(files + 8 + collision));
            if (write_cartesian_bool && write_elliptic_bool){
                  fclose(*(files + 9  + collision));
                  fclose(*(files + 10 + collision));
                  fclose(*(files + 11 + collision));
                  fclose(*(files + 12 + collision));
                  fclose(*(files + 13 + collision));
                  fclose(*(files + 14 + collision));
            }
      }
      free(files);
      files = NULL;
}


void display(struct moonlet * moonlets){

      /******** Outputs the bodies to the output files. If write_cartesian_bool is 1, then 6 files named x.txt, y.txt, etc ... ********/
      /******** up to vz.txt are created and contain the cartesian coordinates of the bodies (in simulation's unit).           ********/
      /******** If write_elliptic_bool is 1, then 6 files named a.txt, lambda.txt, k.txt, h.txt, q.txt and p.txt are created   ********/
      /******** and contain the elliptic elements of the bodies (in simulation's units and radians).                           ********/
      /******** In each file, a timestep is printed on a single line. For example, a line of a.txt contains the semi-major     ********/
      /******** axes of the N bodies existing at that time. The columns of stat.txt are time, total number of bodies, total    ********/
      /******** number of collisions, radius of the largest body, total mass of the bodies, inner fluid disk mass, number      ********/
      /******** of mergers, of super-catastrophic collisions, of half-fragmentations, of full-fragmentations, length of day,   ********/
      /******** value of the J2 and semi-major axis of the evection resonance (assuming everything is planar). If the J2 or    ********/
      /******** perturbations from the star are not taken into account, then 0.0 is displayed for the evection resonance.      ********/
      /******** Files other that stat.txt have a variable number of columns due to the variable number of bodies throughout    ********/
      /******** the simulation. Two files radius.txt and mass.txt that contain the radii and masses are output as well.        ********/
      /******** If central_mass_bool is 1, then the central mass is output first in each file.                                 ********/
      
      
      FILE * filex, * filey, * filez, * filevx, * filevy, * filevz, * filer, * filem, * filea, * filel, * filek, * fileh, * fileq, * filep, * filestat, * filecol;
      typ * alkhqp;
      alkhqp = (typ *)malloc(6*sizeof(typ));
      int collision = collision_bool && write_collisions_bool ? 1 : 0;
      
      /******** Defining the files ********/
      filestat = * files;
      if (collision){
            filecol = *(files + 1);
      }
      if (write_cartesian_bool){
            filem  = *(files + 1 + collision);
            filer  = *(files + 2 + collision);
            filex  = *(files + 3 + collision);
            filey  = *(files + 4 + collision);
            filez  = *(files + 5 + collision);
            filevx = *(files + 6 + collision);
            filevy = *(files + 7 + collision);
            filevz = *(files + 8 + collision);
            if (write_elliptic_bool){
                  filea = *(files + 9  + collision);
                  filel = *(files + 10 + collision);
                  filek = *(files + 11 + collision);
                  fileh = *(files + 12 + collision);
                  fileq = *(files + 13 + collision);
                  filep = *(files + 14 + collision);
            }
      }
      else if (write_elliptic_bool){
            filem = *(files + 1 + collision);
            filer = *(files + 2 + collision);
            filea = *(files + 3 + collision);
            filel = *(files + 4 + collision);
            filek = *(files + 5 + collision);
            fileh = *(files + 6 + collision);
            fileq = *(files + 7 + collision);
            filep = *(files + 8 + collision);
      }
      
      /******** Writing to the files ********/
      int p;
      typ X, Y, Z, vX, vY, vZ, m, R;
      typ total_mass = 0.0;
      typ inner_fluid_disk_mass;
      typ M1, M2, R1, R2;                                 //Masses and radii of the two most massive bodies
      typ a1, a2, l1, l2, K1, K2, H1, H2, Q1, Q2, P1, P2; //Orbital elements of the two most massive bodies
      int i1, i2;                                         //Indexes          of the two most massive bodies
      M1 = 0.; M2 = 0.; i1 = 0; i2 = 0;
      
      if (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
            if (write_cartesian_bool){
                  X  = CM_buffer.x;
                  Y  = CM_buffer.y;
                  Z  = CM_buffer.z;
                  vX = CM_buffer.vx;
                  vY = CM_buffer.vy;
                  vZ = CM_buffer.vz;
                  fprintf(filex,  "%.14lf ",  X);
                  fprintf(filey,  "%.14lf ",  Y);
                  fprintf(filez,  "%.14lf ",  Z);
                  fprintf(filevx, "%.14lf ", vX);
                  fprintf(filevy, "%.14lf ", vY);
                  fprintf(filevz, "%.14lf ", vZ);
            }
            if (write_elliptic_bool){
                  cart2ell(&CM_buffer, 0, alkhqp); //Computing the elliptic elements of the central body.
                  fprintf(filea, "%.14lf ", alkhqp[0]);
                  fprintf(filel, "%.14lf ", alkhqp[1]);
                  fprintf(filek, "%.14lf ", alkhqp[2]);
                  fprintf(fileh, "%.14lf ", alkhqp[3]);
                  fprintf(fileq, "%.14lf ", alkhqp[4]);
                  fprintf(filep, "%.14lf ", alkhqp[5]);
            }
            if (write_cartesian_bool || write_elliptic_bool){
                  m = CM_buffer.mass;
                  R = CM_buffer.radius;
                  fprintf(filer, "%.14lf ", R);
                  fprintf(filem, "%.14lf ", m);
            }
      }
      for (p = 0; p <= largest_id; p ++){
            if(*(exists + p)){
                  if (write_cartesian_bool){
                        X  = (moonlets + p) ->  x;
                        Y  = (moonlets + p) ->  y;
                        Z  = (moonlets + p) ->  z;
                        vX = (moonlets + p) -> vx;
                        vY = (moonlets + p) -> vy;
                        vZ = (moonlets + p) -> vz;
                        fprintf(filex,  "%.14lf ",  X);
                        fprintf(filey,  "%.14lf ",  Y);
                        fprintf(filez,  "%.14lf ",  Z);
                        fprintf(filevx, "%.14lf ", vX);
                        fprintf(filevy, "%.14lf ", vY);
                        fprintf(filevz, "%.14lf ", vZ);
                  }
                  if (write_elliptic_bool){
                        cart2ell(moonlets, p, alkhqp); //Computing the elliptic elements of body p
                        fprintf(filea,   "%.14lf ", alkhqp[0]);
                        fprintf(filel,   "%.14lf ", alkhqp[1]);
                        fprintf(filek,   "%.14lf ", alkhqp[2]);
                        fprintf(fileh,   "%.14lf ", alkhqp[3]);
                        fprintf(fileq,   "%.14lf ", alkhqp[4]);
                        fprintf(filep,   "%.14lf ", alkhqp[5]);
                  }
                  m  = (moonlets + p) ->   mass;
                  R  = (moonlets + p) -> radius;
                  if (write_cartesian_bool || write_elliptic_bool){
                        fprintf(filer, "%.14lf ", R);
                        fprintf(filem, "%.14lf ", m);
                  }
                  total_mass += m;
                  if (m > M1){
                        M2 = M1;
                        i2 = i1;
                        M1 = m;
                        i1 = p;
                  }
                  else if (m > M2){
                        M2 = m;
                        i2 = p;
                  }
            }
      }
      
      /******** Obtaining statistics regarding the two most massive bodies ********/
      if (how_many_moonlets == 0){
            M1 = 0.; M2 = 0.; R1 = 0.; R2 = 0.; a1 = 0.; a2 = 0.; l1 = 0.; l2 = 0.; K1 = 0.; K2 = 0.; H1 = 0.; H2 = 0.; Q1 = 0.; Q2 = 0.; P1 = 0.; P2 = 0.;
      }
      else if(how_many_moonlets == 1){
            M2 = 0.; R2 = 0.; a2 = 0.; l2 = 0.; K2 = 0.; H2 = 0.; Q2 = 0.; P2 = 0.;
            R1 = (moonlets + i1) -> radius;
            cart2ell(moonlets, i1, alkhqp);
            a1 = *alkhqp; l1 = alkhqp[1]; K1 = alkhqp[2]; H1 = alkhqp[3]; Q1 = alkhqp[4]; P1 = alkhqp[5];
      }
      else{
            R1 = (moonlets + i1) -> radius;  R2 = (moonlets + i2) -> radius;
            cart2ell(moonlets, i1, alkhqp);
            a1 = *alkhqp; l1 = alkhqp[1]; K1 = alkhqp[2]; H1 = alkhqp[3]; Q1 = alkhqp[4]; P1 = alkhqp[5];
            cart2ell(moonlets, i2, alkhqp);
            a2 = *alkhqp; l2 = alkhqp[1]; K2 = alkhqp[2]; H2 = alkhqp[3]; Q2 = alkhqp[4]; P2 = alkhqp[5];
      }

      /******** Writing statistics ********/
      fprintf(filestat, "%.14lf %d %d", time_elapsed + t_init, how_many_moonlets, collision_count);
      if (central_mass_bool){
            fprintf(filestat, " %.14lf", CM_buffer.mass);
      }
      else{
            fprintf(filestat, " %.1lf",  0.);
      }
      fprintf(filestat, " %.14lf", total_mass);
      if (inner_fluid_disk_bool){
            inner_fluid_disk_mass = fluid_disk_Sigma*M_PI*(Rroche*Rroche - R_unit*R_unit);
            fprintf(filestat, " %.14lf", inner_fluid_disk_mass);
      }
      else{
            fprintf(filestat, " %.1lf",  0.);
      }
      if (central_tides_bool || (J2_bool && J2_value == 0.)){
            fprintf(filestat, " %.14lf", SideralOmega);
      }
      else{
            fprintf(filestat, " %.1lf", 0.);
      }
      fprintf(filestat, " %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf",
      M1, M2, R1, R2, a1, a2, l1, l2, K1, K2, H1, H2, Q1, Q2, P1, P2);
      
      /******** To be removed ********/
      X       = moonlets -> x  - CM.x;
      Y       = moonlets -> y  - CM.y;
      Z       = moonlets -> z  - CM.z;
      vX      = moonlets -> vx - CM.vx;
      vY      = moonlets -> vy - CM.vy;
      vZ      = moonlets -> vz - CM.vz;
      m       = moonlets -> mass;
      typ v2  = vX*vX + vY*vY + vZ*vZ;
      typ r   = sqrt(X*X + Y*Y + Z*Z);
      typ mu  = G*(M_unit + m);
      typ HK  = 0.5*v2 - mu/r;
      typ r3  = r*r*r;
      typ r4  = r3*r;
      typ r5  = r4*r;
      typ tau = timestep;
      typ rv  = X*vX + Y*vY + Z*vZ;
      typ Hlf = -tau*tau/24.*(mu*v2/r3 - 3.*mu/r5*rv*rv - 2.*mu*mu/r4);
      typ Ht  = -0.5*k2*G*CM.mass*R_unit*R_unit*R_unit*R_unit*R_unit/(r*r5);
      fprintf(filestat, " %.14lf %.14lf %.14lf", HK, Hlf, Ht);

      /******** Terminating the lines ********/
      fprintf(filestat, "\n");
      if (write_cartesian_bool){
            fprintf(filex,  "\n");
            fprintf(filey,  "\n");
            fprintf(filez,  "\n");
            fprintf(filevx, "\n");
            fprintf(filevy, "\n");
            fprintf(filevz, "\n");
      }
      if (write_elliptic_bool){
            fprintf(filea,  "\n");
            fprintf(filel,  "\n");
            fprintf(filek,  "\n");
            fprintf(fileh,  "\n");
            fprintf(fileq,  "\n");
            fprintf(filep,  "\n");
      }
      if (write_cartesian_bool || write_elliptic_bool){
            fprintf(filer,  "\n");
            fprintf(filem,  "\n");
      }
      
      /******** Outputting datas relative to collisions ********/
      if (collision){
            typ t, m1, m2, R1, R2, DeltaV, impact_angle, m_tilde;
            for (p = 0; p < indexCollision; p ++){
                  t            = (collisionDatas + p) -> time;
                  m1           = (collisionDatas + p) -> m1;
                  m2           = (collisionDatas + p) -> m2;
                  R1           = (collisionDatas + p) -> R1;
                  R2           = (collisionDatas + p) -> R2;
                  DeltaV       = (collisionDatas + p) -> DeltaV;
                  impact_angle = (collisionDatas + p) -> impact_angle;
                  m_tilde      = (collisionDatas + p) -> m_tilde;
                  fprintf(filecol, "%.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf\n", t, m1, m2, R1, R2, DeltaV, impact_angle, m_tilde);
            }
            indexCollision = 0;
      }
      
      free(alkhqp);  alkhqp = NULL;
}


void readme(){

      /******** Creates a readme.txt file giving some indications about the content of the files ********/
      
      
      char file_path[800];
      FILE * file;
      
      strcpy(file_path, pth);
      strcat(file_path, "readme.txt");
      file = fopen(file_path, "w");
      if (file == NULL){
            fprintf(stderr, "Error : Cannot create or open file readme.txt.\n");
            abort();
      }
      
      fprintf(file, "Depending on the value that you chose for write_cartesian_bool and write_elliptic_bool in the parameter file, you should find a variable number of files here.\n");
      fprintf(file, "The two files named mass.txt and radius.txt contain the masses and radii of the bodies\n");
      fprintf(file, "Files x.txt, y.txt, z.txt, vx.txt, vy.txt, vz.txt contain the cartesians coordinates of the bodies.\n");
      fprintf(file, "Files a.txt, lambda.txt, k.txt, h.txt, q.txt, p.txt contain the elliptic coordinates of the bodies, where a is the semi-major axis and lambda the mean longitude.\n");
      fprintf(file, "The elements k, h, q and p are defined as (see Laskar & Robutel, 1995) k + ih = e*exp(i*varpi) and q + ip = sin(I/2)*exp(i*Omega),\n");
      fprintf(file, "where e is the eccentricity, I the inclination, varpi the longitude of the periapsis and Omega the longitude of the ascending node.\n");
      fprintf(file, "The eccentricity and longitude of the periapsis can be obtained from e = sqrt(k^2+h^2) and varpi = atan2(h,k).\n");
      fprintf(file, "The inclination and longitude of the ascending node can be obtained from I = 2*asin(sqrt(q^2+p^2)) and Omega = atan2(p,q).\n");
      fprintf(file, "NcorpiON uses these elliptic elements instead of the more traditional ones because they are properly defined even at zero eccentricity and inclination.\n");
      fprintf(file, "\n");
      fprintf(file, "These files contain one line per output and one column per body. For example, the j^th column of the n^th line of x.txt contains the x-coordinate of the j^th\n");
      fprintf(file, "body of the simulation when the n^th output occured. The number of columns can therefore be variable if the number of bodies was variable throughout your\n");
      fprintf(file, "simulation, for example if you used the fragmentation model to resolve collisions, or simply if you lost bodies along the way.\n");
      fprintf(file, "The number of lines of these files depends on the length of the simulation (set with t_end in the parameter file), the timestep (set with time_step)\n");
      fprintf(file, "and the frequency of outputs (set with output_step). For example, if you set output_step to 5, then you have one line for every five timesteps.\n");
      fprintf(file, "If you set central_mass_bool to 1 in the parameter file, then the central body of your simulation is always displayed in the first column.\n");
      fprintf(file, "If you set viscoelastic_bool to 1 in the parameter file, then the perturbing body of your simulation is always displayed in the first column.\n");
      fprintf(file, "\n");
      fprintf(file, "These files are given in the simulation's units, which you defined when you set the values of G, M_unit and R_unit in the parameter file. If you set\n");
      fprintf(file, "random_initial_bool to 0, then your simulation used initial conditions from the file init.txt that you put here, and you should make sure that the units you\n");
      fprintf(file, "used for this file are consistent with your choice of G, M_unit and R_unit. Similarly, if you set random_initial_bool to 1, make sure that the bounds you\n");
      fprintf(file, "set in the parameter file to draw the initial conditions at given in the right units.\n");
      fprintf(file, "\n");
      fprintf(file, "If you set the boolean resume_simulation_bool to 1, then you will find a file init.txt, or the one that was already here has been modified.\n");
      fprintf(file, "You can use it to restart the simulation where it stopped. Just set random_initial_bool to 0 and initial_cartesian_bool to 1. You will also need to\n");
      fprintf(file, "update some parameters like N_0, to restart the simulation exactly like it was left.\n");
      fprintf(file, "\n");
      fprintf(file, "If you set viscoelastic_bool to 1 and resume_simulation_bool to 1, then you will have two files named init.txt and connections.txt that you can use to restart\n");
      fprintf(file, "your simulation. This way, you can run a first simulation with pert_mass to 0.0 and random_initial_bool to 1 to generate a viscoelastic body and let it rest,\n");
      fprintf(file, "and then a second simulation with random_initial_bool to 0 and a non-zero value to pert_mass to simulate the same viscoelastic body being perturbed by the tides\n");
      fprintf(file, "raised by a perturbing body, whose orbit you define in the parameter file.\n");
      fprintf(file, "\n");
      fprintf(file, "You will also find a file stat.txt containing some statistics about the simulation. Like the other files, it contains one line per output,\n");
      fprintf(file, "but unlike the other files, it has a constant number of columns. Its columns are :\n");
      fprintf(file, "Time, N, collisions, central mass, total bodies mass, inner fluid disk mass, sideral rotation, M1, M2, R1, R2, a1, a2, l1, l2, k1, k2, h1, h2, q1, q2, p1, p2.\n");
      fprintf(file, "The last 16 columns are the masses, radii and orbital elements of the most massive and second most massive bodies (central body excluded).\n");
      fprintf(file, "Many of these columns can be 0, depending of the booleans values that you set in the parameter file. The central radius stays at R_unit throughout the simulation\n");
      fprintf(file, "\n");
      fprintf(file, "Finally, the file collisions.txt, that you will find if collision_bool and write_collisions_bool are set to 1, contains statistics about the collisions.\n");
      fprintf(file, "It contains one line per collision, and its columns are (in simulation's units and radians):\n");
      fprintf(file, "Time, impactor mass, target mass, impactor radius, target radius, Delta V, impact angle, mass of the largest remnant.\n");
      fprintf(file, "The time of collision is found assuming that the bodies were moving on a staight line between two timesteps. It can be between two timesteps and smaller than\n");
      fprintf(file, "the initial time if bodies were overlapping at the beginning of the simulation.");
      fclose(file);
}


void rebound(struct moonlet * moonlets){

      /******** Communicates with REBOUND for the real-time visualization of  ********/
      /******** the simulation on port browser_port of the user's web browser ********/

      int j;
      int N = central_mass_bool || (viscoelastic_bool && pert_mass > 0.);
      
      /******** Filling the buffer with the bodies' coordinates ********/
      if(central_mass_bool || (viscoelastic_bool && pert_mass > 0.)){
            sending_buffer[0] = CM.x;
            sending_buffer[1] = CM.y;
            sending_buffer[2] = CM.z;
            sending_buffer[3] = CM.vx;
            sending_buffer[4] = CM.vy;
            sending_buffer[5] = CM.vz;
            sending_buffer[6] = CM.mass;
            sending_buffer[7] = CM.radius;
      }
      for (j = 0; j <= largest_id; j ++){
            if (*(exists + j)){
                  sending_buffer[8*N]     = (moonlets + j) -> x;
                  sending_buffer[8*N + 1] = (moonlets + j) -> y;
                  sending_buffer[8*N + 2] = (moonlets + j) -> z;
                  sending_buffer[8*N + 3] = (moonlets + j) -> vx;
                  sending_buffer[8*N + 4] = (moonlets + j) -> vy;
                  sending_buffer[8*N + 5] = (moonlets + j) -> vz;
                  sending_buffer[8*N + 6] = (moonlets + j) -> mass;
                  sending_buffer[8*N + 7] = (moonlets + j) -> radius;
                  N ++;
            }
      }
      
      /******** Sending the buffer to REBOUND ********/
      #if openGL_bool
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      MPI_Send(&N, 1, MPI_INT, !rank, 0, MPI_COMM_WORLD);
      MPI_Send(sending_buffer, 8*N, MPI_TYP, !rank, 0, MPI_COMM_WORLD);
      #endif
}


void resume(struct moonlet * moonlets){

      /******** Generates a file init.txt at the path indicated in the parameters file.********/
      /******** The file is used to resume the simulation if the user sets the boolean ********/
      /******** resume_simulation_bool to 1. This function is called at the very end   ********/
      /******** of the simulation, before deallocating the array moonlets. If          ********/
      /******** central_mass_bool is 1, the central mass is not put in the file but    ********/
      /******** the bodies' coordinates are written relative to the central mass.      ********/

      
      int j;
      typ X, Y, Z, vX, vY, vZ, m, R, XX, YY, ZZ, vXX, vYY, vZZ;
      char init_path[800];
      strcpy(init_path, pth);
      strcat(init_path, "init.txt");
      FILE * init_file = fopen(init_path, "w");
      if (init_file == NULL){
            fprintf(stderr, "Error : Cannot create or open file init.txt in function resume.\n");
            abort();
      }
      
      /******** Getting the central body's coordinates ********/
      XX  = (central_mass_bool ? CM.x  : 0.);
      YY  = (central_mass_bool ? CM.y  : 0.);
      ZZ  = (central_mass_bool ? CM.z  : 0.);
      vXX = (central_mass_bool ? CM.vx : 0.);
      vYY = (central_mass_bool ? CM.vy : 0.);
      vZZ = (central_mass_bool ? CM.vz : 0.);
      
      for (j = 0; j <= largest_id; j ++){
            if (*(exists + j)){
                  X  = (moonlets + j) -> x  - XX;
                  Y  = (moonlets + j) -> y  - YY;
                  Z  = (moonlets + j) -> z  - ZZ;
                  vX = (moonlets + j) -> vx - vXX;
                  vY = (moonlets + j) -> vy - vYY;
                  vZ = (moonlets + j) -> vz - vZZ;
                  m  = (moonlets + j) -> mass;
                  R  = (moonlets + j) -> radius;                  
                  fprintf(init_file, "%.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf %.14lf\n", X, Y, Z, vX, vY, vZ, m, R);
            }
      }     
      fclose(init_file);
      
      #if viscoelastic_bool
            /******** Printing the connections on the file pth/connections.txt ********/ 
            char connection_path[800];
            int a, b;
            typ rest_length;
            strcpy(connection_path, pth);
            strcat(connection_path, "connections.txt");
            FILE * connection_file = fopen(connection_path, "w");
            if (connection_file == NULL){
                  fprintf(stderr, "Error : Cannot create or open file connections.txt in function resume.\n");
                  abort();
            }
            for (j = 0; j < N_connections; j ++){
                  rest_length = (connections + j) -> rest_length;
                  a           = (connections + j) -> Pair.fst;
                  b           = (connections + j) -> Pair.snd;
                  fprintf(connection_file, "%.1lf %.1lf %.14lf\n", (typ) a, (typ) b, rest_length);
            }
            fclose(connection_file);
      #endif
}
