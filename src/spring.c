/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    spring.c                                                    ********/
/******** @brief   Manages visco-elasticity between the bodies                 ********/
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

struct connection * connections = NULL;
struct chain * first  = NULL;
struct chain * second = NULL;
int N_connections = 0;
typ shapeV;


void precision(struct moonlet * viscoelastic){

      /******** To be removed ********/

      int i,j;
      typ Xi, Yi, Zi, mi, Xj, Yj, Zj, mj, dx, dy, dz, r, r3, aX, aY, aZ;
      struct boxdot * root = NULL;
      
      typ * acc = (typ *)malloc(3*N_0*sizeof(typ));
      if (acc == NULL){
            fprintf(stderr, "Cannot allocate memory for array in function precision\n");
            abort();
      }
      for (i = 0; i < 3*N_0; i ++){
            acc[i] = 0.;
      }
      
      char path[800];
      strcpy(path, pth);
      strcat(path, "acc.txt");
      FILE * file = fopen(path, "r");
      if (file == NULL){
            file = fopen(path, "w");
            if (file == NULL){
                  fprintf(stderr, "Cannot open file in function precision\n");
                  abort();
            }
            clock_t t0, t1;
            t0 = clock();
            for (i = 0; i < N_0; i ++){
                  Xi = (viscoelastic + i) -> x;
                  Yi = (viscoelastic + i) -> y;
                  Zi = (viscoelastic + i) -> z;
                  mi = (viscoelastic + i) -> mass;
                  for (j = 0; j < i; j ++){
                        Xj = (viscoelastic + j) -> x;
                        Yj = (viscoelastic + j) -> y;
                        Zj = (viscoelastic + j) -> z;
                        mj = (viscoelastic + j) -> mass;
                        dx = Xi - Xj;  dy = Yi - Yj;  dz = Zi - Zj;
                        r  = sqrt(dx*dx + dy*dy + dz*dz);
                        r3 = r*r*r;
                        aX = G*dx/r3;
                        aY = G*dy/r3;
                        aZ = G*dz/r3;
                        acc[3*i]     -= aX*mj;
                        acc[3*i + 1] -= aY*mj;
                        acc[3*i + 2] -= aZ*mj;
                        acc[3*j]     += aX*mi;
                        acc[3*j + 1] += aY*mi;
                        acc[3*j + 2] += aZ*mi;
                  }
                  if (!(i%100)){
                        printf("%d/%d\n", i, N_0);
                  }
            }
            t1 = clock();
            printf("Brute-force = %.3lf seconds\n", ((typ) (t1 - t0))/(typ) CLOCKS_PER_SEC);
            for (i = 0; i < N_0; i ++){
                  fprintf(file, "%.14lf %.14lf %.14lf\n", acc[3*i], acc[3*i + 1], acc[3*i + 2]);
            }
            fclose(file);
      }
      else{
            fclose(file);
            readFromFile(path, acc, 3*N_0);

            root     = root_cell(viscoelastic);
            FlatTree = flattree_init(root);
            clear_boxdot(&root);
            com_flattree(FlatTree, viscoelastic);
            Mtot     = FlatTree -> M0;
            rmax_flattree(FlatTree, viscoelastic);
            rcrit_flattree(FlatTree);
            tensor_initialization(FlatTree);
            #if expansion_order >= 3 && mutual_bool
                  multipole_flattree(FlatTree, viscoelastic);
            #endif
            for (j = 0; j < cell_id; j ++){
                  free((FlatTree + j) -> dots);
                  (FlatTree + j) -> dots = NULL;
            }
            free(FlatTree);
            FlatTree       = NULL;
            how_many_cells = 0;
            cell_id        = 0;
            
            int n_timestep = 8;
            typ force[n_timestep];
            typ colli[n_timestep];
            clock_t t0, t1, t2;
            int J;
            for (J = 0; J < n_timestep; J ++){
                  t0       = clock();
                  root     = root_cell(viscoelastic);
                  FlatTree = flattree_init(root);
                  clear_boxdot(&root);
                  com_flattree(FlatTree, viscoelastic);
                  Mtot     = FlatTree -> M0;
                  rmax_flattree(FlatTree, viscoelastic);
                  rcrit_flattree(FlatTree);
                  tensor_initialization(FlatTree);
                  #if expansion_order >= 3 && mutual_bool
                        multipole_flattree(FlatTree, viscoelastic);
                  #endif
                  #if mutual_bool
                  Cm_flattree(FlatTree, viscoelastic);              
                  Cm_downtree(FlatTree, viscoelastic);
                  #endif
                  t1 = clock();
                  
                  center_and_maxR_flattree(FlatTree, viscoelastic);
                  rmax_and_rcrit_flattree (FlatTree, viscoelastic);
                  collision_flattree      (FlatTree, viscoelastic);

                  for (j = 0; j < cell_id; j ++){
                        free((FlatTree + j) -> dots);
                        (FlatTree + j) -> dots = NULL;
                  }
                  free(FlatTree);
                  FlatTree       = NULL;
                  how_many_cells = 0;
                  cell_id        = 0;
                  t2 = clock();
                  
                  force[J] = ((typ) (t1 - t0))/(typ) CLOCKS_PER_SEC;
                  colli[J] = ((typ) (t2 - t1))/(typ) CLOCKS_PER_SEC;
                  
                  printf("%d/%d\n", J + 1, n_timestep);
                  
                  if (J == 0){
                        for (i = 0; i < 3*N_0; i ++){
                              acc[i] = C1Moonlets[i];
                        }
                  }
            }
            typ force_calc   = 0.;
            typ colli_search = 0.;
            for (J = 0; J < n_timestep; J ++){
                  force_calc   += force[J];
                  colli_search += colli[J];
            }
            force_calc   /= (typ) n_timestep;
            colli_search /= (typ) n_timestep;
            printf("Force calculation   = %.3lf seconds\n", force_calc);
            printf("Collision treatment = %.3lf seconds\n", colli_search);
            printf("Total               = %.3lf seconds\n", force_calc + colli_search);

            char path_falcON[800];
            strcpy(path_falcON, pth);
            strcat(path_falcON, "acc_falcON");
            char expansion[20];
            char Theta_min[20];
            sprintf(expansion, "%d", expansion_order);
            sprintf(Theta_min, "%.2lf", theta_min);
            strcat(path_falcON, "_p=");
            strcat(path_falcON, expansion);
            strcat(path_falcON, "_theta_min=");
            strcat(path_falcON, Theta_min);
            strcat(path_falcON, ".txt");
            FILE * file_falcON = fopen(path_falcON, "w");
            if (file_falcON == NULL){
                  fprintf(stderr, "Cannot open file file_falcON in function precision\n");
                  abort();
            }
            for (i = 0; i < N_0; i ++){
                  fprintf(file_falcON, "%.14lf %.14lf %.14lf\n", acc[3*i], acc[3*i + 1], acc[3*i + 2]);
            }
            fclose(file_falcON);
      }
      free(acc);  acc = NULL;
}


struct moonlet * generate_visco_elastic_body(){

      /******** Similar to the function populate in structure.c, but populates    ********/
      /******** the simulation with the particles of a visco-elastic body instead ********/

      int j;
      typ maxD = 0.;
      typ costheta, costheta_max;
      typ distance, distance2, radius;
      typ X, Y, Z;
      int currentN = 0;
      int closest_vertice;
      typ m = M_unit/(typ) N_0;
      
      printf("Generating viscoelastic body\n");

      /******** Generating the array of particles. This will be returned by the function ********/
      struct moonlet * viscoelastic = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //The array of particles. Some will be linked together
      if (viscoelastic == NULL){
            fprintf(stderr, "Error : Cannot allocate array of particles in function generate_visco_elastic_body.\n");
            abort();
      }
      
      if (random_initial_bool){
      
            /******** Trying to retrieve the vertices from the file pth/shape_model.txt ********/
            char fileOfVertices[800]; 
            strcpy(fileOfVertices, pth);
            strcat(fileOfVertices, "shape_model.txt");
            FILE * file = fopen(fileOfVertices, "r");
            if (file == NULL){
                  printf("Could not find the file shape_model.txt in the path 'pth'. NcorpiON will use a sphere of radius R_unit = %lf instead.\n", R_unit);
                  shapeV = 4.0*M_PI*R_unit*R_unit*R_unit/3.0; //The volume of the body is that of a sphere
                  /******** Drawing the points of the viscoelastic from a sphere of radius R_unit ********/
                  while (currentN < N_0){
                        X = rdm(-R_unit, R_unit);
                        Y = rdm(-R_unit, R_unit);
                        Z = rdm(-R_unit, R_unit);
                        if (X*X + Y*Y + Z*Z < R_unit*R_unit){
                              (viscoelastic + currentN) -> x      = X;
                              (viscoelastic + currentN) -> y      = Y;
                              (viscoelastic + currentN) -> z      = Z;
                              (viscoelastic + currentN) -> vx     = 0.;
                              (viscoelastic + currentN) -> vy     = 0.;
                              (viscoelastic + currentN) -> vz     = 0.;
                              (viscoelastic + currentN) -> mass   = m;
                              currentN ++;
                        }
                  }
            }
            else{
                  fclose(file);
                  /******** Generating the array of vertices ********/
                  typ * vertices = NULL;
                  int n_vertices;
                  vertices = readFromFile_withoutConstraint(fileOfVertices, &n_vertices);
                  if (n_vertices % 3 != 0){
                        fprintf(stderr, "Error : The file at path 'pth/shape_model.txt' must have exactly three columns per line.\n");
                        abort();
                  }
                  n_vertices /= 3;
                  if (n_vertices*N_0 > 100000000 || n_vertices*N_0 < 0){
                        printf("Warning : With %d vertices in the shape model and %d nodes, the generation may take a while.\n", n_vertices, N_0);
                  }
      
                  /******** Getting the largest distance between a vertice of the shape-model and the origin ********/
                  for (j = 0; j < n_vertices; j ++){
                        distance = vertices[3*j]*vertices[3*j] + vertices[3*j + 1]*vertices[3*j + 1] + vertices[3*j + 2]*vertices[3*j + 2];
                        if (distance > maxD){
                              maxD = distance;
                        }
                  }
                  maxD = sqrt(maxD);
      
                  /******** Drawing the points of the viscoelastic body from the shape-model ********/
                  typ monteCarloProbability; //Using a Monte-Carlo method to estimate the volume of the shape-model
                  unsigned int N_total = 0;
                  while (currentN < N_0){
                        N_total ++;
                        X = rdm(-maxD, maxD);
                        Y = rdm(-maxD, maxD);
                        Z = rdm(-maxD, maxD);
                        distance        = sqrt(X*X + Y*Y + Z*Z);
                        distance2       = sqrt(vertices[0]*vertices[0] + vertices[1]*vertices[1] + vertices[2]*vertices[2]);
                        costheta_max    = (X*vertices[0] + Y*vertices[1] + Z*vertices[2])/(distance * distance2);
                        closest_vertice = 0;
                        for (j = 1; j < n_vertices; j ++){ //Getting the id of the vertice most collinear with the current point
                              distance2 = sqrt(vertices[3*j]*vertices[3*j] + vertices[3*j + 1]*vertices[3*j + 1] + vertices[3*j + 2]*vertices[3*j + 2]);
                              costheta  = (X*vertices[3*j] + Y*vertices[3*j + 1] + Z*vertices[3*j + 2])/(distance * distance2);
                              if (costheta > costheta_max){
                                    costheta_max = costheta;
                                    closest_vertice = j;
                              }
                        }
                        distance2 = sqrt(vertices[3*closest_vertice]*vertices[3*closest_vertice] 
                        + vertices[3*closest_vertice + 1]*vertices[3*closest_vertice + 1] + vertices[3*closest_vertice + 2]*vertices[3*closest_vertice + 2]);
                        if (distance < distance2){ //The point is inside the viscoelastic body
                              (viscoelastic + currentN) -> x      = X;
                              (viscoelastic + currentN) -> y      = Y;
                              (viscoelastic + currentN) -> z      = Z;
                              (viscoelastic + currentN) -> vx     = 0.;
                              (viscoelastic + currentN) -> vy     = 0.;
                              (viscoelastic + currentN) -> vz     = 0.;
                              (viscoelastic + currentN) -> mass   = m;
                              currentN ++;
                        }
                  }
                  free(vertices);  vertices = NULL;
                  monteCarloProbability = ((typ) N_0) / ((typ) N_total);
                  shapeV                = monteCarloProbability*8.0*maxD*maxD*maxD; //Estimation of the volume of the shape-model
            }
      
            /******** Making sure that particles do not overlap ********/
            radius = 0.5*minimal_distance*pow(shapeV/(typ) N_0, 1.0/3.0);
            for (j = 0; j < N_0; j ++){
                  (viscoelastic + j) -> radius = radius;
            }
            overlap(viscoelastic);
            
            /******** Setting the final radii ********/
            for (j = 0; j < N_0; j ++){
                  (viscoelastic + j) -> radius = 2.0*nodes_radius*radius;
            }
      }
      else{
            /******** Generating particles from files ********/
            char fileOfIC[800]; 
            strcpy(fileOfIC, pth);
            strcat(fileOfIC, "init.txt");
            typ * IC = (typ *)malloc(8*N_0*sizeof(typ));
            if (IC == NULL){
                  fprintf(stderr, "Error : Cannot allocate array for initial conditions in function generate_visco_elastic_body.\n");
                  abort();
            }
            readFromFile(fileOfIC, IC, 8*N_0);
            for (j = 0; j < N_0; j ++){
                  (viscoelastic + j) -> x      = *(IC + 8*j);
                  (viscoelastic + j) -> y      = *(IC + 8*j + 1);
                  (viscoelastic + j) -> z      = *(IC + 8*j + 2);
                  (viscoelastic + j) -> vx     = *(IC + 8*j + 3);
                  (viscoelastic + j) -> vy     = *(IC + 8*j + 4);
                  (viscoelastic + j) -> vz     = *(IC + 8*j + 5);
                  (viscoelastic + j) -> mass   = *(IC + 8*j + 6);
                  (viscoelastic + j) -> radius = *(IC + 8*j + 7);
            }
            free(IC);
            IC = NULL;
      }
      
      /******** Filling the unused cells of the body array with whatever ********/
      for (j = N_0; j < N_max; j ++){
            *(viscoelastic + j) = init(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      }
      
      /******** Reducing to the center of mass ********/
      if (reduce_to_COM_bool){
            typ com[6];
            total_momentum(viscoelastic, com); //Getting the speed and position of the center of mass
            for (j = 0; j < N_0; j ++){        //Substracting the center of mass
                  (viscoelastic + j) -> x  -= com[0];
                  (viscoelastic + j) -> y  -= com[1];
                  (viscoelastic + j) -> z  -= com[2];
                  (viscoelastic + j) -> vx -= com[3];
                  (viscoelastic + j) -> vy -= com[4];
                  (viscoelastic + j) -> vz -= com[5];
            }
      }
      
      /******** Generating the connections ********/
      generate_connections(viscoelastic);
      printf("Number of connections = %d\n", N_connections);
      
      /******** Making the viscoelastic body rotate ********/
      make_rotate(viscoelastic);
      
      /******** Pointing the viscoelastic body's angular momentum towards a particular direction in the inertial frame ********/
      point_angular_momentum(viscoelastic);
      
      /******** Communicating with REBOUND for 3D visualization ********/
      if (openGL_bool){
            rebound(viscoelastic);
      }

      return viscoelastic;
}


void generate_connections(struct moonlet * viscoelastic){

      /******** Called after generate_visco_elastic_body. Generates the connections ********/
      
      
      struct boxdot * root = NULL;
      int i, j, k, index1, index2, index;
      struct chain * fst = NULL;
      struct chain * snd = NULL;
      typ connecting_distance = pow(3.0*connections_per_node*shapeV/(4.0*M_PI*(typ) N_0), 1.0/3.0);
      typ radius = 0.5*minimal_distance*pow(shapeV/(typ) N_0, 1.0/3.0);
      
      if (random_initial_bool){
            /******** Connecting nodes ********/
            if (brute_force_bool){
                  for (i = 0; i < N_0; i ++){
                        for (j = 0; j < i; j ++){
                              if (connects(viscoelastic, i, j)){
                                    N_connections ++;
                                    add(i, &first);
                                    add(j, &second);
                              }
                        }
                  }
            }
            else{ //FalcON
                  for (j = 0; j < N_0; j ++){
                        (viscoelastic + j) -> radius = 0.5*connecting_distance;
                  }
                  root     = root_cell(viscoelastic);
                  FlatTree = flattree_init(root);
                  clear_boxdot(&root);
                  center_and_maxR_flattree(FlatTree, viscoelastic);
                  rmax_and_rcrit_flattree (FlatTree, viscoelastic);
                  viscoelastic_flattree   (FlatTree, viscoelastic, 0);
                  for (j = 0; j < N_0; j ++){
                        (viscoelastic + j) -> radius = 2.0*nodes_radius*radius;
                  }
            }
            if(falcON_bool){ //Resetting the octree and the associated data
                  for (j = 0; j < cell_id; j ++){
                        free((FlatTree + j) -> dots);
                        (FlatTree + j) -> dots = NULL;
                  }
                  free(FlatTree);
                  FlatTree       = NULL;
                  how_many_cells = 0;
                  cell_id        = 0;
            }
            
            /******** Looking for particles not connected enough and connecting them more ********/
            int indexes[3];
            int * connection_per_body = (int *)malloc(N_0*sizeof(int));
            if (connection_per_body == NULL){
                  fprintf(stderr, "Cannot allocate memory for array connection_per_body in function generate_connections\n");
                  abort();
            }
            for (j = 0; j < N_0; j ++){
                  connection_per_body[j] = 0;
            }
            index1 = first  -> how_many - 1;
            index2 = second -> how_many - 1;
            fst    = first;
            snd    = second;
            while (fst != NULL){
                  i = (fst -> ids)[index1];
                  j = (snd -> ids)[index2];
                  connection_per_body[i] ++;
                  connection_per_body[j] ++;
                  index1 --;  index2 --; //Switching to the next pair of connecting particles
                  if (index1 < 0){
                        fst = fst -> queue;
                        index1 = max_ids_per_node - 1;
                  }
                  if (index2 < 0){
                        snd = snd -> queue;
                        index2 = max_ids_per_node - 1;
                  }
            }
            for (j = 0; j < N_0; j ++){
                  if (connection_per_body[j] == 0){      //That node needs three more connections
                        three_closest_nodes(viscoelastic, j, indexes);
                        N_connections += 3;
                        add(j, &first);  add(* indexes,      &second);
                        add(j, &first);  add(*(indexes + 1), &second);
                        add(j, &first);  add(*(indexes + 2), &second);
                  }
                  else if (connection_per_body[j] == 1){ //That node needs two more connections
                        three_closest_nodes(viscoelastic, j, indexes);
                        N_connections += 2;
                        add(j, &first);  add(*(indexes + 1), &second);
                        add(j, &first);  add(*(indexes + 2), &second);
                  }
                  else if (connection_per_body[j] == 2){ //That node needs one more connection
                        three_closest_nodes(viscoelastic, j, indexes);
                        N_connections ++;
                        add(j, &first);  add(*(indexes + 2), &second);
                  }
            }
            free(connection_per_body);  connection_per_body = NULL;
      
            /******** Generating the array of connections ********/
            index       = 0;
            index1      = first  -> how_many - 1;
            index2      = second -> how_many - 1;
            connections = (struct connection *)malloc(N_connections*sizeof(struct connection));
            if (connections == NULL){
                  fprintf(stderr, "Cannot allocate array of %d connections in function generate_connections\n", N_connections);
                  abort();
            }
            fst = first;
            snd = second;
            while (fst != NULL){
                  i = (fst -> ids)[index1];
                  j = (snd -> ids)[index2];
                  *(connections + index) = make_connection(viscoelastic, i, j);
                  index1 --;  index2 --; //Switching to the next pair of connecting particles
                  if (index1 < 0){
                        fst = fst -> queue;
                        index1 = max_ids_per_node - 1;
                  }
                  if (index2 < 0){
                        snd = snd -> queue;
                        index2 = max_ids_per_node - 1;
                  }
                  index ++;
            }
            clear_chain(&first);
            clear_chain(&second);
      }
      else{
            typ rest_length;
            char fileOfIC[800];
            strcpy(fileOfIC, pth);
            strcat(fileOfIC, "connections.txt");
            typ * IC = NULL;
            IC = readFromFile_withoutConstraint(fileOfIC, &N_connections);
            if (N_connections % 3){
                  fprintf(stderr, "There was a problem in reading the file connections.txt in function generate_connections.\n");
                  abort();
            }
            N_connections /= 3;
            connections = (struct connection *)malloc(N_connections*sizeof(struct connection));
            if (connections == NULL){
                  fprintf(stderr, "Cannot allocate array of %d connections in function generate_connections.\n", N_connections);
                  abort();
            }
            for (k = 0; k < N_connections; k ++){
                  i           = (int) (*(IC + 3*k));
                  j           = (int) (*(IC + 3*k + 1));
                  rest_length =        *(IC + 3*k + 2);
                  (connections + k) -> Pair.fst           = i;
                  (connections + k) -> Pair.snd           = j;
                  (connections + k) -> rest_length        = rest_length;
            }
            free(IC);  IC = NULL;
      }
}


void make_rotate(struct moonlet * viscoelastic){

      /******** Makes the viscoelastic body rotate with angular frequency Omega = (OmegaX, OmegaY, OmegaZ) ********/
      /******** in the fixed body reference frame. The vector Omega is given in the parameters file        ********/

      if (OmegaX == 0. && OmegaY == 0. && OmegaZ == 0.){ //No rotation in this case
            return;
      }
      if (!random_initial_bool){
            return;
      }

      int j;
      typ X, Y, Z;
      typ QxOm[3];
      
      for (j = 0; j < N_0; j ++){
            X = (viscoelastic + j) -> x;
            Y = (viscoelastic + j) -> y;
            Z = (viscoelastic + j) -> z;
            cross_product(OmegaX, OmegaY, OmegaZ, X, Y, Z, QxOm);
            (viscoelastic + j) -> vx = QxOm[0];
            (viscoelastic + j) -> vy = QxOm[1];
            (viscoelastic + j) -> vz = QxOm[2];
      }
}


void point_angular_momentum(struct moonlet * viscoelastic){

      /******** Computes the direction of the angular momentum in the body-attached frame and then rotates ********/
      /******** the whole viscoelastic body to make it match the desired direction in the inertial frame   ********/

      if (OmegaX == 0. && OmegaY == 0. && OmegaZ == 0.){ //No rotation in this case
            return;
      }
      if (!random_initial_bool){
            return;
      }

      int j;
      typ X, Y, Z, vX, vY, vZ, m, gx, gy, gz, xr, yr, zr, vxr, vyr, vzr;
      typ g[3] = {0., 0., 0.};
      typ rxv[3];
      struct quaternion q;
      
      /******** Computing the angular momentum ********/
      for (j = 0; j < N_0; j ++){
            X  = (viscoelastic + j) -> x;
            Y  = (viscoelastic + j) -> y;
            Z  = (viscoelastic + j) -> z;
            vX = (viscoelastic + j) -> vx;
            vY = (viscoelastic + j) -> vy;
            vZ = (viscoelastic + j) -> vz;
            m  = (viscoelastic + j) -> mass;
            cross_product(m*X, m*Y, m*Z, vX, vY, vZ, rxv);
            g[0] += rxv[0];  g[1] += rxv[1];  g[2] += rxv[2];
      }
      
      /******** Generating a quaternion that rotates the whole body to make the angular momentum match its desired orientation ********/
      gx = cos(lbd_long)*cos(beta_lat);
      gy = sin(lbd_long)*cos(beta_lat);
      gz = sin(beta_lat);
      q  = get_quaternion(g[0], g[1], g[2], gx, gy, gz);
      
      /******** Rotating the whole viscoelastic body according to the quaternion q ********/
      for (j = 0; j < N_0; j ++){
            X  = (viscoelastic + j) -> x;
            Y  = (viscoelastic + j) -> y;
            Z  = (viscoelastic + j) -> z;
            vX = (viscoelastic + j) -> vx;
            vY = (viscoelastic + j) -> vy;
            vZ = (viscoelastic + j) -> vz;
            rotate_with_quaternion( X,  Y,  Z, q,  &xr,  &yr,  &zr);
            rotate_with_quaternion(vX, vY, vZ, q, &vxr, &vyr, &vzr);
            (viscoelastic + j) -> x  = xr;
            (viscoelastic + j) -> y  = yr;
            (viscoelastic + j) -> z  = zr;
            (viscoelastic + j) -> vx = vxr;
            (viscoelastic + j) -> vy = vyr;
            (viscoelastic + j) -> vz = vzr;
      }
}


void overlap(struct moonlet * viscoelastic){

      /******** Makes sure that no two particle overlap ********/

      int collision_count_before = -1;
      int i, j;
      int count = 0;
      struct boxdot * root = NULL;
      
      while(collision_count > collision_count_before){
            collision_count_before = collision_count;            
            if (brute_force_bool){
                  for (i = 0; i < N_0; i ++){
                        for (j = 0; j < i; j ++){
                              deOverlap(viscoelastic, i, j);
                        }
                  }
            }
            else{ //FalcON
                  root     = root_cell(viscoelastic);
                  FlatTree = flattree_init(root);
                  clear_boxdot(&root);
                  center_and_maxR_flattree(FlatTree, viscoelastic);
                  rmax_and_rcrit_flattree (FlatTree, viscoelastic);                   
                  viscoelastic_flattree   (FlatTree, viscoelastic, 1);   
            }
            if(falcON_bool){ //Resetting the octree and the associated data
                  for (j = 0; j < cell_id; j ++){
                        free((FlatTree + j) -> dots);
                        (FlatTree + j) -> dots = NULL;
                  }
                  free(FlatTree);
                  FlatTree       = NULL;
                  how_many_cells = 0;
                  cell_id        = 0;
            }
            count ++;            
            if (count > 1024){
                  fprintf(stderr, "Error : While loop never exits in function overlap. Maybe minimal_distance is too large in the parameter file?\n");
                  abort();
            }
      }
      collision_count = 0;
}


void deOverlap(struct moonlet * viscoelastic, int a, int b){

      /******** A version adapted to locate overlapping particles in the viscoelastic body ********/

      
      typ xa, ya, za, xb, yb, zb;                            //Cartesian positions of the bodies at the beginning of the timestep
      typ dx, dy, dz;                                        //dx is xa-xb, and so on.
      typ R_a = (viscoelastic + a) -> radius;                //The bodies' radii
      typ R_b = (viscoelastic + b) -> radius;
      typ dr;
      
      /******** Getting the positions ********/
      xa = (viscoelastic + a) -> x;
      ya = (viscoelastic + a) -> y;
      za = (viscoelastic + a) -> z;
      xb = (viscoelastic + b) -> x;
      yb = (viscoelastic + b) -> y;
      zb = (viscoelastic + b) -> z;
      
      /******** Getting the relative positions ********/
      dx  = xa - xb;
      dy  = ya - yb;
      dz  = za - zb;
      
      /******** Checking if they overlap ********/
      dr = sqrt(dx*dx + dy*dy + dz*dz);
      if (dr + 1.0e-13 > R_a + R_b){ //No overlap
            return;
      }
      
      if (dr == 0.){
            dx = R_a + R_b;
            dy = 0.;
            dz = 0.;
      }
      else{
            dx *= (R_a + R_b)/dr;
            dy *= (R_a + R_b)/dr;
            dz *= (R_a + R_b)/dr;
      }
            
      /******** Updating positions. I move the particle closest from the origin so it doesn't get outside the body ********/
      if (xa*xa + ya*ya + za*za < xb*xb + yb*yb + zb*zb){
            (viscoelastic + a) -> x = xb + dx;
            (viscoelastic + a) -> y = yb + dy;
            (viscoelastic + a) -> z = zb + dz;
      }
      else{
            (viscoelastic + b) -> x = xa - dx;
            (viscoelastic + b) -> y = ya - dy;
            (viscoelastic + b) -> z = za - dz;
      }
      collision_count ++;              
}


void viscoelastic_flattree(struct node * FlatTree, struct moonlet * viscoelastic, int generating){

      /******** Tree walk for detecting close particles in time proportional to N ********/
      /******** If generating is 1, then the function is called to make sure that ********/
      /******** the particles inside the viscoelastic body are not too close.     ********/
      /******** Otherwise it is called to determine what particles will connect   ********/
      
      /******** It is hard to tell in advance how many pairs will be treated by the tree walk   ********/
      /******** I expect it will be at most factor * cell_id, but that might have to be changed ********/
      int factor = (int) floor(250.0 * 0.5 / theta_min);
      struct pair * stack = (struct pair *)malloc(factor * cell_id * sizeof(struct pair)); //Stack of pairs of ids of nodes that have to be treated
      if (stack == NULL){
            fprintf(stderr, "Error : Cannot allocate memory for the stack in function viscoelastic_flattree.\n");
            abort();
      }
      
      int i = 0; //Index in the stack of the current pair of ids
      int j = 0; //Index of where to put a pair in the stack
      int a, b;  //Ids of the nodes of the current pair
      int p, q;  //Loop indexes
      int s, u;  //Bodies' indexes
      int * dots_a; //All the bodies in cell a
      int * dots_b; //All the bodies in cell b
      int Na, Nb; //Number of bodies in cells a and b
      int how_many_children_a, how_many_children_b; //Number of children of cells a and b
      int idFirstChild; //Id of first child
      int idLastChild;  //Id of last child
      typ * com_a, * com_b; //Nodes' average body position
      typ R[3];  //com_a - com_b
      typ r_crit_a, r_crit_b; // Critical radii of node a and b
      typ r;
      
      /******** Putting the pair (root_cell, root_cell) in the stack ********/
      (stack + j) -> fst = 0;
      (stack + j) -> snd = 0;
      j ++;
      
      
      /******** I travel the stack of pairs of nodes. At each pair, if NaNb < N_cc_pre, I treat it brute-forcely,  ********/
      /******** otherwise, if the nodes are well-separated, I do nothing, otherwise, if                            ********/
      /******** NaNb < N_cc_post or the pair has no children, I treat it brute-forcely, otherwise, I subdivise the ********/
      /******** largest node of the pair (or the only one that has children). If it is a pair of the same node, I  ********/
      /******** treat it brute-forcely if Na < N_cs or if it has no children, and else I subdivise it              ********/
      while (j > i){
            a = (stack + i) -> fst; //Id of first  node
            b = (stack + i) -> snd; //Id of second node
            Na = (FlatTree + a) -> how_many_dots;
            Nb = (FlatTree + b) -> how_many_dots;
            if (a == b){ //Cell self-interation
                  how_many_children_a = (FlatTree + a) -> how_many_children;
                  if (Na < N_cs_collision || how_many_children_a == 0){ //Direct interaction
                        dots_a = (FlatTree + a) -> dots;
                        for (p = 0; p < Na; p ++){
                              for (q = p + 1; q < Na; q ++){
                                    s = dots_a[p]; //Id of first  body
                                    u = dots_a[q]; //Id of second body
                                    if (generating){
                                          deOverlap(viscoelastic, s, u);
                                    }
                                    else if (connects(viscoelastic, s, u)){
                                          N_connections ++;
                                          add(s, &first);
                                          add(u, &second);
                                    }
                              }
                        }
                  }
                  else{ //Subdivision
                        idFirstChild = (FlatTree + a) -> idFirstChild;
                        idLastChild  = idFirstChild + how_many_children_a;
                        /******** Making sure that the stack is big enough. To be removed when the code is robust ********/
                        if (j + 36 >= factor * cell_id){
                              fprintf(stderr, "Error : The stack is not big enough in function viscoelastic_flattree. Try increasing the value of factor.\n");
                              abort();
                        }
                        for (p = idFirstChild; p < idLastChild; p++){
                              for (q = p; q < idLastChild; q++){                                   
                                    (stack + j) -> fst = p;
                                    (stack + j) -> snd = q;
                                    j++;
                              }
                        }
                  }
            }
            else{ //Interaction between two different cell
                  r_crit_a = (FlatTree + a) -> r_crit;
                  com_a    = (FlatTree + a) -> com;
                  r_crit_b = (FlatTree + b) -> r_crit;
                  com_b    = (FlatTree + b) -> com;
                  R[0]     = com_a[0] - com_b[0];
                  R[1]     = com_a[1] - com_b[1];
                  R[2]     = com_a[2] - com_b[2];
                  r        = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
                  if (r < r_crit_a + r_crit_b){ //The nodes are not well-separated
                        how_many_children_a = (FlatTree + a) -> how_many_children;
                        how_many_children_b = (FlatTree + b) -> how_many_children;
                        if ((!(Na > N_cc_collision && Nb > N_cc_collision) && Na*Nb < N_cc_collision) || (how_many_children_a == 0 && how_many_children_b == 0)){ //Direct interaction
                              dots_a = (FlatTree + a) -> dots;
                              dots_b = (FlatTree + b) -> dots;
                              for (p = 0; p < Na; p ++){
                                    for (q = 0; q < Nb; q ++){
                                          s = dots_a[p]; //Id of first  body
                                          u = dots_b[q]; //Id of second body
                                          if (generating){
                                                deOverlap(viscoelastic, s, u);
                                          }
                                          else if (connects(viscoelastic, s, u)){
                                                N_connections ++;
                                                add(s, &first);
                                                add(u, &second);
                                          }
                                    }
                              }
                        }
                        else{ //Subdivision
                              /******** Making sure that the stack is big enough. To be removed when the code is robust ********/
                              if (j + 8 >= factor * cell_id){
                                    fprintf(stderr, "Error : The stack is not big enough in function viscoelastic_flattree. Try increasing the value of factor.\n");
                                    abort();
                              }
                              if (how_many_children_a == 0){ //Subdivising b
                                    idFirstChild = (FlatTree + b) -> idFirstChild;
                                    idLastChild  = idFirstChild + how_many_children_b;
                                    for (p = idFirstChild; p < idLastChild; p++){                                       
                                          (stack + j) -> fst = a;
                                          (stack + j) -> snd = p;
                                          j++;
                                    }
                              }
                              else if (how_many_children_b == 0){ //Subdivising a
                                    idFirstChild = (FlatTree + a) -> idFirstChild;
                                    idLastChild  = idFirstChild + how_many_children_a;
                                    for (p = idFirstChild; p < idLastChild; p++){
                                          (stack + j) -> fst = b;
                                          (stack + j) -> snd = p;
                                          j++;
                                    }
                              }
                              else{ //Both cells have children, subdivising the largest (in term of critical radius)
                                    if (r_crit_a < r_crit_b){ //Subdivising b
                                          idFirstChild = (FlatTree + b) -> idFirstChild;
                                          idLastChild  = idFirstChild + how_many_children_b;
                                          for (p = idFirstChild; p < idLastChild; p++){
                                                (stack + j) -> fst = a;
                                                (stack + j) -> snd = p;
                                                j++;
                                          }
                                    }
                                    else{ //Subdivising a
                                          idFirstChild = (FlatTree + a) -> idFirstChild;
                                          idLastChild  = idFirstChild + how_many_children_a;
                                          for (p = idFirstChild; p < idLastChild; p++){
                                                (stack + j) -> fst = b;
                                                (stack + j) -> snd = p;
                                                j++;
                                          }
                                    }
                              }
                        }
                  }
            }
            i ++;
      }
      free(stack);
      stack = NULL;
}


int connects(struct moonlet * viscoelastic, int a, int b){

      /******** Returns 1 if the particles a and b are close enough to connect, and 0 else ********/

      typ connecting_distance2 = pow(3.0*connections_per_node*shapeV/(4.0*M_PI*(typ) N_0), 2.0/3.0);
      typ Xa, Ya, Za, Xb, Yb, Zb;
      
      /******** Getting the spring's rest length ********/ 
      Xa = (viscoelastic + a) -> x;
      Ya = (viscoelastic + a) -> y;
      Za = (viscoelastic + a) -> z;
      Xb = (viscoelastic + b) -> x;
      Yb = (viscoelastic + b) -> y;
      Zb = (viscoelastic + b) -> z;
      
      return ((Xa - Xb)*(Xa - Xb) + (Ya - Yb)*(Ya - Yb) + (Za - Zb)*(Za - Zb) < connecting_distance2 ? 1 : 0);
}


struct connection make_connection(struct moonlet * viscoelastic, int a, int b){

      /******** Initializes the connection structure between particles a and b ********/

      typ Xa, Ya, Za, Xb, Yb, Zb;
      struct connection C;
      C.Pair.fst = a;
      C.Pair.snd = b;
      
      /******** Getting the spring's rest length ********/ 
      Xa = (viscoelastic + a) -> x;
      Ya = (viscoelastic + a) -> y;
      Za = (viscoelastic + a) -> z;
      Xb = (viscoelastic + b) -> x;
      Yb = (viscoelastic + b) -> y;
      Zb = (viscoelastic + b) -> z;
      C.rest_length = sqrt((Xa - Xb)*(Xa - Xb) + (Ya - Yb)*(Ya - Yb) + (Za - Zb)*(Za - Zb));
      return C;
}


typ get_perturbing_true_anomaly(typ time){

      /******** Computes the true anomaly of the perturbing body from the differential equation  ********/
      /******** dnu/dt = sqrt(mu)*(1 + e cos nu)^2/(a(1-e^2))^(3/2) using a Runge-Kutta 4 method ********/
      /******** The periodicity is used to prevent any long-term error accumulation.             ********/
      /******** This is equivalent to solving Kepler's equation.                                 ********/

      int j, N_step;
      typ K1, K2, K3, K4;
      typ period, t, dt, n_step, partial_tra, previous_tra, num;
      typ denom = sqrt(pert_sma*(1.0 - pert_ecc*pert_ecc)); denom *= denom*denom;
      typ sq_mu = sqrt(G*(pert_mass + M_unit));

      /******** Establishing the integration time and the timestep to be used ********/
      period  = 2.0*M_PI/sq_mu*sqrt(pert_sma*pert_sma*pert_sma);
      t       = pert_ecc < 1. ? fmod(time, period) : time;
      n_step  = pert_ecc < 1. ? floor(256.0*t/period) + 1.0 : 256.0;
      dt      = t/n_step;
      n_step *= pert_ecc < 1. && pert_ecc > 0.8 ? 4. : 1.; //Decreasing the timestep in case of highly eccentric elliptic trajectory
      dt     /= pert_ecc < 1. && pert_ecc > 0.8 ? 4. : 1.;
      N_step  = (int) n_step;
      if (fabs(t - n_step*dt) > 1.0e-13){
            fprintf(stderr, "Error : Wrong computation of the integration time in function get_perturbing_true_anomaly");
            abort();
      }
      
      /******** Integrating ********/
      previous_tra = pert_tra;
      for (j = 0; j < N_step; j ++){
            partial_tra   = previous_tra;
            num           = 1.0 + pert_ecc * cos(partial_tra);
            K1            = sq_mu*num*num/denom;
            partial_tra   = previous_tra + 0.5*K1*dt;
            num           = 1.0 + pert_ecc * cos(partial_tra);
            K2            = sq_mu*num*num/denom;
            partial_tra   = previous_tra + 0.5*K2*dt;
            num           = 1.0 + pert_ecc * cos(partial_tra);
            K3            = sq_mu*num*num/denom;
            partial_tra   = previous_tra + K3*dt;
            num           = 1.0 + pert_ecc * cos(partial_tra);
            K4            = sq_mu*num*num/denom;
            previous_tra += dt*(K1 + 2.0*K2 + 2.0*K3 + K4)/6.0;
      }
      
      return previous_tra;
}


void quaternion_norm(struct quaternion * q){

      /******** Normalizes the quaternion q ********/

      typ w = q -> w;  typ x = q -> x;  typ y = q -> y;  typ z = q -> z;  
      typ u = sqrt(w*w + x*x + y*y + z*z);
      if (u == 0.){
            return;
      }
      q -> w /= u;
      q -> x /= u;
      q -> y /= u;
      q -> z /= u;
}


struct quaternion get_quaternion(typ ux, typ uy, typ uz, typ vx, typ vy, typ vz){

      /******** Returns a quaternion that describes the shortest arc rotation between vectors u and v                               ********/
      /******** https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another ********/
      /******** Does not manage the exception u x v = 0, in which case any vector orthogonal to u must be put in (q.x, q.y, q.z)    ********/

      typ udotv = ux*vx + uy*vy + uz*vz;
      typ u2    = ux*ux + uy*uy + uz*uz;
      typ v2    = vx*vx + vy*vy + vz*vz;
      struct quaternion q;
      q.w = udotv + sqrt(u2*v2);
      q.x = uy*vz - uz*vy; //u cross v
      q.y = uz*vx - ux*vz;
      q.z = ux*vy - uy*vx;
      quaternion_norm(&q);
      return q;
}


void rotate_with_quaternion(typ x, typ y, typ z, struct quaternion q, typ * xr, typ * yr, typ * zr){
   
      /******** Stores into (xr, yr, zr) the coordinates of the vector (x, y, z) once rotated by the quaternion q ********/
   
      *xr = (1.0 - 2.0*q.y*q.y - 2.0*q.z*q.z)*x +       (2.0*q.x*q.y - 2.0*q.z*q.w)*y +       (2.0*q.x*q.z + 2.0*q.y*q.w)*z;
      *yr =       (2.0*q.x*q.y + 2.0*q.z*q.w)*x + (1.0 - 2.0*q.x*q.x - 2.0*q.z*q.z)*y +       (2.0*q.y*q.z - 2.0*q.x*q.w)*z;
      *zr =       (2.0*q.x*q.z - 2.0*q.y*q.w)*x +       (2.0*q.y*q.z + 2.0*q.x*q.w)*y + (1.0 - 2.0*q.x*q.x - 2.0*q.y*q.y)*z;
}


void three_closest_nodes(struct moonlet * viscoelastic, int k, int * indexes){

      /******** Finds the three nodes closest to node k and stores their indexes.********/
      /******** This fonction is called after the main connection stage on nodes ********/
      /******** connected less than three times. Works brute-forcely             ********/

      typ Xk, Yk, Zk, dX, dY, dZ, d1, d2, d3, dist;
      int i1, i2, i3, j;

      d1 = 1.0e300;  d2 = 1.0e300;  d3 = 1.0e300;
      i1 = 0;        i2 = 0;        i3 = 0;
      Xk = (viscoelastic + k) -> x;
      Yk = (viscoelastic + k) -> y;
      Zk = (viscoelastic + k) -> z;

      for (j = 0; j < N_0; j ++){
            if (j != k){
                  dX   = (viscoelastic + j) -> x - Xk;
                  dY   = (viscoelastic + j) -> y - Yk;
                  dZ   = (viscoelastic + j) -> z - Zk;
                  dist = dX*dX + dY*dY + dZ*dZ;
                  if (dist < d1){
                        i3 = i2; d3 = d2;
                        i2 = i1; d2 = d1;
                        i1 = j;  d1 = dist;
                  }
                  else if (dist < d2){
                        i3 = i2; d3 = d2;
                        i2 = j;  d2 = dist;
                  }
                  else if (dist < d3){
                        i3 = j;  d3 = dist;
                  }
            }
      }
      * indexes      = i1;
      *(indexes + 1) = i2;
      *(indexes + 2) = i3;
}
