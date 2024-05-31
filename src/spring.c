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


struct moonlet * generate_visco_elastic_body(){

      /******** Similar to the function populate in structure.c, but populates    ********/
      /******** the simulation with the particles of a visco-elastic body instead ********/

      int j;
      typ maxD = 0.;
      typ costheta, costheta_max;
      typ distance, distance2;
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
            /******** Generating the array of vertices ********/
            typ * vertices = (typ *)malloc(3*n_vertices*sizeof(typ));
            if (vertices == NULL){
                  fprintf(stderr, "Error : Cannot allocate array of vertices in function generate_visco_elastic_body.\n");
                  abort();
            }
      
            /******** Retrieving all the vertices from the file pth/shape_model.txt ********/
            char fileOfVertices[800]; 
            strcpy(fileOfVertices, pth);
            strcat(fileOfVertices, "shape_model.txt");
            readFromFile(fileOfVertices, vertices, 3*n_vertices);
      
            /******** Getting the largest distance between a vertice and the origin ********/
            for (j = 0; j < n_vertices; j ++){
                  distance = vertices[3*j]*vertices[3*j] + vertices[3*j + 1]*vertices[3*j + 1] + vertices[3*j + 2]*vertices[3*j + 2];
                  if (distance > maxD){
                        maxD = distance;
                  }
            }
            maxD = sqrt(maxD);
      
            /******** Drawing the points of the viscoelastic body ********/
            while (currentN < N_0){
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
                        (viscoelastic + currentN) -> radius = minimal_distance/2.0;
                        currentN ++;
                  }
            }
            free(vertices);  vertices = NULL;
      
            /******** Making sure that particles do not overlap ********/
            overlap(viscoelastic);
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
      
      if (random_initial_bool){
            /******** Getting the total number of connections and the ids of the connecting particles ********/
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
                        (viscoelastic + j) -> radius = connecting_distance/2.0;
                  }
                  root     = root_cell(viscoelastic);
                  FlatTree = flattree_init(root);
                  clear_boxdot(&root);
                  center_and_maxR_flattree(FlatTree, viscoelastic);
                  rmax_and_rcrit_flattree (FlatTree, viscoelastic);                   
                  viscoelastic_flattree   (FlatTree, viscoelastic, 0);   
            }
            if(falcON_bool){ //Resetting the octree and the associated data
                  for (j = 0; j < N_0; j ++){
                        (viscoelastic + j) -> radius = minimal_distance/2.0;
                  }
                  for (j = 0; j < cell_id; j ++){
                        free((FlatTree + j) -> dots);
                        (FlatTree + j) -> dots = NULL;
                  }
                  free(FlatTree);
                  FlatTree       = NULL;
                  how_many_cells = 0;
                  cell_id        = 0;
            }
      
            /******** Generating the array of connections ********/
            index       = 0;
            index1      = first  -> how_many - 1;
            index2      = second -> how_many - 1;
            connections = (struct connection *)malloc(N_connections*sizeof(struct connection));
            if (connections == NULL){
                  fprintf(stderr, "Cannot allocate array of %d connections in function generate_connections\n", N_connections);
                  abort();
            }
            while (first != NULL){
                  i           = (first  -> ids)[index1];
                  j           = (second -> ids)[index2];
                  *(connections + index) = make_connection(viscoelastic, i, j);
                  index1 --;  index2 --; //Switching to the next pair of connecting particles
                  if (index1 < 0){
                        first  = first  -> queue;
                        index1 = max_ids_per_node - 1;
                  }
                  if (index2 < 0){
                        second = second -> queue;
                        index2 = max_ids_per_node - 1;
                  }
                  index ++;
            }
            clear_chain(&first);
            clear_chain(&second);
      }
      else{
            typ Xi, Yi, Zi, Xj, Yj, Zj, rest_length, eq_length;
            char fileOfIC[800]; 
            strcpy(fileOfIC, pth);
            strcat(fileOfIC, "connections.txt");
            typ * IC = (typ *)malloc(150*N_0*sizeof(typ));
            if (IC == NULL){
                  fprintf(stderr, "Error : Cannot allocate array for initial conditions in function generate_visco_elastic_body.\n");
                  abort();
            }
            N_connections = readFromFile_withoutConstraint(fileOfIC, IC, 150*N_0);
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
                  Xi          = (viscoelastic + i) -> x;
                  Yi          = (viscoelastic + i) -> y;
                  Zi          = (viscoelastic + i) -> z;
                  Xj          = (viscoelastic + j) -> x;
                  Yj          = (viscoelastic + j) -> y;
                  Zj          = (viscoelastic + j) -> z;
                  eq_length   = sqrt((Xi - Xj)*(Xi - Xj) + (Yi - Yj)*(Yi - Yj) + (Zi - Zj)*(Zi - Zj));
                  (connections + k) -> Pair.fst           = i;
                  (connections + k) -> Pair.snd           = j;
                  (connections + k) -> rest_length        = rest_length;
                  (connections + k) -> equilibrium_length = eq_length;
            }
            free(IC);
            IC = NULL;
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
            if (count > 512){
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


      typ Xa, Ya, Za, Xb, Yb, Zb;
      
      /******** Getting the spring's rest length ********/ 
      Xa = (viscoelastic + a) -> x;
      Ya = (viscoelastic + a) -> y;
      Za = (viscoelastic + a) -> z;
      Xb = (viscoelastic + b) -> x;
      Yb = (viscoelastic + b) -> y;
      Zb = (viscoelastic + b) -> z;
      
      return ((Xa - Xb)*(Xa - Xb) + (Ya - Yb)*(Ya - Yb) + (Za - Zb)*(Za - Zb) < connecting_distance*connecting_distance ? 1 : 0);
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
      
      C.rest_length        = sqrt((Xa - Xb)*(Xa - Xb) + (Ya - Yb)*(Ya - Yb) + (Za - Zb)*(Za - Zb));
      C.equilibrium_length = 0.; //Will be properly initialized later
      return C;
}



