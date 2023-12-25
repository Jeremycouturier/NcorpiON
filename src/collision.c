/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    collision.c                                                 ********/
/******** @brief   This file manages collision detection                       ********/
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
#include "collision.h"
#include "parameters.h"
#include "structure.h"
#include "physics.h"
#include <errno.h>
#include <math.h>


/*************************************************************************************/
/******** I first implement a mesh-based algorithm to find colliding moonlets ********/
/*************************************************************************************/


typ * closest_approach(struct moonlet * moonlets, int a, int b){

      /******** A much better version than the old version. No more logarithmic search and no more case   ********/
      /******** disjonction according to whether or not the closest approach occurs during or at the end  ********/
      /******** of the timestep. In this version, the time t after the beginning of the timestep when the ********/
      /******** collision occurs is entirely determined by the knowledge of the vectors dr = r_a-r_b and  ********/
      /******** dv = v_a-v_b at the beginning of the timestep. See the PDF draft for details.             ********/
      /******** As in the old version, I return a NULL pointer if no collision occurs during the timestep ********/
      /******** and I approximate the trajectories by straight lines. Fills the array approach with the   ********/
      /******** 6 positions of the two collinding moonlets.                                               ********/
      

      
      /******** One of the moonlet might not exist if it previously collided during this timestep. I check that ********/
      if (!(*(exists + a) && *(exists + b))){
            return NULL;
      }
      
      typ xa, ya, za, xb, yb, zb;                            //Cartesian positions of the moonlets at the beginning of the timestep
      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;                //Cartesian speeds    of the moonlets
      typ dx, dy, dz, dvx, dvy, dvz;                         //dx is xa-xb, and so on.
      typ * to_be_returned = NULL;                           //The array to be returned
      typ R_a = (moonlets+a) -> radius;                      //The moonlets' radii
      typ R_b = (moonlets+b) -> radius;
      
      
      /******** Getting the speeds ********/
      vx_a = (moonlets+a) -> vx;
      vy_a = (moonlets+a) -> vy;
      vz_a = (moonlets+a) -> vz;
      vx_b = (moonlets+b) -> vx;
      vy_b = (moonlets+b) -> vy;
      vz_b = (moonlets+b) -> vz;
      
      
      /******** Getting the positions at the beginning of the timestep ********/
      xa = (moonlets+a) -> x;
      ya = (moonlets+a) -> y;
      za = (moonlets+a) -> z;
      xb = (moonlets+b) -> x;
      yb = (moonlets+b) -> y;
      zb = (moonlets+b) -> z;
      
      /******** Getting the relative positions and velocities ********/
      dx  = xa-xb;
      dy  = ya-yb;
      dz  = za-zb;
      dvx = vx_a-vx_b;
      dvy = vy_a-vy_b;
      dvz = vz_a-vz_b;
      
      
      /******** If the moonlets are getting farther away from each other, no collision will occur. ********/
      /******** This is the case if, and only if, dr.dv > 0                                        ********/
      typ dr_dot_dv = dx*dvx+dy*dvy+dz*dvz;
      if (dr_dot_dv >= 0.0){
            return NULL;
      }
      
      
      /******** If that point is reached, then the moonlets are getting closer to each other at the ********/
      /******** beginning of the timestep and a collision might occur. A collision will occur if,   ********/
      /******** and only if, the discriminant (dr.dv)^2+dv^2(R^2-dr^2) is positive. If that is not  ********/
      /******** the case, then I return a NULL pointer. If that is the case, but the collision      ********/
      /******** will happen after a time larger than the timestep, then I still return a NULL       ********/
      /******** pointer. Otherwise, I return the array [xa, ya, za, xb, yb, zb] at the collision.   ********/
      
      typ dr2 = dx*dx+dy*dy+dz*dz; //Square of the distance between the moonlets at the beginning of the timestep
      typ dv2 = dvx*dvx+dvy*dvy+dvz*dvz; //Square of the relative velocity between the moonlets
      typ R=R_a+R_b; //Sum of moonlets radii
      typ discriminant = dr_dot_dv*dr_dot_dv+dv2*(R*R-dr2);
      
      if (discriminant <= 0.0){ //Then no collision will ever occur because the minimal distance is larger than R
            return NULL;
      }
      
      /******** A collision will occur, but it might be after a time larger than the timestep ********/
      time_until_collision = -(dr_dot_dv+sqrt(discriminant))/dv2; // The time of collision after the beginning of the timestep
      
      if (time_until_collision > timestep){ //The collision won't occur during this timestep
            return NULL;
      }
      
      /******** The collision will occur during this timestep    ********/
      /******** Actualizing the positions to the collision point ********/
      xa += time_until_collision*vx_a;
      ya += time_until_collision*vy_a;
      za += time_until_collision*vz_a;
      xb += time_until_collision*vx_b;
      yb += time_until_collision*vy_b;
      zb += time_until_collision*vz_b;
            
      /******** Filling the array approach ********/                    
      *(approach+0) = xa;
      *(approach+1) = ya;
      *(approach+2) = za;
      *(approach+3) = xb;
      *(approach+4) = yb;
      *(approach+5) = zb;
      return approach;            
          
}      


void hash_table_cell_index(struct moonlet * moonlets, int a){

      /******** Actualizes the array indexes=[i_0,i_1,...,i_26] where i_0 to a_26 are the ********/
      /******** 27 indexes of the cells of the neighbourhood of a in the hash table.      ********/
      /******** If moonlet a is outside the collision cube, then *indexes is set to -1.   ********/
      /******** Since one extra layer was allocated to the collision cube, there are no   ********/
      /******** edge effects to be taken care of.                                         ********/
      

      
      /******** Verifying the existence of the moonlet. To be removed if useless. ********/
      if (*(exists+a)==0){
            fprintf(stderr, "Error : The moonlet doesn't exist.\n");
            abort();
      }
      
      /******** Getting the moonlet's position ********/
      typ X,Y,Z;
      X = (moonlets + a) -> x;
      Y = (moonlets + a) -> y;
      Z = (moonlets + a) -> z;
      
      /******** If the moonlet is outside the collision cube, then *indexes is set to -1 ********/
      if (absolute(X) > collision_cube/2.0 || absolute(Y) > collision_cube/2.0 || absolute(Z) > collision_cube/2.0){
            *indexes = -1;
      }
      else{ //The moonlet is inside the collision cube
            int i,j,k;
            int d  = collision_cube_cells;
            int d2 = d*d;
            i = d/2+1+((int) (integral(X/gam))); //Indexes 0 and collision_cube_cells+1 correspond to the unused outermost layer of the collision cube
            j = d/2+1+((int) (integral(Y/gam))); //That means that the smallest (resp. largest) index for a cell of the collision cube is 1 (resp. collision_cube_cells)
            k = d/2+1+((int) (integral(Z/gam))); //Then, that particular cell is accessed by *(hash+index), where index=i+j*d+k*d^2 and d = collision_cube_cells
            
            
            /******** I compute the indexes by increasing order, so that I subsequently travel along the neighbourhoods in a (supposedly) cache-friendly manner ********/
            
            
            /******** The 9 cells on the bottom of the neighbourhood                                                                          ********/
            /******** Bottom, top, front and back are defined assuming that I are looking at the system from position (x=0, y=-infinity, z=0) ********/
            
            /******** The 3 cells at the front ********/
            *indexes      = (i-1)+(j-1)*d+(k-1)*d2;
            *(indexes+1)  = (i)+(j-1)*d+(k-1)*d2;
            *(indexes+2)  = (i+1)+(j-1)*d+(k-1)*d2;
            /******** The 3 cells in the middle ********/
            *(indexes+3)  = (i-1)+(j)*d+(k-1)*d2;
            *(indexes+4)  = (i)+(j)*d+(k-1)*d2;
            *(indexes+5)  = (i+1)+(j)*d+(k-1)*d2;
            /******** The 3 cells at the back ********/
            *(indexes+6)  = (i-1)+(j+1)*d+(k-1)*d2;
            *(indexes+7)  = (i)+(j+1)*d+(k-1)*d2;
            *(indexes+8)  = (i+1)+(j+1)*d+(k-1)*d2;
            
            
            /******** The 9 cells with same z coordinate as the central cell ********/
            
            /******** The 3 cells at the front ********/
            *(indexes+9)  = (i-1)+(j-1)*d+(k)*d2;
            *(indexes+10) = (i)+(j-1)*d+(k)*d2;
            *(indexes+11) = (i+1)+(j-1)*d+(k)*d2;
            /******** The 3 cells in the middle ********/
            *(indexes+12) = (i-1)+(j)*d+(k)*d2;
            *(indexes+13) = i+j*d+k*d2;
            *(indexes+14) = (i+1)+(j)*d+(k)*d2;
            /******** The 3 cells at the back ********/
            *(indexes+15) = (i-1)+(j+1)*d+(k)*d2;
            *(indexes+16) = (i)+(j+1)*d+(k)*d2;
            *(indexes+17) = (i+1)+(j+1)*d+(k)*d2;

            
            /******** The 9 cells on the top of the neighbourhood ********/
            
            /******** The 3 cells at the front ********/
            *(indexes+18) = (i-1)+(j-1)*d+(k+1)*d2;
            *(indexes+19) = (i)+(j-1)*d+(k+1)*d2;
            *(indexes+20) = (i+1)+(j-1)*d+(k+1)*d2;
            /******** The 3 cells in the middle ********/
            *(indexes+21) = (i-1)+(j)*d+(k+1)*d2;
            *(indexes+22) = (i)+(j)*d+(k+1)*d2;
            *(indexes+23) = (i+1)+(j)*d+(k+1)*d2;
            /******** The 3 cells at the back ********/
            *(indexes+24) = (i-1)+(j+1)*d+(k+1)*d2;
            *(indexes+25) = (i)+(j+1)*d+(k+1)*d2;
            *(indexes+26) = (i+1)+(j+1)*d+(k+1)*d2;
      }
}


void neighbours(struct moonlet * moonlets, int a){

      /******** Fills the chain nghb with the ids of the neighbours of the moonlet a and puts a in the hash table ********/
      /******** If a is outside the collision cube, then nghb is unchanged. If a is in the collision cube but has ********/
      /******** no neighbours, then nghb is also unchanged.                                                       ********/
      
      struct chain * ch = NULL; //The chain of neighbours in a particular cell of the neighbourhood of a
      
      hash_table_cell_index(moonlets, a); //Retrieving the indexes of the cells of the neighbourhood of a in the hash table, and storing them in the array "indexes"
      
      if (*indexes==-1){ //a is outside the collision cube
            return;
      }
      
      /******** a is inside the collision cube. I check its neighbourhood. ********/
      int index;
      int l;
      int node_size;
      int idd;
      
      
      for (l=0; l<27; l++){ // I go over the whole neighbourhood of a
            index = *(indexes+l);    //The index of the neighbouring cell
            ch    = *(hash+index);   //The chain of ids of moonlets contained in this cell
            if (ch!=NULL) { //If there is a moonlet in that cell
            node_size = ch -> how_many;
                  while (ch!=NULL){   //I travel along the chain of moonlets' ids contained in this cell
                        idd = (ch -> ids)[node_size - 1];
                        if (*(exists + idd)){  //idd might not exist anymore if it previously collided during that timestep
                              add(idd, &nghb); //I add the ids of moonlets contained in that cell to the chain nghb
                        }
                        node_size --;
                        if (node_size == 0){
                              ch        = ch -> queue;
                              node_size = max_ids_per_node;
                        }
                  }
            }
      }
      
      /******** I actualize the array modified_cells and the variable how_many_modified to keep track of the modifications to the hash table ********/
      index = *(indexes+13);
      ch    = *(hash+index);
      
      if (ch == NULL){ //If not, then that cell of the hash table was already modified previously in that timestep
            *(modified_cells + how_many_modified) = index;
            how_many_modified ++;
      }
      
      /******** I add a to the hash table ********/
      add(a, hash+index);

}




/******** Implementing the O(N) mesh algorithm for collision detection ********/


void mesh(struct moonlet * moonlets){


      /******** I go over all the moonlets once. For each moonlet, I check whether its radius is smaller  ********/
      /******** or larger than the mesh-size. If it is larger, I look for collisions between that moonlet ********/
      /******** and any other moonlet. If it is smaller, I add the moonlet to the chain of the            ********/
      /******** corresponding cell in the hash table, and I look for collisions between that moonlet and  ********/
      /******** any other moonlet currently in its neighbourhood.                                         ********/
      
      
      int k,p;
      typ * the_approach;
      typ R;
      int current_largest_id = largest_id;
      
      for (k = 0; k <= current_largest_id; k++){
            if(*(exists+k)){ //Checking if there is a moonlet in the k^th cell of the array moonlets
                  R = (moonlets + k) -> radius;

                  /******** If the moonlet is smaller than the mesh-size ********/
                  if (R<=gam){
                        
                        how_many_small+=1;
                        neighbours(moonlets, k); //Adding the moonlet to the hash table, and retrieving its neighbours in the chain nghb
                        
                        while(nghb -> how_many > 0){ //If the moonlet is inside the collision cube and has at least one neighbour
                              p = (nghb -> ids)[nghb -> how_many - 1]; // p is the id of a moonlet neighbour to k
                              the_approach = closest_approach(moonlets, p, k);
                              if (the_approach == NULL || *(did_collide+k) || *(did_collide+p)) { //No collision or one of the moonlet already had a collision during that timestep
                                    if (mutual_bool){ //I register the pair (k,p) to be taken care of for mutual gravitational interactions
                                                      //In case of collision, this is done when treating the collision
                                          (pairs + how_many_pairs) -> fst = k;
                                          (pairs + how_many_pairs) -> snd = p;
                                          how_many_pairs ++;
                                    }
                              }
                              else { //The closest approach leads to a collision and neither k nor p previously had a collision during that timestep
                                    if(elastic_collision_bool){ //All collision are elastic
                                          collision_treatment(moonlets, p, k, 0);
                                    }
                                    else if(inelastic_collision_bool){ //All collision are inelastic
                                          collision_treatment(moonlets, p, k, 1);
                                    }
                                    else if (instant_merger_bool){ //All collisions result in an instant merger
                                          collision_treatment(moonlets, p, k, 2);
                                    }
                                    else if (fragmentation_bool){ //Collision are treated with the fragmentation model described in the PDF draft
                                          collision_treatment(moonlets, p, k, 3);
                                    }
                              }
                               
                              nghb = partial_delete(nghb);
                              total_neighbours ++;
                        }
                  }
                  
                  /******** If the moonlet is larger than the mesh-size ********/
                  else{
                        how_many_big+=1;
                        typ Rp; //Radius of the p^th moonlet
                        for (p=0; p<=largest_id; p++){
                              if(*(exists+p) && p != k){ //If the p^th moonlet exists and k is different from p
                                    Rp = (moonlets+p)->radius;
                                    if (Rp < gam || (Rp >= gam && p > k)){ //If the p^th moonlet is small or if it is big and p > k, so the pair (k,p) is not counted twice
                                          the_approach = closest_approach(moonlets, p, k);
                                          if (the_approach == NULL || *(did_collide+k) || *(did_collide+p)) { //No collision or one of the moonlet already collided during that timestep
                                                if (mutual_bool){ //I register the pair (k,p) to be taken care of for mutual gravitational interactions
                                                                  //In case of collision, this is done when treating the collision
                                                      (pairs+how_many_pairs)->fst = k;
                                                      (pairs+how_many_pairs)->snd = p;
                                                      how_many_pairs ++;
                                                }
                                          }
                                          else { //The closest approach leads to a collision and neither k nor p previously had a collision during that timestep
                                                if(elastic_collision_bool){ //All collision are elastic
                                                      collision_treatment(moonlets, p, k, 0);
                                                }
                                                else if(inelastic_collision_bool){ //All collision are inelastic
                                                      collision_treatment(moonlets, p, k, 1);
                                                }
                                                else if (instant_merger_bool){ //All collisions result in an instant merger
                                                      collision_treatment(moonlets, p, k, 2);
                                                }
                                                else if (fragmentation_bool){ //Collision are treated with the fragmentation model described in the PDF draft
                                                      collision_treatment(moonlets, p, k, 3);
                                                }
                                          }
                                    }
                              }
                        }
                  }    
            }
      }     
}


/*********************************************************************************************/
/******** Implementing the brute-force O(N^2) algorithm for close encounter detection ********/
/*********************************************************************************************/


void brute_force(struct moonlet * moonlets){


      /******** I look for collisions between all pairs of moonlets ********/
      int i,j;
      typ * the_approach;
      int current_largest_id = largest_id;

      for (i=0; i<=current_largest_id; i++){
            if(*(exists+i)){                    //Checking if there is a moonlet in the i^th cell of the array moonlets
                  for (j=i+1; j<=largest_id; j++){
                        if(*(exists+j)){        //Checking if there is a moonlet in the j^th cell of the array moonlets
                              the_approach = closest_approach(moonlets, i, j);
                              if (the_approach != NULL && !(*(did_collide+i)) && !(*(did_collide+j))){ //The closest approach leads to a collision and neither i nor j previously collided
                                    if(elastic_collision_bool){ //All collision are elastic
                                          collision_treatment(moonlets, i, j, 0);
                                    }
                                    else if(inelastic_collision_bool){ //All collision are inelastic
                                          collision_treatment(moonlets, i, j, 1);
                                    }
                                    else if (instant_merger_bool){ //All collisions result in an instant merger
                                          collision_treatment(moonlets, i, j, 2);
                                    }
                                    else if (fragmentation_bool){ //Collision are treated with the fragmentation model described in the PDF draft
                                          collision_treatment(moonlets, i, j, 3);
                                    }
                              }
                        }
                  }
            }
      }
}


/**************************************************************/
/******** Implementing collision detection with falcON ********/
/**************************************************************/


void get_center_and_maxR(struct node * FlatTree, struct moonlet * moonlets, int a){

      /******** Computes the average position of the moonlets and the largest moonlet radius of node a of FlatTree, ********/
      /******** assuming that it has no children. Initializes the fields com and M0 of (FlatTree + a) accordingly   ********/
      
      /******** Checking that the node has indeed no children. To be removed when the code is robust ********/
      if ((FlatTree + a) -> idFirstChild != -1){
            fprintf(stderr, "Error : Node %d has children in function get_center_and_maxR.\n", a);
            abort();
      }
      
      typ com[3] = {0.0, 0.0, 0.0}; //The average position of the moonlets in the cell
      typ maxR = 0.0;               //Largest moonlet radius
      
      typ x, y, z, vx, vy, vz, v, R, crit; //The current moonlet's coordinates
      
      int * dots = (FlatTree + a) -> dots; //All the moonlets in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int j; //Id of the current moonlet
      
      for (i = 0; i < how_many_dots; i++){
            j  = dots[i];
            R  = (moonlets + j) -> radius;
            x  = (moonlets + j) -> x;
            y  = (moonlets + j) -> y;
            z  = (moonlets + j) -> z;
            vx = (moonlets + j) -> vx;
            vy = (moonlets + j) -> vy;
            vz = (moonlets + j) -> vz;
            v  = sqrt(vx*vx + vy*vy +vz*vz);
            crit = R + v*timestep;
            
            if (crit > maxR){
                  maxR = crit;
            }
            com[0] += x;
            com[1] += y;
            com[2] += z;
      }
      
      /******** Initializing the relevant fields ********/
      (FlatTree + a) -> M0 = maxR; //When treating collision, the field M0 of the FlatTree gives the largest moonlet radius, not the cell's mass
      ((FlatTree + a) -> com)[0] = com[0]/((typ) how_many_dots); //Similarly, the field com gives the average position of the moonlets, not the center of mass
      ((FlatTree + a) -> com)[1] = com[1]/((typ) how_many_dots);
      ((FlatTree + a) -> com)[2] = com[2]/((typ) how_many_dots);
}


void get_rmax_and_rcrit(struct node * FlatTree, struct moonlet * moonlets, int a){

      /******** Initializes the fields rmax and rcrit of node a ********/
      /******** of FlatTree, assuming that it has no children.  ********/

      /******** Checking that the node has indeed no children. To be removed when the code is robust ********/
      if ((FlatTree + a) -> idFirstChild != -1){
            fprintf(stderr, "Error : Node %d has children in function get_rmax_and_rcrit.\n", a);
            abort();
      }
      
      typ rmax = 0.0;
      typ distance;
      
      /******** Retrieving the center of mass of the cell ********/
      typ com_x, com_y, com_z;
      com_x = ((FlatTree +a) -> com)[0];  com_y = ((FlatTree +a) -> com)[1];  com_z = ((FlatTree +a) -> com)[2];
      
      typ dx, dy, dz; //Distance with a dot along each axis
      
      int * dots = (FlatTree + a) -> dots; //All the moonlets in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int j; //Id of the moonlet whose distance to the center of mass is to be computed
      
      for (i = 0; i < how_many_dots; i++){
            j = dots[i];
            dx = (moonlets + j) -> x - com_x;
            dy = (moonlets + j) -> y - com_y;
            dz = (moonlets + j) -> z - com_z;
            distance = sqrt(dx*dx + dy*dy + dz*dz);
            if (distance > rmax){
                  rmax = distance;            
            }
      }
      
      /******** Initializing the relevant field ********/
      (FlatTree + a) -> r_max  = rmax;
      (FlatTree + a) -> r_crit = rmax + (FlatTree + a) -> M0;
}


void get_center_and_maxR_from_children(struct node * FlatTree, int a){

      /******** Computes the average position of the moonlets and the largest moonlet radius of node a of FlatTree ********/
      /******** from the children. Initializes the fields com and M0 of (FlatTree + a) accordingly ********/

      int idFirstChild = (FlatTree + a) -> idFirstChild;
      int idLastChild  = idFirstChild + (FlatTree + a) -> how_many_children;
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      
      
      typ com[3] = {0.0, 0.0, 0.0}; //The average position of the moonlets in the cell
      typ maxR = 0.0;               //Largest moonlet radius
      typ * com_child;              //The average position of the moonlets in the child
      typ maxR_child;               //Largest moonlet radius of the child
      int how_many_dots_in_child;   //Number of moonlets in the child
      
      typ x, y, z; //The current moonlet's coordinates
      typ R;       //The current moonlet's radius
      
      for (int i = idFirstChild; i < idLastChild; i++){
            maxR_child = (FlatTree + i) -> M0;
            com_child  = (FlatTree + i) -> com;
            how_many_dots_in_child = (FlatTree + i) -> how_many_dots;
            
            if (maxR_child > maxR){
                  maxR = maxR_child;
            }
            com[0] += (typ) how_many_dots_in_child * com_child[0];
            com[1] += (typ) how_many_dots_in_child * com_child[1];
            com[2] += (typ) how_many_dots_in_child * com_child[2];
            
      }

      /******** Initializing the relevant fields ********/
      (FlatTree + a) -> M0 = maxR; //When treating collision, the field M0 of the FlatTree gives the largest moonlet radius, not the cell's mass
      ((FlatTree + a) -> com)[0] = com[0]/((typ) how_many_dots); //Similarly, the field com gives the average position of the moonlets, not the center of mass
      ((FlatTree + a) -> com)[1] = com[1]/((typ) how_many_dots);
      ((FlatTree + a) -> com)[2] = com[2]/((typ) how_many_dots);
}


void get_rmax_and_rcrit_from_children(struct node * FlatTree, int a){

      /******** Initializes the fields rmax and rcrit of node a ********/
      /******** of FlatTree from its children.                  ********/

      int idFirstChild = (FlatTree + a) -> idFirstChild;
      int idLastChild  = idFirstChild + (FlatTree + a) -> how_many_children;
      int i;

      typ * center = (FlatTree + a) -> center; //The center of the node
      typ * com = (FlatTree + a) -> com;       //The average position of the moonlets in the node
      typ D = (FlatTree + a) -> sidelength;    //The sidelength of the node
      typ corner[8][3] = {{center[0] - D - com[0], center[1] - D - com[1], center[2] - D - com[2]},
                          {center[0] - D - com[0], center[1] - D - com[1], center[2] + D - com[2]},
                          {center[0] - D - com[0], center[1] + D - com[1], center[2] - D - com[2]},
                          {center[0] - D - com[0], center[1] + D - com[1], center[2] + D - com[2]},
                          {center[0] + D - com[0], center[1] - D - com[1], center[2] - D - com[2]},
                          {center[0] + D - com[0], center[1] - D - com[1], center[2] + D - com[2]},
                          {center[0] + D - com[0], center[1] + D - com[1], center[2] - D - com[2]},
                          {center[0] + D - com[0], center[1] + D - com[1], center[2] + D - com[2]}}; //corner - average position of the moonlets
      typ distance_to_farthest_corner = 0.0;
      typ distance;

      /******** Computing the distance between the center of mass and the most distant corner ********/
      for (i = 0; i < 8; i++){
            distance = sqrt(corner[i][0]*corner[i][0] + corner[i][1]*corner[i][1] + corner[i][2]*corner[i][2]);
            if (distance > distance_to_farthest_corner){
                  distance_to_farthest_corner = distance;
            }
      }
      
      typ * com_child; //Center of mass of a child
      typ rmax_child;  //Maximal radius of a child
      typ com_difference[3];
      typ max_ri = 0.0;

      /******** Eq. (9) of Dehnen (2002) ********/
      for (int i = idFirstChild; i < idLastChild; i++){
            com_child = (FlatTree + i) -> com;
            rmax_child = (FlatTree + i) -> r_max;
            com_difference[0] = com[0] - com_child[0];
            com_difference[1] = com[1] - com_child[1];
            com_difference[2] = com[2] - com_child[2];
            distance = rmax_child + sqrt(com_difference[0]*com_difference[0] + com_difference[1]*com_difference[1] + com_difference[2]*com_difference[2]);
            if (distance > max_ri){
                  max_ri = distance;
            }
      }
      
      /******** Initializing the relevant field ********/
      if (max_ri < distance_to_farthest_corner){
            (FlatTree + a) -> r_max  = max_ri;
            (FlatTree + a) -> r_crit = max_ri + (FlatTree + a) -> M0;
      }
      else{
            (FlatTree + a) -> r_max  = distance_to_farthest_corner;
            (FlatTree + a) -> r_crit = distance_to_farthest_corner + (FlatTree + a) -> M0;
      }
}


void center_and_maxR_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Computes the average moonlet position and the largest moonlet radius of all the nodes of FlatTree ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that could not be treated due to their child not having been treated yet
      int i; //Id of the current node
      int j = 0; //Index of where to put a node in the stack
      int how_many_children;
      
      /******** I travel the flattree. If a node is childless, I compute the average moonlet position and the largest moonlet radius directly ********/
      /******** Otherwise, I store it in the stack for future treatment                                                                       ********/
      for (i = 0; i < cell_id; i++){
            how_many_children = (FlatTree + i) -> how_many_children;
            if (how_many_children == 0){ //If the node has no children
                  get_center_and_maxR(FlatTree, moonlets, i);
            }
            else{ //If the node has children
                  stack[j] = i;
                  j++;
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function center_and_maxR_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** I now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j--;
            i = stack[j];
            get_center_and_maxR_from_children(FlatTree, i);
      }
      
      free(stack);
      stack = NULL;
}


void rmax_and_rcrit_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Computes and initializes the radii of convergence of all the nodes of FlatTree ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that could not be treated due to their child not having been treated yet
      int i; //Id of the current node
      int j = 0; //Index of where to put a node in the stack
      int how_many_children;
      
      /******** I travel the flattree. If a node is childless, I compute its convergence radius ********/
      /******** Otherwise, I store it in the stack for future treatment                         ********/
      for (i = 0; i < cell_id; i++){
            how_many_children = (FlatTree + i) -> how_many_children;
            if (how_many_children == 0){ //If the node has no children
                  get_rmax_and_rcrit(FlatTree, moonlets, i);
            }
            else{ //If the node has children
                  stack[j] = i;
                  j++;
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function rmax_and_rcrit_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** I now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j--;
            i = stack[j];
            get_rmax_and_rcrit_from_children(FlatTree, i);
      }
      
      free(stack);
      stack = NULL;
}


void collision_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Tree walk for collision detection ********/
      
      /******** It is hard to tell in advance how many pairs will be treated by the tree walk   ********/
      /******** I expect it will be at most factor * cell_id, but that might have to be changed ********/
      int factor = (int) integral(250.0 * 0.5 / theta_min);
      struct pair * stack = (struct pair *)malloc(factor * cell_id * sizeof(struct pair)); //Stack of pairs of ids of nodes that have to be treated
      if (stack == NULL){
            fprintf(stderr, "Error : Cannot allocate memory for the stack in function collision_flattree.\n");
            abort();
      }
      
      int i = 0; //Index in the stack of the current pair of ids
      int j = 0; //Index of where to put a pair in the stack
      int a, b;  //Ids of the nodes of the current pair
      int p, q;  //Loop indexes
      int s, u;  //Moonlet indexes
      int * dots_a; //All the moonlets in cell a
      int * dots_b; //All the moonlets in cell b
      int Na, Nb; //Number of moonlets in cells a and b
      int how_many_children_a, how_many_children_b; //Number of children of cells a and b
      int idFirstChild; //Id of first child
      int idLastChild;  //Id of last child
      typ * com_a, * com_b; //Nodes' average moonlet position
      typ R[3];  //com_a - com_b
      typ r_crit_a, r_crit_b; // Critical radii of node a and b
      typ * the_approach;
      typ r;
      
      /******** Putting the pair (root_cell, root_cell) in the stack ********/
      (stack + j) -> fst = 0;
      (stack + j) -> snd = 0;
      j++;
      
      
      /******** I travel the stack of pairs of nodes. At each pair, if NaNb < N_cc_pre, I treat it brute-forcely,  ********/
      /******** otherwise, if the nodes are well-separated, I do nothing, otherwise, if                            ********/
      /******** NaNb < N_cc_post or the pair has no children, I treat it brute-forcely, otherwise, I subdivise the ********/
      /******** largest node of the pair (or the only one that has children). If it is a pair of the same node, I  ********/
      /******** treat it brute-forcely if Na < N_cs or if it has no children, and I subdivise it else              ********/
      while (j > i){
            a = (stack + i) -> fst; //Id of first  node
            b = (stack + i) -> snd; //Id of second node
            Na = (FlatTree + a) -> how_many_dots;
            Nb = (FlatTree + b) -> how_many_dots;
            if (a == b){ //Cell self-interation
                  how_many_children_a = (FlatTree + a) -> how_many_children;
                  if (Na < N_cs_collision || how_many_children_a == 0){ //Direct interaction
                        dots_a = (FlatTree + a) -> dots;
                        for (p = 0; p < Na; p++){
                              for (q = p + 1; q < Na; q++){
                                    s = dots_a[p]; //Id of first  moonlet
                                    u = dots_a[q]; //Id of second moonlet
                                    if (!did_collide[s] && !did_collide[u]){
                                          the_approach = closest_approach(moonlets, s, u);
                                          if (the_approach != NULL) { //The closest approach leads to a collision and neither s nor u previously had a collision during that timestep
                                                if(elastic_collision_bool){ //All collisions are elastic
                                                      collision_treatment(moonlets, s, u, 0);
                                                }
                                                else if(inelastic_collision_bool){ //All collisions are inelastic
                                                      collision_treatment(moonlets, s, u, 1);
                                                }
                                                else if (instant_merger_bool){ //All collisions result in an instant merger
                                                      collision_treatment(moonlets, s, u, 2);
                                                }
                                                else if (fragmentation_bool){ //Collisions are treated with the fragmentation model described in the PDF draft
                                                      collision_treatment(moonlets, s, u, 3);
                                                }
                                          }
                                    } 
                              }
                        }
                  }
                  else{ //Subdivision
                        idFirstChild = (FlatTree + a) -> idFirstChild;
                        idLastChild  = idFirstChild + how_many_children_a;
                        /******** Making sure that the stack is big enough. To be removed when the code is robust ********/
                        if (j + 36 >= factor * cell_id){
                              fprintf(stderr, "Error : The stack is not big enough in function collision_flattree. Try increasing the value of factor.\n");
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
                              for (p = 0; p < Na; p++){
                                    for (q = 0; q < Nb; q++){
                                          s = dots_a[p]; //Id of first  moonlet
                                          u = dots_b[q]; //Id of second moonlet
                                          if (!did_collide[s] && !did_collide[u]){
                                                the_approach = closest_approach(moonlets, s, u);
                                                if (the_approach != NULL) { //The closest approach leads to a collision and neither s nor u previously collided during that timestep
                                                      if(elastic_collision_bool){ //All collisions are elastic
                                                            collision_treatment(moonlets, s, u, 0);
                                                      }
                                                      else if(inelastic_collision_bool){ //All collisions are inelastic
                                                            collision_treatment(moonlets, s, u, 1);
                                                      }
                                                      else if (instant_merger_bool){ //All collisions result in an instant merger
                                                            collision_treatment(moonlets, s, u, 2);
                                                      }
                                                      else if (fragmentation_bool){ //Collisions are treated with the fragmentation model described in the PDF draft
                                                            collision_treatment(moonlets, s, u, 3);
                                                      }
                                                }
                                          }
                                    }
                              }
                        }
                        else{ //Subdivision
                              /******** Making sure that the stack is big enough. To be removed when the code is robust ********/
                              if (j + 8 >= factor * cell_id){
                                    fprintf(stderr, "Error : The stack is not big enough in function collision_flattree. Try increasing the value of factor.\n");
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
            i++;
      }
      free(stack);
      stack = NULL;
}


/******************************************************************************/
/******** Implementing collision detection with the standard tree code ********/
/******************************************************************************/


void standard_tree_collision(struct node * FlatTree, struct moonlet * moonlets, int b){

      /******** Standard tree code to detect collisions involving moonlet b ********/

      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that have to be considered
      
      int i = 0; //Index in the stack of the current node's id
      int j = 0; //Index of where to put an id in the stack
      int a;     //Ids of the current node
      int s;     //Loop indexes
      int k;     //Moonlet indexes
      int * dots;//All the moonlets in the current cell
      int how_many_dots;  //Number of moonlets in the current cell
      int how_many_children; //Number of children of the current cell
      int idFirstChild; //Id of first child
      int idLastChild;  //Id of last child
      typ * com; //Current node's average moonlet position
      typ R[3];  //com - r_b
      typ r_crit; // Critical radii of the current node
      typ * the_approach;
      typ r;
      typ Rad= (moonlets + b) -> radius;
      typ vx = (moonlets + b) -> vx;
      typ vy = (moonlets + b) -> vy;
      typ vz = (moonlets + b) -> vz;
      typ v  = sqrt(vx*vx + vy*vy + vz*vz);
      
      stack[j] = 0;
      j++;
      

      while (j > i){
      
            if (did_collide[b]){
                  free(stack);
                  stack = NULL;
                  return;
            }
      
            a = stack[i];
            how_many_dots = (FlatTree + a) -> how_many_dots;
            r_crit = (FlatTree + a) -> r_crit;
            com    = (FlatTree + a) -> com;
            R[0]   = (moonlets + b) -> x - com[0];
            R[1]   = (moonlets + b) -> y - com[1];
            R[2]   = (moonlets + b) -> z - com[2];
            r = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
            if (r < r_crit + Rad+v*timestep){
                  how_many_children = (FlatTree + a) -> how_many_children;
                  if (how_many_dots < N_cb_collision || how_many_children == 0){ //Direct summation
                        dots = (FlatTree + a) -> dots;
                        for (s = 0; s < how_many_dots; s++){
                              k = dots[s];
                              if (k < b && !(did_collide[k]) && !(did_collide[b])){
                                    the_approach = closest_approach(moonlets, k, b);
                                    if (the_approach != NULL){ //The closest approach leads to a collision and neither k nor b previously had a collision during that timestep
                                          if(elastic_collision_bool){ //All collisions are elastic
                                                collision_treatment(moonlets, k, b, 0);
                                          }
                                          else if(inelastic_collision_bool){ //All collisions are inelastic
                                                collision_treatment(moonlets, k, b, 1);
                                          }
                                          else if (instant_merger_bool){ //All collisions result in an instant merger
                                                collision_treatment(moonlets, k, b, 2);
                                          }
                                          else if (fragmentation_bool){ //Collisions are treated with the fragmentation model described in the PDF draft
                                                collision_treatment(moonlets, k, b, 3);
                                          }
                                    }      
                              }
                        }   
                  }
                  else{ //Subdivision
                        idFirstChild = (FlatTree + a) -> idFirstChild;
                        idLastChild = idFirstChild + how_many_children;
                        for (s = idFirstChild; s < idLastChild; s++){
                              stack[j] = s;
                              j++; 
                        }
                  }
            }
            i++;
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function standard_tree_collision. Aborting before segmentation fault.\n");
            abort();
      }
      free(stack);
      stack = NULL;
}



