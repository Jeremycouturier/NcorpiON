/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    physics.c                                                   ********/
/******** @brief   This file computes the vector field and resolves collisions ********/
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
#include <errno.h>
#include <math.h>



void vector_field(struct moonlet * moonlets){

      /******** The vector field passed as argument to the integrator                 ********/
      /******** Stores the accelerations into the velocity field of moonlet structure ********/
      
      
      int k, p;
      typ rk, r2, r3, r5;
      typ X, Y, Z, XX, YY, ZZ, K, M;
      typ mu, mk;
      typ Xp, Yp, Zp, mp, DX, DY, DZ, D, D3, softening, Rk, Rp;
      typ X_sun, Y_sun, Z_sun, KK;
      typ aX = 0.0, aY = 0.0, aZ = 0.0;
      
      if (central_mass_bool){
            XX        = CM.x;  YY = CM.y;  ZZ = CM.z;
            M         = (inner_fluid_disk_bool ? CM.mass + fluid_disk_Sigma*M_PI*(Rroche*Rroche - R_unit*R_unit) : CM.mass);
            CM_acc[0] = 0.;  CM_acc[1] = 0.;  CM_acc[2] = 0.;
      }
      
      /******** Getting the coordinates of the three largest bodies. Useful only when treating mutual gravitational interactions with the O(N) mesh algorithm ********/
      #if mesh_bool
      typ X0, Y0, Z0, m0, R0, X1, Y1, Z1, m1, R1, X2, Y2, Z2, m2, R2;
      if (mutual_bool && mesh_bool && !force_naive_bool){
            if (exists[0]){
                  X0 =  moonlets       -> x;
                  Y0 =  moonlets       -> y;
                  Z0 =  moonlets       -> z;
                  m0 =  moonlets       -> mass;
                  R0 =  moonlets       -> radius;
            }
            if (exists[1]){
                  X1 =  (moonlets + 1) -> x;
                  Y1 =  (moonlets + 1) -> y;
                  Z1 =  (moonlets + 1) -> z;
                  m1 =  (moonlets + 1) -> mass;
                  R1 =  (moonlets + 1) -> radius;
            }
            if (exists[2]){
                  X2 =  (moonlets + 2) -> x;
                  Y2 =  (moonlets + 2) -> y;
                  Z2 =  (moonlets + 2) -> z;
                  m2 =  (moonlets + 2) -> mass;
                  R2 =  (moonlets + 2) -> radius;
            }
      }
      #endif
      
      for (k = 0; k <= largest_id; k++){
            if (*(exists + k)){ //Checking whether or not there is a body in the kth cell of the body array
            
                  X  = (moonlets + k) -> x;
                  Y  = (moonlets + k) -> y;
                  Z  = (moonlets + k) -> z;
                  mk = (moonlets + k) -> mass;
                  DX = X - XX;
                  DY = Y - YY;
                  DZ = Z - ZZ;
                  rk = sqrt(DX*DX + DY*DY + DZ*DZ);
                  r3 = rk*rk*rk;
            
                  /******** Contribution from the central mass ********/
                  mu = G*M;
                  
                  if (central_mass_bool){
                        aX         =  -mu*DX/r3;
                        aY         =  -mu*DY/r3;
                        aZ         =  -mu*DZ/r3;
                        CM_acc[0] += G*mk*DX/r3;
                        CM_acc[1] += G*mk*DY/r3;
                        CM_acc[2] += G*mk*DZ/r3;
                  }
                  
                  
                  /******** Contribution from the Earth symmetrical equatorial bulge ********/
                  if (J2_bool){
                        r2 = rk*rk;
                        r5 = r3*r2;
                        K  = mu*R_unit*R_unit*J2/r5;
                        /******** Updating the accelerations ********/
                        aX        += K*(-1.5*DX + 7.5/r2*DZ*DZ*DX);
                        aY        += K*(-1.5*DY + 7.5/r2*DZ*DZ*DY);
                        aZ        += K*(-4.5*DZ + 7.5/r2*DZ*DZ*DZ);
                        CM_acc[0] += K*( 1.5*DX - 7.5/r2*DZ*DZ*DX)*mk/M;
                        CM_acc[1] += K*( 1.5*DY - 7.5/r2*DZ*DZ*DY)*mk/M;
                        CM_acc[2] += K*( 4.5*DZ - 7.5/r2*DZ*DZ*DZ)*mk/M;
                  }
                  
                  
                  /******** Contribution from the star or companion star ********/
                  if (Sun_bool){
                        X_sun      = * sun_vector;
                        Y_sun      = *(sun_vector + 1);
                        Z_sun      = *(sun_vector + 2);
                        K          = G*star_mass/(star_semi_major*star_semi_major*star_semi_major);
                        KK         = 3.0*(X*X_sun + Y*Y_sun + Z*Z_sun)/(star_semi_major*star_semi_major);
                        /******** Updating the acceleration of body k ********/
                        aX -= K*(X - KK*X_sun);
                        aY -= K*(Y - KK*Y_sun);
                        aZ -= K*(Z - KK*Z_sun);
                  }
                  
                  
                  /******** Mutual gravitational interactions with the brute-force O(N^2) algorithm ********/
                  if (mutual_bool && (brute_force_bool || force_naive_bool)){
                        for (p = 0; p < k; p ++){
                              if (*(exists + p)){
                              
                                    /******** Getting the positions, masses and radii ********/
                                    mp = (moonlets + p) -> mass;
                                    Xp = (moonlets + p) -> x;
                                    Yp = (moonlets + p) -> y;
                                    Zp = (moonlets + p) -> z;
                                    Rp = (moonlets + p) -> radius;
                                    Rk = (moonlets + k) -> radius;
                                    
                                    DX = X - Xp;  DY = Y - Yp;  DZ = Z - Zp;
                                    softening = softening_parameter*(Rk + Rp);
                                    D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                                    D3 = D*D*D;

                                    /******** Updating the acceleration of body k ********/
                                    aX -= G*mp*DX/D3;
                                    aY -= G*mp*DY/D3;
                                    aZ -= G*mp*DZ/D3;
                                    
                                    /******** Updating the acceleration of body p ********/
                                    (moonlets + p) -> vx  += G*mk*DX/D3; //dV/dt=A
                                    (moonlets + p) -> vy  += G*mk*DY/D3;
                                    (moonlets + p) -> vz  += G*mk*DZ/D3;
                              }
                        }
                  }
                  
                  /******** Mutual gravitational interactions with the mesh O(N) algorithm ********/
                  #if mesh_bool
                  if (mutual_bool && mesh_bool && !force_naive_bool){
                  
                        /******** For now, I only consider gravitational interactions between pairs containing one of the three largest bodies and another body ********/
                        if (k > 2){
                              Rk = (moonlets + k) -> radius;
                              
                              /******** pair (0,k) ********/
                              if (exists[0]){
                                    DX = X - X0;  DY = Y - Y0;  DZ = Z - Z0;
                                    softening = softening_parameter*(Rk + R0);
                                    D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                                    D3 = D*D*D;
                                    /******** Updating the acceleration of body k ********/
                                    aX -= G*m0*DX/D3;
                                    aY -= G*m0*DY/D3;
                                    aZ -= G*m0*DZ/D3;
                                    /******** Updating the acceleration of body 0 ********/
                                    moonlets -> vx += G*mk*DX/D3;
                                    moonlets -> vy += G*mk*DY/D3;
                                    moonlets -> vz += G*mk*DZ/D3;
                              }
                              
                              /******** pair (1,k) ********/
                              if (exists[1]){
                                    DX = X - X1;  DY = Y - Y1;  DZ = Z - Z1;
                                    softening = softening_parameter*(Rk + R1);
                                    D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                                    D3 = D*D*D;
                                    /******** Updating the acceleration of body k ********/
                                    aX -= G*m1*DX/D3;
                                    aY -= G*m1*DY/D3;
                                    aZ -= G*m1*DZ/D3;
                                    /******** Updating the acceleration of body 1 ********/
                                    (moonlets + 1) -> vx += G*mk*DX/D3;
                                    (moonlets + 1) -> vy += G*mk*DY/D3;
                                    (moonlets + 1) -> vz += G*mk*DZ/D3;
                              }
                              
                              /******** pair (2,k) ********/
                              if (exists[2]){
                                    DX = X - X2;  DY = Y - Y2;  DZ = Z - Z2;
                                    softening = softening_parameter*(Rk + R2);
                                    D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                                    D3 = D*D*D;
                                    /******** Updating the acceleration of body k ********/
                                    aX -= G*m2*DX/D3;
                                    aY -= G*m2*DY/D3;
                                    aZ -= G*m2*DZ/D3;
                                    /******** Updating the acceleration of body 2 ********/
                                    (moonlets + 2) -> vx += G*mk*DX/D3;
                                    (moonlets + 2) -> vy += G*mk*DY/D3;
                                    (moonlets + 2) -> vz += G*mk*DZ/D3;
                              }
                        }      
                  }
                  #endif
                  
                  /******** Actualizing the acceleration of body k ********/
                  (moonlets + k) -> vx = aX; //dV/dt=A
                  (moonlets + k) -> vy = aY;
                  (moonlets + k) -> vz = aZ;
                  
            }     
      }
      
      
      /******** Mutual gravitational interactions with a tree algorithm (falcON or the standard tree code) ********/
      if (mutual_bool && (falcON_bool || standard_tree_bool) && !force_naive_bool){
            if (falcON_bool){
                  Cm_flattree(FlatTree, moonlets);              
                  Cm_downtree(FlatTree, moonlets);
            }
            else if (standard_tree_bool){
                  for (k = 0; k <= largest_id; k++){
                        if (exists[k]){
                              standard_tree_acceleration(FlatTree, moonlets, k);
                        }
                  }
            }
            for (k = 0; k <= largest_id; k ++){
                  if (exists[k]){                  
                        (moonlets + k) -> vx += C1Moonlets[3*k]  ;
                        (moonlets + k) -> vy += C1Moonlets[3*k + 1];
                        (moonlets + k) -> vz += C1Moonlets[3*k + 2];
                  }
            }
      }

      
      /******** Mutual gravitational interactions with the mesh O(N) algorithm ********/
      #if mesh_bool
      if (mutual_bool && mesh_bool && !brute_force_bool && !force_naive_bool){
      
            /******** Gravitational interactions between the three largest bodies ********/
            /******** pair (0,1) ********/
            if (exists[0] && exists[1]){
                  DX = X0 - X1;  DY = Y0 - Y1;  DZ = Z0 - Z1;
                  softening = softening_parameter*(R0 + R1);
                  D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                  D3 = D*D*D;
                  moonlets       -> vx  -= G*m1*DX/D3;
                  moonlets       -> vy  -= G*m1*DY/D3;
                  moonlets       -> vz  -= G*m1*DZ/D3;
                  (moonlets + 1) -> vx  += G*m0*DX/D3;
                  (moonlets + 1) -> vy  += G*m0*DY/D3;
                  (moonlets + 1) -> vz  += G*m0*DZ/D3;
            }
            
            /******** pair (0,2) ********/
            if (exists[0] && exists[2]){
                  DX = X0 - X2;  DY = Y0 - Y2;  DZ = Z0 - Z2;
                  softening = softening_parameter*(R0 + R2);
                  D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                  D3 = D*D*D;
                  moonlets       -> vx  -= G*m2*DX/D3;
                  moonlets       -> vy  -= G*m2*DY/D3;
                  moonlets       -> vz  -= G*m2*DZ/D3;
                  (moonlets + 2) -> vx  += G*m0*DX/D3;
                  (moonlets + 2) -> vy  += G*m0*DY/D3;
                  (moonlets + 2) -> vz  += G*m0*DZ/D3;
            }
            
            /******** pair (1,2) ********/
            if (exists[1] && exists[2]){
                  DX = X1 - X2;  DY = Y1 - Y2;  DZ = Z1 - Z2;
                  softening = softening_parameter*(R1 + R2);
                  D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                  D3 = D*D*D;
                  (moonlets + 1) -> vx  -= G*m2*DX/D3;
                  (moonlets + 1) -> vy  -= G*m2*DY/D3;
                  (moonlets + 1) -> vz  -= G*m2*DZ/D3;
                  (moonlets + 2) -> vx  += G*m1*DX/D3;
                  (moonlets + 2) -> vy  += G*m1*DY/D3;
                  (moonlets + 2) -> vz  += G*m1*DZ/D3;
            }
      
            /******** I now consider gravitational interactions between pairs in the same neighbourhood, excluding pairs containing a big body. ********/
            int j;
            
            for (j = 0; j < how_many_pairs; j ++){ //I go over all such pairs. The expected value of how_many_pairs is 0.5*N*how_many_neighbours = O(N)
                  
                  k = (pairs + j) -> fst; //The array "pairs" was updated by the function mesh in collision.c
                  p = (pairs + j) -> snd;
                  
                  if (*(exists + k) && *(exists + p) && k > 2 && p > 2){
                  
                        /******** Getting the positions, masses and radii ********/
                        X  = (moonlets + k) -> x;
                        Y  = (moonlets + k) -> y;
                        Z  = (moonlets + k) -> z;
                        mk = (moonlets + k) -> mass;
                        Rk = (moonlets + k) -> radius;
                        Xp = (moonlets + p) -> x;
                        Yp = (moonlets + p) -> y;
                        Zp = (moonlets + p) -> z;
                        mp = (moonlets + p) -> mass;
                        Rp = (moonlets + p) -> radius;
                  
                        DX = X - Xp;  DY = Y - Yp;  DZ = Z - Zp;
                        softening = softening_parameter*(Rk + Rp);
                        D  = sqrt(DX*DX + DY*DY + DZ*DZ + softening*softening);
                        D3 = D*D*D;
                        /******** Updading the acceleration of body k ********/
                        (moonlets + k) -> vx -= G*mp*DX/D3;
                        (moonlets + k) -> vy -= G*mp*DY/D3;
                        (moonlets + k) -> vz -= G*mp*DZ/D3;
                        /******** Updading the acceleration of body p ********/
                        (moonlets + p) -> vx += G*mk*DX/D3;
                        (moonlets + p) -> vy += G*mk*DY/D3;
                        (moonlets + p) -> vz += G*mk*DZ/D3;           
                  }
            }       
      }
      #endif
}


void tides(struct moonlet * X){

      /******** When this function is called, the body array xx contains the acceleration without tides    ********/
      /******** This function uses the knowledge of the acceleration without tides in order to compute the ********/
      /******** speed at the middle of the kick phase. Then, the acceleration due to tides can be computed ********/
      /******** and is added to the acceleration without tides in the array xx                             ********/
      
      three_largest_moonlets(X); //Only the three largest bodies raise tides. I retrieve their indexes
      typ largest_positions [9]; //Positions  of the three largest bodies at t - Delta t
      typ largest_velocities[9]; //Velocities of the three largest bodies at t
      int j, k;
      
      /******** Estimating the velocities in the middle of the kick phase ********/
      for (k = 0; k <= 2; k ++){
            if (three_largest_indexes[k] != -1 && *(exists + three_largest_indexes[k])){
                  largest_velocities[3*k]     = (X + three_largest_indexes[k]) -> vx + 0.5*timestep * (xx + three_largest_indexes[k]) -> vx;
                  largest_velocities[3*k + 1] = (X + three_largest_indexes[k]) -> vy + 0.5*timestep * (xx + three_largest_indexes[k]) -> vy;
                  largest_velocities[3*k + 2] = (X + three_largest_indexes[k]) -> vz + 0.5*timestep * (xx + three_largest_indexes[k]) -> vz;
            }
      }
      
      /******** Computing the r_k^\star at t - Delta t ********/
      for (k = 0; k <= 2; k ++){
            if (three_largest_indexes[k] != -1 && *(exists + three_largest_indexes[k])){
                  largest_positions [3*k]     = (X + three_largest_indexes[k]) -> x - Delta_t * (largest_velocities[3*k]     + SideralOmega*(X + three_largest_indexes[k]) -> y);
                  largest_positions [3*k + 1] = (X + three_largest_indexes[k]) -> y - Delta_t * (largest_velocities[3*k + 1] - SideralOmega*(X + three_largest_indexes[k]) -> x);
                  largest_positions [3*k + 2] = (X + three_largest_indexes[k]) -> z - Delta_t * (largest_velocities[3*k + 2]);
            }
      }
      
      /******** Computing the tidal acceleration ********/
      typ R_unit3 = R_unit *R_unit*R_unit;
      typ R_unit5 = R_unit3*R_unit*R_unit;
      typ scalar_product;
      typ xj, yj, zj, xk, yk, zk;
      typ rj2, rk2, rj5rk5, mj, mk, factor;
      for (j = 0; j <= largest_id; j ++){ //Looping over all bodies
            if (*(exists + j)){
                  for (k = 0; k <= 2; k ++){ //Looping over all three perturbing bodies
                        if (three_largest_indexes[k] != -1 && *(exists + three_largest_indexes[k])){
                              xj = (X + j) -> x;
                              yj = (X + j) -> y;
                              zj = (X + j) -> z;
                              xk = largest_positions[3*k];
                              yk = largest_positions[3*k + 1];
                              zk = largest_positions[3*k + 2];
                              scalar_product = xj*xk + yj*yk + yj*yk;
                              rj2 = xj*xj + yj*yj + zj*zj;
                              rk2 = xk*xk + yk*yk + zk*zk;
                              rj5rk5 = rj2 * rj2 * rk2 * rk2 * sqrt(rj2*rk2);
                              mj  = (X + j) -> mass;
                              mk  = (X + three_largest_indexes[k]) -> mass;
                              
                              /******** Adding the contributions from tides to the acceleration ********/
                              factor = 1.5*k2*G*mk*R_unit5/rj5rk5;
                              (xx + j) -> vx += factor*(rk2*xj - 5.0/rj2*scalar_product*scalar_product*xj + 2.0*scalar_product*xk);
                              (xx + j) -> vy += factor*(rk2*yj - 5.0/rj2*scalar_product*scalar_product*yj + 2.0*scalar_product*yk);
                              (xx + j) -> vz += factor*(rk2*zj - 5.0/rj2*scalar_product*scalar_product*zj + 2.0*scalar_product*zk);
                              
                              /******** Computing the new sideral rotation of the central body ********/
                              SideralOmega -= 3.0*k2*timestep*G*mj*mk*R_unit3/(dimensionless_moi*M_unit*rj5rk5)*scalar_product*(xj*yk - yj*xk);
                              
                              /******** If j is not one of the three largest bodies, I must consider the reciprocal interaction to conserve angular momentum ********/
                              /******** I limit myself to elastic tides for the reciprocal interaction                                                       ********/
                              if (j != *three_largest_indexes && j != three_largest_indexes[1] && j != three_largest_indexes[2]){
                                    xk = (X + three_largest_indexes[k]) -> x;
                                    yk = (X + three_largest_indexes[k]) -> y;
                                    zk = (X + three_largest_indexes[k]) -> z;
                                    rk2 = xk*xk + yk*yk + zk*zk;
                                    rj5rk5 = rj2 * rj2 * rk2 * rk2 * sqrt(rj2*rk2);
                                    factor = 1.5*k2*G*mj*R_unit5/rj5rk5;
                                    scalar_product = xj*xk + yj*yk + yj*yk;
                                    (xx + three_largest_indexes[k]) -> vx += factor*(rj2*xk - 5.0/rk2*scalar_product*scalar_product*xk + 2.0*scalar_product*xj);
                                    (xx + three_largest_indexes[k]) -> vy += factor*(rj2*yk - 5.0/rk2*scalar_product*scalar_product*yk + 2.0*scalar_product*yj);
                                    (xx + three_largest_indexes[k]) -> vz += factor*(rj2*zk - 5.0/rk2*scalar_product*scalar_product*zk + 2.0*scalar_product*zj);
                                    SideralOmega -= 3.0*k2*timestep*G*mj*mk*R_unit3/(dimensionless_moi*M_unit*rj5rk5)*scalar_product*(xk*yj - yk*xj);
                              }
                        }
                  }
            }
      }
}


void collision(struct moonlet * moonlets, int a, int b, typ f){

      /******** Treats the collision between bodies a and b. The collision's elasticity is determined       ********/
      /******** by the parameter 1 <= collision_parameter <= 2 defined in the parameter file                ********/
      /******** The array approach returned by the function closest_approach contains                       ********/
      /******** the positions of the bodies at the collision. See top of function closest_approach          ********/

      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;     //Cartesian speeds of the bodies.
      typ xa, ya, za, xb, yb, zb;                 //Cartesian positions at the collision
      typ R_a = (moonlets + a) -> radius;         //The bodies' radii
      typ R_b = (moonlets + b) -> radius;
      typ R   = R_a + R_b;                        //Sum of the radii;
      typ m_a = (moonlets + a) -> mass;           //The bodies's masses
      typ m_b = (moonlets + b) -> mass;              
      
      
      /******** Getting the speeds ********/
      vx_a = (moonlets + a) -> vx;
      vy_a = (moonlets + a) -> vy;
      vz_a = (moonlets + a) -> vz;
      vx_b = (moonlets + b) -> vx;
      vy_b = (moonlets + b) -> vy;
      vz_b = (moonlets + b) -> vz;
      
      
      /******** Getting the positions ********/
      xa = *(approach);
      ya = *(approach + 1);
      za = *(approach + 2);
      xb = *(approach + 3);
      yb = *(approach + 4);
      zb = *(approach + 5);
      
      
      /******** Calculating the new speeds after impact ********/
      typ dr_dot_dv = (xa-xb)*(vx_a-vx_b)+(ya-yb)*(vy_a-vy_b)+(za-zb)*(vz_a-vz_b); //Scalar product dv.dr where dr=r_b-r_a and dv=v_b-v_a are the relative position and speed
      typ alpha = f*m_a*m_b/((m_a+m_b)*R*R);                                       //Factor such that J=alpha*(dv.dr)*dr
      typ Jx = alpha*dr_dot_dv*(xb-xa);                                            //x component of J
      typ Jy = alpha*dr_dot_dv*(yb-ya);                                            //y component of J
      typ Jz = alpha*dr_dot_dv*(zb-za);                                            //z component of J
      
      
      /******** Calculating the speeds after impact ********/
      vx_a += Jx/m_a;  //v_a=v_a+J/m_a
      vy_a += Jy/m_a;
      vz_a += Jz/m_a;
      vx_b -= Jx/m_b;  //v_b=v_b-J/m_b
      vy_b -= Jy/m_b;
      vz_b -= Jz/m_b;
      
      
      /******** Actualizing the speeds and positions in the array moonlets                 ********/
      /******** Bodies are put back where they were at the beginning of the timestep, but  ********/
      /******** with their post-impact velocity, then they drift with all the other bodies ********/
      (moonlets + a) -> x  = xa - time_until_collision * vx_a;  //Actualizing body a
      (moonlets + a) -> y  = ya - time_until_collision * vy_a;
      (moonlets + a) -> z  = za - time_until_collision * vz_a;
      (moonlets + a) -> vx = vx_a;
      (moonlets + a) -> vy = vy_a;
      (moonlets + a) -> vz = vz_a;
      (moonlets + b) -> x  = xb - time_until_collision * vx_b;  //Actualizing body b
      (moonlets + b) -> y  = yb - time_until_collision * vy_b;
      (moonlets + b) -> z  = zb - time_until_collision * vz_b;
      (moonlets + b) -> vx = vx_b;
      (moonlets + b) -> vy = vy_b;
      (moonlets + b) -> vz = vz_b;
      

}


void merger(struct moonlet * moonlets, int a, int b){

      /******** Merges body a and body b together ********/


      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;     //Cartesian speeds of the bodies.
      typ xa, ya, za, xb, yb, zb;                 //Cartesian positions at the collision
      typ m_a = (moonlets + a) -> mass;           //The bodies's masses
      typ m_b = (moonlets + b) -> mass;           
      typ R_a = (moonlets + a) -> radius;         //The bodies's radii
      typ R_b = (moonlets + b) -> radius;
      typ rho_a = 3.0*m_a/(4.0*M_PI*R_a*R_a*R_a); //Density of body a
      typ rho_b = 3.0*m_b/(4.0*M_PI*R_b*R_b*R_b); //Density of body b
      typ m = m_a + m_b;                          //Sum of the masses
      typ r_tilde[3]; //Position of the merger
      typ v_tilde[3]; //Speed of the merger
      typ average_density = (m_a*rho_a + m_b*rho_b)/m; //The average density of the moonlets, weighted by mass
      
      /******** Getting the speeds ********/
      vx_a = (moonlets + a) -> vx;
      vy_a = (moonlets + a) -> vy;
      vz_a = (moonlets + a) -> vz;
      vx_b = (moonlets + b) -> vx;
      vy_b = (moonlets + b) -> vy;
      vz_b = (moonlets + b) -> vz;

      /******** Getting the positions ********/
      xa = * approach;
      ya = *(approach + 1);
      za = *(approach + 2);
      xb = *(approach + 3);
      yb = *(approach + 4);
      zb = *(approach + 5);
      
      /******** The speed after merging ********/ 
      v_tilde[0] = (m_a*vx_a + m_b*vx_b)/m;
      v_tilde[1] = (m_a*vy_a + m_b*vy_b)/m;
      v_tilde[2] = (m_a*vz_a + m_b*vz_b)/m; 
                  
      /******** The position after merging ********/
      r_tilde[0] = (m_a*xa + m_b*xb)/m;
      r_tilde[1] = (m_a*ya + m_b*yb)/m;
      r_tilde[2] = (m_a*za + m_b*zb)/m;     
      
      /******** Actualizing the position ********/
      (moonlets + a) -> x  = r_tilde[0] - time_until_collision*v_tilde[0];
      (moonlets + a) -> y  = r_tilde[1] - time_until_collision*v_tilde[1];
      (moonlets + a) -> z  = r_tilde[2] - time_until_collision*v_tilde[2];
      
      /******** Actualizing the position ********/
      (moonlets + a) -> vx = v_tilde[0];
      (moonlets + a) -> vy = v_tilde[1];
      (moonlets + a) -> vz = v_tilde[2];
      
      /******** Actualizing the mass and radius ********/
      (moonlets + a) -> mass   = m;
      (moonlets + a) -> radius = pow(3.0*m/(4.0*M_PI*average_density),1.0/3.0);
      
      /******** Body b does not exist anymore and I disallow body a to collide again for that timestep ********/        
      lose_moonlet(b);
      *(did_collide + a) = one_collision_only_bool;
}

void fragmentation(struct moonlet * moonlets, int a, int b){

      /******** Treats the fragmentation due to the collision between bodies a and b ********/
      
      typ C1_3mu = pow(C1_parameter, 3.0*mu_parameter);
      typ stigma = (3.0*mu_parameter-1.0)/(3.0*mu_parameter);
      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;                  //Cartesian speeds of the bodies.
      typ xa, ya, za, xb, yb, zb;                              //Cartesian positions at the collision
      typ dx, dy, dz, dvx, dvy, dvz;                           //dx is xa-xb, and so on.
      typ R_a = (moonlets + a) -> radius;                      //The bodies' radii
      typ R_b = (moonlets + b) -> radius;
      typ R   = R_a + R_b;                                     //Sum of the radii;
      typ m_a = (moonlets + a) -> mass;                        //The bodies's masses
      typ m_b = (moonlets + b) -> mass;
      typ M   = m_a + m_b;                                     //Sum of the masses
      typ m_1;                                                 //Mass of the impactor
      typ rho_12_3nu1;                                         //Ratio between the density of the impactor and that of the target at the power 3*nu - 1
      typ rho_a = 3.0*m_a/(4.0*M_PI*R_a*R_a*R_a);              //Density of body a
      typ rho_b = 3.0*m_b/(4.0*M_PI*R_b*R_b*R_b);              //Density of body b
      typ average_density = (m_a*rho_a + m_b*rho_b)/M;         //The average density of the moonlets, weighted by mass
      
      /******** Defining the impactor and the target. ********/
      if (R_a > R_b){ // a is the target and b the impactor
            m_1         = m_b;
            rho_12_3nu1 = pow(rho_b/rho_a, 3.0*nu_parameter - 1.0);
      }
      else { // a is the impactor and b the target
            m_1         = m_a;
            rho_12_3nu1 = pow(rho_a/rho_b, 3.0*nu_parameter - 1.0);
      }            
      
      /******** Getting the speeds at the impact ********/
      vx_a = (moonlets + a) -> vx;
      vy_a = (moonlets + a) -> vy;
      vz_a = (moonlets + a) -> vz;
      vx_b = (moonlets + b) -> vx;
      vy_b = (moonlets + b) -> vy;
      vz_b = (moonlets + b) -> vz;
      
      /******** Getting the positions at the impact ********/
      xa = * approach;
      ya = *(approach + 1);
      za = *(approach + 2);
      xb = *(approach + 3);
      yb = *(approach + 4);
      zb = *(approach + 5);
      
      /******** Getting the relative positions and velocities ********/
      dx  = xa   - xb;
      dy  = ya   - yb;
      dz  = za   - zb;
      dvx = vx_a - vx_b;
      dvy = vy_a - vy_b;
      dvz = vz_a - vz_b;
      
      typ dr_dot_dv = dx*dvx + dy*dvy + dz*dvz;          //  (r_1-r_2).(v_1-v_2)
      typ dv_norm   = sqrt(dvx*dvx + dvy*dvy + dvz*dvz); //  ||v_1-v_2||
      typ costheta  = -dr_dot_dv/(R*dv_norm);            //Cosine of impact angle;
      typ R_eq = pow(3.0*M/(4.0*M_PI*average_density),1.0/3.0);    //Hypothetical radius if the impactor were to merge on the target
      typ vesc = sqrt(2.0*G*M/R_eq);                     //Escape velocity at the surface if the impactor were to merge on the target
      typ m_check = 3.0*k_parameter*C1_3mu*m_1/(4.0*M_PI)*pow(dv_norm*costheta/vesc, 3.0*mu_parameter)*rho_12_3nu1; //Mass of the tail
      typ m_tilde = M - m_check;                         //Mass of the largest fragment
      typ m_tilde_2 = m_check / (typ) N_tilde;           //Mass of the tail's fragments
      
      
      /******** Merger case. Bodies a and b merge together ********/
 
      if (m_check < frag_threshold && m_tilde >= 0.1 * M){
            merger(moonlets, a ,b);
            merger_count ++;
            return;
      }
      
      
      /******** Some data used for the fragmentation case ********/
      typ v_cm[3], r_cm[3]; //Velocity and position of the center of mass before impact
      typ r_tilde[3]; //Position of the largest fragment
      typ v_tilde[3]; //Speed of the largest fragment
      v_cm[0] = (m_a*vx_a + m_b*vx_b)/M;  v_cm[1] = (m_a*vy_a + m_b*vy_b)/M;  v_cm[2] = (m_a*vz_a + m_b*vz_b)/M;
      r_cm[0] = (m_a*xa   + m_b*  xb)/M;  r_cm[1] = (m_a*ya   + m_b*  yb)/M;  r_cm[2] = (m_a*za   + m_b*  zb)/M;
      
      
      /******** Super-catastrophic fragmentation. The mass of the largest fragment is less than 10 % of the total mass. The ejecta is discarded ********/
      
      if (m_tilde < 0.1*M){ //m_tilde is not proportionnal to Qr/Qr* in this regime
      
            m_tilde = 0.1*M*pow(2.0/1.8*m_check/M, -1.5); //Eq. (44) of Leinhardt and Stewart (2012)
            
            /******** Only the largest fragment remains in the super-catastrophic regime ********/
            /******** Actualizing its speed ********/
            (moonlets + a) -> vx = v_cm[0];
            (moonlets + a) -> vy = v_cm[1];
            (moonlets + a) -> vz = v_cm[2];
            /******** Actualizing its position ********/
            (moonlets + a) -> x  = r_cm[0] - time_until_collision*v_cm[0];
            (moonlets + a) -> y  = r_cm[1] - time_until_collision*v_cm[1];
            (moonlets + a) -> z  = r_cm[2] - time_until_collision*v_cm[2];
            /******** Actualizing its mass and radius ********/
            (moonlets + a) -> mass   = m_tilde;
            (moonlets + a) -> radius = pow(3.0*m_tilde/(4.0*M_PI*average_density), 1.0/3.0);
            /******** Reducing the center of mass to compensate for the lost momentum ********/
            need_to_reduce_COM_bool = 1;
            /******** Body b does not exist anymore and I disallow body a to collide again for that timestep ********/
            lose_moonlet(b);
            *(did_collide + a) = one_collision_only_bool;
            super_catastrophic_count ++;
            if (m_tilde < frag_threshold){ //Body a is discarded as well
                  lose_moonlet(a);
            }
            return;
      }
      
      
      /******** Partial fragmentation. The tail is reunited into a single body. ********/
      
      if (m_check >= frag_threshold && m_tilde_2 < frag_threshold){
            typ r_k[3]; //Position of the tail with respect to the largest fragment
            typ v_k[3]; //Velocity of the tail with respect to the largest fragment
            typ v_k_scalar = vesc/stigma; //Scalar velocity of the tail with respect to the largest fragment
            typ in_front_of = v_k_scalar/R;
            if (R_a > R_b){
                  r_k[0] = -dx;  r_k[1] = -dy;  r_k[2] = -dz;
            }
            else {
                  r_k[0] = dx;  r_k[1] = dy;  r_k[2] = dz;
            }
            v_k[0] = in_front_of * r_k[0];  v_k[1] = in_front_of * r_k[1];  v_k[2] = in_front_of * r_k[2];           

            /******** Total momentum is conserved ********/
            r_tilde[0] = r_cm[0] - m_check/M*r_k[0];
            r_tilde[1] = r_cm[1] - m_check/M*r_k[1];
            r_tilde[2] = r_cm[2] - m_check/M*r_k[2];
            v_tilde[0] = v_cm[0] - m_check/M*v_k[0];
            v_tilde[1] = v_cm[1] - m_check/M*v_k[1];
            v_tilde[2] = v_cm[2] - m_check/M*v_k[2];  
            
            /******** Actualizing the bodies ********/
            typ R_tilde = pow(3.0*m_tilde/(4.0*M_PI*average_density),1.0/3.0); //Radius of the largest fragment
            typ R_check = pow(3.0*m_check/(4.0*M_PI*average_density),1.0/3.0); //Radius of the tail fragment
            (moonlets + a) -> x      = r_tilde[0] - time_until_collision*v_tilde[0];
            (moonlets + a) -> y      = r_tilde[1] - time_until_collision*v_tilde[1];
            (moonlets + a) -> z      = r_tilde[2] - time_until_collision*v_tilde[2];
            (moonlets + a) -> vx     = v_tilde[0];
            (moonlets + a) -> vy     = v_tilde[1];
            (moonlets + a) -> vz     = v_tilde[2];
            (moonlets + a) -> mass   = m_tilde;
            (moonlets + a) -> radius = R_tilde;
            (moonlets + b) -> x      = r_tilde[0] + r_k[0] - time_until_collision*(v_tilde[0] + v_k[0]);
            (moonlets + b) -> y      = r_tilde[1] + r_k[1] - time_until_collision*(v_tilde[1] + v_k[1]);
            (moonlets + b) -> z      = r_tilde[2] + r_k[2] - time_until_collision*(v_tilde[2] + v_k[2]);
            (moonlets + b) -> vx     = v_tilde[0] + v_k[0];
            (moonlets + b) -> vy     = v_tilde[1] + v_k[1];
            (moonlets + b) -> vz     = v_tilde[2] + v_k[2];
            (moonlets + b) -> mass   = m_check;
            (moonlets + b) -> radius = R_check;
            
            if (mutual_bool && mesh_bool && !force_naive_bool){ //If the mutual gravitational interactions are considered, I register the pair (a,b) to be taken care of
                  (pairs + how_many_pairs) -> fst = a;
                  (pairs + how_many_pairs) -> snd = b;
                  how_many_pairs ++;
            }
            
            /******** I disallow bodies a and be to collide again during that timestep ********/
            *(did_collide + a) = one_collision_only_bool;
            *(did_collide + b) = one_collision_only_bool;
            half_fragmentation_count ++;
            return;
      }
      
      
      /******** Full fragmentation. The tail is made up of N_tilde bodies ********/

      if (m_tilde_2 >= frag_threshold){
            /******** Determination of the r_k' and v_k' ********/ 
            typ r_k[3*N_tilde]; //Position of the fragments of the tail with respect to the largest fragment
            typ v_k[3*N_tilde]; //Velocity of the fragments of the tail with respect to the largest fragment
            typ dr[3]; //r_1-r_2
            typ dv[3]; //v_1-v_2
            typ R_tilde_2 = pow(3.0*m_tilde_2/(4.0*M_PI*average_density), 1.0/3.0); //Radius of the fragments of the tail
            typ u[3]; //Vectors used to locate the fragments of the tail
            typ v[3];
            if (R_a > R_b){
                  dr[0] = -dx;   dr[1] = -dy;   dr[2] = -dz;
                  dv[0] = -dvx;  dv[1] = -dvy;  dv[2] = -dvz;
            }
            else {
                  dr[0] = dx;   dr[1] = dy;   dr[2] = dz;
                  dv[0] = dvx;  dv[1] = dvy;  dv[2] = dvz;
            }
            typ dr_x_dv[3]; // dr x dv
            typ dr_x_dv_norm; // ||dr x dv||
            int pq[4] = pq_min_max; // {p_k_min, p_k_max, q_k_min, q_k_max}
            int p,q;
            cross_product(dr[0], dr[1], dr[2], dv[0], dv[1], dv[2], dr_x_dv);
            dr_x_dv_norm = sqrt(dr_x_dv[0]*dr_x_dv[0] + dr_x_dv[1]*dr_x_dv[1] + dr_x_dv[2]*dr_x_dv[2]);
            if (dr_x_dv_norm > 1.0e-5){ //Oblique collision
                  /******** v = (dr x dv) / ||dr x dv|| ********/
                  v[0] = dr_x_dv[0]/dr_x_dv_norm;  v[1] = dr_x_dv[1]/dr_x_dv_norm;  v[2] = dr_x_dv[2]/dr_x_dv_norm; //Defining the unit vector v
            }
            else{ //Nearly Frontal collision
                  /******** v is any unit vector orthogonal to dr ********/
                  int m = maximum(fabs(dr[0]), fabs(dr[1]), fabs(dr[2]));
                  typ numerator;
                  if (m == 0){
                        numerator = dr[1] + dr[2];
                        v[1] = 1.0/sqrt(2.0 + numerator*numerator/(dr[0]*dr[0])); //Defining the unit vector v
                        v[2] =  v[1];
                        v[0] = -v[1]*numerator/dr[0];
                  }
                  else if (m == 1){
                        numerator = dr[0] + dr[2];
                        v[0] = 1.0/sqrt(2.0 + numerator*numerator/(dr[1]*dr[1])); //Defining the unit vector v
                        v[2] =  v[0];
                        v[1] = -v[0]*numerator/dr[1];
                  }
                  else {
                        numerator = dr[0] + dr[1];
                        v[0] = 1.0/sqrt(2.0 + numerator*numerator/(dr[2]*dr[2])); //Defining the unit vector v
                        v[1] =  v[0];
                        v[2] = -v[0]*numerator/dr[2];
                  }
            }
            typ v_x_dr[3]; // v x dr
            cross_product(v[0], v[1], v[2], dr[0], dr[1], dr[2], v_x_dr);
            u[0] = v_x_dr[0]/R;  u[1] = v_x_dr[1]/R;  u[2] = v_x_dr[2]/R; //Defining the unit vector u
            int n = 0;
            typ v_k_scalar, two_p_R_tilde_2, two_q_R_tilde_2, in_front_of;
            for (p = pq[0]; p <= pq[1]; p ++){ //I travel along the rectangle of integer coordinate points to define the position and speeds of the fragments of the tail
                  for (q = pq[2]; q <= pq[3]; q ++){
                        two_p_R_tilde_2 = ((typ) (2*p))*R_tilde_2;  two_q_R_tilde_2 = ((typ) (2*q))*R_tilde_2;
                        r_k[3*n]    = dr[0] + two_p_R_tilde_2*u[0] + two_q_R_tilde_2*v[0]; // x-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        r_k[3*n+1]  = dr[1] + two_p_R_tilde_2*u[1] + two_q_R_tilde_2*v[1]; // y-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        r_k[3*n+2]  = dr[2] + two_p_R_tilde_2*u[2] + two_q_R_tilde_2*v[2]; // z-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        v_k_scalar  = vesc * N_tilde/stigma *(pow(1.0-((typ) n)/N_tilde,stigma)-pow(1.0-((typ) (n+1))/N_tilde,stigma)); //Scalar velocity of the tail's (n+1)^th fragments
                        in_front_of = v_k_scalar/sqrt(R*R+4.0*R_tilde_2*R_tilde_2*(p*p+q*q));
                        v_k[3*n]    = in_front_of * r_k[3*n];  v_k[3*n + 1] = in_front_of * r_k[3*n + 1];  v_k[3*n + 2] = in_front_of * r_k[3*n + 2];
                        n++;
                  }
            }
            
            /******** Total momentum is conserved ********/
            typ corr[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            for (n = 0; n < N_tilde; n ++){
                  corr[0] += r_k[3*n];
                  corr[1] += r_k[3*n + 1];
                  corr[2] += r_k[3*n + 2];
                  corr[3] += v_k[3*n];
                  corr[4] += v_k[3*n + 1];
                  corr[5] += v_k[3*n + 2];
            }
            r_tilde[0] = r_cm[0] - m_tilde_2/M*corr[0];
            r_tilde[1] = r_cm[1] - m_tilde_2/M*corr[1];
            r_tilde[2] = r_cm[2] - m_tilde_2/M*corr[2];
            v_tilde[0] = v_cm[0] - m_tilde_2/M*corr[3];
            v_tilde[1] = v_cm[1] - m_tilde_2/M*corr[4];
            v_tilde[2] = v_cm[2] - m_tilde_2/M*corr[5];
            
            /******** Managing the indexes of all the bodies ********/
            int id[N_tilde + 1]; //The N_tilde + 1 indexes where the N_tilde + 1 bodies will be stored in the array moonlets
            id[0] = a;           //Putting the largest fragment there
            for (n = 1; n <= N_tilde; n ++){
                  if (mesh_bool && !force_naive_bool){ //By construction of the mesh algorithm, the fragments must be put at the end 
                        id[n] = get_free_index(1); //Putting the remaining N_tilde fragments of the tail there
                  }
                  else{ //No need to put the fragments at the end
                        id[n] = get_free_index(0); //Putting the remaining N_tilde fragments of the tail there
                  }
            }
            if (mutual_bool && mesh_bool && !force_naive_bool){ //Adding the N_tilde*(N_tilde+1)/2 pairs to be taken into account for gravitational interactions
                  for (p = 0; p <= N_tilde; p++){
                        for (q = 0; q < p; q++){
                              (pairs + how_many_pairs) -> fst = id[p];
                              (pairs + how_many_pairs) -> snd = id[q];
                              how_many_pairs ++;
                        }
                  }
            }
            
            lose_moonlet(b); //Body b does not exist anymore
            for (n = 0; n <= N_tilde; n ++){
                  *(did_collide + id[n]) = one_collision_only_bool; //Bodies won't be able to collide during that timestep if one_collision_only_bool
            }

            /******** Actualizing the largest fragment ********/
            (moonlets + id[0]) -> x      = r_tilde[0] - time_until_collision*v_tilde[0];
            (moonlets + id[0]) -> y      = r_tilde[1] - time_until_collision*v_tilde[1];
            (moonlets + id[0]) -> z      = r_tilde[2] - time_until_collision*v_tilde[2];
            (moonlets + id[0]) -> vx     = v_tilde[0];
            (moonlets + id[0]) -> vy     = v_tilde[1];
            (moonlets + id[0]) -> vz     = v_tilde[2];
            (moonlets + id[0]) -> mass   = m_tilde;
            (moonlets + id[0]) -> radius = pow(3.0*m_tilde/(4.0*M_PI*average_density),1.0/3.0);

            /******** Actualizing the fragments of the tail ********/
            for (n = 0; n < N_tilde; n++){
                  (moonlets + id[n + 1]) -> x      = r_tilde[0] + r_k[3*n]     - time_until_collision*(v_tilde[0] + v_k[3*n]);
                  (moonlets + id[n + 1]) -> y      = r_tilde[1] + r_k[3*n + 1] - time_until_collision*(v_tilde[1] + v_k[3*n + 1]);
                  (moonlets + id[n + 1]) -> z      = r_tilde[2] + r_k[3*n + 2] - time_until_collision*(v_tilde[2] + v_k[3*n + 2]);
                  (moonlets + id[n + 1]) -> vx     = v_tilde[0] + v_k[3*n];
                  (moonlets + id[n + 1]) -> vy     = v_tilde[1] + v_k[3*n + 1];
                  (moonlets + id[n + 1]) -> vz     = v_tilde[2] + v_k[3*n + 2];
                  (moonlets + id[n + 1]) -> mass   = m_tilde_2;
                  (moonlets + id[n + 1]) -> radius = R_tilde_2;
            }
            full_fragmentation_count ++;
            return;
      }
}


void collision_treatment(struct moonlet * moonlets, int a, int b, int type_of_collision){

      /******** Treats the collision between body a and body b                  ********/
      /******** The integer type_of_collision determines the type of collision. ********/
      /******** 0 --> Elastic collision                                         ********/
      /******** 1 --> Inelastic collision                                       ********/
      /******** 2 --> Instant merger                                            ********/
      /******** 3 --> Fragmentation                                             ********/
      
      typ m_a = (moonlets + a) -> mass;
      typ m_b = (moonlets + b) -> mass;
      int A   = (m_a > m_b ? a : b);
      int B   = (m_a > m_b ? b : a);
      
      if (!(*(exists + a) && *(exists + b))){
            return;
      }
      
      if (type_of_collision == 0){
            collision(moonlets, A, B, 2.0);
            *(did_collide + a) = one_collision_only_bool;
            *(did_collide + b) = one_collision_only_bool;
            if (mutual_bool && mesh_bool && !force_naive_bool){ //If the mutual gravity is considered, I register the pair (a,b) to be taken care of
                  (pairs + how_many_pairs) -> fst = a;
                  (pairs + how_many_pairs) -> snd = b;
                  how_many_pairs ++;
            }
      }
      
      else if (type_of_collision == 1){
            collision(moonlets, A, B, collision_parameter);
            *(did_collide + a) = one_collision_only_bool;
            *(did_collide + b) = one_collision_only_bool;
            if (mutual_bool && mesh_bool && !force_naive_bool){ //If the mutual gravity is considered, I register the pair (a,b) to be taken care of
                  (pairs + how_many_pairs) -> fst = a;
                  (pairs + how_many_pairs) -> snd = b;
                  how_many_pairs ++;
            }
      }
      
      else if (type_of_collision == 2){
            merger(moonlets, A, B);
      }
      
      else if (type_of_collision == 3){
            fragmentation(moonlets, A, B);
      }
      collision_count ++;
}


void get_neighbours_mesh(struct moonlet * moonlets){


      /******** When mutual interactions are treated with the mesh algorithm and neighbouring pairs ********/
      /******** were not previously found by the function mesh (e.g. if collision_bool is 0)        ********/


      int k,p;
      int current_largest_id = largest_id;
      
      for (k = 0; k <= current_largest_id; k++){
            if(*(exists + k)){ //Checking if there is a body in the k^th cell of the array moonlets
                  neighbours(moonlets, k); //Adding the body to the hash table, and retrieving its neighbours in the chain nghb
                        
                  while(nghb -> how_many > 0){ //If the body is inside the collision cube and has at least one neighbour
                        p = (nghb -> ids)[nghb -> how_many - 1]; // p is the id of a body neighbour to k
                        (pairs + how_many_pairs) -> fst = k;
                        (pairs + how_many_pairs) -> snd = p;
                        how_many_pairs ++;
                        nghb = partial_delete(nghb);
                        total_neighbours ++;
                  }
            }
      }
}





