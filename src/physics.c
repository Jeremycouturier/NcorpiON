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
#include "spring.h"
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
      typ aX, aY, aZ;
      
      if (central_mass_bool){
            XX        = CM.x;  YY = CM.y;  ZZ = CM.z;
            M         = (inner_fluid_disk_bool ? CM.mass + fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) : CM.mass);
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
      
      for (k = 0; k <= largest_id; k ++){
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
                  else{
                        aX = 0.;  aY = 0.;  aZ = 0.;
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
                  if (mutual_bool && (brute_force_bool || force_naive_bool) && 0){
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
                  for (k = 0; k <= largest_id; k ++){
                        if (exists[k]){
                              standard_tree_acceleration(FlatTree, moonlets, k);
                        }
                  }
            }
            for (k = 0; k <= largest_id; k ++){
                  if (exists[k]){                  
                        (moonlets + k) -> vx += C1Moonlets[3*k];
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
                        /******** Updating the acceleration of body k ********/
                        (moonlets + k) -> vx -= G*mp*DX/D3;
                        (moonlets + k) -> vy -= G*mp*DY/D3;
                        (moonlets + k) -> vz -= G*mp*DZ/D3;
                        /******** Updating the acceleration of body p ********/
                        (moonlets + p) -> vx += G*mk*DX/D3;
                        (moonlets + p) -> vy += G*mk*DY/D3;
                        (moonlets + p) -> vz += G*mk*DZ/D3;           
                  }
            }       
      }
      #endif
      
      /******** Accelerations of the particles of the viscoelastic body ********/
      if (viscoelastic_bool){
            int j;
            typ xk, yk, zk, xp, yp, zp, L, r, dL;
            /******** Acceleration due to the springs ********/
            for (j = 0; j < N_connections; j ++){ //Travelling through all the connections
                  k  = (connections + j) -> Pair.fst;
                  p  = (connections + j) -> Pair.snd;
                  L  = (connections + j) -> rest_length;
                  if (exists[k] && exists[p]){
                        xk = (moonlets    + k) -> x;
                        yk = (moonlets    + k) -> y;
                        zk = (moonlets    + k) -> z;
                        mk = (moonlets    + k) -> mass;
                        xp = (moonlets    + p) -> x;
                        yp = (moonlets    + p) -> y;
                        zp = (moonlets    + p) -> z;
                        mp = (moonlets    + p) -> mass;
                        DX = xk - xp;  DY = yk - yp;  DZ = zk - zp;
                        r  = sqrt(DX*DX + DY*DY + DZ*DZ);
                        dL = r - L;
                        K  = spring_modulus*L*dL/r;
                        /******** Updating the accelerations ********/
                        (moonlets + k) -> vx -= K*DX/mk;
                        (moonlets + k) -> vy -= K*DY/mk;
                        (moonlets + k) -> vz -= K*DZ/mk;
                        (moonlets + p) -> vx += K*DX/mp;
                        (moonlets + p) -> vy += K*DY/mp;
                        (moonlets + p) -> vz += K*DZ/mp;
                  }
            }
            /******** Tidal acceleration due to the perturbator ********/
            if (pert_mass > 0.){ //The perturbing body affects the center of mass of the viscoelastic body
                  need_to_reduce_COM_bool = 1;  
                  typ cart[6];
                  typ nu = get_perturbing_true_anomaly(time_elapsed + t_init);              //Retrieving the true anomaly of the perturbing body
                  mu = G*(pert_mass + M_unit);
                  ell2cart(pert_sma, pert_ecc, pert_inc, nu, pert_aop, pert_lan, mu, cart); //Retrieving the coordinates  of the perturbing body
                  if (openGL_bool || (viscoelastic_bool && pert_mass > 0.)){                //For visualization in webGL or writing to files
                        CM.x  = cart[0] + 0.5*timestep*cart[3];  CM.y  = cart[1] + 0.5*timestep*cart[4];  CM.z  = cart[2] + 0.5*timestep*cart[5];
                        CM.vx = cart[3];                         CM.vy = cart[4];                         CM.vz = cart[5];
                        
                  }
                  mu = G*pert_mass;
                  XX = cart[0] + 0.5*timestep*cart[3];  YY = cart[1] + 0.5*timestep*cart[4];  ZZ = cart[2] + 0.5*timestep*cart[5];
                  Rk = sqrt(XX*XX + YY*YY + ZZ*ZZ);
                  aX = mu*XX/(Rk*Rk*Rk);  aY = mu*YY/(Rk*Rk*Rk);  aZ = mu*ZZ/(Rk*Rk*Rk); //Acceleration at the center of mass
                  for (j = 0; j <= largest_id; j ++){
                        if (exists[j]){
                              /******** Vector going from the current point in the viscoelastic body towards the perturbing body ********/
                              X  = XX - (moonlets + j) -> x;
                              Y  = YY - (moonlets + j) -> y;
                              Z  = ZZ - (moonlets + j) -> z;
                              Rp = sqrt(X*X + Y*Y + Z*Z);
                              /******** Updating the accelerations ********/
                              (moonlets + j) -> vx += mu*X/(Rp*Rp*Rp) - aX;
                              (moonlets + j) -> vy += mu*Y/(Rp*Rp*Rp) - aY;
                              (moonlets + j) -> vz += mu*Z/(Rp*Rp*Rp) - aZ;
                        }
                  }
            }
      }
}


void KelvinVoigtDamping(struct moonlet * X){

      /******** Computes the acceleration due to the dampers in the Kelvin-Voigt models of the viscoelastic body ********/
      /******** The array xx is updated. This function is only called (by "kick") when viscoelastic_bool is 1    ********/
      
      
      int j, k, p;
      typ xk, yk, zk, vxk, vyk, vzk, mk, xp, yp, zp, vxp, vyp, vzp, mp, L, drdv, K, dX, dY, dZ, dvX, dvY, dvZ, dr;

      for (j = 0; j < N_connections; j ++){ //Travelling through all the connections
            k    = (connections + j) -> Pair.fst;
            p    = (connections + j) -> Pair.snd;
            L    = (connections + j) -> rest_length;
            if (exists[k] && exists[p]){
                  xk   = (X + k) -> x;
                  yk   = (X + k) -> y;
                  zk   = (X + k) -> z;
                  vxk  = (X + k) -> vx;
                  vyk  = (X + k) -> vy;
                  vzk  = (X + k) -> vz;
                  mk   = (X + k) -> mass;
                  xp   = (X + p) -> x;
                  yp   = (X + p) -> y;
                  zp   = (X + p) -> z;
                  vxp  = (X + p) -> vx;
                  vyp  = (X + p) -> vy;
                  vzp  = (X + p) -> vz;
                  mp   = (X + p) -> mass;
                  /******** Updating the speeds ********/
                  vxk += 0.5*timestep*(xx + k) -> vx;
                  vyk += 0.5*timestep*(xx + k) -> vy;
                  vzk += 0.5*timestep*(xx + k) -> vz;
                  vxp += 0.5*timestep*(xx + p) -> vx;
                  vyp += 0.5*timestep*(xx + p) -> vy;
                  vzp += 0.5*timestep*(xx + p) -> vz;
                  dX   =  xk -  xp;  dY  =  yk -  yp;  dZ  =  zk -  zp;
                  dvX  = vxk - vxp;  dvY = vyk - vyp;  dvZ = vzk - vzp;
                  drdv = dX*dvX + dY*dvY + dZ*dvZ;
                  dr   = sqrt(dX*dX + dY*dY + dZ*dZ);
                  K    = damping_coefficient*L*drdv/(dr*dr);
                  /******** Updating the accelerations ********/
                  (xx + k) -> vx -= K*dX/mk;
                  (xx + k) -> vy -= K*dY/mk;
                  (xx + k) -> vz -= K*dZ/mk;
                  (xx + p) -> vx += K*dX/mp;
                  (xx + p) -> vy += K*dY/mp;
                  (xx + p) -> vz += K*dZ/mp;
            }
      }
}


void tides(struct moonlet * bodies){

      /******** When this function is called (by "kick"), the body array xx contains the acceleration without tides ********/
      /******** This function uses the knowledge of the acceleration without tides in order to compute the speed at ********/
      /******** the middle of the kick phase. Then, the acceleration due to tides can be computed and is added to   ********/
      /******** the acceleration without tides in the array xx                                                      ********/
      
      
      int j;
      typ X, Y, Z, vX, vY, vZ, m, K, r2, r10, rv, r_x_OmX, r_x_OmY, aX, aY, aZ;
      typ R5 = R_unit*R_unit*R_unit*R_unit*R_unit;
      typ A  = 3.0*k2*G*R5;
      typ M  = (inner_fluid_disk_bool ? CM.mass + fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) : CM.mass);
      
      for (j = 0; j <= largest_id; j ++){ //Looping over all bodies
            if (*(exists + j)){
                  X       = (bodies + j) -> x  - CM.x;
                  Y       = (bodies + j) -> y  - CM.y;
                  Z       = (bodies + j) -> z  - CM.z;
                  vX      = (bodies + j) -> vx - CM.vx + 0.5*timestep*((xx + j) -> vx - CM_acc[0]);
                  vY      = (bodies + j) -> vy - CM.vy + 0.5*timestep*((xx + j) -> vy - CM_acc[1]);
                  vZ      = (bodies + j) -> vz - CM.vz + 0.5*timestep*((xx + j) -> vz - CM_acc[2]);
                  m       = (bodies + j) -> mass;
                  r2      = X*X + Y*Y + Z*Z;
                  r10     = r2*r2*r2*r2*r2;
                  rv      = X*vX + Y*vY + Z*vZ;
                  r_x_OmX =  Y*SideralOmega; //r_cross_Omega
                  r_x_OmY = -X*SideralOmega;
                  K       = A*m/r10;
                  aX      = K*(r2*X + Delta_t*(2.0*rv*X + r2*(vX + r_x_OmX)));
                  aY      = K*(r2*Y + Delta_t*(2.0*rv*Y + r2*(vY + r_x_OmY)));
                  aZ      = K*(r2*Z + Delta_t*(2.0*rv*Z + r2*vZ));
                  /******** Updating the acceleration of body n° j ********/
                  (xx + j) -> vx -= aX;
                  (xx + j) -> vy -= aY;
                  (xx + j) -> vz -= aZ;
                  /******** Updating the acceleration of the central body ********/
                  CM_acc[0] += aX*m/M;
                  CM_acc[1] += aY*m/M;
                  CM_acc[2] += aZ*m/M;
                  /******** Updating the sideral rotation of the central body (to conserve angular momentum along Z) ********/
                  SideralOmega += timestep*m/(dimensionless_moi*M*R_unit*R_unit)*(X*aY - Y*aX);
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
      typ dr_dot_dv = (xa - xb)*(vx_a - vx_b) + (ya - yb)*(vy_a - vy_b) + (za - zb)*(vz_a - vz_b); //Dot product dv.dr where dr=r_b-r_a and dv=v_b-v_a are the relative position and speed
      typ alpha     = f*m_a*m_b/((m_a + m_b)*R*R);                                                 //Factor such that J=alpha*(dv.dr)*dr
      typ Jx        = alpha*dr_dot_dv*(xb - xa);                                                   //Vector J
      typ Jy        = alpha*dr_dot_dv*(yb - ya);
      typ Jz        = alpha*dr_dot_dv*(zb - za);     
      
      /******** Calculating the speeds after impact ********/
      vx_a += Jx/m_a;  //v_a = v_a + J/m_a
      vy_a += Jy/m_a;
      vz_a += Jz/m_a;
      vx_b -= Jx/m_b;  //v_b = v_b - J/m_b
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
      
      #if write_collisions_bool //Writing collision's data
      typ DeltaV = sqrt((vx_a - vx_b)*(vx_a - vx_b) + (vy_a - vy_b)*(vy_a - vy_b) + (vz_a - vz_b)*(vz_a - vz_b));
      if (write_to_files_bool){
            (collisionDatas + indexCollision) -> time         = t_init + time_elapsed + time_until_collision;
            (collisionDatas + indexCollision) -> m1           = m_b;
            (collisionDatas + indexCollision) -> m2           = m_a;
            (collisionDatas + indexCollision) -> R1           = R_b;
            (collisionDatas + indexCollision) -> R2           = R_a;
            (collisionDatas + indexCollision) -> DeltaV       = DeltaV;
            (collisionDatas + indexCollision) -> impact_angle = acos(-dr_dot_dv/(R*DeltaV));
            (collisionDatas + indexCollision) -> m_tilde      = m_a;
            indexCollision ++;
            if (indexCollision == 2*output_step*N_max){
                  fprintf(stderr, "Error: The array collisionDatas is not large enough.\n");
                  abort();
            }
      }
      #endif
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
      
      #if write_collisions_bool //Writing collision's data
      typ DeltaV = sqrt((vx_a - vx_b)*(vx_a - vx_b) + (vy_a - vy_b)*(vy_a - vy_b) + (vz_a - vz_b)*(vz_a - vz_b));
      if (write_to_files_bool && !fragmentation_bool){
            (collisionDatas + indexCollision) -> time         = t_init + time_elapsed + time_until_collision;
            (collisionDatas + indexCollision) -> m1           = m_b;
            (collisionDatas + indexCollision) -> m2           = m_a;
            (collisionDatas + indexCollision) -> R1           = R_b;
            (collisionDatas + indexCollision) -> R2           = R_a;
            (collisionDatas + indexCollision) -> DeltaV       = DeltaV;
            (collisionDatas + indexCollision) -> impact_angle = acos(-((xa - xb)*(vx_a - vx_b) + (ya - yb)*(vy_a - vy_b) + (za - zb)*(vz_a - vz_b))/((R_a + R_b)*DeltaV));
            (collisionDatas + indexCollision) -> m_tilde      = m; 
            indexCollision ++;
            if (indexCollision == 2*output_step*N_max){
                  fprintf(stderr, "Error: The array collisionDatas is not large enough.\n");
                  abort();
            }
      }
      #endif
}

void fragmentation(struct moonlet * moonlets, int a, int b){

      /******** Treats the fragmentation due to the collision between bodies a and b ********/
      
      typ stigma = (3.0*mu_parameter-1.0)/(3.0*mu_parameter);
      typ vx_a, vy_a, vz_a, vx_b, vy_b, vz_b;                  //Cartesian speeds of the bodies.
      typ xa, ya, za, xb, yb, zb;                              //Cartesian positions at the collision
      typ dx, dy, dz, dvx, dvy, dvz;                           //dx is xa-xb, and so on.
      typ R_2 = (moonlets + a) -> radius;                      //The bodies' radii
      typ R_1 = (moonlets + b) -> radius;
      typ R   = R_1 + R_2;                                     //Sum of the radii;
      typ m_2 = (moonlets + a) -> mass;                        //The bodies's masses
      typ m_1 = (moonlets + b) -> mass;
      typ M   = m_1 + m_2;                                     //Sum of the masses
      typ rho_1 = 3.0*m_1/(4.0*M_PI*R_1*R_1*R_1);              //Densities
      typ rho_2 = 3.0*m_2/(4.0*M_PI*R_2*R_2*R_2);
      typ average_density = (m_1*rho_1 + m_2*rho_2)/M;         //The average density of the moonlets, weighted by mass
      typ rho_12 = rho_1/rho_2;
      typ rho_12_3nu = pow(rho_12, 3.0*nu_parameter);          //Ratio between the density of the impactor and that of the target at the power 3*nu           
      
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
      
      typ dr_dot_dv   = dx*dvx + dy*dvy + dz*dvz;                          //  (r_1 - r_2).(v_1 - v_2)
      typ dv_norm     = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);                 //  ||v_1 - v_2||
      typ costheta    = -dr_dot_dv/(R*dv_norm);                            //Cosine of impact angle;
      typ costheta3mu = pow(costheta, 3.0*mu_parameter);                   //Cosine of impact angle at the power 3*mu;
      typ R_eq        = pow(3.0*M/(4.0*M_PI*average_density), 1.0/3.0);    //Hypothetical radius if the impactor were to merge on the target
      typ vesc        = sqrt(2.0*G*M/R_eq);                                //Escape velocity at the surface if the impactor were to merge on the target
      typ K           = 3.0*k_parameter*m_1/(4.0*M_PI)/rho_12*costheta3mu;
      typ m_check     = K*(pow(C1_parameter*dv_norm/vesc, 3.0*mu_parameter)*rho_12_3nu - 1.0); //Mass of the tail : Ejected mass
      typ m_tilde     = M - m_check;                                       //Mass of the largest fragment
      typ m_tilde_2   = m_check / (typ) N_tilde;                           //Mass of the tail's fragments
           
      /******** Some data used for the fragmentation case ********/
      typ v_cm[3], r_cm[3]; //Velocity and position of the center of mass before impact
      typ r_tilde[3];       //Position of the largest fragment
      typ v_tilde[3];       //Speed of the largest fragment
      v_cm[0] = (m_2*vx_a + m_1*vx_b)/M;  v_cm[1] = (m_2*vy_a + m_1*vy_b)/M;  v_cm[2] = (m_2*vz_a + m_1*vz_b)/M;
      r_cm[0] = (m_2*xa   + m_1*  xb)/M;  r_cm[1] = (m_2*ya   + m_1*  yb)/M;  r_cm[2] = (m_2*za   + m_1*  zb)/M;
      
      #if write_collisions_bool //Writing collision's data
      if (write_to_files_bool){
            (collisionDatas + indexCollision) -> time         = t_init + time_elapsed + time_until_collision;
            (collisionDatas + indexCollision) -> m1           = m_1;
            (collisionDatas + indexCollision) -> m2           = m_2;
            (collisionDatas + indexCollision) -> R1           = R_1;
            (collisionDatas + indexCollision) -> R2           = R_2;
            (collisionDatas + indexCollision) -> DeltaV       = dv_norm;
            (collisionDatas + indexCollision) -> impact_angle = acos(costheta);
            (collisionDatas + indexCollision) -> m_tilde      = m_tilde;
            indexCollision ++;
            if (indexCollision == 2*output_step*N_max){
                  fprintf(stderr, "Error: The array collisionDatas is not large enough.\n");
                  abort();
            }
      }
      #endif
      
      /******** Merger case. If there is no ejecta (v_max <= vesc) then bodies a and b merge together ********/
      if (m_check <= merging_threshold*M){ //No ejecta or the ejected mass is less than merging_threshold*(m1 + m2)
            merger(moonlets, a, b);
            merger_count ++;
            return;
      }
      
      /******** Super-catastrophic fragmentation. The mass of the largest fragment is less than 10 % of the total mass. The ejecta is discarded ********/
      if (m_tilde < 0.1*M){ //m_tilde is not proportionnal to Qr/Qr* in this regime but is proportionnal to (Qr/Qr*)^(-3/2)
      
            m_tilde = 0.1*M*pow(1.0/0.9*m_check/M, -1.5); //Eq. (44) of Leinhardt and Stewart (2012)
            
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
            lose_moonlet(b);
            *(did_collide + a) = one_collision_only_bool;
            super_catastrophic_count ++;
            /******** The discarded mass is put in the inner fluid disk to prevent mass from just vanishing, or in the central body ********/
            typ discarded_mass = M - m_tilde;
            if (m_tilde < fragment_threshold/4.0){ //Body a is discarded as well
                  lose_moonlet(a);
                  discarded_mass = M;
            }
            if (inner_fluid_disk_bool){
                  fluid_disk_Sigma += discarded_mass/(M_PI*(Rout*Rout - R_unit*R_unit));
            }
            else if (central_mass_bool){
                  CM.mass += discarded_mass;
            }
            
            #if write_collisions_bool //Writing collision's data
            if (write_to_files_bool){
                  (collisionDatas + indexCollision - 1) -> m_tilde = m_tilde;
            }
            #endif
            
            return;
      }
      
      
      /******** Partial fragmentation. The tail is reunited into a single body. ********/ 
      if (m_tilde_2 < fragment_threshold && m_check <= m_tilde){
            typ r_k[3]; //Position of the tail with respect to the largest fragment
            typ v_k[3]; //Velocity of the tail with respect to the largest fragment
            typ v_k_scalar = vesc/stigma*(m_check + K)/m_check*(1.0 - pow(K/(K + m_check), stigma)); //Scalar velocity of the tail with respect to the largest fragment
            if (v_k_scalar < vesc){
                  fprintf(stderr, "Error: In function fragmentation, the ejection velocity cannot be less than the escape velocity.\n");
                  abort();
            }
            typ in_front_of = v_k_scalar/R;
            r_k[0] = -dx;  r_k[1] = -dy;  r_k[2] = -dz;
            v_k[0] = in_front_of*r_k[0];  v_k[1] = in_front_of*r_k[1];  v_k[2] = in_front_of*r_k[2];
            typ R_tilde   = pow(3.0*m_tilde/(4.0*M_PI*average_density), 1.0/3.0); //Radius of the largest fragment
            typ R_check   = pow(3.0*m_check/(4.0*M_PI*average_density), 1.0/3.0); //Radius of the fragment of the tail
            typ R_tilde12 = R_tilde + R_check;
            typ dr_norm   = sqrt(r_k[0]*r_k[0] + r_k[1]*r_k[1] + r_k[2]*r_k[2]);
            r_k[0] *= R_tilde12/dr_norm;  r_k[1] *= R_tilde12/dr_norm;  r_k[2] *= R_tilde12/dr_norm;

            /******** Total momentum is conserved ********/
            r_tilde[0] = r_cm[0] - m_check/M*r_k[0];
            r_tilde[1] = r_cm[1] - m_check/M*r_k[1];
            r_tilde[2] = r_cm[2] - m_check/M*r_k[2];
            v_tilde[0] = v_cm[0] - m_check/M*v_k[0];
            v_tilde[1] = v_cm[1] - m_check/M*v_k[1];
            v_tilde[2] = v_cm[2] - m_check/M*v_k[2];  
            
            /******** Actualizing the bodies ********/
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
      }
      
      
      /******** Full fragmentation. The tail is made up of N_tilde bodies ********/
      else{
            /******** Determination of the r_k' and v_k' ********/ 
            typ r_k[3*N_tilde]; //Position of the fragments of the tail with respect to the largest fragment
            typ v_k[3*N_tilde]; //Velocity of the fragments of the tail with respect to the largest fragment
            typ dr[3];          //r_1-r_2
            typ dv[3];          //v_1-v_2
            typ R_tilde   = pow(3.0*m_tilde  /(4.0*M_PI*average_density), 1.0/3.0); //Radius of the largest fragment
            typ R_tilde_2 = pow(3.0*m_tilde_2/(4.0*M_PI*average_density), 1.0/3.0); //Radius of the fragments of the tail
            typ R_tilde12 = R_tilde + R_tilde_2;
            typ u[3]; //Vectors used to locate the fragments of the tail
            typ v[3];
            dr[0] = -dx;   dr[1] = -dy;   dr[2] = -dz;
            dv[0] = -dvx;  dv[1] = -dvy;  dv[2] = -dvz;
            typ dr_norm = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
            dr[0] *= R_tilde12/dr_norm;  dr[1] *= R_tilde12/dr_norm;  dr[2] *= R_tilde12/dr_norm;
            typ dr_x_dv[3]; // dr x dv
            typ dr_x_dv_norm; // ||dr x dv||
            int pq[4] = pq_min_max; // {p_k_min, p_k_max, q_k_min, q_k_max}
            int p, q;
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
            u[0] = v_x_dr[0]/R_tilde12;  u[1] = v_x_dr[1]/R_tilde12;  u[2] = v_x_dr[2]/R_tilde12; //Defining the unit vector u
            int n = 0;
            typ v_k_scalar, two_p_R_tilde_2, two_q_R_tilde_2, in_front_of, zk_stg, zk1_stg;
            for (p = pq[0]; p <= pq[1]; p ++){ //I travel along the rectangle of integer coordinate points to define the position and speeds of the fragments of the tail
                  for (q = pq[2]; q <= pq[3]; q ++){
                        two_p_R_tilde_2 = ((typ) (2*p))*R_tilde_2;  two_q_R_tilde_2 = ((typ) (2*q))*R_tilde_2;
                        r_k[3*n]     = dr[0] + two_p_R_tilde_2*u[0] + two_q_R_tilde_2*v[0]; // x-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        r_k[3*n + 1] = dr[1] + two_p_R_tilde_2*u[1] + two_q_R_tilde_2*v[1]; // y-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        r_k[3*n + 2] = dr[2] + two_p_R_tilde_2*u[2] + two_q_R_tilde_2*v[2]; // z-coordinate of the position of the (n+1)^th fragment of the tail wrt the largest fragment
                        zk_stg       = pow(1.0 - ((typ) n)     *m_tilde_2/(K + m_check), stigma);
                        zk1_stg      = pow(1.0 - ((typ) n + 1.)*m_tilde_2/(K + m_check), stigma);
                        v_k_scalar   = vesc/stigma*(m_check + K)/m_tilde_2*(zk_stg - zk1_stg); //Scalar velocity of the tail's fragment
                        if (v_k_scalar < vesc){
                              fprintf(stderr, "Error: In function fragmentation, the velocity of a fragment cannot be less than the escape velocity.\n");
                              abort();
                        }
                        in_front_of  = v_k_scalar/sqrt(1.0 + ((typ) p)*((typ) p) + ((typ) q)*((typ) q));
                        v_k[3*n]     = in_front_of*(dr[0]/R_tilde12 + ((typ) p)*u[0] + ((typ) q)*v[0]);
                        v_k[3*n + 1] = in_front_of*(dr[1]/R_tilde12 + ((typ) p)*u[1] + ((typ) q)*v[1]);
                        v_k[3*n + 2] = in_front_of*(dr[2]/R_tilde12 + ((typ) p)*u[2] + ((typ) q)*v[2]);
                        n ++;
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
                        id[n] = get_free_index(1);     //Putting the remaining N_tilde fragments of the tail there
                  }
                  else{                                //No need to put the fragments at the end
                        id[n] = get_free_index(0);     //Putting the remaining N_tilde fragments of the tail there
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
            (moonlets + id[0]) -> radius = R_tilde;

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


int willMergeWithDisk(struct moonlet * bodies, int id){

      /******** This function is called if body n° id goes below the disruption_threshold.                ********/
      /******** It determines if it will collide with the central body or merge with the inner fluid disk.********/
      /******** If the periapsis is above the surface or if the body will cross the xy plane before       ********/
      /******** hitting the surface, then 1 is returned. Otherwise, 0 is returned.                        ********/
      
      typ alkhqp[6];
      typ cart[6];
      typ a, k, h, q, p, e, i, nu, omega, Omega, mu, Z, z, disk_angular_momentum;
      
      /******** Retrieving the orbital elements of the bodies ********/
      mu = G*(CM.mass + fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) + (bodies + id) -> mass);
      cart2ell(bodies, id, alkhqp, mu);
      a     = *alkhqp;
      k     = alkhqp[2];
      h     = alkhqp[3];
      e     = sqrt(k*k + h*h);
      q     = alkhqp[4];
      p     = alkhqp[5];
      i     = 2.0*asin(sqrt(q*q + p*p));
      if (a*(1.0 - e) > R_unit){ //The periapsis is above the surface
            disk_angular_momentum = 0.8*fluid_disk_Sigma*M_PI*(Rout*Rout*sqrt(Rout) - R_unit*R_unit*sqrt(R_unit));
            innerFluidDiskAngularMomentum((bodies + id) -> mass, a, e, cos(i), disk_angular_momentum, fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit));
            return 1;
      }
      Omega = atan2(p, q);
      omega = atan2(h, k) - Omega;
      
      /******** Getting the true anomaly of the intersection of the orbit with the surface ********/
      nu = -acos((a*(1 - e*e) - R_unit)/(e*R_unit));

      /******** Getting the cartesian coordinates of the intersection of the orbit with the surface ********/
      ell2cart(a, e, i, nu, omega, Omega, mu, cart);

      /******** The body will cross the xy plane before hitting the surface if, and only if, its Z coordinate ********/
      /******** has a different sign than the z coordinate of the intersection of the orbit with the surface  ********/
      Z = (bodies + id) -> z - CM.z;
      z = cart[2];
      if (z*Z < 0.){
            disk_angular_momentum = 0.8*fluid_disk_Sigma*M_PI*(Rout*Rout*sqrt(Rout) - R_unit*R_unit*sqrt(R_unit));
            innerFluidDiskAngularMomentum((bodies + id) -> mass, a, e, cos(i), disk_angular_momentum, fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit));
            return 1;
      }
      return 0;
}


void innerFluidDiskAngularMomentum(typ m, typ a, typ e, typ cosi, typ g1, typ m1){

      /******** Updates Rout to the new outer radius of the inner fluid disk so that the angular ********/
      /******** momentum is conserved. Called when a body spawns or merges with the inner disk   ********/
      /******** The body's mass m is given positive (resp. negative) for a merge (resp. a spawn) ********/
      /******** There is no analytical expression for Rout so a Newton-Raphson method is used.   ********/
      /******** The mass per unit area of the inner fluid disk is also updated                   ********/
      
      typ g         = m*cosi*sqrt(a*(1. - e*e)); //z-component of the angular momentum of the body that spawned/merged 
      typ X         = sqrt(Rout);                //We take the current outer edge as initial condition
      typ precision = 1.;
      int n_steps   = 0;
      typ mm1       = 0.8*(m + m1);
      typ gg1       = g + g1;
      typ X0        = sqrt(R_unit); typ X02 = X0*X0; typ X03 = X02*X0; typ X04 = X02*X02;
      typ fX, dfX, X2, X3, X4, dX;
      
      while(precision > 1.e-6){
            n_steps ++;
            X2        = X*X;  X3 = X2*X;  X4 = X2*X2;
            fX        = mm1*(X4 + X0*X3 + X02*X2 + X03*X + X04) - gg1*(X + X0)*(X2 + X02);
            dfX       = mm1*(4.*X3 + 3.*X0*X2 + 2.*X02*X + X03) - gg1*(3.*X2 + X02 + 2.*X*X0);
            dX        = -fX/dfX;
            X        += dX;
            precision = fabs(dX/X);
            
            /******** If the method gets lost, I try to put it back on track ********/
            if (X < 0. || X != X){
                  X = rdm(0., 50.*sqrt(R_unit));
            }
            
            /******** Checking the convergence of the Newton-Raphson method. To be removed when the code is robust ********/
            if (n_steps >= 150){
                  printf("Warning : The Newton-Raphson method to find the inner fluid disk's outer edge does not converge.\n");
                  fluid_disk_Sigma = (m + m1)/(M_PI*(Rout*Rout - R_unit*R_unit));
                  return;
            }
      }
      
      Rout = X*X;
      Rout = Rout < 1.1*disruption_threshold ? 1.1*disruption_threshold : Rout;
      fluid_disk_Sigma = (m + m1)/(M_PI*(Rout*Rout - R_unit*R_unit));
}
