/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    structure.c                                                 ********/
/******** @brief   Miscellaneous structural implementations                    ********/
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
#include "structure.h"
#include "parameters.h"
#include "ffm.h"
#include "display.h"
#include "spring.h"
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>


/******** Declaring external array ********/
struct moonlet * xx;
int * exists;
struct chain ** hash;
int * modified_cells;
int * indexes;
struct chain * nghb;
int * free_indexes;
struct pair * pairs;
struct moonlet * three_largest;
int * three_largest_indexes;
typ * tam_loss;
typ * approach;
int * did_collide;
int * already_in_tree;
typ CM_acc[3] = {0., 0., 0.};
typ * sending_buffer;


/******** Declaring external variables ********/
struct moonlet CM        = {0., 0., 0., 0., 0., 0., M_unit, R_unit};
struct moonlet CM_buffer = {0., 0., 0., 0., 0., 0., M_unit, R_unit};
typ J2;
typ timestep;
typ gam;
typ gam_min;
int how_many_big;
int how_many_small;
int largest_id;
int how_many_modified;
typ collision_cube;
int total_neighbours;
typ average_neighbours;
int first_passage;
typ time_elapsed;
int how_many_free;
int how_many_moonlets;
int force_naive_bool;
typ time_until_collision;
int how_many_pairs;
int super_catastrophic_count;
int half_fragmentation_count;
int full_fragmentation_count;
int merger_count;
int collision_count;
typ fluid_disk_Sigma;
typ flowed_since_last_spawn;
typ Rout;
typ SideralOmega;
typ star_mean_motion;
typ evection_resonance;
int need_to_reduce_COM_bool;
typ previous_tra;
typ pert_M;


void ell2cart(typ a, typ e, typ i, typ nu, typ omega, typ Omega, typ mu, typ * cart){

      /******** Returns the array [X,Y,Z,vX,vY,vZ] of the cartesian coordinates                               ********/
      /******** a is the semi-major axis, e is the eccentricity, i is the inclination, nu is the true anomaly,********/
      /******** omega is the argument of periapsis and Omega is the longitude of the ascending node           ********/
      
      typ X,Y,Z,vX,vY,vZ;                //Cartesian coordinates
      typ X_buff,Y_buff,vX_buff,vY_buff; //Buffer for cartesian coordinates
      typ r;                             //Body's distance to Earth's center
      typ g;                             //Angular momentum per unit mass
      typ dnudt;
      typ drdt;
      typ cosnu    = cos(nu);
      typ sinnu    = sin(nu);
      typ cosvarpi = cos(omega + Omega);
      typ sinvarpi = sin(omega + Omega);
      typ q        = sin(i/2.)*cos(Omega);
      typ p        = sin(i/2.)*sin(Omega);
      typ chi      = cos(i/2.);
      typ pp       = 1. - 2.*p*p;
      typ qq       = 1. - 2.*q*q;
      typ dpq      = 2.*p*q;

      /******** In the orbital plane (see e.g. Laskar & Robutel 1995) ********/
      r = a*(1. - e*e)/(1. + e*cosnu);
      if (J2_bool && central_mass_bool){ //Using the geometric elliptical elements instead of the osculating ones in case of an oblate Earth (See Greenberg, 1981)
            mu = mu*(1. + 1.5*J2*R_unit*R_unit/(a*a));
      }
      g       = sqrt(mu*a*(1. - e*e));
      dnudt   = g/(r*r);
      drdt    = a*e*dnudt*sinnu*(1. - e*e)/((1. + e*cosnu)*(1. + e*cosnu));
      X_buff  = r*cosnu;
      Y_buff  = r*sinnu;
      vX_buff = drdt*cosnu - r*dnudt*sinnu;
      vY_buff = drdt*sinnu + r*dnudt*cosnu;
      
      /******** Rotations to convert to reference plane (see e.g. Laskar & Robutel 1995) ********/      
      X  =  X_buff*(pp*cosvarpi + dpq*sinvarpi) +  Y_buff*(dpq*cosvarpi - pp*sinvarpi);
      vX = vX_buff*(pp*cosvarpi + dpq*sinvarpi) + vY_buff*(dpq*cosvarpi - pp*sinvarpi);
      Y  =  X_buff*(qq*sinvarpi + dpq*cosvarpi) +  Y_buff*(qq*cosvarpi  - dpq*sinvarpi);
      vY = vX_buff*(qq*sinvarpi + dpq*cosvarpi) + vY_buff*(qq*cosvarpi  - dpq*sinvarpi);
      Z  =  X_buff*(2.*q*chi*sinvarpi - 2.*p*chi*cosvarpi) +  Y_buff*(2.*p*chi*sinvarpi + 2.*q*chi*cosvarpi);
      vZ = vX_buff*(2.*q*chi*sinvarpi - 2.*p*chi*cosvarpi) + vY_buff*(2.*p*chi*sinvarpi + 2.*q*chi*cosvarpi);
      
      /******** Writing the cartesian coordinates ********/
      * cart      = X;
      *(cart + 1) = Y;
      *(cart + 2) = Z;
      *(cart + 3) = vX;
      *(cart + 4) = vY;
      *(cart + 5) = vZ;
}


void cart2ell(struct moonlet * moonlets, int id, typ * alkhqp, typ mu){

      /******** Computes the elliptic elements a, l, k, h, q & p where a is the semi-major axis and l is the mean      ********/
      /******** longitude. k, h, q & p are defined as k + ih = e*exp(i*varpi) and q + ip = sin(I/2)*exp(i*Omega). The  ********/
      /******** eccentricity e, inclination I, longitude of the pericentre varpi and longitude of the ascending node   ********/
      /******** Omega are straightforward to obtain from k, h, q & p, but these variables are better since they are    ********/
      /******** regular even et zero eccentricity and inclination. We have e = sqrt(k^2+h^2), sin(I/2) = sqrt(q^2+p^2),********/
      /******** varpi = atan2(h,k) and Omega = atan2(p,q). The cartesian coordinates are those of body n° id           ********/
      /******** The elliptic elements a, l, k, h, q & p are written in the vector alkhqp. This vector is overwritten   ********/
      /******** This function originally comes from the IMCCE lab (Paris Observatory) and was adapted for NcorpiON     ********/

      typ X, Y, Z, vX, vY, vZ, R, V2, RV, C1, C2, C3, DC, CC, AA, aux0;
      typ a11, a12, a21, a22, c11, c12, c21, c22, K1, H1, K2, H2, K, H;
      typ SMU, FAC1, FAC2, b12, b22, sinF, cosF, F, USQA, aux1;

      /******** Getting the cartesian coordinates ********/
      X  = (moonlets + id) -> x;
      Y  = (moonlets + id) -> y;
      Z  = (moonlets + id) -> z;
      vX = (moonlets + id) -> vx;
      vY = (moonlets + id) -> vy;
      vZ = (moonlets + id) -> vz;
      X  -= central_mass_bool ? CM_buffer.x  : 0.;
      Y  -= central_mass_bool ? CM_buffer.y  : 0.;
      Z  -= central_mass_bool ? CM_buffer.z  : 0.;
      vX -= central_mass_bool ? CM_buffer.vx : 0.;
      vY -= central_mass_bool ? CM_buffer.vy : 0.;
      vZ -= central_mass_bool ? CM_buffer.vz : 0.;

      /******** Computing the semi-major axis ********/
      R   = sqrt(X*X + Y*Y + Z*Z);
      V2  = vX*vX + vY*vY + vZ*vZ;
      AA  = R*mu/(2.0*mu - R*V2); //Division by zero if the trajectory is perfectly parabolic.
      if (J2_bool && AA > 0.0 && central_mass_bool){  //If the Earth is oblate and the trajectory is elliptic, correcting Kepler third law (See Greenberg, 1981)
            mu *= 1.0 + 1.5*J2*R_unit*R_unit/(AA*AA);
            AA  = R*mu/(2.0*mu - R*V2);
      }
      *alkhqp = AA;

      /******** Normalizing the velocities (Adopting the convention of J. Laskar's 2004 lectures notes) ********/
      SMU = sqrt(mu);
      vX /= SMU;
      vY /= SMU;
      vZ /= SMU;

      /******** Computing the angular momentum ********/
      V2 = vX*vX + vY*vY + vZ*vZ;
      RV = X*vX  + Y*vY  + Z*vZ;
      C1 = Y*vZ  - Z*vY;
      C2 = Z*vX  - X*vZ;
      C3 = X*vY  - Y*vX;
      CC = C1*C1 + C2*C2 + C3*C3;
      DC = sqrt(CC);

      /******** Computing (q, p) ********/
      aux0          = sqrt(2.0*(CC + DC*C3));
      *(alkhqp + 4) = (aux0 == 0. ? 0. : -C2/aux0);
      *(alkhqp + 5) = (aux0 == 0. ? 0. :  C1/aux0);

      if (R == 0. || V2 == 0.){
            K = 0.;
            H = 0.;
      }
      else{
            /******** Computing the matrix coefficients needed for (k, h) ********/
            a11 = V2 - 1.0/R;
            a12 = RV/(R*DC);
            a21 = -RV;
            a22 = DC - R/DC;

            /******** Computing (k, h) ********/
            c11  =  X*a11  + vX*a21;
            c12  =  X*a12  + vX*a22;
            c21  =  Y*a11  + vY*a21;
            c22  =  Y*a12  + vY*a22;
            FAC1 =  C1/(DC + C3);
            FAC2 =  C2/(DC + C3);
            K1   =  c11 - FAC1*(Z*a11 + vZ*a21);
            H1   = -c12 + FAC1*(Z*a12 + vZ*a22);
            H2   =  c21 - FAC2*(Z*a11 + vZ*a21); //Should be equal to H1
            K2   =  c22 - FAC2*(Z*a12 + vZ*a22); //Should be equal to K1
            if (fabs(H1 - H2) + fabs(K1 - K2) > 1.0e-6){
                  printf("Warning : Bad computation of (k,h) in function cart2ell. (K1 - K2, H1 - H2) = (%.13lf, %.13lf)\n", K1 - K2, H1 - H2);
            }
            K    =  0.5*(K1 + K2);
            H    =  0.5*(H1 + H2);
      }
      *(alkhqp + 2) = K;
      *(alkhqp + 3) = H;

      if (R == 0. || V2 == 0.){
            *(alkhqp + 1) = 0.;
      }
      else{
            /******** Computing the mean longitude l = M + varpi ********/
            USQA = sqrt(2.0/R - V2);
            if ((USQA) >= 0.0){ //Elliptic case
                  b12  = vX - FAC1*vZ;
                  b22  = vY - FAC2*vZ;
                  aux1 = (R*V2 - 1.0)/(1.0 + DC*USQA);
                  sinF = -b12*R*USQA + H*aux1;
                  cosF =  b22*R*USQA + K*aux1;
                  F    =  atan2(sinF, cosF);
                  *(alkhqp + 1) = F - RV*USQA;
            }
            else{ //Hyperbolic case
                  USQA  = sqrt(-(2.0/R - V2));
                  typ E = atanh(RV*USQA/(R*V2 - 1.0));
                  typ M = sqrt(K*K + H*H) * sinh(E) - E;
                  *(alkhqp + 1) = M + atan2(H, K);
            }
      }
}


void cart2aei(struct moonlet * moonlets, int id, typ * aei){

      /******** Modifies the array [a,e,i] of the semimajor axis, the eccentricity and the inclination ********/
      /******** of the body whose id is id. Stores these three quantities in the array aei.            ********/
      /******** This function is currently not used by NcorpiON. To be removed eventually              ********/
      
      
      typ r,v2,r_cross_v,r_cross_v_square,X,Y,Z,vX,vY,vZ,a,e,i,cosi,r_cross_v_1,r_cross_v_2,r_cross_v_3,mu;
      
      /******** Getting the cartesian coordinates ********/
      X  = (moonlets + id) -> x;
      Y  = (moonlets + id) -> y;
      Z  = (moonlets + id) -> z;
      vX = (moonlets + id) -> vx;
      vY = (moonlets + id) -> vy;
      vZ = (moonlets + id) -> vz;
      
      /******** Getting the semimajor axis from the orbital energy ********/
      r  = sqrt(X*X + Y*Y + Z*Z);
      v2 = vX*vX + vY*vY + vZ*vZ;
      mu = central_mass_bool ? G*CM.mass : G*M_unit;
      a  = mu*r/(2.0*mu - r*v2);
      
      /******** If the Earth is oblate, correcting Kepler third law (See Greenberg, 1981) ********/
      if (J2_bool && central_mass_bool){
            mu = mu*(1.0 + 1.5*J2*R_unit*R_unit/(a*a));
            a  = mu*r/(2.0*mu - r*v2);
      }
      
      /******** Getting the eccentricity from the momentum ********/
      r_cross_v_1 = Y*vZ - vY*Z;
      r_cross_v_2 = vX*Z - X*vZ;
      r_cross_v_3 = X*vY - vX*Y;
      r_cross_v_square = r_cross_v_1*r_cross_v_1 + r_cross_v_2*r_cross_v_2 + r_cross_v_3*r_cross_v_3;
      e = sqrt(1.0 - r_cross_v_square/(mu*a));
      
      /******** Getting the inclination from the third component of the angular momentum ********/
      r_cross_v = sqrt(r_cross_v_square);
      cosi = r_cross_v_3/r_cross_v;
      i = acos(cosi);
      
      /******** Filling the array aei ********/
      * aei      = a;
      *(aei + 1) = e;
      *(aei + 2) = i; 

}


/******************************************/
/******** Initialization functions ********/
/******************************************/


struct moonlet init(typ a, typ e, typ i, typ nu, typ omega, typ Omega, typ density, typ rad){


      /******** a is the semi-major axis, e is the eccentricity, i is the inclination, nu is the true anomaly,     ********/
      /******** omega is the argument of periapsis, Omega is the longitude of the ascending node and m is the mass ********/

      struct moonlet mlt;
      mlt.radius = rad;                                 //Initializing the radius
      typ m      = 4.0/3.0*M_PI*density*rad*rad*rad;    //Defining the mass
      mlt.mass   = m;                                   //Initializing the mass
      typ mu     = central_mass_bool ? CM.mass : M_unit;  mu += inner_fluid_disk_bool ? fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) : 0.;
      mu        += m;  mu *= G;
      
      typ cart[6] = {a, e, i, nu, omega, Omega};
      ell2cart(a, e, i, nu, omega, Omega, mu , cart);   //Computing the cartesian coordinates from the orbital elements
      mlt.x  = * cart;                                  //Initializing the cartesian coordinates
      mlt.y  = *(cart + 1);
      mlt.z  = *(cart + 2);
      mlt.vx = *(cart + 3);
      mlt.vy = *(cart + 4);
      mlt.vz = *(cart + 5);
      
      return mlt;
}


typ rdm(typ min, typ max){
      
      /******** Returns a random number between min and max according to a uniform distribution ********/

      /******** Generating the random number ********/
      typ MyRand = ((typ) rand())/((typ) RAND_MAX); //between 0 and 1
      MyRand = min + (max - min)*MyRand;            //between min and max

      return MyRand;
}


struct moonlet * populate(){

      /******** Populates the simulation with N_0 bodies. If random_initial_bool is 1, then the initial conditions are chosen ********/
      /******** uniformly at random between bounds defined in the parameter file. Otherwise, the initial conditions are read  ********/
      /******** in the file pth/init.txt, where the path pth is defined in the parameter file.                                ********/

      struct moonlet * moonlets = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //The array of bodies
      if (moonlets == NULL){
            fprintf(stderr, "Error : Cannot allocate array of bodies in function populate.\n");
            abort();
      }
      
      /******** Populating the body array ********/
      typ a, e, i, nu, omega, Omega, rad, density, m;
      int k;
      
      if (random_initial_bool){ //The initial conditions are chosen at random
            for (k = 0; k < N_0; k ++){
                  a               = rdm(sma_min, sma_max);
                  e               = rdm(eccentricity_min, eccentricity_max);
                  i               = rdm(inclination_min, inclination_max);
                  nu              = rdm(0.0, 2.0*M_PI);
                  omega           = rdm(0.0, 2.0*M_PI);
                  Omega           = rdm(0.0, 2.0*M_PI);
                  rad             = rdm(radius_min, radius_max);
                  density         = rdm(density_min, density_max);
                  *(moonlets + k) = init(a, e, i, nu, omega, Omega, density, rad);
            }
      }
      else{                     //The initial conditions are read from the file pth/init.txt where pth is defined in the parameter file
            char fileOfIC[800]; 
            strcpy(fileOfIC, pth);
            strcat(fileOfIC, "init.txt");
            typ * IC = (typ *)malloc(8*N_0*sizeof(typ));
            if (IC == NULL){
                  fprintf(stderr, "Error : Cannot allocate array for initial conditions in function populate.\n");
                  abort();
            }
            readFromFile(fileOfIC, IC, 8*N_0);
            for (k = 0; k < N_0; k ++){
                  if (initial_cartesian_bool){ //Initial conditions are given in cartesian coordinates
                        (moonlets + k) -> x      = *(IC + 8*k);
                        (moonlets + k) -> y      = *(IC + 8*k + 1);
                        (moonlets + k) -> z      = *(IC + 8*k + 2);
                        (moonlets + k) -> vx     = *(IC + 8*k + 3);
                        (moonlets + k) -> vy     = *(IC + 8*k + 4);
                        (moonlets + k) -> vz     = *(IC + 8*k + 5);
                        (moonlets + k) -> mass   = *(IC + 8*k + 6);
                        (moonlets + k) -> radius = *(IC + 8*k + 7);
                  }
                  else{                       //Initial conditions are given in elliptic elements
                        a               = *(IC + 8*k);
                        e               = *(IC + 8*k + 1);
                        i               = *(IC + 8*k + 2);
                        nu              = *(IC + 8*k + 3);
                        omega           = *(IC + 8*k + 4);
                        Omega           = *(IC + 8*k + 5);
                        m               = *(IC + 8*k + 6);
                        rad             = *(IC + 8*k + 7);
                        density         = 3.0*m/(4.0*M_PI*rad*rad*rad);
                        *(moonlets + k) = init(a, e, i, nu, omega, Omega, density, rad);
                        if (*(IC + 8*k + 1) < 0.){
                              printf("Warning : A body with negative eccentricity was given in the initial conditions. Discarding the body.\n");
                              *(exists + k) = 0;
                        }
                  }
                  if (*(IC + 8*k + 6) == 0.){
                        printf("Warning : NcorpiON does not support bodies with zero mass at the moment. Discarding the body.\n");
                        *(exists + k) = 0;
                  }
                  if (*(IC + 8*k + 6) < 0.){
                        printf("Warning : A body with negative mass was given in the initial conditions.\n");
                  }
                  if (*(IC + 8*k + 7) < 0.){
                        printf("Warning : A body with negative radius was given in the initial conditions. Discarding the body.\n");
                        *(exists + k) = 0;
                  }
            }
            free(IC);
            IC = NULL;
      }
      
      /******** Filling the unused cells of the body array with whatever ********/
      for (k = N_0; k < N_max; k ++){
            *(moonlets + k) = init(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      }
      
      /******** Reducing to the center of mass ********/
      if (reduce_to_COM_bool){
            typ com[6];
            total_momentum(moonlets, com); //Getting the speed and position of the center of mass
            for (k = 0; k < N_0; k ++){    //Substracting the center of mass
                  (moonlets + k) -> x  -= com[0];
                  (moonlets + k) -> y  -= com[1];
                  (moonlets + k) -> z  -= com[2];
                  (moonlets + k) -> vx -= com[3];
                  (moonlets + k) -> vy -= com[4];
                  (moonlets + k) -> vz -= com[5];
            }
            if (central_mass_bool){
                  CM.x -= com[0];  CM.y -= com[1];  CM.z -= com[2];  CM.vx -= com[3];  CM.vy -= com[4];  CM.vz -= com[5];
            }
      }
      
      /******** Communicating with REBOUND for 3D visualization ********/
      if (openGL_bool){
            rebound(moonlets);
      }
      
      return moonlets;
}


void variable_initialization(){

      /******** Defines and initializes external variables ********/

      timestep                = time_step;
      largest_id              = N_0 - 1;
      first_passage           = 1;
      time_elapsed            = 0.;
      how_many_free           = 0;
      how_many_moonlets       = N_0;
      how_many_cells          = 0;
      cell_id                 = 0;
      force_naive_bool        = N_0 < switch_to_brute_force ? 1 : 0;
      IndexPeanoHilbertOrder  = N_0;
      SideralOmega            = 2.*M_PI/Tearth;
      J2                      = J2_value == 0. ? 0.5*SideralOmega*SideralOmega*R_unit*R_unit*R_unit/(G*CM.mass) : J2_value;
      evection_resonance      = pow(1.5*sqrt(M_unit/pert_mass)*J2, 2./7.)*pow(pert_sma/R_unit, 3./7.)*R_unit;
      need_to_reduce_COM_bool = 0;
      previous_tra            = pert_tra;
      indexCollision          = 0;
      Rout                    = R_roche;
      fluid_disk_Sigma        = inner_mass/(M_PI*(Rout*Rout - R_unit*R_unit));
      flowed_since_last_spawn = 0.;
      if(!brute_force_bool){
            typ sinsigma = sin(inclination_max);
            gam = pow(sma_max*sma_max*sma_max-sma_min*sma_min*sma_min, 1./3.)*pow(4.*M_PI*how_many_neighbours*sinsigma/(((typ) N_0)*81.), 1./3.); //The mesh-size for the O(N) 
                                                                                                                                                  //collision detection algorithm
            gam_min = collision_cube_min/((typ) collision_cube_cells);
            if (gam < gam_min){
                  gam = gam_min;
            }

            how_many_big       = 0;
            how_many_small     = 0;
            how_many_modified  = 0;
            total_neighbours   = 0;
            collision_cube     = gam*((typ) collision_cube_cells);
            average_neighbours = 0.;
      }
      collision_count = 0;
      if (collision_bool && fragmentation_bool){
            super_catastrophic_count = 0;
            half_fragmentation_count = 0;
            full_fragmentation_count = 0;
            merger_count             = 0;
      }
      /******** Obtaining the mean anomaly of the point-mass perturbator at time t=0 ********/
      if (pert_ecc > 1.e-9){
            typ cart[6];
            typ alkhqp[6];
            ell2cart(pert_sma, pert_ecc, pert_inc, pert_tra, pert_aop, pert_lan, 1., cart);
            struct moonlet mlt = {cart[0], cart[1], cart[2], cart[3], cart[4], cart[5], 1., 1.};
            cart2ell(&mlt, 0, alkhqp, 1.);
            pert_M = alkhqp[1] - atan2(alkhqp[3], alkhqp[2]);
      }
      else{
            pert_M = pert_tra;
      }
}


void array_initialization(){

      /******** Defines and initializes external arrays ********/     
      
      
      /******** Array of bodies used as buffer ********/
      xx = (struct moonlet *)malloc(sizeof(struct moonlet)*N_max);

      /******** The kth cell of this array contains 1 if the kth cell of the array moonlets contains a body ********/
      exists = (int *)malloc(sizeof(int)*N_max);

      int p;

      /******** Initializing the array exists with 1 if the kth cell of the array moonlets contains a body, 0 otherwise ********/
      for (p = 0; p < N_0; p ++){
            *(exists + p) = 1;
      }
      for (p = N_0; p < N_max; p ++){
            *(exists + p) = 0;
      }
      
      free_indexes = (int *)malloc(N_max*sizeof(int));
      
      /******** If the O(N) algorithm is used for collision detection, initializing the hash table and its dependencies ********/
      if (mesh_bool && (mutual_bool || collision_bool)){
      
            bigint n = collision_cube_cells+2; //Adding one layer to the collision cube so I don't have to worry about edges.
            typ hash_table_size = ((typ) (n*n*n*(bigint) (sizeof(struct chain *))))/(1024.0*1024.0*1024.0); //Size of the hash table in GiB
            printf("Allocating %.1lf GiB of RAM to the hash table.\n", hash_table_size);
            hash = (struct chain **)malloc(sizeof(struct chain *)*n*n*n); //A pointer weighs 8 bytes, so this command should allocate 8GiB of RAM if collision_cube_cells = 1024
            if (hash == NULL){
                  fprintf(stderr, "Error : Can't allocate memory for the hash table in function array_initialization.\n");
                  abort();
            }
            for (p = 0; p < n*n*n; p ++){
                  *(hash + p) = NULL;
            }       
            modified_cells = (int *)malloc(sizeof(int)*N_max); //An array that contains the indexes of modified cells of the hash table. Initialized to {-1, -1, ..., -1}
            for (p = 0; p < N_max; p ++){
                  *(modified_cells + p) = -1;
            }            
            indexes = (int *)malloc(27*sizeof(int)); //An array that stores the indexes of the cells of the neighbourhood of a body in the hash table
            for (p = 0; p < 27; p ++){
                  *(indexes + p) = -1;
            }            
            nghb = NULL;
            add(0, &nghb);
            nghb -> how_many = 0;
            
            /******** If there are mutual gravitational interactions, I allocate the array of pairs that need to be treated for mutual gravitational interactions ********/
            pairs = (struct pair *)malloc(sizeof(struct pair) * 5 * N_max * (int) how_many_neighbours); //Array of pairs of bodies being in the same neighbourhood
                                                                                                        //The expected number of such pairs is N*how_many_neighbours/2
            how_many_pairs = 0;
            if (pairs == NULL){
                  fprintf(stderr, "Error : Can't allocate memory for the array pairs in function array_initialization.\n");
                  abort();
            }
      }     
      three_largest_indexes      = (int *)malloc(3*sizeof(int)); //Array of the indexes of the three largest bodies 
      if (collision_bool){
            approach    = (typ *)malloc(6    *sizeof(typ));
            did_collide = (int *)malloc(N_max*sizeof(int));
            for (p = 0; p < N_max; p ++){
                  *(did_collide + p) = 0;
            }
      }
      if ((mutual_bool || collision_bool) && !brute_force_bool && (falcON_bool || standard_tree_bool)){
            PeanoHilbertOrder = (int *)malloc(  N_max*sizeof(int));
            already_in_tree   = (int *)malloc(  N_max*sizeof(int));
            C1Moonlets        = (typ *)malloc(3*N_max*sizeof(typ));
            for (p = 0; p <= largest_id; p ++){
                  *(PeanoHilbertOrder + p) = p;
            }
      }
      if (openGL_bool){
            sending_buffer = (typ *)malloc(8*(N_max + (central_mass_bool || (viscoelastic_bool && pert_mass > 0.)))*sizeof(typ));
            for (p = 0; p <= largest_id; p ++){
                  *(sending_buffer + p) = 0.;
            }
      }
      if (viscoelastic_bool){
            typ cart[6];
            typ nu = get_perturbing_true_anomaly(t_init);
            ell2cart(pert_sma, pert_ecc, pert_inc, nu, pert_aop, pert_lan, G*(pert_mass + M_unit), cart);
            CM.x = cart[0];  CM.y = cart[1];  CM.z = cart[2];  CM.vx = cart[3];  CM.vy = cart[4];  CM.vz = cart[5];  CM.radius = pert_radius;
      }
      if (collision_bool && write_to_files_bool && write_collisions_bool){
            int size = output_step*N_max/5 + 10;
            collisionDatas = (struct collisionData *)malloc(size*sizeof(struct collisionData));
            if (collisionDatas == NULL){
                  fprintf(stderr, "Error: Cannot allocate array collisionDatas in function array_initialization.\n");
                  abort();
            }
      }
}


/*************************************/
/******** Memory deallocation ********/
/*************************************/


void deallocation(){

      /******** Deallocates the memory used globally ********/
      
      free(exists);
      exists = NULL;
      free(free_indexes);
      free_indexes = NULL;
      free(xx);
      xx = NULL;
      if(mesh_bool && (mutual_bool || collision_bool)){
            free(hash);
            hash = NULL;
            free(indexes);
            indexes = NULL;
            free(nghb);
            nghb = NULL;
            free(modified_cells);
            modified_cells = NULL;
            free(pairs);
            pairs = NULL;
      }
      free(three_largest_indexes);
      three_largest_indexes = NULL;
      if (collision_bool){
            free(approach);
            approach = NULL;
            free(did_collide);
            did_collide = NULL;
      }
      if ((mutual_bool || collision_bool) && !brute_force_bool && (falcON_bool || standard_tree_bool)){
            free(PeanoHilbertOrder);
            PeanoHilbertOrder = NULL;
            free(already_in_tree);
            already_in_tree = NULL;
            free(C1Moonlets);
            C1Moonlets = NULL;
      }
      if (openGL_bool){
            free(sending_buffer);
            sending_buffer = NULL;
      }
      if (viscoelastic_bool){
            free(connections);
            connections = NULL;
      }
      if (collision_bool && write_to_files_bool && write_collisions_bool){
            free(collisionDatas);
            collisionDatas = NULL;
      }
}

/**********************************************************************************/
/******** Functions relative to chains (unrolled linked list) manipulation ********/
/**********************************************************************************/


void add(int head, struct chain ** ch){

      /******** Adds the new element head to the chain *ch ********/
      
      if (*ch == NULL){
            struct chain * to_be_returned = (struct chain *)malloc(sizeof(struct chain));
            (to_be_returned -> ids)[0] = head;
            to_be_returned -> how_many = 1;
            to_be_returned -> queue = *ch;
            *ch = to_be_returned;
      }
      else{
            int hwmn = (*(ch)) -> how_many;
            if (hwmn < max_ids_per_node){
                  ((*ch) -> ids)[hwmn] = head;
                  (*ch) -> how_many += 1;
            }
            else{
                  struct chain * to_be_returned = (struct chain *)malloc(sizeof(struct chain));
                  (to_be_returned -> ids)[0] = head;
                  to_be_returned -> how_many = 1;
                  to_be_returned -> queue = *ch;
                  *ch = to_be_returned;
            }
      }  
}


struct chain * Add(int head, struct chain * ch){

      /******** Same as above but with different data types ********/
      
      if (ch == NULL){
            struct chain * to_be_returned = (struct chain *)malloc(sizeof(struct chain));
            (to_be_returned -> ids)[0] = head;
            to_be_returned -> how_many = 1;
            to_be_returned -> queue = ch;
            return to_be_returned;
      }
      else{
            int hwmn = ch -> how_many;
            if (hwmn < max_ids_per_node){
                  (ch -> ids)[hwmn] = head;
                  ch -> how_many += 1;
                  return ch;
            }
            else{
                  struct chain * to_be_returned = (struct chain *)malloc(sizeof(struct chain));
                  (to_be_returned -> ids)[0] = head;
                  to_be_returned -> how_many = 1;
                  to_be_returned -> queue = ch;
                  return to_be_returned;
            }
      }
}


struct chain * delete_chain(struct chain * ch){

      /******** Deletes the first element of the chain ch ********/
      
      if (ch != NULL){
            int hwmn = ch -> how_many;
            if (hwmn > 1){
                  ch -> how_many -= 1;
                  return ch;
            }
            else{
                  struct chain * to_be_returned = NULL;
                  to_be_returned = ch -> queue;
                  free(ch);
                  ch = NULL;
                  return to_be_returned;
            }
      }
      return NULL;
}


struct chain * partial_delete(struct chain * ch){

      /******** Deletes the first element of the chain ch. Does not unallocate the chain if empty ********/
      
      if (ch != NULL){
            int hwmn = ch -> how_many;
            if (hwmn > 1){
                  ch -> how_many -= 1;
                  return ch;
            }
            else{
                  if (ch -> queue != NULL){
                        struct chain * to_be_returned = NULL;
                        to_be_returned = ch -> queue;
                        free(ch);
                        ch = NULL;
                        return to_be_returned;
                  }
                  else{
                        ch -> how_many = 0;
                        return ch;
                  }
            }
      }
      return NULL;
}


struct chain * delete_node(struct chain * ch){

      /******** Deletes a full node of the chain ch ********/      
      
      struct chain * to_be_returned = NULL;
      if (ch != NULL){
            to_be_returned = ch -> queue;
            free(ch);
            ch = NULL;
            return to_be_returned;
      }
      return to_be_returned;
}


void clear_chain(struct chain ** ch){

      /******** Clears and deallocates the chain ch ********/
      
      while ((*ch) != NULL){
            *ch = delete_node(*ch);
      } 
}


/******************************************************************************/
/******** Miscellaneous functions relative to the numerical simulation ********/
/******************************************************************************/


int get_free_index(int should_be_put_at_the_end){

      /******** Returns the index of where to store a new body in the simulation ********/
      /******** The new body will be stored at *(moonlets + index)               ********/
      
      
      int index;

      if (how_many_free == 0 || should_be_put_at_the_end){
            largest_id += 1;
            index = largest_id;
            if (index == N_max){
                  fprintf(stderr, "Error : Cannot add a body to the simulation. Try increasing the value of N_max.\n");
                  abort();
            }           
      }
      else {
            how_many_free -= 1;
            index = *(free_indexes + how_many_free);
            if (index > largest_id){
                  largest_id = index;
            }  
      }
      if (*(exists + index)){ //Checking that the supposedly free index is indeed free. To be removed when the code is robust
            fprintf(stderr, "Error : This index is not free.\n");
            abort();
      }
      *(exists + index) = 1;
      return index;
}


void tidy_up(struct moonlet * moonlets){

      /******** Reorders the array moonlets so that the N bodies occupy the indexes 0 through N-1 ********/
      /******** This function is called at the end of the timestep, if N/largest_id < 0.9         ********/

      int i = 0;
      int j;
      
      for (j = 0; j <= largest_id; j ++){
            if (*(exists + j)){
                  if (j > i){
                        /******** Body i becomes body j ********/
                        *(moonlets + i) = *(moonlets + j);
                  }
                  i ++;
            }
      }
      
      /******** Updating the array exists ********/
      for (j = 0; j <= largest_id; j ++){
            if (j < i){
                  *(exists + j) = 1;
            }
            else{
                  *(exists + j) = 0;
            }
      }
      
      largest_id    = i - 1;
      how_many_free = 0;    
}


void three_largest_moonlets(struct moonlet * moonlets){

      /******** Actualizes the array three_largest_indexes with the indexes of the ********/
      /******** currently three largest bodies of the simulation. After a call to  ********/
      /******** this function, *three_largest_indexes, *(three_largest_indexes+1)  ********/
      /******** and *(three_largest_indexes+2) are the indexes of the largest,     ********/
      /******** second largest and third largest body of the simulation            ********/
      /******** This function is only called if mesh_bool is 1                     ********/

      int i = 0;
      int j;
      typ R_0 = 0., R_1 = 0., R_2 = 0., R = 0.;
      
      three_largest_indexes[0] = -1;
      three_largest_indexes[1] = -1;
      three_largest_indexes[2] = -1;

      for (j = 0; j <= largest_id; j ++){
            if (*(exists + j)){
            
                  if (i == 0){
                        R_0 = (moonlets + j) -> radius;
                        *three_largest_indexes = j;
                  }
                  
                  else if (i == 1){
                        R_1 = (moonlets + j) -> radius;
                        if (R_0 > R_1){
                              *(three_largest_indexes + 1) = j;
                        }
                        else {
                              *(three_largest_indexes + 1) = *three_largest_indexes;
                              *three_largest_indexes = j;
                              R = R_1;
                              R_1 = R_0;
                              R_0 = R;
                        }
                  }
                  
                  else if (i == 2){
                        R_2 = (moonlets + j) -> radius;
                        if (R_1 >= R_2){
                              *(three_largest_indexes + 2) = j;
                        }
                        else if (R_2 > R_1 && R_2 <= R_0){
                              *(three_largest_indexes + 2) = *(three_largest_indexes + 1);
                              *(three_largest_indexes + 1) = j;
                              R   = R_2;
                              R_2 = R_1;
                              R_1 = R;
                        }
                        else {
                              *(three_largest_indexes + 2) = *(three_largest_indexes + 1);
                              *(three_largest_indexes + 1) = *three_largest_indexes;
                              *three_largest_indexes = j;
                              R   = R_2;
                              R_2 = R_1;
                              R_1 = R_0;
                              R_0 = R;
                        }
                  }
                  
                  /******** At this stage, the three first body indexes occupy three_largest_indexes by decreasing radius ********/
                  
                  else { // If i > 2
                        R = (moonlets + j) -> radius;
                        if (R > R_0){
                              *(three_largest_indexes + 2) = *(three_largest_indexes + 1);
                              *(three_largest_indexes + 1) = *three_largest_indexes;
                              *three_largest_indexes = j;
                              R_2 = R_1;
                              R_1 = R_0;
                              R_0 = R;
                        }
                        else if (R > R_1){
                              *(three_largest_indexes + 2) = *(three_largest_indexes + 1);
                              *(three_largest_indexes + 1) = j;
                              R_2 = R_1;
                              R_1 = R;
                        }
                        else if (R > R_2){
                              *(three_largest_indexes + 2) = j;
                              R_2 = R;
                        }
                  }
                  i++;
            }
      }
}


void three_largest_three_first(struct moonlet * moonlets){

      /******** Makes sure that the three largest bodies are the three first bodies ********/
      /******** This function is only called if mesh_bool is 1                     ********/
      
      struct moonlet buffer;    // A buffer used to temporarily store a body
      int first, second, third; //The indexes of the first, second and third largest body
      
      first  = *three_largest_indexes;
      second = *(three_largest_indexes + 1);
      third  = *(three_largest_indexes + 2);
      
      if (first != 0){
            if (*exists){ //There is a body at index zero, but it is not the largest body. I exchange it with the largest body
                  buffer = *moonlets;
                  *moonlets = *(moonlets + first);
                  *(moonlets+first) = buffer;
            }
            else { //Index zero is not occupied. I occupy it with the largest body
                  *moonlets = *(moonlets + first);
                  lose_moonlet(first);
                  *exists = 1;
            }
            if (second == 0){
                  second = first;
            }
            else if (third == 0){
                  third = first;
            }
      }
      if (second != 1){
            if (*(exists + 1)){ //There is a body at index one, but it is not the second largest body. I exchange it with the second largest body
                  buffer = *(moonlets + 1);
                  *(moonlets + 1) = *(moonlets + second);
                  *(moonlets + second) = buffer;
            }
            else { //Index one is not occupied. I occupy it with the second largest body
                  *(moonlets + 1) = *(moonlets + second);
                  lose_moonlet(second);
                  *(exists + 1) = 1;
            }
            if (third == 1){
                  third = second;
            }
      }
      if (third != 2){
            if (*(exists + 2)){ //There is a body at index two, but it is not the third largest body. I exchange it with the third largest body
                  buffer = *(moonlets + 2);
                  *(moonlets + 2) = *(moonlets + third);
                  *(moonlets + third) = buffer;
            }
            else { //Index two is not occupied. I occupy it with the third largest body
                  *(moonlets + 2) = *(moonlets + third);
                  lose_moonlet(third);
                  *(exists + 2) = 1;
            }
      }

}


void cross_product(typ u_1, typ u_2, typ u_3, typ v_1, typ v_2, typ v_3, typ * uxv){

      /******** Fills uxv with the coordinates of the cross product u x v ********/
      
      *uxv       = u_2*v_3 - u_3*v_2;
      *(uxv + 1) = u_3*v_1 - u_1*v_3;
      *(uxv + 2) = u_1*v_2 - u_2*v_1;
      
}


void lose_moonlet(int a){

      /******** Treats the loss of body a from the simulation ********/
      
      *(free_indexes + how_many_free) = a;
      how_many_free ++;
      *(exists + a) = 0;
      if (a <= 2 && mesh_bool){
            how_many_free --;
      }
}


int maximum(typ i, typ j, typ k){

      /******** Returns the position of the maximum of (i,j,k) ********/
      /******** For example, if max(i,j,k) = k, then returns 2 ********/
      
      int m;
      typ mx;
      if (i >= j){
            m  = 0;
            mx = i;
      }
      else {
            m  = 1;
            mx = j;
      }
      
      if (mx >= k){
            return m;
      }
      else {
            return 2;
      }
}


void readFromFile(char * file_name, typ * storage, int n_data){

      /******** Reads floating point data from the file file_name and stores them in storage ********/
      /******** It is assumed that storage will not overflow. Expects to read n_data values  ********/

      FILE * file = fopen(file_name, "r");
      if (file == NULL){
            fprintf(stderr, "Error : Could not open file in function readFromFile. Did you specify the path 'pth' in the parameter file ?\n");
            abort();
      }
      typ i = 0.0;
      int j = 0;
      int returnValue = 1;
 
      while (returnValue == 1){
            if (j > n_data){
                  fprintf(stderr, "Error : There are too many data in path 'pth/init.txt' given the value of N_0. Update N_0 in src/parameters.h\n");
                  abort();
            }
            returnValue = fscanf(file, "%lf", &i);
            if (returnValue == 1){
                  *(storage + j) = i;
            }
            j ++;
      }
      fclose(file);
      if (j != n_data + 1){
            fprintf(stderr, "Error : There are not enough data in path 'pth/init.txt' given the value of N_0. Update N_0 in src/parameters.h\n");
            abort();
      }   
}


typ * readFromFile_withoutConstraint(char * file_name, int * size){

      /******** Same as above but the number of data to be read does not need to be specified. ********/
      /******** Writes the number of data read in *size and returns a buffer                   ********/

      FILE * file = fopen(file_name, "r");
      if (file == NULL){
            fprintf(stderr, "Error : Could not open file in function readFromFile_withoutConstraint. Did you specify the path 'pth' in the parameter file ?\n");
            abort();
      }
      typ i = 0.0;
      int j = 0;
      int returnValue = 1;
 
      while (returnValue == 1){ //Obtaining the number of data to be read
            returnValue = fscanf(file, "%lf", &i);
            j ++;
      }
      j --;
      *size = j;
      typ * buffer = (typ *)malloc(j*sizeof(typ)); //Allocating memory for the buffer
      if (buffer == NULL){
            fprintf(stderr, "Could not allocate memory for the buffer in function readFromFile_withoutConstraint\n");
            abort();
      }
      
      rewind(file); //Going back to the beginning of the file
      returnValue = 1;
      j           = 0;
      while (returnValue == 1){ //Reading the file a second time to store its data
            if (j > *size){
                  fprintf(stderr, "Error : The buffer is not big enough in function readFromFile_withoutConstraint.\n");
                  abort();
            }
            returnValue = fscanf(file, "%lf", &i);
            if (returnValue == 1){
                  *(buffer + j) = i;
            }
            j ++;
      }
      fclose(file);
      return buffer;
}


void total_momentum(struct moonlet * moonlets, typ * momentum){

      /******** Computes the total momentum of the system and stores it in the array momentum ********/


      typ X, Y, Z, vX, vY, vZ, m, M;
      int j;
      
      /******** Contribution from the central mass and inner fluid disk ********/
      X   = (central_mass_bool ? CM.x    : 0.);
      Y   = (central_mass_bool ? CM.y    : 0.);
      Z   = (central_mass_bool ? CM.z    : 0.);
      vX  = (central_mass_bool ? CM.vx   : 0.);
      vY  = (central_mass_bool ? CM.vy   : 0.);
      vZ  = (central_mass_bool ? CM.vz   : 0.);
      M   = (central_mass_bool ? CM.mass : 0.);
      M  += (inner_fluid_disk_bool ? fluid_disk_Sigma*M_PI*(Rout*Rout - R_unit*R_unit) : 0.);
      momentum[0] = M*X;  momentum[1] = M*Y;  momentum[2] = M*Z;  momentum[3] = M*vX;  momentum[4] = M*vY;  momentum[5] = M*vZ;
      
      /******** Contribution from the other bodies ********/
      for (j = 0; j <= largest_id; j ++){
            if (*(exists + j)){
                  X  = (moonlets + j) -> x;
                  Y  = (moonlets + j) -> y;
                  Z  = (moonlets + j) -> z;
                  vX = (moonlets + j) -> vx;
                  vY = (moonlets + j) -> vy;
                  vZ = (moonlets + j) -> vz;
                  m  = (moonlets + j) -> mass;
                  M += m;
                  momentum[0] += m*X;  momentum[1] += m*Y;  momentum[2] += m*Z;  momentum[3] += m*vX;  momentum[4] += m*vY;  momentum[5] += m*vZ;
            }
      }
      momentum[0] /= M;  momentum[1] /= M;  momentum[2] /= M;  momentum[3] /= M;  momentum[4] /= M;  momentum[5] /= M;  
}


void verify(){

      /******** Checks if the parameter file is properly setup ********/
      
      
      /******** I first check that booleans are given as integers ********/
      if(!type_check(typeof(write_to_files_bool),      int)){fprintf(stderr, "Error : write_to_files_bool must be given as an integer (0 or 1).\n");     abort();}
      if(!type_check(typeof(make_animation_bool),      int)){fprintf(stderr, "Error : make_animation_bool must be given as an integer (0 or 1).\n");     abort();}
      if(!type_check(typeof(write_cartesian_bool),     int)){fprintf(stderr, "Error : write_cartesian_bool must be given as an integer (0 or 1).\n");    abort();}
      if(!type_check(typeof(write_elliptic_bool),      int)){fprintf(stderr, "Error : write_elliptic_bool must be given as an integer (0 or 1).\n");     abort();}
      if(!type_check(typeof(central_mass_bool),        int)){fprintf(stderr, "Error : central_mass_bool must be given as an integer (0 or 1).\n");       abort();}
      if(!type_check(typeof(reduce_to_COM_bool),       int)){fprintf(stderr, "Error : reduce_to_COM_bool must be given as an integer (0 or 1).\n");      abort();}
      if(!type_check(typeof(random_initial_bool),      int)){fprintf(stderr, "Error : random_initial_bool must be given as an integer (0 or 1).\n");     abort();}
      if(!type_check(typeof(initial_cartesian_bool),   int)){fprintf(stderr, "Error : initial_cartesian_bool must be given as an integer (0 or 1).\n");  abort();}
      if(!type_check(typeof(seed_bool),                int)){fprintf(stderr, "Error : seed_bool must be given as an integer (0 or 1).\n");               abort();}
      if(!type_check(typeof(one_collision_only_bool),  int)){fprintf(stderr, "Error : one_collision_only_bool must be given as an integer (0 or 1).\n"); abort();}
      if(!type_check(typeof(openGL_bool),              int)){fprintf(stderr, "Error : openGL_bool must be given as an integer (0 or 1).\n");             abort();}
      if(!type_check(typeof(resume_simulation_bool),   int)){fprintf(stderr, "Error : resume_simulation_bool must be given as an integer (0 or 1).\n");  abort();}
      if(!type_check(typeof(viscoelastic_bool),        int)){fprintf(stderr, "Error : viscoelastic_bool must be given as an integer (0 or 1).\n");       abort();}
      if(!type_check(typeof(J2_bool),                  int)){fprintf(stderr, "Error : J2_bool must be given as an integer (0 or 1).\n");                 abort();}
      if(!type_check(typeof(inner_fluid_disk_bool),    int)){fprintf(stderr, "Error : inner_fluid_disk_bool must be given as an integer (0 or 1).\n");   abort();}
      if(!type_check(typeof(central_tides_bool),       int)){fprintf(stderr, "Error : central_tides_bool must be given as an integer (0 or 1).\n");      abort();}
      if(!type_check(typeof(collision_bool),           int)){fprintf(stderr, "Error : collision_bool must be given as an integer (0 or 1).\n");          abort();}
      if(!type_check(typeof(mutual_bool),              int)){fprintf(stderr, "Error : mutual_bool must be given as an integer (0 or 1).\n");             abort();}
      if(!type_check(typeof(brute_force_bool),         int)){fprintf(stderr, "Error : brute_force_bool must be given as an integer (0 or 1).\n");        abort();}
      if(!type_check(typeof(falcON_bool),              int)){fprintf(stderr, "Error : falcON_bool must be given as an integer (0 or 1).\n");             abort();}
      if(!type_check(typeof(standard_tree_bool),       int)){fprintf(stderr, "Error : standard_tree_bool must be given as an integer (0 or 1).\n");      abort();}
      if(!type_check(typeof(mesh_bool),                int)){fprintf(stderr, "Error : mesh_bool must be given as an integer (0 or 1).\n");               abort();}
      if(!type_check(typeof(elastic_collision_bool),   int)){fprintf(stderr, "Error : elastic_collision_bool must be given as an integer (0 or 1).\n");  abort();}
      if(!type_check(typeof(inelastic_collision_bool), int)){fprintf(stderr, "Error : inelastic_collision_bool must be given as an integer (0 or 1).\n");abort();}
      if(!type_check(typeof(instant_merger_bool),      int)){fprintf(stderr, "Error : instant_merger_bool must be given as an integer (0 or 1).\n");     abort();}
      if(!type_check(typeof(fragmentation_bool),       int)){fprintf(stderr, "Error : fragmentation_bool must be given as an integer (0 or 1).\n");      abort();}
      
      /******** I now check if all the booleans are either 0 or 1 ********/
      int OK = 1;
      if (write_to_files_bool      != 0 && write_to_files_bool      != 1){OK = 0;}
      if (make_animation_bool      != 0 && make_animation_bool      != 1){OK = 0;}
      if (write_cartesian_bool     != 0 && write_cartesian_bool     != 1){OK = 0;}
      if (write_elliptic_bool      != 0 && write_elliptic_bool      != 1){OK = 0;}
      if (central_mass_bool        != 0 && central_mass_bool        != 1){OK = 0;}
      if (reduce_to_COM_bool       != 0 && reduce_to_COM_bool       != 1){OK = 0;}
      if (random_initial_bool      != 0 && random_initial_bool      != 1){OK = 0;}
      if (initial_cartesian_bool   != 0 && initial_cartesian_bool   != 1){OK = 0;}
      if (seed_bool                != 0 && seed_bool                != 1){OK = 0;}
      if (one_collision_only_bool  != 0 && one_collision_only_bool  != 1){OK = 0;}
      if (openGL_bool              != 0 && openGL_bool              != 1){OK = 0;}
      if (resume_simulation_bool   != 0 && resume_simulation_bool   != 1){OK = 0;}
      if (viscoelastic_bool        != 0 && viscoelastic_bool        != 1){OK = 0;}
      if (J2_bool                  != 0 && J2_bool                  != 1){OK = 0;}
      if (inner_fluid_disk_bool    != 0 && inner_fluid_disk_bool    != 1){OK = 0;}
      if (central_tides_bool       != 0 && central_tides_bool       != 1){OK = 0;}
      if (collision_bool           != 0 && collision_bool           != 1){OK = 0;}
      if (mutual_bool              != 0 && mutual_bool              != 1){OK = 0;}
      if (brute_force_bool         != 0 && brute_force_bool         != 1){OK = 0;}
      if (falcON_bool              != 0 && falcON_bool              != 1){OK = 0;}
      if (standard_tree_bool       != 0 && standard_tree_bool       != 1){OK = 0;}
      if (mesh_bool                != 0 && mesh_bool                != 1){OK = 0;}
      if (elastic_collision_bool   != 0 && elastic_collision_bool   != 1){OK = 0;}
      if (inelastic_collision_bool != 0 && inelastic_collision_bool != 1){OK = 0;}
      if (instant_merger_bool      != 0 && instant_merger_bool      != 1){OK = 0;}
      if (fragmentation_bool       != 0 && fragmentation_bool       != 1){OK = 0;}
      if (!OK){
            fprintf(stderr, "Error : All booleans must be either 0 or 1 in the parameter file.\n");
            abort();
      }
      
      /******** I now check the booleans relative to collision resolution and mutual interactions ********/
      if (collision_bool && elastic_collision_bool + inelastic_collision_bool + instant_merger_bool + fragmentation_bool != 1){
            fprintf(stderr, "Error : Exactly one of these booleans must be 1 : (elastic_collision_bool, inelastic_collision_bool, instant_merger_bool, fragmentation_bool).\n");
            abort();
      }
      if (mutual_bool && brute_force_bool + falcON_bool + standard_tree_bool + mesh_bool != 1){
            fprintf(stderr, "Error : Exactly one of these booleans must be 1 : (brute_force_bool, falcON_bool, standard_tree_bool, mesh_bool).\n");
            abort();
      }
      
      /******** I now verify that floating point numbers stayed that way ********/
      if(!type_check(typeof(R_unit),                    typ)){fprintf(stderr, "Error : R_unit must be given as a floating-point number.\n");                abort();}
      if(!type_check(typeof(M_unit),                    typ)){fprintf(stderr, "Error : M_unit must be given as a floating-point number.\n");                abort();}
      if(!type_check(typeof(G)     ,                    typ)){fprintf(stderr, "Error : G must be given as a floating-point number.\n");                     abort();}
      if(!type_check(typeof(Tearth),                    typ)){fprintf(stderr, "Error : Tearth must be given as a floating-point number.\n");                abort();}
      if(!type_check(typeof(J2_value),                  typ)){fprintf(stderr, "Error : J2_value must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(k2),                        typ)){fprintf(stderr, "Error : k2 must be given as a floating-point number.\n");                    abort();}
      if(!type_check(typeof(Delta_t),                   typ)){fprintf(stderr, "Error : Delta_t must be given as a floating-point number.\n");               abort();}
      if(!type_check(typeof(dimensionless_moi),         typ)){fprintf(stderr, "Error : dimensionless_moi must be given as a floating-point number.\n");     abort();}
      if(!type_check(typeof(inner_mass),                typ)){fprintf(stderr, "Error : inner_mass must be given as a floating-point number.\n");            abort();}
      if(!type_check(typeof(spawned_density),           typ)){fprintf(stderr, "Error : spawned_density must be given as a floating-point number.\n");       abort();}
      if(!type_check(typeof(f_tilde),                   typ)){fprintf(stderr, "Error : f_tilde must be given as a floating-point number.\n");               abort();}
      if(!type_check(typeof(R_roche),                   typ)){fprintf(stderr, "Error : R_roche must be given as a floating-point number.\n");               abort();}
      if(!type_check(typeof(disruption_threshold),      typ)){fprintf(stderr, "Error : disruption_threshold must be given as a floating-point number.\n");  abort();}
      if(!type_check(typeof(pert_sma),                  typ)){fprintf(stderr, "Error : pert_sma must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(pert_ecc),                  typ)){fprintf(stderr, "Error : pert_ecc must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(pert_inc),                  typ)){fprintf(stderr, "Error : pert_inc must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(pert_tra),                  typ)){fprintf(stderr, "Error : pert_tra must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(pert_aop),                  typ)){fprintf(stderr, "Error : pert_aop must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(pert_lan),                  typ)){fprintf(stderr, "Error : pert_lan must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(pert_mass),                 typ)){fprintf(stderr, "Error : pert_mass must be given as a floating-point number.\n");             abort();}
      if(!type_check(typeof(pert_radius),               typ)){fprintf(stderr, "Error : pert_radius must be given as a floating-point number.\n");           abort();}
      if(!type_check(typeof(t_end),                     typ)){fprintf(stderr, "Error : t_end must be given as a floating-point number.\n");                 abort();}
      if(!type_check(typeof(time_step),                 typ)){fprintf(stderr, "Error : time_step must be given as a floating-point number.\n");             abort();}
      if(!type_check(typeof(high_dumping_threshold),    typ)){fprintf(stderr, "Error : high_dumping_threshold must be given as a floating-point number.\n");abort();}
      if(!type_check(typeof(softening_parameter),       typ)){fprintf(stderr, "Error : softening_parameter must be given as a floating-point number.\n");   abort();}
      if(!type_check(typeof(radius_min),                typ)){fprintf(stderr, "Error : radius_min must be given as a floating-point number.\n");            abort();}
      if(!type_check(typeof(radius_max),                typ)){fprintf(stderr, "Error : radius_max must be given as a floating-point number.\n");            abort();}
      if(!type_check(typeof(density_min),               typ)){fprintf(stderr, "Error : density_min must be given as a floating-point number.\n");           abort();}
      if(!type_check(typeof(density_max),               typ)){fprintf(stderr, "Error : density_max must be given as a floating-point number.\n");           abort();}
      if(!type_check(typeof(eccentricity_min),          typ)){fprintf(stderr, "Error : eccentricity_min must be given as a floating-point number.\n");      abort();}
      if(!type_check(typeof(eccentricity_max),          typ)){fprintf(stderr, "Error : eccentricity_max must be given as a floating-point number.\n");      abort();}
      if(!type_check(typeof(sma_min),                   typ)){fprintf(stderr, "Error : sma_min must be given as a floating-point number.\n");               abort();}
      if(!type_check(typeof(sma_max),                   typ)){fprintf(stderr, "Error : sma_max must be given as a floating-point number.\n");               abort();}
      if(!type_check(typeof(inclination_min),           typ)){fprintf(stderr, "Error : inclination_min must be given as a floating-point number.\n");       abort();}
      if(!type_check(typeof(inclination_max),           typ)){fprintf(stderr, "Error : inclination_max must be given as a floating-point number.\n");       abort();}
      if(!type_check(typeof(radius_blow_up_factor),     typ)){fprintf(stderr, "Error : radius_blow_up_factor must be given as a floating-point number.\n"); abort();}
      if(!type_check(typeof(spring_modulus),            typ)){fprintf(stderr, "Error : spring_modulus must be given as a floating-point number.\n");        abort();}
      if(!type_check(typeof(damping_coefficient),       typ)){fprintf(stderr, "Error : damping_coefficient must be given as a floating-point number.\n");   abort();}
      if(!type_check(typeof(connections_per_node),      typ)){fprintf(stderr, "Error : connections_per_node must be given as a floating-point number.\n");  abort();}
      if(!type_check(typeof(nodes_radius),              typ)){fprintf(stderr, "Error : nodes_radius must be given as a floating-point number.\n");          abort();}
      if(!type_check(typeof(minimal_distance),          typ)){fprintf(stderr, "Error : minimal_distance must be given as a floating-point number.\n");      abort();}
      if(!type_check(typeof(OmegaX),                    typ)){fprintf(stderr, "Error : OmegaX must be given as a floating-point number.\n");                abort();}
      if(!type_check(typeof(OmegaY),                    typ)){fprintf(stderr, "Error : OmegaY must be given as a floating-point number.\n");                abort();}
      if(!type_check(typeof(OmegaZ),                    typ)){fprintf(stderr, "Error : OmegaZ must be given as a floating-point number.\n");                abort();}
      if(!type_check(typeof(lbd_long),                  typ)){fprintf(stderr, "Error : lbd_long must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(beta_lat),                  typ)){fprintf(stderr, "Error : beta_lat must be given as a floating-point number.\n");              abort();}
      if(!type_check(typeof(theta_min),                 typ)){fprintf(stderr, "Error : theta_min must be given as a floating-point number.\n");             abort();}
      if(!type_check(typeof(root_sidelength),           typ)){fprintf(stderr, "Error : root_sidelength must be given as a floating-point number.\n");       abort();}
      if(!type_check(typeof(collision_cube_min),        typ)){fprintf(stderr, "Error : collision_cube_min must be given as a floating-point number.\n");    abort();}
      if(!type_check(typeof(how_many_neighbours),       typ)){fprintf(stderr, "Error : how_many_neighbours must be given as a floating-point number.\n");   abort();}
      if(!type_check(typeof(collision_parameter),       typ)){fprintf(stderr, "Error : collision_parameter must be given as a floating-point number.\n");   abort();}
      if(!type_check(typeof(mu_parameter),              typ)){fprintf(stderr, "Error : mu_parameter must be given as a floating-point number.\n");          abort();}
      if(!type_check(typeof(nu_parameter),              typ)){fprintf(stderr, "Error : nu_parameter must be given as a floating-point number.\n");          abort();}
      if(!type_check(typeof(C1_parameter),              typ)){fprintf(stderr, "Error : C1_parameter must be given as a floating-point number.\n");          abort();}
      if(!type_check(typeof(k_parameter),               typ)){fprintf(stderr, "Error : k_parameter must be given as a floating-point number.\n");           abort();}
      if(!type_check(typeof(merging_threshold),         typ)){fprintf(stderr, "Error : merging_threshold must be given as a floating-point number.\n");     abort();}
      if(!type_check(typeof(fragment_threshold),        typ)){fprintf(stderr, "Error : fragment_threshold must be given as a floating-point number.\n");    abort();}

      /******** I now verify that integer numbers stayed that way ********/
      if(!type_check(typeof(N_max),                     int)){fprintf(stderr, "Error : N_max must be given as an integer.\n");                              abort();}
      if(!type_check(typeof(N_0),                       int)){fprintf(stderr, "Error : N_0 must be given as an integer.\n");                                abort();}
      if(!type_check(typeof(output_step),               int)){fprintf(stderr, "Error : output_step must be given as an integer.\n");                        abort();}
      if(!type_check(typeof(max_ids_per_node),          int)){fprintf(stderr, "Error : max_ids_per_node must be given as an integer.\n");                   abort();}
      if(!type_check(typeof(seed),                      int)){fprintf(stderr, "Error : seed must be given as an integer.\n");                               abort();}
      if(!type_check(typeof(switch_to_brute_force),     int)){fprintf(stderr, "Error : switch_to_brute_force must be given as an integer.\n");              abort();}
      if(!type_check(typeof(browser_port),              int)){fprintf(stderr, "Error : browser_port must be given as an integer.\n");                       abort();}
      if(!type_check(typeof(expansion_order),           int)){fprintf(stderr, "Error : expansion_order must be given as an integer.\n");                    abort();}
      if(!type_check(typeof(subdivision_threshold),     int)){fprintf(stderr, "Error : subdivision_threshold must be given as an integer.\n");              abort();}
      if(!type_check(typeof(level_max),                 int)){fprintf(stderr, "Error : level_max must be given as an integer.\n");                          abort();}
      if(!type_check(typeof(N_cc_pre),                  int)){fprintf(stderr, "Error : N_cc_pre must be given as an integer.\n");                           abort();}
      if(!type_check(typeof(N_cc_post),                 int)){fprintf(stderr, "Error : N_cc_post must be given as an integer.\n");                          abort();}
      if(!type_check(typeof(N_cs),                      int)){fprintf(stderr, "Error : N_cs must be given as an integer.\n");                               abort();}
      if(!type_check(typeof(N_cb_pre),                  int)){fprintf(stderr, "Error : N_cb_pre must be given as an integer.\n");                           abort();}
      if(!type_check(typeof(N_cb_post),                 int)){fprintf(stderr, "Error : N_cb_post must be given as an integer.\n");                          abort();}
      if(!type_check(typeof(N_cc_collision),            int)){fprintf(stderr, "Error : N_cc_collision must be given as an integer.\n");                     abort();}
      if(!type_check(typeof(N_cs_collision),            int)){fprintf(stderr, "Error : N_cs_collision must be given as an integer.\n");                     abort();}
      if(!type_check(typeof(N_cb_collision),            int)){fprintf(stderr, "Error : N_cb_collision must be given as an integer.\n");                     abort();}
      if(!type_check(typeof(collision_cube_cells),      int)){fprintf(stderr, "Error : collision_cube_cells must be given as an integer.\n");               abort();}
      if(!type_check(typeof(N_tilde),                   int)){fprintf(stderr, "Error : N_tilde must be given as an integer.\n");                            abort();}
      
      /******** Miscellaneous verifications ********/
      if (N_max < N_0){
            fprintf(stderr, "Error : N_max cannot be smaller than N_0.\n");
            abort();
      }
      if (mutual_bool && (falcON_bool || standard_tree_bool) && (expansion_order < 1 || expansion_order > 8)){
            fprintf(stderr, "Error : The expansion order in the multipole expansion must be between 1 and 8 included.\n");
            abort();
      }
      if (collision_bool && fragmentation_bool && N_max < N_tilde*N_0/2){
            printf("Warning : Your choice of N_max might be too small for the fragmentation model given your choice of N_tilde.\n");
      }
      if (mutual_bool && falcON_bool && N_cc_post < N_cc_pre){
            fprintf(stderr, "Error : N_cc_post cannot be smaller than N_cc_pre.\n");
            abort();
      }
      if (mutual_bool && standard_tree_bool && N_cb_post < N_cb_pre){
            fprintf(stderr, "Error : N_cb_post cannot be smaller than N_cb_pre.\n");
            abort();
      }
      if (!central_mass_bool && (J2_bool || inner_fluid_disk_bool || central_tides_bool)){
            fprintf(stderr, "Error : If central_mass_bool is 0, then none of these booleans can be 1 : (J2_bool, inner_fluid_disk_bool, central_tides_bool).\n");
            abort();
      }
      if (pert_mass > 0. && !reduce_to_COM_bool){
            fprintf(stderr, "Error : Perturbation from a point-mass perturbator can only be taken into account if reduce_to_COM_bool is 1.\n");
            abort();
      }
      if (make_animation_bool && !write_to_files_bool){
            fprintf(stderr, "Error : write_to_files_bool must be set to 1 when make_animation_bool is set to 1.\n");
            abort();
      }
      if (make_animation_bool && !write_elliptic_bool){
            fprintf(stderr, "Error : write_elliptic_bool must be set to 1 when make_animation_bool is set to 1.\n");
            abort();
      }
      if (viscoelastic_bool && central_mass_bool){
            fprintf(stderr, "Error : central_mass_bool must be set to 0 if viscoelastic_bool is set to 1.\n");
            abort();
      }
      if (viscoelastic_bool && (mesh_bool || standard_tree_bool)){
            fprintf(stderr, "Error : Only falcON algorithm and the brute-force method are supported when NcorpiON is used to simulate a viscoelastic body.\n");
            abort();
      }
      if (viscoelastic_bool && collision_bool && (instant_merger_bool || fragmentation_bool)){
            fprintf(stderr, "Error : Collisions can only be resolved elastically and inelastically when NcorpiON is used to simulate a viscoelastic body.\n");
            abort();
      }
      if (viscoelastic_bool && !reduce_to_COM_bool){
            printf("Warning : The boolean reduce_to_COM_bool should be set to 1 so the viscoelastic body can be properly rotated.\n");
      }
}


void Lyapunov(struct moonlet * moonlets){

      /******** To be removed eventually ********/
      
      time_t t;
      time(&t);
      srand((unsigned) t);
      typ D = 1.0e-11;
      int j;
      typ dx, dy, dz, dvx, dvy, dvz;
      
      for (j = 0; j < N_0; j ++){
            dx  = rdm(-D, D);
            dy  = rdm(-D, D);
            dz  = rdm(-D, D);
            dvx = rdm(-D, D);
            dvy = rdm(-D, D);
            dvz = rdm(-D, D);
            (moonlets + j) -> x  += dx;
            (moonlets + j) -> y  += dy;
            (moonlets + j) -> z  += dz;
            (moonlets + j) -> vx += dvx;
            (moonlets + j) -> vy += dvy;
            (moonlets + j) -> vz += dvz;
      }
}
