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
/******** along with rebound.  If not, see <http://www.gnu.org/licenses/>.     ********/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structure.h"
#include "parameters.h"
#include "ffm.h"
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/sysinfo.h>


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
typ * sun_vector;


/******** Declaring external variables ********/
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
typ time_since_last_spawn;
typ time_between_spawn;
typ fluid_disk_Sigma;
typ SideralOmega;
typ star_mean_motion;


typ * ell2cart(typ a, typ e, typ i, typ nu, typ omega, typ Omega){

      /******** Returns the array [X,Y,Z,vX,vY,vZ] or the cartesian coordinates                                ********/
      /******** a is the semi-major axis, e is the eccentricity, i is the inclination, nu is the true anomaly, ********/
      /******** omega is the argument of periapsis and Omega is the longitude of the ascending node            ********/
      
      typ X,Y,Z,vX,vY,vZ; //Cartesian coordinates
      typ X_buff,Y_buff,Z_buff,vX_buff,vY_buff,vZ_buff; //Buffer for cartesian coordinates
      typ r; //Moonlet's distance to Earth's center
      typ n; //Moonlet's mean motion
      typ mu; // Gravitational parameter
      typ cosnu=cos(nu);
      typ sinnu=sin(nu);
      typ cosomega=cos(omega);
      typ sinomega=sin(omega);
      typ cosi=cos(i);
      typ sini=sin(i);
      typ cosOmega=cos(Omega);
      typ sinOmega=sin(Omega);

      /********  In the orbital plane (see Laskar & Robutel 1995)  ********/
      r=a*(1.0-e*e)/(1.0+e*cosnu);
      mu=G*Mearth;
      if (J2_bool){ //Using the geometric elliptical elements instead of the osculating ones in case of an oblate Earth. See Greenberg (1981)
            mu = mu*(1.0+1.5*J2*Rearth*Rearth/(a*a));
      }
      n=sqrt(mu/(a*a*a));
      X=r*cosnu;
      Y=r*sinnu;
      Z=0.0;
      vX=-n*a/sqrt(1.0-e*e)*sinnu;
      vY=n*a/sqrt(1.0-e*e)*(e+cosnu);
      vZ=0.0;
      
      /********  Rotations to convert to reference plane (see Laskar & Robutel 1995)  ********/
      X_buff=X; Y_buff=Y; vX_buff=vX; vY_buff=vY;
      X=cosomega*X_buff-sinomega*Y_buff;
      vX=cosomega*vX_buff-sinomega*vY_buff;
      Y=sinomega*X_buff+cosomega*Y_buff;
      vY=sinomega*vX_buff+cosomega*vY_buff;
      
      Y_buff=Y; Z_buff=Z; vY_buff=vY; vZ_buff=vZ;
      Y=cosi*Y_buff-sini*Z_buff;
      vY=cosi*vY_buff-sini*vZ_buff;
      Z=sini*Y_buff+cosi*Z_buff;
      vZ=sini*vY_buff+cosi*vZ_buff;
      
      X_buff=X; Y_buff=Y; vX_buff=vX; vY_buff=vY;
      X=cosOmega*X_buff-sinOmega*Y_buff;
      vX=cosOmega*vX_buff-sinOmega*vY_buff;
      Y=sinOmega*X_buff+cosOmega*Y_buff;
      vY=sinOmega*vX_buff+cosOmega*vY_buff;
      
      /********  Returning the cartesian coordinates  ********/
      typ * cart=(typ *)malloc(sizeof(typ)*6);
      *cart=X;
      *(cart+1)=Y;
      *(cart+2)=Z;
      *(cart+3)=vX;
      *(cart+4)=vY;
      *(cart+5)=vZ;
      return cart;

}


void cart2aei(struct moonlet * moonlets, int id, typ * aei){

      /******** Modifies the array [a,e,i] of the semimajor axis, the eccentricity and the inclination ********/
      /******** of the moonlet whose id is id. Stores these three quantities in the array aei.         ********/
      
      
      typ r,v2,r_cross_v,r_cross_v_square,X,Y,Z,vX,vY,vZ,a,e,i,cosi,r_cross_v_1,r_cross_v_2,r_cross_v_3,mu;
      
      /******** Getting the cartesian coordinates ********/
      X  = (moonlets + id) -> x;
      Y  = (moonlets + id) -> y;
      Z  = (moonlets + id) -> z;
      vX = (moonlets + id) -> vx;
      vY = (moonlets + id) -> vy;
      vZ = (moonlets + id) -> vz;
      
      /******** Getting the semimajor axis from the orbital energy ********/
      r  = sqrt(X*X+Y*Y+Z*Z);
      v2 = vX*vX+vY*vY+vZ*vZ;
      mu = G*Mearth;
      a  = mu*r/(2.0*mu-r*v2);
      
      /******** If the Earth is oblate, correcting Kepler third law (See Greenberg, 1981) ********/
      if (J2_bool){
            mu = mu*(1.0+1.5*J2*Rearth*Rearth/(a*a));
            a = mu*r/(2.0*mu-r*v2);
      }
      
      /******** Getting the eccentricity from the momentum ********/
      r_cross_v_1 = Y*vZ-vY*Z;
      r_cross_v_2 = vX*Z-X*vZ;
      r_cross_v_3 = X*vY-vX*Y;
      r_cross_v_square = r_cross_v_1*r_cross_v_1 + r_cross_v_2*r_cross_v_2 + r_cross_v_3*r_cross_v_3;
      e = sqrt(1.0-r_cross_v_square/(mu*a));
      
      /******** Getting the inclination from the third component of the angular momentum ********/
      r_cross_v = sqrt(r_cross_v_square);
      cosi = r_cross_v_3/r_cross_v;
      i = acos(cosi);
      
      /******** Filling the array aei ********/
      *aei = a;
      *(aei+1) = e;
      *(aei+2) = i; 

}


/******************************************/
/******** Initialization functions ********/
/******************************************/


struct moonlet init(typ a, typ e, typ i, typ nu, typ omega, typ Omega, typ rad){      //Initializes a moonlet

      /******** a is the semi-major axis, e is the eccentricity, i is the inclination, nu is the true anomaly,     ********/
      /******** omega is the argument of periapsis, Omega is the longitude of the ascending node and m is the mass ********/

      struct moonlet mlt;
      mlt.radius = rad;                                 //Initializing the radius
      typ m = 4.0/3.0*M_PI*density*rad*rad*rad;         //Defining the mass
      mlt.mass = m;                                     //Initializing the mass
      
      typ * cart = ell2cart(a, e, i, nu, omega, Omega); //Computing the cartesian coordinates from the orbital elements
      mlt.x  = * cart;                                      //Initializing the cartesian coordinates
      mlt.y  = *(cart+1);
      mlt.z  = *(cart+2);
      mlt.vx = *(cart+3);
      mlt.vy = *(cart+4);
      mlt.vz = *(cart+5);
      free(cart);
      cart=NULL;
      
      return mlt;

}


typ rdm(typ min, typ max){      //Returns a random number between min and max according to a (a priori) uniform distribution
      
      /******** Generating the random number ********/
      typ MyRand = ((typ) rand())/((typ) RAND_MAX);  //between 0 and 1
      MyRand = min+(max-min)*MyRand;                 //between min and max
      return MyRand;

}


struct moonlet * populate(typ M, typ dR){

      /******** Populates the simulation with N_0 moonlets. The semi-major axes are chosen randomly and uniformly in [sma_min, sma_max], and        ********/
      /******** similarly for the eccentricities and inclinations, which are defined by constants of the file structure.h                           ********/
      /******** The angles nu, omega and Omega are chosen randomly and uniformly in [0,2*pi]. N_max should be chosen larger than N_0 to handle the  ********/
      /******** variation in moonlet number. M is the approximate total mass and dR<1 is such that the moonlets have a radius chosen randomly and   ********/
      /******** uniformly in [(1-dR)*mean_rad, (1+dR)*mean_rad] where mean_rad is the radius corresponding to a mass M/N_0                          ********/


      struct moonlet * moonlets = (struct moonlet *)malloc(N_max*sizeof(struct moonlet)); //The array of moonlets
      if (moonlets == NULL){
            fprintf(stderr, "Error : Can't allocate array of moonlets.\n");
            abort();
      }
      if (radius_stddev >= 1.0/sqrt(3.0)){
            fprintf(stderr, "Error : radius_stddev must be less than 1/sqrt(3).\n");
            abort();
      }
      
      /******** Populating the moonlet array ********/
      typ a, e, i, nu, omega, Omega, rad, m, mean_rad;
      int k;
      m = M/((typ) N_0);
      mean_rad = pow(3.0*m/(4.0*M_PI*density*(1.0+3.0*radius_stddev*radius_stddev)), 1.0/3.0);
      for (k=0; k<N_0; k++){
            a               = rdm(sma_min, sma_max);
            e               = rdm(eccentricity_min, eccentricity_max);
            i               = rdm(inclination_min, inclination_max);
            nu              = rdm(0.0, 2.0*M_PI);
            omega           = rdm(0.0, 2.0*M_PI);
            Omega           = rdm(0.0, 2.0*M_PI);
            rad             = rdm((1.0-sqrt(3.0)*radius_stddev)*mean_rad, (1.0+sqrt(3.0)*radius_stddev)*mean_rad);
            *(moonlets + k) = init(a, e, i, nu, omega, Omega, rad);
      }
      for (k=N_0; k<N_max; k++){ //Filling the unused cells of the moonlet array with whatever
            *(moonlets + k) = init(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
      }
      
      return moonlets;

}


void variable_initialization(){

      /******** Defines and initializes external variables ********/
      
      if (J2_value == 0.0){
            J2 = 0.5/(Tearth*Tearth); //Earth's equatorial bulge parameter. This equation is valid only if Earth is fluid and Tearth^2>>1
      }
      else{
            J2 = J2_value;
      }
      timestep               = time_step;
      largest_id             = N_0-1;
      first_passage          = 1;
      time_elapsed           = 0.0;
      how_many_free          = 0;
      how_many_moonlets      = N_0;
      how_many_cells         = 0;
      cell_id                = 0;
      force_naive_bool       = 0;
      IndexPeanoHilbertOrder = N_0;
      SideralOmega           = 2.0*M_PI/Tearth;
      star_mean_motion       = sqrt(G*(Mearth + star_mass)/(star_semi_major*star_semi_major*star_semi_major));
      if(!brute_force_bool){
            typ sinsigma = sin(inclination_max);
            gam = pow(sma_max*sma_max*sma_max-sma_min*sma_min*sma_min,1.0/3.0)*pow(4.0*M_PI*how_many_neighbours*sinsigma/(((typ) N_0)*81.0),1.0/3.0); //The mesh-size for the O(N) 
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
            average_neighbours = 0.0;
      }
      collision_count = 0;
      if (collision_bool && fragmentation_bool){
            super_catastrophic_count = 0;
            half_fragmentation_count = 0;
            full_fragmentation_count = 0;
            merger_count             = 0;
      }
      if (inner_fluid_disk_bool){
            time_since_last_spawn = 0.0;
            fluid_disk_Sigma      = inner_mass/(M_PI*(Rroche*Rroche-Rearth*Rearth));
            typ x                 = Rroche/Rearth;
            typ x52               = fast_pow(sqrt(x), 5) - 1.0;
            typ x32               = fast_pow(sqrt(x), 3);
            typ x2                = x*x - 1.0;
            typ gx                = (4.0*x*x52-5.0*x2*x32)/(4.0*x52-5.0*x2);
            typ Rroche3           = Rroche*Rroche*Rroche;
            typ Omega             = sqrt(G*Mearth/Rroche3);
            time_between_spawn    = 16.0*M_PI*f_tilde*f_tilde*Rroche3*Rroche3*Omega*Omega*Omega*gx*(1.0-x)*(1.0-gx)/(Mearth*Mearth*x*G*G); //From Salmon & Canup 2012
            printf("Characteristic timescale between spawn = %.8lf\n", time_between_spawn);
      }
}


void array_initialization(){

      /******** Defines and initializes external arrays ********/
      
      
      
      /******** Array of moonlets used as buffer ********/
      xx = (struct moonlet *)malloc(sizeof(struct moonlet)*N_max);

      /******** The kth cell of this array contains 1 if the kth cell of the array moonlets contains a moonlet ********/
      exists = (int *)malloc(sizeof(int)*N_max);


      int p;

      /******** Initializing the array exists with 1 if the kth cell of the array moonlets contains a moonlet, 0 otherwise ********/
      for (p = 0; p < N_0; p++){
            *(exists+p) = 1;
      }
      for (p = N_0; p < N_max; p++){
            *(exists+p) = 0;
      }
      
      free_indexes = (int *)malloc(N_max*sizeof(int));
      
      /******** If the O(N) algorithm is used for collision detection, initializing the hash table and its dependencies ********/
      if (mesh_bool && (mutual_bool || collision_bool)){
      
            bigint n = collision_cube_cells+2; //Adding one layer to the collision cube so I don't have to worry about edges.
            typ hash_table_size=((typ) (n*n*n*(bigint) (sizeof(struct chain *))))/(1024.0*1024.0*1024.0); //Size of the hash table in GiB
            printf("Allocating %.1lf GiB of RAM to the hash table.\n", hash_table_size);
            hash=(struct chain **)malloc(sizeof(struct chain *)*n*n*n); //A pointer weighs 8 bytes, so this command should allocate 8GiB of RAM if collision_cube_cells = 1024
            if (hash==NULL){
                  fprintf(stderr, "Error : Can't allocate memory for the hash table in function array_initialization.\n");
                  abort();
            }
            for (p = 0; p < n*n*n; p++){
                  *(hash+p) = NULL;
            }
            
            modified_cells = (int *)malloc(sizeof(int)*N_max); //An array that contains the indexes of modified cells of the hash table. Initialized to {-1, -1, ..., -1}
            for (p = 0; p < N_max; p++){
                  *(modified_cells+p) = -1;
            }
            
            indexes = (int *)malloc(27*sizeof(int)); //An array that stores the indexes of the cells of the neighbourhood of a moonlet in the hash table
            for (p = 0; p < 27; p++){
                  *(indexes+p) = -1;
            }
            
            nghb = NULL;
            add(0,&nghb);
            nghb -> how_many = 0;
            
            /******** If there are mutual gravitational interactions, I allocate the array of pairs that need to be treated for mutual gravitational interactions ********/
            pairs = (struct pair *)malloc(sizeof(struct pair) * 5 * N_max * (int) how_many_neighbours); //Array of pairs of moonlets being in the same neighbourhood
                                                                                                        //The expected number of such pairs is N*how_many_neighbours/2
            how_many_pairs = 0;
            if (pairs == NULL){
                  fprintf(stderr, "Error : Can't allocate memory for the array pairs in function array_initialization.\n");
                  abort();
            }
      }
      
      three_largest_indexes      = (int *)malloc(3 * sizeof(int)); //Array of the indexes of the three largest moonlets

      if (Sun_bool){
            sun_vector = (typ *)malloc(3 * sizeof(typ));
            *sun_vector       = star_semi_major;
            *(sun_vector + 1) = 0.0;
            *(sun_vector + 2) = 0.0;
      }
      
      if (collision_bool){
            approach = (typ *)malloc(6 * sizeof(typ));
            did_collide = (int *)malloc(N_max * sizeof(int));
            for (p=0; p<N_max; p++){
                  *(did_collide+p) = 0;
            }
      }
      
      if ((mutual_bool || collision_bool) && !brute_force_bool && (falcON_bool || standard_tree_bool)){
            PeanoHilbertOrder = (int *)malloc(N_max * sizeof(int));
            already_in_tree   = (int *)malloc(N_max * sizeof(int));
            C1Moonlets = (typ *)malloc(3 * N_max * sizeof(typ));
            for (p=0; p <= largest_id; p++){
                  *(PeanoHilbertOrder + p) = p;
            }
      }
}


/*************************************/
/******** Memory deallocation ********/
/*************************************/


void deallocation(){

      /******** Deallocates the memory used globally ********/

      int p;
      
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
      if (Sun_bool){
            free(sun_vector);
            sun_vector = NULL;
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
}

/**********************************************************************************/
/******** Functions relative to chains (unrolled linked list) manipulation ********/
/**********************************************************************************/


void add(int head, struct chain ** ch){

      /******** Adds the new element head to the chain *ch ********/
      
      if (*ch==NULL){
            struct chain * to_be_returned=(struct chain *)malloc(sizeof(struct chain));
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
                  struct chain * to_be_returned=(struct chain *)malloc(sizeof(struct chain));
                  (to_be_returned -> ids)[0] = head;
                  to_be_returned -> how_many = 1;
                  to_be_returned -> queue = *ch;
                  *ch = to_be_returned;
            }
      }  
}


struct chain * Add(int head, struct chain * ch){

      /******** Same as above but with different data types ********/
      
      if (ch==NULL){
            struct chain * to_be_returned=(struct chain *)malloc(sizeof(struct chain));
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
                  struct chain * to_be_returned=(struct chain *)malloc(sizeof(struct chain));
                  (to_be_returned -> ids)[0] = head;
                  to_be_returned -> how_many = 1;
                  to_be_returned -> queue = ch;
                  return to_be_returned;
            }
      }
}


struct chain * delete_chain(struct chain * ch){

      /******** Deletes the first element of the chain ch ********/
      
      if (ch!=NULL){
            int hwmn = ch -> how_many;
            if (hwmn > 1){
                  ch -> how_many-=1;
                  return ch;
            }
            else{
                  struct chain * to_be_returned=NULL;
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
      
      if (ch!=NULL){
            int hwmn = ch->how_many;
            if (hwmn>1){
                  ch->how_many-=1;
                  return ch;
            }
            else{
                  if (ch->queue!=NULL){
                        struct chain * to_be_returned=NULL;
                        to_be_returned=ch->queue;
                        free(ch);
                        ch=NULL;
                        return to_be_returned;
                  }
                  else{
                        ch->how_many=0;
                        return ch;
                  }
            }
      }
      return NULL;
}


struct chain * delete_node(struct chain * ch){

      /******** Deletes a full node of the chain ch ********/
      
      
      struct chain * to_be_returned=NULL;
      if (ch!=NULL){
            to_be_returned=ch->queue;
            free(ch);
            ch=NULL;
            return to_be_returned;
      }
      return to_be_returned;

}


void clear_chain(struct chain ** ch){

      /******** Clears and deallocates the chain ch ********/
 
      
      while ((*ch)!=NULL){
            *ch=delete_node(*ch);
      }
      
}


/******************************************************************************/
/******** Miscellaneous functions relative to the numerical simulation ********/
/******************************************************************************/


void new_moonlet(struct moonlet * moonlets, int index, typ rad){

      /******** Creates a new moonlets of radius rad and orbital elements chosen at random like in the  ********/
      /******** function "populate". Stores it at the index index in the array moonlets. Updates the    ********/
      /******** array exists and the variable largest_id accordingly.                                   ********/
      
      typ a,e,i,nu,omega,Omega;
      a = rdm(sma_min, sma_max);
      e = rdm(eccentricity_min, eccentricity_max);
      i = rdm(inclination_min, inclination_max);
      nu = rdm(0.0, 2.0*M_PI);
      omega = rdm(0.0, 2.0*M_PI);
      Omega = rdm(0.0, 2.0*M_PI);

      *(moonlets+index) = init(a, e, i, nu, omega, Omega, rad);
      *(exists+index) = 1;
      if (index > largest_id){
            largest_id = index;
      }
      
}


int get_free_index(int should_be_put_at_the_end){

      /******** Returns the index of where to store a new moonlet in the simulation ********/
      /******** The new moonlet will be stored at *(moonlets+index)                 ********/
      
      
      int index;

      if (how_many_free == 0 || should_be_put_at_the_end){
            largest_id += 1;
            index = largest_id;
            if (index == N_max){
                  fprintf(stderr, "Error : Can't add a moonlet to the simulation. Try increasing the value of N_max.\n");
                  abort();
            }           
      }
      else {
            how_many_free -= 1;
            index = *(free_indexes+how_many_free);
            if (index > largest_id){
                  largest_id = index;
            }  
      }

      if (*(exists+index)){ //Checking that the supposedly free index is indeed free. To be removed when the code is robust
            fprintf(stderr, "Error : This index is not free.\n");
            abort();
      }

      *(exists+index) = 1;
      return index;

}


void tidy_up(struct moonlet * moonlets){

      /******** Reorders the array moonlets so that the N moonlets occupy the indexes 0 through N-1 ********/
      /******** This function is called at the end of the timestep, if N/largest_id < 0.9           ********/

      int i = 0;
      int j;
      
      for (j=0; j<=largest_id; j++){
            if (*(exists+j)){
                  if (j>i){
                        /******** Moonlet i becomes moonlet j ********/
                        *(moonlets+i) = *(moonlets+j);
                  }
                  i++;
            }
      }
      
      /******** Updating the array exists ********/
      for (j=0; j<=largest_id; j++){
            if (j<i){
                  *(exists+j)=1;
            }
            else{
                  *(exists+j)=0;
            }
      }
      
      largest_id = i-1;
      how_many_free = 0;
      
}


void three_largest_moonlets(struct moonlet * moonlets){

      /******** Actualizes the array three_largest_indexes with the indexes of the  ********/
      /******** currently three largest moonlets of the simulation. After a call to ********/
      /******** this function, *three_largest_indexes, *(three_largest_indexes+1)   ********/
      /******** and *(three_largest_indexes+2) are the indexes of the largest,      ********/
      /******** second largest and third largest moonlet of the simulation          ********/


      int i = 0;
      int j;
      typ R_0, R_1, R_2, R;
      
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
                              *(three_largest_indexes+1) = j;
                        }
                        else {
                              *(three_largest_indexes+1) = *three_largest_indexes;
                              *three_largest_indexes = j;
                              R = R_1;
                              R_1 = R_0;
                              R_0 = R;
                        }
                  }
                  
                  else if (i == 2){
                        R_2 = (moonlets + j) -> radius;
                        if (R_1 >= R_2){
                              *(three_largest_indexes+2) = j;
                        }
                        else if (R_2 > R_1 && R_2 <= R_0){
                              *(three_largest_indexes+2) = *(three_largest_indexes+1);
                              *(three_largest_indexes+1) = j;
                              R = R_2;
                              R_2 = R_1;
                              R_1 = R;
                        }
                        else {
                              *(three_largest_indexes+2) = *(three_largest_indexes+1);
                              *(three_largest_indexes+1) = *three_largest_indexes;
                              *three_largest_indexes = j;
                              R = R_2;
                              R_2 = R_1;
                              R_1 = R_0;
                              R_0 = R;
                        }
                  }
                  
                  /******** At this stage, the three first moonlets indexes occupy three_largest_indexes by decreasing radius ********/
                  
                  else { // If i > 2
                        R = (moonlets + j) -> radius;
                        if (R > R_0){
                              *(three_largest_indexes+2) = *(three_largest_indexes+1);
                              *(three_largest_indexes+1) = *three_largest_indexes;
                              *three_largest_indexes = j;
                              R_2 = R_1;
                              R_1 = R_0;
                              R_0 = R;
                        }
                        else if (R > R_1){
                              *(three_largest_indexes+2) = *(three_largest_indexes+1);
                              *(three_largest_indexes+1) = j;
                              R_2 = R_1;
                              R_1 = R;
                        }
                        else if (R > R_2){
                              *(three_largest_indexes+2) = j;
                              R_2 = R;
                        }
                  }
                  i++;
            }
      }
}


void three_largest_three_first(struct moonlet * moonlets){

      /******** Makes sure that the three largest moonlets are the three first moonlets ********/
      
      struct moonlet buffer; // A buffer used to temporarily store a moonlet
      int first, second, third; //The indexes of the first, second and third largest moonlet
      
      first = *three_largest_indexes;
      second = *(three_largest_indexes+1);
      third = *(three_largest_indexes+2);
      
      if (first != 0){
            if (*exists){ //There is a moonlet at index zero, but it is not the largest moonlet. I exchange it with the largest moonlet
                  buffer = *moonlets;
                  *moonlets = *(moonlets+first);
                  *(moonlets+first) = buffer;
            }
            else { //Index zero is not occupied. I occupy it with the largest moonlet
                  *moonlets = *(moonlets+first);
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
            if (*(exists+1)){ //There is a moonlet at index one, but it is not the second largest moonlet. I exchange it with the second largest moonlet
                  buffer = *(moonlets+1);
                  *(moonlets+1) = *(moonlets+second);
                  *(moonlets+second) = buffer;
            }
            else { //Index one is not occupied. I occupy it with the second largest moonlet
                  *(moonlets+1) = *(moonlets+second);
                  lose_moonlet(second);
                  *(exists+1) = 1;
            }
            if (third == 1){
                  third = second;
            }
      }
      if (third != 2){
            if (*(exists+2)){ //There is a moonlet at index two, but it is not the third largest moonlet. I exchange it with the third largest moonlet
                  buffer = *(moonlets+2);
                  *(moonlets+2) = *(moonlets+third);
                  *(moonlets+third) = buffer;
            }
            else { //Index two is not occupied. I occupy it with the third largest moonlet
                  *(moonlets+2) = *(moonlets+third);
                  lose_moonlet(third);
                  *(exists+2) = 1;
            }
      }

}


void cross_product(typ u_1, typ u_2, typ u_3, typ v_1, typ v_2, typ v_3, typ * uxv){

      /******** Fills uxv with the coordinates of the cross product u x v ********/
      
      *uxv     = u_2*v_3-u_3*v_2;
      *(uxv+1) = u_3*v_1-u_1*v_3;
      *(uxv+2) = u_1*v_2-u_2*v_1;
      
}


void lose_moonlet(int a){

      /******** Treats the loss of moonlet a from the simulation ********/
      
      *(free_indexes + how_many_free) = a;
      how_many_free ++;
      *(exists + a) = 0;
      if (a <= 2 && mesh_bool && !force_naive_bool){
            how_many_free --;
      }
}


int maximum(typ i, typ j, typ k){

      /******** Returns the position of the maximum of (i,j,k) ********/
      /******** For example, if max(i,j,k)=k, then returns 2   ********/
      
      int m;
      typ mx;
      if (i >= j){
            m = 0;
            mx = i;
      }
      else {
            m = 1;
            mx = j;
      }
      
      if (mx >= k){
            return m;
      }
      else {
            return 2;
      }
}




