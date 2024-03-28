/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    structure.h                                                 ********/
/******** @brief   Header file to structure.c                                  ********/
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


#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include "parameters.h"


/******** Defining some external variables ********/
extern typ J2;                       //Earth's equatorial bulge flattening parameter
extern typ timestep;                 //The timestep
extern int largest_id;               //The largest if of a moonlet. This variable allows us not to travel the whole moonlet array
extern int first_passage;            //To be removed later. Just for the benchmark
extern typ time_elapsed;             //The physical time elapsed since the beginning of the simulation (in units of the surface orbital period)
extern int how_many_free;            //Total number of indexes used in priority when a new moonlet is created
extern int how_many_moonlets;        //Total number of moonlets in the simulation
extern int force_naive_bool;         //If 1, then forces the program to use the brute-force algorithm for collision detection. Automatically set to 1 by the program when N becomes small
extern typ time_until_collision;     //Time until collision
extern typ time_since_last_spawn;    //The time elapsed since a moonlet last spawned from the inner fluid disk
extern typ time_between_spawn;       //The characteristic timescale between moonlet spawning
extern typ fluid_disk_Sigma;         //The mass of the inner fluid disk per unit area
extern typ SideralOmega;             //The sideral rotation frequency of the central body
extern typ star_mean_motion;         //The mean motion of the central body around its star or companion star



/******** Defining the moonlet structure ********/
struct moonlet {
      typ x;              //Moonlet geocentric position (or velocity if derived with respect to time)
      typ y;
      typ z;
      typ vx;             //Moonlet geocentric velocity (or acceleration if derived with respect to time)
      typ vy;
      typ vz;
      typ mass;           //Moonlet mass.
      typ radius;         //Moonlet radius. While the radius can be computed from the mass, storing the radius reduces computation time, at the cost of more memory usage.
                          //This does not even yield more memory usage due to padding if the field radius is removed. The structure moonlet weighs 64 bytes if typ is double.
};


/******** Global arrays of moonlets used as buffer ********/
extern struct moonlet * xx;


/******** Global arrays ********/
extern int * exists;                          //The kth cell contains 1 if *(moonlets+k) contains a moonlet, and zero else
extern int * free_indexes;                    //Arrays of ids of non-existant moonlets


/****************************************************************************************/
/******** Defining a cache-friendly chain structure that holds ids of moonlets   ********/
/****************************************************************************************/


/******** The regular definition of a chain (or linked list) is                  ********/
/********                                                                        ********/
/******** struct chain {                                                         ********/
/********       struct chain * queue;                                            ********/
/********       int id;                                                          ********/
/******** };                                                                     ********/
/********                                                                        ********/
/******** but this definition is such that two consecutives ids are far away in  ********/
/******** memory. In a cache friendly chain, most ids consecutive in the chain   ********/
/******** are consecutive in memory too, which favours cache hits and reduces    ********/
/******** calls to the functions malloc and free                                 ********/ 
struct chain {
      int ids[max_ids_per_node]; //The array of ids of moonlets in that node of the chain. Choosing max_ids_per_node = 13 ensures sizeof(chain) = 64 = one cache line
      struct chain * queue;      //A pointer towards the queue of the chain
      int how_many;              //How many ids are in this node
};


/******** Datas for time-efficient collision detection ********/
extern struct chain ** hash;            //Hash table
extern int * modified_cells;            //Array containing the cells of the hash table that were modified during the current timestep.
extern int how_many_modified;           //The number of cells in the hash tables that were modified during the current timestep
extern typ gam;                         //The mesh-size updated every timestep to ensure an average number of neighbours of how_many_neighbours
extern typ gam_min;                     //The minimum possible value of the mesh-size, defined as collision_cube_min/collision_cube_cells
extern int how_many_big;                //The number of moonlets larger  than the mesh-size
extern int how_many_small;              //The number of moonlets smaller than the mesh-size
extern int * indexes;                   //Array containing the indexes of the cells of the neighbourhood of a moonlet in the hash table
extern typ collision_cube;              //The sidelength of the collision cube, defined as collision_cube_cells*gam
extern int total_neighbours;            //The total number of neighbours visited during the current timestep
extern typ average_neighbours;          //total_neighbours/how_many_small
extern struct chain * nghb;             //The neighbours of the current moonlet
extern typ * approach;                  //Array containing the position of two colliding moonlets

struct pair {                           //Simple data structure defining a pair of moonlets. Contains a pair of close moonlets when mutual_bool is 1. 
      int fst;
      int snd;
};
extern struct pair * pairs;             //Array of pairs of moonlets undergoing mutual gravitational interactions (being in each-other neighbourhood). Used only if mutual_bool is 1.
extern int how_many_pairs;              //Number of pairs contained in the array pairs.
extern int * three_largest_indexes;     //The indexes of the three largest moonlets of the simulation.


/******** Datas for fragmentation ********/
extern int super_catastrophic_count;      //Number of collisions having resulted in a super-catastrophic event so far
extern int half_fragmentation_count;      //Number of collisions such that ejected mass > frag_threshold > ejected mass / N_tilde
extern int full_fragmentation_count;      //Number of collisions such that ejected mass / N_tilde > frag_threshold
extern int merger_count;                  //Number of collisions having resulted in a merger, that is, such that ejected mass < frag_threshold
extern int collision_count;               //Number of collisions so far
extern int * did_collide;                 //The k^th cell contains 1 if moonlet k did collide during that timestep, and 0 otherwise


/******** array relative to pertubations from the star or companion star ********/
extern typ sun_vector[3];                 //The position of the star in the geocentric reference frame. Unused if Sun_bool is 0


/******** Global array of output files ********/
extern FILE ** files;




typ * ell2cart(typ a, typ e, typ i, typ nu, typ omega, typ Omega);


void cart2aei(struct moonlet * moonlets, int id, typ * aei);


struct moonlet init(typ a, typ e, typ i, typ nu, typ omega, typ Omega, typ rad);


typ rdm(typ min, typ max);


struct moonlet * populate(typ M, typ dR);


void variable_initialization();


void array_initialization();


void deallocation();


void add(int head, struct chain ** ch);


struct chain * Add(int head, struct chain * ch);


struct chain * delete_chain(struct chain * ch);


struct chain * partial_delete(struct chain * ch);


struct chain * delete_node(struct chain * ch);


void clear_chain(struct chain ** ch);


void new_moonlet(struct moonlet * moonlets, int index, typ rad);


int get_free_index(int should_be_put_at_the_end);


void tidy_up(struct moonlet * moonlets);


void three_largest_moonlets(struct moonlet * moonlets);


void three_largest_three_first(struct moonlet * moonlets);


void cross_product(typ u_1, typ u_2, typ u_3, typ v_1, typ v_2, typ v_3, typ * uxv);


void lose_moonlet(int a);


int maximum(typ i, typ j, typ k);


#endif
