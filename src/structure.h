#ifndef _STRUCTURE_H_
#define _STRUCTURE_H_

#include "parameters.h"


/******** Defining some external variables ********/
extern typ J2;                       //Earth's equatorial bulge flattening parameter
extern typ timestep;                 //The timestep
extern int largest_id;               //The largest if of a moonlet. This variable allows us not to travel the whole moonlet array
extern int first_passage;            //To be removed later. Just for the benchmark
extern typ time_elapsed;             //The physical time elapsed since the beginning of the simulation (in units of the surface orbital period)
extern typ time_since_last_spawn;    //The physical time elapsed since a moonlet last spawned from the inner fluid disk
extern typ inner_fluid_mass;         //Mass of the inner fluid disk
extern int how_many_free;            //Total number of indexes used in priority when a new moonlet is created
extern typ spawned_mass_factor;      //Factor used to compute the mass   of a spawned moonlet
extern typ spawned_radius_factor;    //Factor used to compute the radius of a spawned moonlet
extern typ inner_disk_mass_factor;   //Factor used to compute the variation in the inner fluid disk's mass
extern int how_many_moonlets;        //Total number of moonlets in the simulation
extern int force_naive_bool;         //If 1, then forces the program to use the brute-force algorithm for collision detection. Automatically set to 1 by the program when N becomes small
extern typ time_until_collision;     //Time until collision



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
                          //This does not even yield more memory usage due to padding if the field radius is removed. The structure moonlet weighs 64 bytes in both cases.
};


/******** Global arrays of moonlets used as buffer in the function iteration and leapfrog ********/
extern struct moonlet * xx;
extern struct moonlet * k1;
extern struct moonlet * k2;
extern struct moonlet * k3;
extern struct moonlet * k4;
extern struct moonlet * mid_point_speed;


/******** Global arrays ********/
extern int * exists;                          //The kth cell contains 1 if *(moonlets+k) contains a moonlet, and zero else
extern struct chain * to_be_added_fluid_disk; //Array of ids of moonlets that need to be put in the inner fluid disk at the end of the timestep
extern int * free_indexes;                    //Arrays of ids of non-existant moonlets


/******** Defining a cache-friendly chain structure that holds ids of moonlets   ********/


/******** The regular definition of a chain (or linked list) is                  ********/

/******** struct chain {                                                         ********/
/********       struct chain * queue;                                            ********/
/********       int id; };                                                       ********/

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
extern int * catastrophic_pdf;            //Probability density function of the catastrophicity of collisions
extern typ * tam_loss;                    //Total angular momentum loss due to dumped moonlets
extern int * did_collide;                 //The k^th cell contains 1 if moonlet k did collide during that timestep, and 0 otherwise

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


void catastrophicity(typ m_tilde, typ M);



#endif