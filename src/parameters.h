#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_



/*****************************************************************************************/
/******** Redefining keywords. Do not touch if you are fine with double precision ********/
/*****************************************************************************************/

#define typ double                   //Renaming double as typ. Define typ as long double (resp. float) for quadruple precision (resp. simple precision)
#define bigint long int              //Renaming long int as bigint. If sizeof(int) == sizeof(long int) == 4 on your system, then define bigint as long long int instead
                                     //Open the file /usr/share/gtksourceview/language-specs/c.lang and then find the field
                                     //<context id="types" style-ref="type">. Add the line <keyword>typ</keyword> and <keyword>bigint</keyword> to it
                                     //and color gedit with C. The keywords typ and bigint should now be colored after a reboot
#define absolute fabs                //The function returning the absolute value of a typ. Choose fabs (resp. fabsf or fabsl) if typ is double (resp. float or long double)
#define integral floor               //The function returning the floor of a typ. Choose floor (resp. floorf or floorl) if typ is double (resp. float or long double)



/*********************************************************************************************************/
/******** Defining the path where data have to be output. Unimportant if write_to_files_bool is 0 ********/
/******** Put / at the end of the path                                                            ********/
/*********************************************************************************************************/

#define path "/run/media/jeremy/Windows/Users/miste/linux/Moon/dump2/"



/**************************************************/
/******** Defining some physical constants ********/
/**************************************************/

#define Rearth 1.0                   //Radius of the Earth (or central mass) is the unit of length
#define Mearth 1.0                   //Mass of the Earth (or central mass) is the unit of mass
#define G 39.47841760435743          //Gravitational constant is set to 4*pi^2, so that the unit of time is the surface orbital period
#define density 0.1448               //The density of the moonlets (3344 kg/m^3) in Mearth/Rearth^3.
#define density_earth 0.238764       //The density of the Earth    (5514 kg/m^3) in Mearth/Rearth^3.
#define Rroche 2.9                   //The Roche radius
#define Tearth 3.4076                //Earth's sideral period in units of the surface orbital period. Must be larger than 1. Today value is 17.038



/*******************************************************/
/******** Parameters relative to the simulation ********/
/*******************************************************/

#define N_max 10000                  //The maximum number of moonlets the simulation can handle
#define N_0 1000                     //The initial number of moonlets. Must be less than or equal to N_max.
#define M_0 0.02214                  //Total (expected) moonlet mass at t=0
#define t_end 512.0                  //Total simulation time (in surface orbital period)
#define time_step 0.0078125          //Timestep of the simulation.
#define output_step 128              //Output occurs every output_step timestep. Matters only if write_to_files_bool is 1
#define radius_stddev 0.57           //Standard deviation of moonlet's radii at t=0 (drawn uniformly), in units of the mean radius. Must be less than 1/sqrt(3) to prevent negative radius
#define eccentricity_min 0.0         //Minimal eccentricity             for a moonlet at t=0
#define eccentricity_max 0.2         //Maximal eccentricity             for a moonlet at t=0
#define sma_min 2.9                  //Minimal semi-major axis          for a moonlet at t=0
#define sma_max 14.0                 //Maximal semi-major axis          for a moonlet at t=0
#define inclination_min 0.0          //Minimal inclination (in radians) for a moonlet at t=0
#define inclination_max 0.174533     //Maximal inclination (in radians) for a moonlet at t=0
#define low_dumping_threshold 2.0    //Moonlets falling below this threshold (in Earth radii) are dumped from the simulation (collision with the Earth or disruption by tidal forces)
#define high_dumping_threshold 200.0 //Moonlets going  beyond this threshold (in Earth radii) are dumped from the simulation (assumed unbounded)
#define max_ids_per_node 173         //The maximum number of ids in each node of the unrolled linked lists (chains). Choose such that sizeof(struct chain) be a multiple of the cache line
#define softening_parameter 0.00     //The softening parameter for mutual gravitational interations, in units of the sum of the radii.
#define seed 778345128               //The seed used for random number generation. Does not matter if seed_bool is 0.
#define switch_to_brute_force 512    //Threshold for N below which the program switches to the brute-force method for mutual interactions. Does not matter if brute_force_bool is 1



/****************************************************************/
/******** Parameters relative to the mesh-grid algorithm ********/
/****************************************************************/

/******** Parameters relative to collision detection with the mesh algorithm ********/
#define collision_cube_min 50.0      //Minimal sidelength of the cube centered on the Earth where collisions are looked for, for the minimal size we allow for the mesh-size
#define collision_cube_cells 1024    //The number of mesh cells per dimension of the collision cube. Must be even. For 16+ GiB of RAM (resp. 8 or 4 GiB), choose ~1000 (resp. ~800 or ~500)
#define how_many_neighbours 16.0     //The desired expected number of neighbours for a moonlet



/**************************************************************************************************/
/******** Parameters relative to tree-based algorithms (falcON and the standard tree code) ********/
/**************************************************************************************************/

#define expansion_order 3            //The order p of the Taylor expansions. This is p = 3 in Dehnen (2002). NcorpiON allows up to p = 6. Minimum is 1 since order 0 yields no acceleration
#define theta_min 0.5                //Minimal value of the tolerance parameter theta. Must be strictly less than 1. Sensible values are 0.2 < theta_min < 0.7
                                     //The precision (and computational time) of the mutual gravity computed by Ncorpion increases with increasing p and decreasing theta_min.
#define subdivision_threshold 17     //A cubic cell is not divided as long as it contains at most that many moonlets. Called s in Dehnen (2002). Must be > 0
                                     //The computational time depends a lot on this threshold. Suggested values are s ~ 26 if p = 3 and s ~ 80 if p = 6, but should be tweaked by the user.
#define root_sidelength 80.0         //Sidelength of the root cell. Must be machine representable. Particles outside of the root cell only feel Earth's gravity and do not affect others
#define level_max 25                 //The maximum allowed number of levels in the tree. Root is at level 0 and a cell at level level_max - 1 is never divided.
#define child_multipole_threshold 1  //If number of particules/number of children is at most this threshold then the multipole moments of a cell are computed directly from the particules.
                                     //Otherwise, they are computed from that of the children

/******** Parameters relative to mutual gravity computation with falcON algorithm ********/
#define N_cc_pre 8                   //For two cells with N1 and N2 moonlets, if N1N2 < N_cc_pre, then the interaction is computed brute-forcely, no matter if they are well-separated or not
#define N_cc_post 64                 //For two cells with N1 and N2 moonlets, if N1N2 < N_cc_post and they are not well-separated, then the interaction is computed brute-forcely
#define N_cs 64                      //If N1 < N_cs, the self-interaction of a cell containing N1 moonlets is computed brute-forcely
                                     //See TreeWalk algorithm in the draft for details about these thresholds. Since the computational cost of the multipole method increases greatly with the
                                     //expansion order, the values of N_cc_pre, N_cc_post and N_cs should be adapted to the expansion order p. Note, however, that the computational cost
                                     //depends more on the subdivision threshold s that on these thresholds. Suggested values are (s, N_cc_pre, N_cc_post, N_cs) = (26, 8, 64, 64) for p = 3
                                     //and (s, N_cc_pre, N_cc_post, N_cs) = (80, 256, 1024, 128) for p = 6 but this probably depends on the architecture and should be tweaked by the user.
                                     
/******** Parameters relative to mutual gravity computation with the standard tree code ********/
#define N_cb_pre 6                   //If N1 < N_cb_pre, the interaction between a moonlet and a cell with N1 moonlets is performed brute-forcely, no matter if they are well-separated or not
#define N_cb_post 16                 //If N1 < N_cb_post and they are not well-separated, the interaction between a moonlet and a cell with N1 moonlets is performed brute-forcely

/******** Parameters relative to collision detection with falcON algorithm ********/
#define N_cc_collision 16            //This is N_cc_post. For two non well-separated nodes with Na and Nb moonlets, if NaNb < N_cc_collision, then collisions are searched brute-forcely,
                                     //otherwise, the largest node is subdivised. N_cc_pre is 0 when falcON is used for collision detection.
#define N_cs_collision 12            //For a node with Na moonlets, if Na < N_cs_collision, then collisions are searched brute-forcely, otherwise, the node is subdivised

/******** Parameters relative to collision detection with the standard tree code ********/
#define N_cb_collision 16            //This is N_cb_post. For a moonlet not well-separated from a node with Na moonlets, if Na < N_cb_collision, then collisions are searched brute-forcely,
                                     //otherwise, the node is subdivised. N_cb_pre is 0 when the standard tree code is used for collision detection.



/*****************************************************************************************/
/******** Parameters relative to collision resolution and the fragmentation model ********/
/*****************************************************************************************/

/******** Collision resolution when all collisions are inelastic ********/
#define collision_parameter 1.1      //The collision parameter f to model inelastic collision. Must be between 1 (completely inelastic) and 2 (elastic)

/******** Collision resolution with the fragmentation model described in the PDF draft ********/
#define beta_slope 2.6470588235294118//Slope of the power law for fragment size distribution in Leinhardt and Stewart (2012). Must be such that N_tilde is integer.
#define N_tilde 15                   //2*beta_slope/(3 - beta_slope). This is the ratio between the ejected mass and the mass of the second largest fragment -> N° of fragments in the tail
#define mu_parameter 0.55            //The exponent of the impact velocity in the coupling parameter. See Table 3 of Housen & Holsapple (2011)
#define C1_parameter 1.5             //A dimensionless parameter of impact theories. See Table 3 of Housen & Holsapple (2011)
#define k_parameter 0.2              //A dimensionless parameter of impact theories. See Table 3 of Housen & Holsapple (2011)
#define frag_threshold 0.000000005   //Ejected mass threshold below which the collision results in a merger. If the ejected mass is above that threshold but the mass of the second
                                     //largest fragment is below, then the tail is reunited into a single moonlet. If the mass of the second largest fragment is above that threshold
                                     //then the collision results in a full fragmentation.                             
#define pq_min_max {-1,3,-1,1}       //Extremal integer values for p_k and q_k to determine the position of the tail fragments with respect to the largest fragment.
                                     //Must define a rectangle containing exactly N_tilde points with integer coordinates. More precisely, if pq_min_max = {a,b,c,d},
                                     //then we must have N_tilde = (b-a+1)(d-c+1). See PDF draft for details                                 
                                     


/****************************************/
/******** Defining some booleans ********/
/****************************************/

#define write_to_files_bool      1   //Determines if the simulation writes to output files. Set to 0 to run speed tests, or if you are satisfied with what is displayed in the terminal   
#define seed_bool                0   //Determines if the seed for random number generation is chosen by the user. If seed_bool is 0, then the seed is the number of seconds since 01/01/1970
#define J2_bool                  1   //Determines if the contribution from the symmetrical equatorial bulge is taken into account in the simulation
#define Sun_bool                 0   //Determines if the perturbations from the Sun are taken into account in the simulation
#define disk_bool                0   //Determines if the perturbations from the inner fluid disk are taken into account in the simulation
#define collision_bool           1   //Determines if the moonlets are able to collide
#define mutual_bool              1   //Determines if there are mutual gravitational interactions between the moonlets.
                                     
/******** Boolean relative to conservation of the total momentum or total angular momentum ********/
#define tam_bool                 0   //Determines if the total angular momentum should be conserved upon merging or fragmenting impact. If tam_bool is 0 then the total momentum is conserved

/******** Booleans relative to collision resolution. Exactly one of them must be 1 ********/
#define elastic_collision_bool   0   //Determines if the collisions are all elastic.
#define inelastic_collision_bool 0   //Determines if the collisions are all inelastic with inelasticity parameter collision_parameter.
#define instant_merger_bool      0   //Determines if the collisions all result in a merger.
#define fragmentation_bool       1   //Determines if the outcome of the collision should be either a merge, a partial fragmentation, a full fragmentation, or a catastrophic disruption,
                                     //depending on the relative velocity and according to the fragmentation model of NcorpiON. 

/******** Booleans relative to the treatment of mutual interactions (collisions and self-gravity). Exactly one of them must be 1. ********/
#define brute_force_bool         0   //Determines if a brute-force method should be used for mutual interactions.
#define falcON_bool              1   //Determines if falcON algorithm     should be used for mutual interactions. (Often the best speed/accuracy compromise for large N)
#define standard_tree_bool       0   //Determines if a standard tree code should be used for mutual interactions. (Significantly slower than falcON for the same accuracy).
#define mesh_bool                0   //Determines if the  mesh algorithm  should be used for mutual interactions. (Fastest but gravity with neighbours and three largest moonlets only).



#endif
