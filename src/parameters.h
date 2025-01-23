/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    parameters.h                                                ********/
/******** @brief   Main header file. The user who does not modify NcorpiON     ********/
/********          shall only modify this file                                 ********/
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


#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_



/***************************************************************************************************/
/******** Defining the input/output path where data are output if write_to_files_bool is 1  ********/
/******** and where initial conditions are read if random_initial_bool is 0. Put / at the   ********/
/******** end of the path. The path must exist and be absolute. This is also where the file ********/
/******** init.txt used to resume a simulation is written when resume_simulation_bool is 1  ********/
/***************************************************************************************************/

#define pth "/path/towards/input/output/location/" //The path must be absolute, end with / and already exist



/**************************************************************************/
/******** Defining booleans to determine the behaviour of NcorpiON ********/
/**************************************************************************/

/******** General booleans relative to the simulation ********/
#define write_to_files_bool      0   //Determines if the simulation writes to output files. Set to 0 to run speed tests, or if you are satisfied with what is displayed in the terminal
                                     //You can also set this boolean to 0 if you only want to 3D visualize the simulation in your browser.
#define make_animation_bool      0   //Determines if animations of the simulation are produced. write_to_files_bool and write_elliptic_bool must both be set to 1
#define write_cartesian_bool     0   //Determines if the cartesian elements x, y, z, vx, vy, vz   should be output. Unimportant if write_to_files_bool is 0. Output in simulation's units
#define write_elliptic_bool      0   //Determines if the elliptic  elements a, lambda, k, h, q, p should be output. Unimportant if write_to_files_bool is 0. Output in simulation's units
#define write_collisions_bool    1   //Determines if statistics regarding the collisions are output. Unimportant if not both write_to_files_bool and collision_bool are 1
#define central_mass_bool        1   //Determines if there is a central mass. If 0, none of the bodies in the simulation play any particular role. If 1, the central body is initially
                                     //at (x,y,z,vx,vy,vz) = {0.0} and treated independently. Should be set to 1 if one body is very massive and if brute_force_bool is 0, so gravity
                                     //with that body is computed directly and not from a tree or a mesh. If 1, some physical effects can be considered as well(see below).
#define reduce_to_COM_bool       1   //Determines if the position and speed of the center of mass of the system are cancelled by a translation before simulating. If 0, then the center
                                     //of mass will drift linearly with time, which could be a problem due to the limit of floating-point representation. Setting to 1 is recommended.
                                     //Note that the position and speed of the central body (if any) will slightly deviate from 0 after the reduction of the center of mass.
#define random_initial_bool      1   //Determines if the initial conditions are drawn randomly between bounds defined below. If 0, the initial conditions are read from the file init.txt
                                     //located at the path defined above, with one line per body and 8 columns, the first six for the initial conditions and the last two for mass and
                                     //radius (in simulation unit). If central_mass_bool is 1 and random_initial_bool is 0, do not include the central body in init.txt.
#define initial_cartesian_bool   1   //Determines if the initial conditions are given in cartesian coordinates (X,Y,Z,vX,vY,vZ). If 0, then the initial conditions are (semi-major axis,
                                     //eccentricity, inclination, true anomaly, argument of periapsis, longitude of ascending node). Unimportant if random_initial_bool is 1.
                                     //The initial conditions in the file init.txt must be given in simulation's units and radians.
#define seed_bool                1   //Determines if the seed for random number generation is chosen by the user. If seed_bool is 0, the seed is the number of seconds since 01/01/1970
#define one_collision_only_bool  0   //Determines if bodies are only allowed to collide once per timestep. If 0, there is no restriction on the number of collisions a body can experience
                                     //during a timestep. Setting first to 1 and then to 0 is a good way to know if the timestep is adapted to the bodies' mean free path.
#define openGL_bool              0   //MATCH TO THE SAME VALUE IN THE MAKEFILE. Determines if a 3D real-time visualization of the simulation with WebGL in REBOUND is enabled.
#define resume_simulation_bool   0   //Determines if, at the end of the simulation, NcorpiON generates a file named init.txt that can be used to resume the simulation. The file init.txt
                                     //is stored at the path indicated above. To resume the simulation, you need to set random_initial_bool to 0, initial_cartesian_bool to 1, and N_0 to
                                     //the number of lines of init.txt. If init.txt already exists in path pth, it will be overwritten. Simulation's variables should be updated.
#define viscoelastic_bool        0   //Determines if NcorpiON is used to simulate a viscoelastic body. The body is discretized by N_0 points linked by Kelvin-Voigt models. If random_
                                     //initial_bool is 1, then a file pth/shape_model.txt should be provided to define the boundary of of the viscoelastic body. If it is not provided, a
                                     //sphere of radius R_unit will be taken. If random_initial_bool is 0, the body is retrieved from the files init.txt and connections.txt that can have
                                     //been created by a previous simulation where resume_simulation_bool was set to 1.

/******** Booleans relative to interactions with the central mass. Set all to 0 if central_mass_bool is 0 ********/
#define J2_bool                  0   //Determines if the contribution from the J2 is taken into account in the simulation. The (x,y) plane of the simulation must be the equatorial plane
#define central_tides_bool       0   //Determines if orbiting bodies raise tides on the central body. The tidal model used by NcorpiON is the constant timelag model
#define inner_fluid_disk_bool    0   //Determines if there is an inner fluid disk (disk of liquid material below the Roche radius from which bodies spawn). See Salmon & Canup 2012
                                     //If set to 1, its mass is added to that of the central body when computing gravitational interactions and when preserving the total momemtum

/******** Booleans relative to mutual interactions between the bodies ********/
#define collision_bool           1   //Determines if the bodies are able to collide. If 0, the bodies pass through each other and you should set a non-zero value for softening_parameter
#define mutual_bool              1   //Determines if there is mutual gravity between the bodies. If central_mass_bool is 1, the bodies interact with the central mass even if set to 0.

/******** Booleans relative to the treatment of mutual interactions (collisions and self-gravity). Exactly one of them must be 1. Unimportant if mutual_bool is 0 ********/
#define brute_force_bool         0   //Determines if a brute-force method should be used for mutual interactions.
#define falcON_bool              1   //Determines if falcON algorithm     should be used for mutual interactions. (Best speed/accuracy compromise for large N + preserve total momentum)
#define standard_tree_bool       0   //Determines if a standard tree code should be used for mutual interactions. (Much slower than falcON for the same accuracy, momentum not preserved).
#define mesh_bool                0   //Determines if the  mesh algorithm  should be used for mutual interactions. (Fastest but gravity with neighbours and three largest bodies only).

/******** Booleans relative to collision resolution. Exactly one of them must be 1. Unimportant if collision_bool is 0 ********/
#define elastic_collision_bool   0   //Determines if the collisions are all elastic.
#define inelastic_collision_bool 0   //Determines if the collisions are all inelastic with inelasticity parameter collision_parameter.
#define instant_merger_bool      0   //Determines if the collisions all result in a merger.
#define fragmentation_bool       1   //Determines if the outcome of the collision should be either a merge, a partial fragmentation, a full fragmentation, or a catastrophic disruption,
                                     //depending on the relative velocity. Set to 1 to use the fragmentation model of NcorpiON.



/**************************************************/
/******** Defining some physical constants ********/
/**************************************************/

/******** Defining a system of units for the simulation. If random_initial_bool is 0, then the units in the file init.txt must be the simulation's units ********/
#define R_unit 1.                    //If central_mass_bool is 1, then radius of the central body. If viscoelastic_bool is 1, then radius of the viscoelastic body
#define M_unit 1.                    //If central_mass_bool is 1, then   mass of the central body. If viscoelastic_bool is 1, then mass   of the viscoelastic body
#define G 39.47841760435743447533796 //The gravitational constant. If you set it to 4*pi^2, then a particle of semi-major axis 1 orbiting a mass 1 has a period 1.
                                     //If viscoelastic_bool is 1, M_unit and R_unit are the mass and mean radius of the viscoelastic body to be simulated.
                                     //Note that R_unit and M_unit do not necessarily have to be 1, they should just be equal to the mass and radius of the central body or viscoelastic
                                     //body in the system of units that you use. It is generally advised to use a system of units such that the simulation will not have to manipulate
                                     //absurdely large or small floating points numbers. If both central_mass_bool and viscoelastic_bool are 0, then you don't have to define R_unit and
                                     //you should define M_unit as the mass, in your system of units, that is used to convert cartesian coordinates into elliptic coordinates.
                                     //Regardless of central_mass_bool and viscoelastic_bool, you have to define G as the value of the gravitational constant in your system of units.

/******** Physical constants relative to interactions with the central body (J2, inner disk, central tides) or a distant object. Unimportant if central_mass_bool is 0 ********/
#define Tearth 3.556721809           //Central body's sideral period in units of the surface orbital period. Must be > 1.
                                     //In case of tides, the sideral period changes and this is the value at initial time.
#define J2_value 0.                  //The J2 of the central body. If you choose J2_value = 0.0, then J2 is obtained from J2 = 1/2*Omega^2/Omega_crit^2 (fluid body) where Omega is the
                                     //sideral frequency and Omega_crit = sqrt(G*M_unit/R_unit^3). In that case, J2 is variable throughout the simulation.
#define k2 1.5                       //Second Love number of the central body. Here the value is for a fluid body (zero shear modulus). The constant timelag model is used
#define Delta_t 0.0015807            //The timelag between tidal stress and response. In simulation's units
#define dimensionless_moi 0.3307     //The moment of inertia of the central body, in simulation units (in units of its mass times its radius squared).
#define inner_mass 0.                //Mass of the inner fluid disk at initial time.
#define spawned_density 0.1448       //Density of the bodies that spawn from the inner fluid disk, in simulation's units.
#define f_tilde 0.3                  //A parameter controlling the mass of bodies spawned from the inner fluid disk. Must be < 1. Salmon & Canup (2012) choose 0.3
#define R_roche 2.9                  //The initial outer radius of the inner fluid disk where bodies spawn (in simulation's units). Must be larger than disruption_threshold and R_unit.
#define disruption_threshold 1.8     //Threshold (in simulation's units) below which bodies are tidally disrupted by the central mass. If inner_fluid_disk_bool is 0, then the mass of
                                     //the dumped body is added to the central body. Otherwise, the mass of the dumped body is added to the inner fluid disk if the body's periapsis is
                                     //above the surface or if it will cross the xy plane before hitting the surface, and to the central mass else.
                                     
/******** Orbit of a point-mass perturbator, in an inertial reference frame. Set pert_mass to 0.0 if you do not want a perturbator ********/
#define pert_sma 23481.0897818238895 //The perturbator is on a Keplerian trajectory defined by the six elements (semi-major axis, eccentricity, inclination, true anomaly, argument of
#define pert_ecc 0.                  //periapsis, longitude of ascending node), given by these six parameters (in radians and simulation's units). The gravitational parameter used to
#define pert_inc 0.                  //convert these elliptic elements into cartesian coordinates is mu = G*(pert_mass + M_unit). The Keplerian orbit can be hyperbolic
#define pert_tra 0.                  //(pert_ecc can exceed 1), but then, the semi-major axis pert_sma must be negative and the true anomaly must verify |pert_tra| < acos(-1/e).
#define pert_aop 0.                  //The eccentricity pert_ecc cannot be exactly equal to 1. The true anomaly is given at the time t = 0, not at the initial time t_init. This means
#define pert_lan 0.                  //that if you set pert_tra to 0., then the periapsis happens at time t = 0.
#define pert_mass 0.//332946.04346   //Mass of the perturbator (in simulation's units).
#define pert_radius 109.1979         //Radius of the perturbator (for visualization purposes only. Does not matter if openGL_bool is 0)



/*******************************************************/
/******** Parameters relative to the simulation ********/
/*******************************************************/

/******** General parameters ********/
#define N_max 5000                   //Maximum number of bodies that the simulation can handle. The simulation will stop if the number of bodies ever exceeds N_max.
#define N_0 1000                     //Initial number of bodies, central body excluded (if any). Must be less than N_max. If random_initial_bool is 0, number of lines of init.txt
#define t_init 0.                    //Time at the beginning of the simulation (in simulation's units)
#define t_end 512.                   //Time at the end       of the simulation (in simulation's units). The actual final time will be larger if (t_end - t_init)/time_step is not integer
#define time_step 0.015625           //Timestep of the simulation (in simulation's units)
#define output_step 32               //Output occurs every output_step timestep. Unimportant if write_to_files_bool is 0
#define high_dumping_threshold 230.95//Threshold (in simulation's units) beyond which bodies are dumped from the simulation (assumed unbounded)

/******** Specific parameters ********/
#define max_ids_per_node 173         //The maximum number of ids in each node of the unrolled linked lists (chains). Choose such that sizeof(struct chain) be a multiple of the cache line
#define softening_parameter 0.       //The softening parameter for mutual gravitational interations, in units of the sum of the radii.
#define seed 129425373               //The seed used for random number generation. Unimportant if seed_bool is 0.
#define switch_to_brute_force 110    //Threshold for N below which NcorpiON switches to the brute-force method for mutual interactions. Unimportant if brute_force_bool is already 1

/******** Bounds for initial conditions. Unimportant if random_initial_bool is 0. The initial conditions are drawn uniformly at random between these bounds ********/
#define radius_min 0.004             //Minimal radius                   of a body at initial time
#define radius_max 0.05              //Maximal radius                   of a body at initial time
#define density_min 0.1048           //Minimal density                  of a body at initial time
#define density_max 0.1848           //Maximal density                  of a body at initial time
#define eccentricity_min 0.          //Minimal eccentricity             of a body at initial time
#define eccentricity_max 0.2         //Maximal eccentricity             of a body at initial time
#define sma_min 2.9                  //Minimal semi-major axis          of a body at initial time
#define sma_max 21.                  //Maximal semi-major axis          of a body at initial time
#define inclination_min 0.           //Minimal inclination (in radians) of a body at initial time
#define inclination_max 0.174533     //Maximal inclination (in radians) of a body at initial time
                                     //The true longitude, argument of pericenter and longitude of the ascending node are drawn uniformly at random between 0 and 2*M_PI
                                     //These bounds must be defined in the simulation's units

/******** Parameters relative to 3D visualization with REBOUND. Unimportant if openGL_bool is 0 ********/
#define browser_port 1234            //The http port where your browser will communicate with REBOUND. You can visualize several simulations at the same time if you change the port
#define radius_blow_up_factor 4.0    //All the bodies, except the central mass (if any), are displayed with a radius that much larger than their true radius. Can enhance visualization



/**********************************************************************************************************************************************************/
/******** Parameters relative to viscoelasticity. Unimportant if viscoelastic_bool is 0. NcorpiON allows for a visco-elastic body to be simulated. ********/
/******** The body is discretized by N_0 nodes connected by Kelvin-Voigt models (a spring and a damper in parallel). If the user sets              ********/
/******** random_initial_bool to 1, then the nodes of the discretization are drawn at random inside a shape model (if file pth/shape_model.txt is  ********/
/******** provided) or inside a sphere of radius R_unit (if no such file is provided). Once the nodes have been generated, the initial distance    ********/
/******** between them is the rest_length of the springs. At the beginning of the simulation, the body will slightly collapse due to gravity but   ********/
/******** will then reach equilibrium due to the damping and the spring's reaction (as long as pert_mass is set to 0.0, as to remove tides).       ********/
/******** The idea is to run a first simulation with random_initial_bool set to 1 to allow the body to rest. If resume_simulation_bool is also set ********/
/******** to 1, then files init.txt and connections.txt are generated and allow for a second simulation to start from a body at rest (by setting   ********/
/******** random_initial_bool to 0 this time). The perturbing mass can then be set to a non-zero value to generate tidal forces in the second      ********/
/******** simulation. The file shape_model.txt has 3 columns (x,y,z) giving the coordinates of vertices on the surface of the viscoelastic body.   ********/
/******** NcorpiON assumes that the shape model is in the principal axis frame (X,Y,Z) with Z towards the shortest axis, but this is not required  ********/
/**********************************************************************************************************************************************************/

#define spring_modulus 500.0         //The expected bulk modulus of the body. Spring stiffness is k = spring_modulus*L. Force is -spring_modulus*L*dL with L the rest_length
#define damping_coefficient 0.125    //The damping coefficient of the dampers in the Kelvin-Voigt models. Force is -damping_coefficient*L*dL/dt with L the rest_length
#define minimal_distance 0.6         //Minimal initial distance between two particles in units of (V/N_0)^(1/3) where V is the volume of the body. 0.3 < minimal_distance < 0.7 is best
#define connections_per_node 35.     //Expected value of the number of connections per node. Values larger than 12.0 are advised for structural integrity. Some nodes will be connected
                                     //less than that as this is just the expected value. NcorpiON makes sure that no node is connected less than three times to prevent wandering.
#define nodes_radius 0.25            //Nodes' radii in units of the minimal initial distance. Can be used to check the structural integrity. In the resting simulation, if collision_bool
                                     //is set to 1 and this parameter is set to a small value (e.g. 0.1), then no collisions should occur if spring_modulus is large enough

/******** Rotation of the viscoelastic body, in the fixed body frame (same reference frame as the shape model). Set all three to 0.0 for no rotation  ********/
#define OmegaX 0.152752              //X-component of the rotation vector. The rotation vector is Omega = (OmegaX, OmegaY, OmegaZ)
#define OmegaY 0.080406              //Y-component of the rotation vector. To be given in radians per unit of time of the simulation
#define OmegaZ 0.419375              //Z-component of the rotation vector. To be given in radians per unit of time of the simulation

/******** These two parameters define the direction of the angular momentum of the viscoelastic body in the inertial reference frame (the perturbator's orbit frame) ********/
#define lbd_long 4.31096             //The corresponding unit vector is (X,Y,Z) = (cos lbd_long*cos beta_lat, sin lbd_long*cos beta_lat, sin beta_lat). After generating the rotation
#define beta_lat -1.029744           //in the fixed-body frame, NcorpiON computes the direction of the angular momentum and rotates the whole body in order to make it match (X,Y,Z)
                                     //To be given in radians. See https://ssp.imcce.fr/forms/ssocard for these two values in the case of small solar system bodies.



/*************************************************************************************************************************/
/******** Parameters relative to tree-based algorithms (falcON and the standard tree code). If either falcON_bool ********/
/******** or standard_tree_bool is 1, then Dehnen's falcON algorithm or Barnes & Hut's standard tree code is used ********/
/******** to compute mutual gravity and to detect collisions. In theses cases, these parameters must be defined   ********/
/*************************************************************************************************************************/

#define expansion_order 3            //The order p of the multipole expansions. NcorpiON allows up to p = 8. Minimum is 1 as order 0 yields no acceleration
#define theta_min 0.45               //Minimal value of the opening angle theta. Must be strictly less than 1. Advised values are 0.25 < theta_min < 0.75
                                     //Larger expansion orders p or smaller opening angles theta_min yield a better precision on the gravity computation.
                                     //The precision (and computational time) of the mutual gravity computed by Ncorpion increases with increasing p and decreasing theta_min.
#define subdivision_threshold 17     //A cubic cell is not divided as long as it contains at most that many bodies. Called s in Dehnen (2002). Must be > 0. The precision does not depend
                                     //on this threshold, but the computational time does. Suggested values are s = (10, 10, 15, 30, 50, 50, 75, 110) for p = (1, 2, 3, 4, 5, 6, 7, 8)
#define root_sidelength 462.0        //Sidelength of the root cell (in simulation's units). Should be machine representable. Particles outside of the root cell don't feel others.
#define level_max 25                 //The maximum allowed number of levels in the tree. Root is at level 0 and a cell at level level_max - 1 is never divided.

/******** Parameters specifically relative to mutual gravity computation with falcON algorithm ********/
#define N_cc_pre 8                   //For two cells with N1 and N2 bodies, if N1N2 < N_cc_pre, then the interaction is computed brute-forcely, regardless of their well-separation
#define N_cc_post 64                 //For two cells with N1 and N2 bodies, if N1N2 < N_cc_post and they are not well-separated, then the interaction is computed brute-forcely
#define N_cs 64                      //If N1 < N_cs, then the self-interaction of a cell containing N1 bodies is computed brute-forcely
                                     //See TreeWalk algorithm in the paper for details about these thresholds. The precision does not depend on these thresholds, but the computational
                                     //time slightly does. However, note that the computational cost depends much more on the subdivision threshold s that on these thresholds. Possible
                                     //values are (N_cc_pre, N_cc_post, N_cs) = (8, 64, 64) for p <= 4 and (N_cc_pre, N_cc_post, N_cs) = (256, 1024, 128) for p >= 5 but this depends
                                     //on the architecture and should be tweaked. N_cc_post must be larger that N_cc_pre.

/******** Parameters specifically relative to mutual gravity computation with the standard tree code ********/
#define N_cb_pre 6                   //If N1 < N_cb_pre, the interaction between a body and a cell with N1 bodies is performed brute-forcely, regardless of their well-separation
#define N_cb_post 16                 //If N1 < N_cb_post and they are not well-separated, the interaction between a body and a cell with N1 bodies is performed brute-forcely.

/******** Parameters specifically relative to collision detection with falcON algorithm ********/
#define N_cc_collision 16            //This is N_cc_post. For two nodes not well-separated with Na and Nb bodies, if NaNb < N_cc_collision, then collisions are searched brute-forcely,
                                     //otherwise, the largest node is subdivised. N_cc_pre is always 0 when falcON is used for collision detection and is not a parameter
#define N_cs_collision 12            //For a node with Na bodies, if Na < N_cs_collision, then collisions are searched brute-forcely, otherwise, the node is subdivised

/******** Parameters specifically relative to collision detection with the standard tree code ********/
#define N_cb_collision 16            //This is N_cb_post. For a body not well-separated from a node with Na bodies, if Na < N_cb_collision, then collisions are searched brute-forcely,
                                     //otherwise, the node is subdivised. N_cb_pre is always 0 when the standard tree code is used for collision detection and is not a parameter



/*****************************************************************************************************************/
/******** Parameters relative to the mesh-grid algorithm. If mesh_bool is 1 then a mesh algorithm is used ********/
/******** to compute mutual gravity and detect collisions. In that case, these parameters must be defined ********/
/*****************************************************************************************************************/

#define collision_cube_min 80.0      //Minimal sidelength of the mesh-grid (in simulation's units). The mesh-size will never be less than collision_cube_min/collision_cube_cells
#define collision_cube_cells 1024    //Number of mesh cells per dimension of the collision cube. Must be even. For 16+ GiB of RAM (resp. 8 or 4 GiB), choose ~1000 (resp. ~800 or ~500)
#define how_many_neighbours 16.0     //Expected value of the number of neighbours for a body



/*****************************************************************************************/
/******** Parameters relative to collision resolution and the fragmentation model ********/
/*****************************************************************************************/

/******** Collision resolution when all collisions are inelastic. Important only if both collision_bool and inelastic_collision_bool are 1 ********/
#define collision_parameter 2.       //The collision parameter f to model inelastic collision. Must be between 1 (completely inelastic) and 2 (elastic)

/******** Collision resolution with the fragmentation model described in NcorpiON's paper. Important only if both collision_bool and fragmentation_bool are 1 ********/
#define N_tilde 15                   //Ratio between ejected mass and mass of the second largest fragment : N° of fragments in the tail. 2*beta/(3-beta) in Leinhardt & Stewart(2012)
#define mu_parameter 0.55            //The exponent of the impact velocity in the coupling parameter. See Table 3 of Housen & Holsapple (2011)
#define nu_parameter 0.4             //The exponent of the density         in the coupling parameter. See Table 3 of Housen & Holsapple (2011)
#define C1_parameter 1.5             //A dimensionless parameter of impact theories. See Table 3 of Housen & Holsapple (2011)
#define k_parameter 0.2              //A dimensionless parameter of impact theories. See Table 3 of Housen & Holsapple (2011)
#define merging_threshold 0.0002     //Threshold on the total ejected mass. If (total ejected mass)/(m1 + m2) is less than that, then the collision results in a merger.
#define fragment_threshold 1.5e-8    //Threshold on the mass of the fragments. If the fragments of the ejecta tail have a mass smaller than that and if the tail is less massive than
                                     //the largest fragment and if the collision is not a merger, then the tail is reunited into a single body instead of being fragmented into
                                     //N_tilde fragments. These two threshold prevent the number of bodies to grow uncontrollably

#define pq_min_max {-1,3,-1,1}       //Extremal integer values for p_k and q_k to determine the position of the tail fragments with respect to the largest fragment.
                                     //Must define a rectangle containing exactly N_tilde points with integer coordinates. More precisely, if pq_min_max = {a, b, c, d},
                                     //then we must have N_tilde = (b - a + 1)*(d - c + 1). See NcorpiON's paper for details



/*********************************************************/
/******** Redefining keywords. Not to be modified ********/
/*********************************************************/

#define typ double                   //Renaming double as typ
#define MPI_TYP MPI_DOUBLE           //Renaming double as typ for mpi
#define bigint long int              //Renaming long int as bigint. If sizeof(int) == sizeof(long int) == 4 on your system, then define bigint as long long int instead
                                     //Open the file /usr/share/gtksourceview/language-specs/c.lang and then find the field
                                     //<context id="types" style-ref="type">. Add the lines <keyword>typ</keyword> and <keyword>bigint</keyword> to it
                                     //and color gedit with C. The keywords typ and bigint should now be colored after a reboot
#define type_check __builtin_types_compatible_p



#endif
