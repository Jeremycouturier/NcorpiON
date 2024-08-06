/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    ffm.h                                                       ********/
/******** @brief   Header file to ffm.c                                        ********/
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


#ifndef _FFM_H_
#define _FFM_H_

#include "parameters.h"
#include "structure.h"


extern typ Mtot;                            //Total body mass inside the root cell. Mass of the root cell
extern int how_many_cells;                  //Total number of cells in the boxdot tree
extern int cell_id;                         //The current unique id of a cell
extern int * already_in_tree;               //Takes care of bodies no yet put in the tree in Hilbert order
extern typ * C1Moonlets;                    //The array of acceleration of the bodies due to their mutual gravity. This is the output of the ffm algorithm


/******** A tree structure employed prior to the three stages of Dehnen's algorithm ********/
struct boxdot {
      struct boxdot * oct[8];    //The eight octants, or children of that cell.
      struct chain * dots;       //A chain containing the ids of the bodies in that cell.
      typ corner[3];             //Top-left-front corner of the node.
      typ sidelength;            //Sidelength of the node (half that of its parent).
      int id;                    //The unique id of that cell, determined by its Hilbert-Peano order : index in the flattree
      int how_many;              //The number of bodies in that cell.
      int rotation;              //The rotation id such that digit = DigitFromOctant[rotParent][oct] and rotation = RotationFromOctant[rotParent][oct]
      int level;                 //The level of that node in the boxdot tree
};

/******** The three phases of Dehnen's algorithm are performed on an array of nodes, defined as follow ********/
struct node {
      #if expansion_order >= 8 && mutual_bool
      typ C8[45];            //Eighth  order tensor field of interactions
      typ M7[36];            //Seventh order multipole moment
      #endif
      #if expansion_order >= 7 && mutual_bool
      typ C7[36];            //Seventh order tensor field of interactions
      typ M6[28];            //Sixth   order multipole moment
      #endif
      #if expansion_order >= 6 && mutual_bool
      typ C6[28];            //Sixth order tensor field of interactions
      typ M5[21];            //Fifth order multipole moment
      #endif
      #if expansion_order >= 5 && mutual_bool
      typ C5[21];            //Fifth  order tensor field of interactions
      typ M4[15];            //Fourth order multipole moment
      #endif
      #if expansion_order >= 4 && mutual_bool
      typ C4[15];            //Fourth order tensor field of interactions
      typ M3[10];            //Octupole
      #endif
      #if expansion_order >= 3 && mutual_bool
      typ C3[10];            //Third order tensor field of interactions
      typ M2[6];             //Quadrupole
      #endif
      #if expansion_order >= 2 && mutual_bool
      typ C2[6];             //Second order tensor field of interactions (tidal tensor)
      #endif
      typ com[3];            //When treating self-gravity : Center of mass and expansion center. When treating collisions : Average position of the bodies in the node
      typ C1[3];             //First order tensor field of interactions (acceleration)
      typ center[3];         //Center of the node
      typ sidelength;        //Sidelength of the node
      typ M0;                //When treating self-gravity : Mass of the cell or zeroth-order multipole moment. When treating collisions : max(R_i + v_i*timestep)
      typ r_max;             //Convergence radius of the Taylor expansion
      typ r_crit;            //When treating self-gravity : r_max/theta. When treating collisions : r_max + M0
      int * dots;            //Array of ids of bodies contained in that node
      int idParent;          //Index of the parent node
      int idFirstChild;      //Index of the first child (Child with smallest index)
      int how_many_children; //The number of children that node has
      int how_many_dots;     //The number of bodies that this node contains
};


typ fast_pow(typ x, int power);


void subdivision(struct boxdot * BoxDot, struct moonlet * moonlets);


void add_boxdot(struct boxdot * BoxDot, struct moonlet * moonlets, int a);


struct boxdot * root_cell(struct moonlet * moonlets);


void fill_boxdot_int(struct boxdot * BoxDot, struct boxdot * Parent, int octantChild);


struct node * flattree_init(struct boxdot * BoxDot);


void tensor_initialization(struct node * FlatTree);


void get_s1_s2_s3(int k, int n, int * s1, int * s2, int * s3);


void get_s2_s3(int k, int * s2, int * s3);


void s1s2s3_from_kn_init();


void s2s3_from_k_init();


void get_k(int s2, int s3, int * k);


void k_from_s2s3_init();


void indexes_from_kn(int k, int n, int * ijklmn);


void ijklmn_from_kn_init();


void k_from_indexes(int * k, int * ijklmn);


void k_from_ijklmn_init();


void permutation_from_kn(int k, int n, int * perm);


void perm_from_kn_init();


void q1_from_q2q3(int * q1, int q2, int q3);


void q1fromq2q3_init();


void get_com(struct node * FlatTree, struct moonlet * moonlets, int a);


void get_rmax(struct node * FlatTree, struct moonlet * moonlets, int a);


void get_tolerance_parameter(struct node * FlatTree, int a, typ precision);


void get_Xn(int k, int n, typ * X, typ m, typ * M);


void get_Xn_overwrite(int k, int n, typ * X, typ m, typ * M);


#if expansion_order >= 3
void get_Mn(struct node * FlatTree, struct moonlet * moonlets, int a);
#endif


void get_com_from_children(struct node * FlatTree, int a);


void get_rmax_from_children(struct node * FlatTree, int a);


#if expansion_order >= 3
void get_Mn_from_children(struct node * FlatTree, int a);
#endif


void com_flattree(struct node * FlatTree, struct moonlet * moonlets);


void rmax_flattree(struct node * FlatTree, struct moonlet * moonlets);


void rcrit_flattree(struct node * FlatTree);


#if expansion_order >= 3
void multipole_flattree(struct node * FlatTree, struct moonlet * moonlets);
#endif


void gradR(typ * R, typ * grad);


void inner_product(typ * T1, typ * T2, typ * T3, int p, int q, typ factor);


void double_IP(typ * T1, typ * T2, typ * TT2, typ * T3, typ * TT3, int p, int q, typ factor1, typ factor2);


#if mutual_bool
void Cm_flattree(struct node * FlatTree, struct moonlet * moonlets);


void Cm_downtree(struct node * FlatTree, struct moonlet * moonlets);


void standard_tree_acceleration(struct node * FlatTree, struct moonlet * moonlets, int a);
#endif


void create_boxdot(struct boxdot ** BoxDot, typ * corner_coordinates, typ D);


void clear_boxdot(struct boxdot ** BoxDot);


int get_octant(typ x, typ y, typ z, typ xa, typ ya, typ za, typ D);


void get_corner_coordinates(typ X, typ Y, typ Z, typ D, int i, typ * corner);


void print_boxdot(struct boxdot * BoxDot);


/******** Material for Hilbert-Peano order ********/
extern int RotationFromOctant[48][8];
extern int DigitFromOctant[48][8];
extern int OctantFromDigit[48][8];
extern int * PeanoHilbertOrder;
extern int IndexPeanoHilbertOrder;

/******** Material to take advantage of the symmetry of the manipulated tensors ********/
extern int k_from_s2s3[expansion_order + 1][expansion_order + 1];
extern int s1s2s3_from_kn[((expansion_order + 1)*(expansion_order + 2))/2][expansion_order + 1][3];
extern int s2s3_from_k[((expansion_order + 1)*(expansion_order + 2))/2][2];
extern int ijklmn_from_kn[(expansion_order*(expansion_order + 1))/2][expansion_order][expansion_order - 1];
extern int perm_from_kn[((expansion_order + 1)*(expansion_order + 2))/2][expansion_order + 1];
extern int q1fromq2q3[((expansion_order + 1)*(expansion_order + 2))/2][((expansion_order + 1)*(expansion_order + 2))/2];
extern int k_from_ijklmn[4][4][4][4][4][4];


#endif
