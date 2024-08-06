/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    ffm.c                                                       ********/
/******** @brief   This file manages gravitation with falcON and Barnes & Hut  ********/
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
#include "ffm.h"
#include <errno.h>
#include <math.h>

#include "parameters.h"

/*******************************************************************************************************************/
/******** In this file, I implement the fast multipole algorithm falcON for O(N) mutual gravity evaluation ********/
/*******************************************************************************************************************/


typ Mtot;
int how_many_cells;
int cell_id;
int * PeanoHilbertOrder;
int IndexPeanoHilbertOrder;
typ * C1Moonlets;

int k_from_s2s3[expansion_order + 1][expansion_order + 1];
int s1s2s3_from_kn[((expansion_order + 1)*(expansion_order + 2))/2][expansion_order + 1][3];
int s2s3_from_k[((expansion_order + 1)*(expansion_order + 2))/2][2];
int ijklmn_from_kn[(expansion_order*(expansion_order + 1))/2][expansion_order][expansion_order - 1];
int perm_from_kn[((expansion_order + 1)*(expansion_order + 2))/2][expansion_order + 1];
int q1fromq2q3[((expansion_order + 1)*(expansion_order + 2))/2][((expansion_order + 1)*(expansion_order + 2))/2];
int k_from_ijklmn[4][4][4][4][4][4];


/******** Arrays relative to Hilbert-Peano order      ********/
/******** These arrays were provided by Walter Dehnen ********/
int RotationFromOctant[48][8] = {
      {36,10,25,25,28,10,27,27},
      {29,37,24,26,11,11,24,26},
      { 8,30,25,25, 8,38,27,27},
      { 9, 9,24,26,39,31,24,26},
      {40,40,44,44,24, 6,32, 6},
      {25,41,33,45, 7,41, 7,45},
      { 4,26, 4,34,42,42,46,46},
      {43, 5,47, 5,43,27,47,35},
      {33,33,36, 2,35,35,28, 2},
      {32,34,29,37,32,34, 3, 3},
      {33,33, 0,30,35,35, 0,38},
      {32,34, 1, 1,32,34,39,31},
      {24,14,32,14,42,42,46,46},
      {43,25,47,33,43,15,47,15},
      {40,40,44,44,12,26,12,34},
      {13,41,13,45,27,41,35,45},
      {28,38,28,38,41,43,22,22},
      {42,29,23,29,40,39,23,39},
      {41,43,20,20,36,30,36,30},
      {37,42,37,21,31,40,31,21},
      {28,38,28,38,18,18,45,47},
      {19,29,46,29,19,39,44,39},
      {16,16,45,47,36,30,36,30},
      {37,17,37,46,31,17,31,44},
      {12,34, 1, 1, 4,34, 3, 3},
      { 5,13, 0, 2,35,35, 0, 2},
      {32, 6, 1, 1,32,14, 3, 3},
      {33,33, 0, 2,15, 7, 0, 2},
      {16,16,20,20, 0,30, 8,30},
      { 1,17, 9,21,31,17,31,21},
      {28, 2,28,10,18,18,22,22},
      {19,29,23,29,19, 3,23,11},
      { 9, 9,12,26,11,11, 4,26},
      { 8,10, 5,13, 8,10,27,27},
      { 9, 9,24, 6,11,11,24,14},
      { 8,10,25,25, 8,10,15, 7},
      { 0,38, 8,38,18,18,22,22},
      {19, 1,23, 9,19,39,23,39},
      {16,16,20,20,36, 2,36,10},
      {37,17,37,21, 3,17,11,21},
      { 4,14, 4,14,17,19,46,46},
      {18, 5,47, 5,16,15,47,15},
      {17,19,44,44,12, 6,12, 6},
      {13,18,13,45, 7,16, 7,45},
      { 4,14, 4,14,42,42,21,23},
      {43, 5,22, 5,43,15,20,15},
      {40,40,21,23,12, 6,12, 6},
      {13,41,13,22, 7,41, 7,20}
};


int OctantFromDigit[48][8] = {
      { 0, 2, 3, 1, 5, 7, 6, 4},
      { 1, 3, 7, 5, 4, 6, 2, 0},
      { 5, 7, 6, 4, 0, 2, 3, 1},
      { 4, 6, 2, 0, 1, 3, 7, 5},
      { 4, 0, 1, 5, 7, 3, 2, 6},
      { 0, 1, 5, 4, 6, 7, 3, 2},
      { 1, 5, 4, 0, 2, 6, 7, 3},
      { 5, 4, 0, 1, 3, 2, 6, 7},
      { 6, 4, 5, 7, 3, 1, 0, 2},
      { 2, 0, 4, 6, 7, 5, 1, 3},
      { 3, 1, 0, 2, 6, 4, 5, 7},
      { 7, 5, 1, 3, 2, 0, 4, 6},
      { 2, 6, 7, 3, 1, 5, 4, 0},
      { 3, 2, 6, 7, 5, 4, 0, 1},
      { 7, 3, 2, 6, 4, 0, 1, 5},
      { 6, 7, 3, 2, 0, 1, 5, 4},
      { 5, 1, 3, 7, 6, 2, 0, 4},
      { 4, 5, 7, 6, 2, 3, 1, 0},
      { 0, 4, 6, 2, 3, 7, 5, 1},
      { 1, 0, 2, 3, 7, 6, 4, 5},
      { 6, 2, 0, 4, 5, 1, 3, 7},
      { 2, 3, 1, 0, 4, 5, 7, 6},
      { 3, 7, 5, 1, 0, 4, 6, 2},
      { 7, 6, 4, 5, 1, 0, 2, 3},
      { 4, 6, 7, 5, 1, 3, 2, 0},
      { 0, 2, 6, 4, 5, 7, 3, 1},
      { 1, 3, 2, 0, 4, 6, 7, 5},
      { 5, 7, 3, 1, 0, 2, 6, 4},
      { 6, 2, 3, 7, 5, 1, 0, 4},
      { 2, 3, 7, 6, 4, 5, 1, 0},
      { 3, 7, 6, 2, 0, 4, 5, 1},
      { 7, 6, 2, 3, 1, 0, 4, 5},
      { 2, 0, 1, 3, 7, 5, 4, 6},
      { 3, 1, 5, 7, 6, 4, 0, 2},
      { 7, 5, 4, 6, 2, 0, 1, 3},
      { 6, 4, 0, 2, 3, 1, 5, 7},
      { 0, 4, 5, 1, 3, 7, 6, 2},
      { 1, 0, 4, 5, 7, 6, 2, 3},
      { 5, 1, 0, 4, 6, 2, 3, 7},
      { 4, 5, 1, 0, 2, 3, 7, 6},
      { 4, 0, 2, 6, 7, 3, 1, 5},
      { 0, 1, 3, 2, 6, 7, 5, 4},
      { 1, 5, 7, 3, 2, 6, 4, 0},
      { 5, 4, 6, 7, 3, 2, 0, 1},
      { 7, 3, 1, 5, 4, 0, 2, 6},
      { 6, 7, 5, 4, 0, 1, 3, 2},
      { 2, 6, 4, 0, 1, 5, 7, 3},
      { 3, 2, 0, 1, 5, 4, 6, 7}
};

/************************************************************************************************************/
/******** First, I implement functions allowing to put bodies in a tree-like structure called boxdot ********/
/************************************************************************************************************/

typ fast_pow(typ x, int power){

      /******** Raises x to the strictly positive integer power power ********/
      
      int current_power = power;
      typ to_be_returned = x;
      typ to_be_multiplied_by = 1.0;
      while (current_power >= 2){
            if (current_power % 2 == 0){
                  current_power /= 2;
            }
            else {
                  current_power = (current_power - 1)/2;
                  to_be_multiplied_by *= to_be_returned;
            }
            to_be_returned *= to_be_returned;
      }
      return to_be_multiplied_by*to_be_returned;
}


void subdivision(struct boxdot * BoxDot, struct moonlet * moonlets){

      /******** Auxiliary function to add_boxdot.          ********/
      /******** Adds the bodies of BoxDot to its children. ********/

      int a; //The current body
      typ x, y, z, xa, ya, za, D; //Coordinates
      int octant; //The octant where the body must go
      struct chain * dots = BoxDot -> dots; //Chain of ids of bodies inside BoxDot
      int index = dots -> how_many - 1; //The index of a in BoxDot -> dots
      
      /******** Getting the box coordinates ********/
      x = (BoxDot -> corner)[0];
      y = (BoxDot -> corner)[1];
      z = (BoxDot -> corner)[2];
      D =  BoxDot -> sidelength/2.0;

      while (dots != NULL){
            a = (dots -> ids)[index]; //The body's id
            
            /******** Getting the body's coordinates ********/
            xa = (moonlets + a) -> x;
            ya = (moonlets + a) -> y;
            za = (moonlets + a) -> z;
            
            octant = get_octant(x, y, z, xa, ya, za, D); //Retrieving the octant where the body belongs
            
            /******** Initializing the child if necessary ********/
            if ((BoxDot -> oct)[octant] == NULL){
                  typ corner[3];
                  get_corner_coordinates(x, y, z, D, octant, corner);
                  create_boxdot(&((BoxDot -> oct)[octant]), corner, D);
                  fill_boxdot_int((BoxDot -> oct)[octant], BoxDot, octant);
                  how_many_cells ++;
            }
            
            add_boxdot((BoxDot -> oct)[octant], moonlets, a); //The while loop of add_boxdot will never be reached.
                                                              //subdivision will be called again by add_boxdot only if all the bodies of BoxDot go into the same octant
            /******** Switching to the next body ********/                                                  
            index --;
            if (index < 0){
                  dots = dots -> queue;
                  index = max_ids_per_node - 1;
            }
      }
}


void add_boxdot(struct boxdot * BoxDot, struct moonlet * moonlets, int a){

      /******** Adds body a to the boxdot BoxDot, and to its descendants.    ********/
      /******** It is assumed that it has been checked, prior to calling     ********/
      /******** this function, that the body a belongs to the boxdot BoxDot. ********/
      /******** It is also assumed that BoxDot was previously initialized.   ********/
      
      int n = BoxDot -> how_many; //Number of bodies currently in BoxDot
      BoxDot -> dots = Add(a, BoxDot -> dots);
      
      if (n < subdivision_threshold || BoxDot -> level == level_max - 1){ //That box has no children and don't need any or it is too deep to be divided
            BoxDot -> how_many ++;
            return;
      }
      
      if (n == subdivision_threshold){ //That box has no children but needs some.
            BoxDot -> how_many ++;
            subdivision(BoxDot, moonlets);
            return;
      }
      BoxDot -> how_many ++;
      
      /******** If that point is reached, BoxDot has descendants and is divisible. ********/
      
      /******** Retrieving the box coordinates ********/
      typ x, y, z, D; //xyz are the coordinates of the corner of octant 0 and D is the sidelength
      D =  BoxDot -> sidelength;
      x = (BoxDot -> corner)[0];
      y = (BoxDot -> corner)[1];
      z = (BoxDot -> corner)[2];
      
      /******** Retrieving the body's coordinates ********/
      typ xa, ya, za;
      xa = (moonlets + a) -> x;
      ya = (moonlets + a) -> y;
      za = (moonlets + a) -> z;
      
      int octant;
      
      while (1){ //While the box has children, I keep getting deeper into the tree. The while loop will be escaped thanks to the return instructions.
            D /= 2.0;
            octant = get_octant(x, y, z, xa, ya, za, D); //Retrieving the corresponding octant
            
            if ((BoxDot -> oct)[octant] == NULL){ //If the child is still NULL, I initialize it 
                  typ corner[3];
                  get_corner_coordinates(x, y, z, D, octant, corner);
                  create_boxdot(&((BoxDot -> oct)[octant]), corner, D);
                  fill_boxdot_int((BoxDot -> oct)[octant], BoxDot, octant);
                  how_many_cells ++;
            }
            
            BoxDot = (BoxDot -> oct)[octant]; //The new considered box is now the child
            n = BoxDot -> how_many;           //Retrieving the number of bodies it contains
            BoxDot -> dots = Add(a, BoxDot -> dots);
            
            if (n < subdivision_threshold || BoxDot -> level == level_max - 1){ //That box has no children and don't need any or it is too deep to be divided
                  BoxDot -> how_many ++;
                  return;
            }
            if (n == subdivision_threshold){ //That box has no children yet but needs some.
                  BoxDot -> how_many ++;
                  subdivision(BoxDot, moonlets);
                  return;
            }
            
            BoxDot -> how_many ++;
            /******** Actualizing the box coordinates ********/
            x = (BoxDot -> corner)[0];
            y = (BoxDot -> corner)[1];
            z = (BoxDot -> corner)[2];
      }
}


struct boxdot * root_cell(struct moonlet * moonlets){

      /******** Creates and initializes the root cell ********/
      /******** Adds all bodies to it                 ********/
      
      int i, a;
      typ x, y, z; //Body's coordinates
      for (i = 0; i <= largest_id; i ++){
            already_in_tree[i] = 0;
      }
      
      /******** Creating and initializing the root cell ********/
      struct boxdot * root = NULL;
      typ corner[3] = {-root_sidelength/2.0, -root_sidelength/2.0, root_sidelength/2.0};
      create_boxdot(&root, corner, root_sidelength);
      root -> rotation = 0;
      root -> level    = 0; //The root cell is at level 0
      how_many_cells ++;
      
      /******** Adding all moonlets to it ********/
      for (i = 0; i < IndexPeanoHilbertOrder; i ++){ //Adding the bodies in Hilbert order
            a = PeanoHilbertOrder[i];  
            if (*(exists + a)){
                  x = (moonlets + a) -> x;
                  y = (moonlets + a) -> y;
                  z = (moonlets + a) -> z;
                  if (fabs(x) <= root_sidelength/2.0 && fabs(y) <= root_sidelength/2.0 && fabs(z) <= root_sidelength/2.0){ //If body a belongs to the root cell
                        add_boxdot(root, moonlets, a);
                  }
                  already_in_tree[a] = 1;
            }
      }
      for (i = 0; i <= largest_id; i ++){ //Adding the few remaining bodies in random order
            if (exists[i] && !already_in_tree[i]){
                  x = (moonlets + i) -> x;
                  y = (moonlets + i) -> y;
                  z = (moonlets + i) -> z;
                  if (fabs(x) <= root_sidelength/2.0 && fabs(y) <= root_sidelength/2.0 && fabs(z) <= root_sidelength/2.0){ //If body i belongs to the root cell
                        add_boxdot(root, moonlets, i);
                  }
            }
      }
      
      return root;
}


/***********************************************************************************************************************/
/******** I now link the boxdot tree to a simple array containing the multipole moments and field tensor C^(n)  ********/
/******** Each node of the boxdot tree is attributed a unique id that respects the Hilbert-Peano order and that ********/
/******** corresponds to its index position into the aforementioned array.                                      ********/
/******** The three phases of Dehnen's algorithm are performed on that array instead of on the boxdot tree      ********/
/***********************************************************************************************************************/


void fill_boxdot_int(struct boxdot * BoxDot, struct boxdot * Parent, int octantChild){

      /******** Initializes the integer fields rotation and level of the boxdot BoxDot        ********/
      /******** using those fields of the parent cell Parent. BoxDot is in octant octantChild ********/
      
      int rotationParent = Parent -> rotation;
      int levelParent = Parent -> level;
      
      /******** Computing the child's integer fields ********/
      int rotationChild = RotationFromOctant[rotationParent][octantChild];
      int levelChild = levelParent + 1;
      
      /******** Initializing the integer fields of BoxDot ********/
      BoxDot -> rotation = rotationChild;
      BoxDot -> level    = levelChild;
}


struct node * flattree_init(struct boxdot * BoxDot){

      /******** This function, when called on the root cell, initializes its unique id and that of all its descendants          ********/
      /******** It is assumed that the global variable cell_id is 0 prior to calling this function and that the global variable ********/
      /******** how_many_cells is the total number of cells. This function also returns an array FlatTree on which the three    ********/
      /******** stages of Dehnen's FalcON algorithm are performed, and initializes the fields idFirstChild, how_many_children,  ********/
      /******** idParent, how_many_dots, dots, sidelength and corner of each node of FlatTree                                   ********/
      
      int p;
      int octant;
      int rotation;
      int how_many_child = 0;
      int isFirstChild   = 1;
      typ D;
      int how_many_dots;
      int * dots;
      struct chain * ch;
      int index;
      IndexPeanoHilbertOrder = 0;
      
      /******** Allocation of memory to FlatTree ********/
      struct node * FlatTree = (struct node *)malloc(how_many_cells * sizeof(struct node));
      if (FlatTree == NULL){ //Checking that the memory was properly allocated
            fprintf(stderr, "Error : Cannot allocate %d cells to the array FlatTree. Try increasing subdivision_threshold or decreasing level_max.\n", how_many_cells);
            abort();
      }
      (FlatTree + cell_id) -> idParent = -1; //Arbitrarily setting to -1 the unique id of the root cell's parent
      
      /******** Stack of nodes still to be treated ********/
      struct boxdot ** stack = (struct boxdot **)malloc(how_many_cells * sizeof(struct boxdot *));

      int j = 0; // j is the index of where to store a node to be treated
      
      stack[j] = BoxDot;
      j ++;
      
      struct boxdot * to_be_treated;
      struct boxdot * child;
      
      while (j > cell_id){ //While there are still nodes to be treated
            to_be_treated       = stack[cell_id];
            to_be_treated -> id = cell_id; //Initializing the unique id of that cell
            rotation            = to_be_treated -> rotation;
            D                   = to_be_treated -> sidelength;
            
            for (p = 0; p < 8; p ++){ // p is the Hilbert-Peano digit
                  octant = OctantFromDigit[rotation][p];
                  child  = (to_be_treated -> oct)[octant];
                  if (child != NULL){
                        stack[j] = child;
                        (FlatTree + j) -> idParent = cell_id; //Initializing the field idParent of all the children
                        if (isFirstChild){                    //Initializing the field idFirstChild 
                              (FlatTree + cell_id) -> idFirstChild = j;
                              isFirstChild = 0;
                        }
                        j ++;
                        how_many_child ++;
                  }
            }
            (FlatTree + cell_id) -> how_many_children = how_many_child;
            how_many_dots                             = to_be_treated -> how_many;
            (FlatTree + cell_id) -> how_many_dots     = how_many_dots;
            
            /******** Initializing the field array dots ********/
            dots                         = (int *)malloc(how_many_dots * sizeof(int));
            (FlatTree + cell_id) -> dots = dots;
            if (dots == NULL){
                  fprintf(stderr, "Error : Could not allocate memory for %d dots in function flattree_init.\n", how_many_dots);
                  abort();
            }
            ch    = to_be_treated -> dots;                   
            index = ch -> how_many - 1;
            for (p = 0; p < how_many_dots; p ++){
                  *(dots + p) = (ch -> ids)[index];
                  if (how_many_child == 0){
                        PeanoHilbertOrder[IndexPeanoHilbertOrder] = (ch -> ids)[index];
                        IndexPeanoHilbertOrder ++;
                  }
                  /******** Switching to the next body ********/                                                  
                  index --;
                  if (index < 0){
                        ch = ch -> queue;
                        index = max_ids_per_node - 1;
                  }
            }
            
            /******** Some initialization ********/           
            (FlatTree  + cell_id) -> sidelength = D;
            ((FlatTree + cell_id) -> center)[0] = (to_be_treated -> corner)[0] + D/2.0;
            ((FlatTree + cell_id) -> center)[1] = (to_be_treated -> corner)[1] + D/2.0;
            ((FlatTree + cell_id) -> center)[2] = (to_be_treated -> corner)[2] - D/2.0;
            if (how_many_child == 0){ //If to_be_treated has no children, then I arbitrarily set the unique id of its first child to -1 
                  (FlatTree + cell_id) -> idFirstChild = -1;
            }
            how_many_child = 0;
            isFirstChild   = 1;
            cell_id ++;
      }

      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > how_many_cells){
            fprintf(stderr, "Error : The stack is not big enough in function unique_id. Aborting before segmentation fault.\n");
            abort();
      }
      free(stack);
      stack = NULL;
      
      return FlatTree;
}


void tensor_initialization(struct node * FlatTree){

      int p, k;

      /******** Initializing to zero the multipole and interaction tensors ********/
      for (p = 0; p <= 3*largest_id + 2; p ++){
            *(C1Moonlets + p) = 0.;
      }
      for (p = 0; p < how_many_cells; p ++){
            for (k = 0; k < 3; k ++){
                  ((FlatTree + p) -> C1)[k] = 0.;
            }
            #if expansion_order >= 2 && mutual_bool
                  for (k = 0; k < 6; k ++){
                        ((FlatTree + p) -> C2)[k] = 0.;
                  }
                  #if expansion_order >= 3
                        for (k = 0; k < 10; k ++){
                              ((FlatTree + p) -> C3)[k] = 0.;
                        }
                        for (k = 0; k < 6; k ++){
                              ((FlatTree + p) -> M2)[k] = 0.;
                        }
                        #if expansion_order >= 4
                              for (k = 0; k < 15; k ++){
                                    ((FlatTree + p) -> C4)[k] = 0.;
                              }
                              for (k = 0; k < 10; k ++){
                                    ((FlatTree + p) -> M3)[k] = 0.;
                              }
                              #if expansion_order >= 5
                                    for (k = 0; k < 21; k ++){
                                          ((FlatTree + p) -> C5)[k] = 0.;
                                    }
                                    for (k = 0; k < 15; k ++){
                                          ((FlatTree + p) -> M4)[k] = 0.;
                                    }
                                    #if expansion_order >= 6
                                          for (k = 0; k < 28; k ++){
                                                ((FlatTree + p) -> C6)[k] = 0.;
                                          }
                                          for (k = 0; k < 21; k ++){
                                                ((FlatTree + p) -> M5)[k] = 0.;
                                          }
                                          #if expansion_order >= 7
                                                for (k = 0; k < 36; k ++){
                                                      ((FlatTree + p) -> C7)[k] = 0.;
                                                }
                                                for (k = 0; k < 28; k ++){
                                                      ((FlatTree + p) -> M6)[k] = 0.;
                                                }
                                                #if expansion_order >= 8
                                                      for (k = 0; k < 45; k ++){
                                                            ((FlatTree + p) -> C8)[k] = 0.;
                                                      }
                                                      for (k = 0; k < 36; k ++){
                                                            ((FlatTree + p) -> M7)[k] = 0.;
                                                      }
                                                #endif
                                          #endif
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      }
}


/****************************************************************************/
/******** Some functions allowing symmetrical tensors of any order n ********/
/******** to be stored into a simple array of (n+1)(n+2)/2 elements  ********/
/****************************************************************************/


void get_s1_s2_s3(int k, int n, int * s1, int * s2, int * s3){

      /******** Updates the number of 1, 2 and 3 (denoted by s1, s2 and s3) in the indexes of the k^th distinct component      ********/
      /******** of symmetrical tensor T^(n). For example, if k = 12 and n = 5, then s1 = 1, s2 = 2 and s3 = 2. Therefore, when ********/
      /******** the tensor T^(5) is stored in an array T5[21], then T5[12] stores T_12233 and all the permutations             ********/
      
      int n1 = (int) floor(0.5*(sqrt(1.0 + 8.0*(typ) k) - 1.0));
      int n2 = k - (n1*(n1 + 1))/2;
      *s1 = n - n1;
      *s2 = n1 - n2;
      *s3 = n2;
}


void get_s2_s3(int k, int * s2, int * s3){

      /******** Same as above but depends only on k and computes s2 and s3 only ********/
      
      int n1 = (int) floor(0.5*(sqrt(1.0 + 8.0*(typ) k) - 1.0));
      int n2 = k - (n1*(n1 + 1))/2;
      *s2 = n1 - n2;
      *s3 = n2;
}


void s1s2s3_from_kn_init(){

      /******** In order not to call function get_s1_s2_s3 too many times, I store its return values in a table s1s2s3_from_kn ********/
      
      int k, n, s1, s2, s3;
      for (k = 0; k < ((expansion_order + 1)*(expansion_order + 2))/2; k ++){
            for (n = 0; n <= expansion_order; n ++){
                  get_s1_s2_s3(k, n, &s1, &s2, &s3);
                  s1s2s3_from_kn[k][n][0] = s1;
                  s1s2s3_from_kn[k][n][1] = s2;
                  s1s2s3_from_kn[k][n][2] = s3;
            }
      }
}


void s2s3_from_k_init(){

      /******** Same as above but for the function get_s2_s3 ********/
      
      int k, s2, s3;
      for (k = 0; k < ((expansion_order + 1)*(expansion_order + 2))/2; k ++){
            get_s2_s3(k, &s2, &s3);
            s2s3_from_k[k][0] = s2;
            s2s3_from_k[k][1] = s3;
      }
}


void get_k(int s2, int s3, int * k){

      /******** Inverse of get_s1_s2_s3. ********/
      
      int n1 = s2 + s3;
      int n2 = s3;
      *k = (n1*(n1 + 1))/2 + n2;

}

void k_from_s2s3_init(){

      /******** In order not to call function get_k too many times, I store its return values in a table k_from_s2s3 ********/
      int k;
      int i,j;
      
      for (i = 0; i <= expansion_order; i ++){       // i is the number of 2's
            for (j = 0; j <= expansion_order; j ++){ // j is the number of 3's
                  get_k(i, j, &k);
                  k_from_s2s3[i][j] = k;
            }
      }
}


void indexes_from_kn(int k, int n, int * ijklmn){

      /******** Stores the indexes i,j,k,l,m,n of the k^th component of a symmetrical tensor of order n <= 6 ********/
      /******** If n < 6, then the last indexes are set to 0                                                 ********/

      int s1, s2, s3;
      s1 = s1s2s3_from_kn[k][n][0];
      s2 = s1s2s3_from_kn[k][n][1];
      s3 = s1s2s3_from_kn[k][n][2];
      int p = 0;
      
      while (s1 > 0){
            ijklmn[p] = 1;
            p ++;
            s1 --;
      }
      while (s2 > 0){
            ijklmn[p] = 2;
            p ++;
            s2 --;
      }
      while (s3 > 0){
            ijklmn[p] = 3;
            p ++;
            s3 --;
      }
      while (p < expansion_order - 1){
            ijklmn[p] = 0;
            p ++;
      }
}


void ijklmn_from_kn_init(){

      /******** Stores the return value of indexes_from_kn into ijklmn_from_kn ********/
      
      int ijklmn[expansion_order - 1];
      int k, n, l;
      
      for (k = 0; k < (expansion_order*(expansion_order + 1))/2; k ++){
            for (n = 0; n < expansion_order; n ++){
                  indexes_from_kn(k, n, ijklmn);
                  for (l = 0; l < expansion_order - 1; l ++){
                        ijklmn_from_kn[k][n][l] = ijklmn[l];
                  }
            }
      }
}


void k_from_indexes(int * k, int * ijklmn){

      /******** Inverse of function indexes_from_kn. Retrieves k from the indexes ********/
      
      
      int i;
      int s2 = 0, s3 = 0;
      for (i = 0; i < 6; i ++){
            if (ijklmn[i] == 2){
                  s2 ++;
            }
            else if (ijklmn[i] == 3){
                  s3 ++;
            }
      }
      *k = k_from_s2s3[s2][s3];
}


void k_from_ijklmn_init(){

      /******** Stores the return value of k_from_indexes into k_from_ijklmn ********/
      
      int p;
      int ijklmn[6];
      int i, j, k, l, m , n;
      
      for (i = 0; i < 4; i ++){
            for (j = 0; j < 4; j ++){
                  for (k = 0; k < 4; k ++){
                        for (l = 0; l < 4; l ++){
                              for (m = 0; m < 4; m ++){
                                    for (n = 0; n < 4; n ++){
                                          ijklmn[0] = i; ijklmn[1] = j; ijklmn[2] = k; ijklmn[3] = l; ijklmn[4] = m; ijklmn[5] = n;
                                          k_from_indexes(&p, ijklmn);
                                          k_from_ijklmn[i][j][k][l][m][n] = p;
                                    }
                              }
                        }
                  }
            }
      }
}


void permutation_from_kn(int k, int n, int * perm){

      /******** Computes how many permutations can be formed from the k^th independant ********/
      /******** component of an order n symmetrical tensor                             ********/
      
      int factorial[9] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320};
      
      int s1, s2, s3;
      s1 = s1s2s3_from_kn[k][n][0];
      s2 = s1s2s3_from_kn[k][n][1];
      s3 = s1s2s3_from_kn[k][n][2];
      if (s1 >= 0 && s2 >= 0 && s3 >= 0){
            *perm = factorial[n]/(factorial[s1]*factorial[s2]*factorial[s3]);
      }
}


void perm_from_kn_init(){

      /******** Stores the return value of permutation_from_kn into the array perm_from_kn ********/

      int perm = 0;
      int k, n;
      
      for (k = 0; k < ((expansion_order + 1)*(expansion_order + 2))/2; k ++){
            for (n = 0; n <= expansion_order; n ++){
                  permutation_from_kn(k, n, &perm);
                  perm_from_kn[k][n] = perm;
            }
      }
}

void q1_from_q2q3(int * q1, int q2, int q3){

      /******** Function auxiliary to inner_product. Gives q1 as a function of q2 and q3 ********/

      int n1_2, n1_3, n2_2, n2_3, n1, n2, s2_2, s2_3, s3_2, s3_3, s2, s3;
      n1_2 = (int) floor(0.5*(sqrt(1.0 + 8.0*(typ) q2) - 1.0));
      n1_3 = (int) floor(0.5*(sqrt(1.0 + 8.0*(typ) q3) - 1.0));
      n2_2 = q2 - (n1_2*(n1_2 + 1))/2;
      n2_3 = q3 - (n1_3*(n1_3 + 1))/2;
      s2_2 = n1_2 - n2_2;
      s3_2 = n2_2;
      s2_3 = n1_3 - n2_3;
      s3_3 = n2_3;
      s2   = s2_2 + s2_3;
      s3   = s3_2 + s3_3;
      n1   = s2 + s3;
      n2   = s3;
      *q1  = (n1*(n1 + 1))/2 + n2;
}


void q1fromq2q3_init(){

      /******** Stores the return value of q1_from_q2q3 into the array q1fromq2q3 ********/

      int q1, q2, q3;
      for (q2 = 0; q2 < ((expansion_order + 1)*(expansion_order + 2))/2; q2 ++){
            for (q3 = 0; q3 < ((expansion_order + 1)*(expansion_order + 2))/2; q3 ++){
                  q1_from_q2q3(&q1, q2, q3);
                  q1fromq2q3[q2][q3] = q1;
            }
      }
}


/******************************************************************************************/
/******** I now write functions relative to the first stage of Dehnen's algorithm ********/
/******************************************************************************************/


void get_com(struct node * FlatTree, struct moonlet * moonlets, int a){

      /******** Computes the mass and center of mass of node a of FlatTree, assuming that it has ********/
      /******** no children. Initializes the corresponding fields of (FlatTree + a)              ********/
      
      /******** Checking that the node has indeed no children. To be removed when the code is robust ********/
      if ((FlatTree + a) -> idFirstChild != -1){
            fprintf(stderr, "Error : Node %d has children in function get_com.\n", a);
            abort();
      }
      
      typ M0 = 0.0; //The cell's mass
      typ com[3] = {0.0, 0.0, 0.0}; //The cell's center of mass
      
      typ x, y, z;  //The current body's coordinates
      typ m;        //The current body's mass
      
      int * dots        = (FlatTree + a) -> dots; //All the bodies in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int j; //Id of the current body
      
      for (i = 0; i < how_many_dots; i++){
            j = dots[i];
            m = (moonlets + j) -> mass;
            x = (moonlets + j) -> x;
            y = (moonlets + j) -> y;
            z = (moonlets + j) -> z;
            
            M0     += m;
            com[0] += m*x;
            com[1] += m*y;
            com[2] += m*z;
      }
      
      /******** Initializing the relevant fields ********/
      (FlatTree  + a) -> M0      = M0;
      ((FlatTree + a) -> com)[0] = com[0]/M0;
      ((FlatTree + a) -> com)[1] = com[1]/M0;
      ((FlatTree + a) -> com)[2] = com[2]/M0;
}


void get_rmax(struct node * FlatTree, struct moonlet * moonlets, int a){

      /******** Initializes the convergence radius of node a of FlatTree ********/
      /******** assuming that it has no children.                        ********/

      /******** Checking that the node has indeed no children. To be removed when the code is robust ********/
      if ((FlatTree + a) -> idFirstChild != -1){
            fprintf(stderr, "Error : Node %d has children in function get_rmax.\n", a);
            abort();
      }
      
      typ rmax = 0.0;
      typ distance;
      
      /******** Retrieving the center of mass of the cell ********/
      typ com_x, com_y, com_z;
      com_x = ((FlatTree +a) -> com)[0];  com_y = ((FlatTree +a) -> com)[1];  com_z = ((FlatTree +a) -> com)[2];
      
      typ dx, dy, dz; //Distance with a dot along each axis
      
      int * dots = (FlatTree + a) -> dots; //All the bodies in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int j; //Id of the body whose distance to the center of mass is to be computed
      
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
      (FlatTree + a) -> r_max = rmax;
}


void get_tolerance_parameter(struct node * FlatTree, int a, typ precision){

      /******** Computes the tolerance parameters theta of node a of FlatTree ********/
      /******** Solves Eq. (13) of Dehnen (2002) with a Newton-Raphson method ********/
      /******** Initializes the field r_crit of (FlatTree + a)                ********/
      
      typ M0                = (FlatTree + a) -> M0;
      typ one_theta_min     = 1. - theta_min;
      typ power             = ((typ) expansion_order) + 2.;
      typ rhs               = fast_pow(theta_min, power)/one_theta_min/one_theta_min*pow(M0/Mtot, -1.0/3.0);
      typ current_precision = 1.0;
      typ theta             = 0.7;
      
      if (M0/Mtot < 0.00000000001){
            (FlatTree + a) -> r_crit = (FlatTree + a) -> r_max; //theta = 1
            return;
      }
      else if (M0 < 0.0000001*Mtot){
            theta = 0.95;
      }
      else if (M0 < 0.0001*Mtot){
            theta = 0.8;
      }
      typ fX0, dfX0, one_theta, one_theta_2, one_theta_3, theta_power, dtheta;
      int nstep = 0;
      while(current_precision > precision){
            one_theta         = 1. - theta;
            one_theta_2       = one_theta*one_theta;
            one_theta_3       = one_theta_2*one_theta;
            theta_power       = fast_pow(theta, power);
            fX0               = theta_power/one_theta_2-rhs;
            dfX0              = (2.*theta_power + power*one_theta*theta_power/theta)/one_theta_3;
            dtheta            = -fX0/dfX0;
            theta            += dtheta;
            current_precision = fabs(dtheta/theta);
            nstep ++;
            
            /******** If the method gets lost, I try to put it back on track ********/
            if (theta >= 1.0 || theta <= theta_min || theta != theta){
                  theta = rdm(theta_min, 1.0);
            }
            
            /******** Checking the convergence of the Newton-Raphson method. To be removed when the code is robust ********/
            if (nstep >= 150){
                  fprintf(stderr, "Error : The Newton-Raphson method to find the tolerance parameter does not converge.\n");
                  abort();
            }
      }
      
      /******** Initializing the relevant field ********/
      (FlatTree + a) -> r_crit = (FlatTree + a) -> r_max / theta;
}


void get_Xn(int k, int n, typ * X, typ m, typ * M){

      /******** Accumulates in tensor M the k^th component of tensor m * X^(n) ********/
      
      int s1, s2, s3;
      int * s2s3 = s2s3_from_k[k];
      s2 = s2s3[0];
      s3 = s2s3[1];
      s1 = n - s2 - s3;
      while (s1 > 0){
            m *= X[0];
            s1--;
      }
      while (s2 > 0){
            m *= X[1];
            s2--;
      }
      while (s3 > 0){
            m *= X[2];
            s3--;
      }
      *(M + k) += m;
}


void get_Xn_overwrite(int k, int n, typ * X, typ m, typ * M){

      /******** Stores in tensor M the k^th component of tensor m * X^(n) ********/
      /******** Tensor M is overwritten                                   ********/
      
      int s1, s2, s3;
      int * s2s3 = s2s3_from_k[k];
      s2 = s2s3[0];
      s3 = s2s3[1];
      s1 = n - s2 - s3;
      while (s1 > 0){
            m *= X[0];
            s1--;
      }
      while (s2 > 0){
            m *= X[1];
            s2--;
      }
      while (s3 > 0){
            m *= X[2];
            s3--;
      }
      *(M + k) = m;
}


#if expansion_order >= 3 && mutual_bool
void get_Mn(struct node * FlatTree, struct moonlet * moonlets, int a){

      /******** Computes the multipole moments of node a of FlatTree, assuming that it has ********/
      /******** no children. Initializes the corresponding fields of the FlatTree          ********/
      /******** This function is called only if the expansion order is at least 3          ********/

      /******** Retrieving the mass and center of mass of the cell ********/
      typ com[3];
      com[0] = ((FlatTree + a) -> com)[0];  com[1] = ((FlatTree + a) -> com)[1];  com[2] = ((FlatTree + a) -> com)[2];  
      typ M2[6]  = {0., 0., 0., 0., 0., 0.};                                                             //The 6  distinct components of the quadrupole
      #if expansion_order >= 4
      typ M3[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};                                             //The 10 distinct components of the octupole
      #endif
      #if expansion_order >= 5
      typ M4[15] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};                         //The 15 distinct components of the fourth multipole moment
      #endif
      #if expansion_order >= 6
      typ M5[21] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 21 distinct components of the fifth multipole moment
      #endif
      #if expansion_order >= 7
      typ M6[28] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      #endif                                                                                             //The 28 distinct components of the sixth multipole moment
      #if expansion_order >= 8
      typ M7[36] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; 
      #endif                                                                                             //The 36 distinct components of the seventh multipole moment

      typ m; //The body's mass   
      int * dots        = (FlatTree + a) -> dots; //All the bodies in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int p; //Id of the current body
      int k; //Array index
      typ dX[3]; // x_i - com
      
      for (i = 0; i < how_many_dots; i ++){
            p     = dots[i];
            m     = (moonlets + p) -> mass;
            dX[0] = (moonlets + p) -> x - com[0];
            dX[1] = (moonlets + p) -> y - com[1];
            dX[2] = (moonlets + p) -> z - com[2];
            
            /******** Accumulating the multipole moments ********/
            #if expansion_order >= 3 //M2 needs to be computed
                  M2[0] += m * dX[0] * dX[0]; M2[1] += m * dX[0] * dX[1]; M2[2] += m * dX[0] * dX[2];
                  M2[3] += m * dX[1] * dX[1]; M2[4] += m * dX[1] * dX[2]; M2[5] += m * dX[2] * dX[2];
                  #if expansion_order >= 4 //M3 needs to be computed
                        for (k = 0; k < 10; k ++){
                              get_Xn(k, 3, dX, m, M3);
                        }
                        #if expansion_order >= 5 //M4 needs to be computed
                              for (k = 0; k < 15; k ++){
                                    get_Xn(k, 4, dX, m, M4);
                              }
                              #if expansion_order >= 6 //M5 needs to be computed
                                    for (k = 0; k < 21; k ++){
                                          get_Xn(k, 5, dX, m, M5);
                                    }
                                    #if expansion_order >= 7 //M6 needs to be computed
                                          for (k = 0; k < 28; k ++){
                                                get_Xn(k, 6, dX, m, M6);
                                          }
                                          #if expansion_order >= 8 //M7 needs to be computed
                                                for (k = 0; k < 36; k ++){
                                                      get_Xn(k, 7, dX, m, M7);
                                                }
                                          #endif
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      }
      
      /******** Initializing the fields M2 to M5 of the FlatTree ********/
      #if expansion_order >= 3
            for (k = 0; k < 6; k ++){
                  ((FlatTree + a) -> M2)[k] = M2[k];
            }
            #if expansion_order >= 4
                  for (k = 0; k < 10; k ++){
                        ((FlatTree + a) -> M3)[k] = M3[k];
                  }
                  #if expansion_order >= 5
                        for (k = 0; k < 15; k ++){
                              ((FlatTree + a) -> M4)[k] = M4[k];
                        }
                        #if expansion_order >= 6
                              for (k = 0; k < 21; k ++){
                                    ((FlatTree + a) -> M5)[k] = M5[k];
                              }
                              #if expansion_order >= 7
                                    for (k = 0; k < 28; k ++){
                                          ((FlatTree + a) -> M6)[k] = M6[k];
                                    }
                                    #if expansion_order >= 8
                                          for (k = 0; k < 36; k ++){
                                                ((FlatTree + a) -> M7)[k] = M7[k];
                                          }
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      #endif
}
#endif


void get_com_from_children(struct node * FlatTree, int a){

      /******** Computes the mass and center of mass of node a of FlatTree from that of ********/
      /******** its children. Initializes the corresponding fields of (FlatTree + a)    ********/

      int idFirstChild = (FlatTree + a) -> idFirstChild;
      int idLastChild  = idFirstChild + (FlatTree + a) -> how_many_children;
      int i;
      
      typ M0 = 0.0;
      typ com[3] = {0.0, 0.0, 0.0};
      typ M0_child;
      typ * com_child; //Center of mass of a child
      
      for (i = idFirstChild; i < idLastChild; i++){
            M0_child  = (FlatTree + i) -> M0;
            com_child = (FlatTree + i) -> com;
            
            M0 += M0_child;
            com[0] += M0_child * com_child[0];
            com[1] += M0_child * com_child[1];
            com[2] += M0_child * com_child[2];
            
      }

      /******** Initializing the relevant fields ********/
      (FlatTree  + a) -> M0 = M0;
      ((FlatTree + a) -> com)[0] = com[0]/M0;
      ((FlatTree + a) -> com)[1] = com[1]/M0;
      ((FlatTree + a) -> com)[2] = com[2]/M0;
}


void get_rmax_from_children(struct node * FlatTree, int a){

      /******** Initializes the convergence radius of node a of FlatTree ********/
      /******** from that of its children.                               ********/

      int idFirstChild = (FlatTree + a) -> idFirstChild;
      int idLastChild  = idFirstChild + (FlatTree + a) -> how_many_children;
      int i;

      typ * center = (FlatTree + a) -> center;     //The center of the node
      typ * com    = (FlatTree + a) -> com;        //The center of mass of the node
      typ D        = (FlatTree + a) -> sidelength; //The sidelength of the node
      typ corner[8][3] = {{center[0] - D - com[0], center[1] - D - com[1], center[2] - D - com[2]},
                          {center[0] - D - com[0], center[1] - D - com[1], center[2] + D - com[2]},
                          {center[0] - D - com[0], center[1] + D - com[1], center[2] - D - com[2]},
                          {center[0] - D - com[0], center[1] + D - com[1], center[2] + D - com[2]},
                          {center[0] + D - com[0], center[1] - D - com[1], center[2] - D - com[2]},
                          {center[0] + D - com[0], center[1] - D - com[1], center[2] + D - com[2]},
                          {center[0] + D - com[0], center[1] + D - com[1], center[2] - D - com[2]},
                          {center[0] + D - com[0], center[1] + D - com[1], center[2] + D - com[2]}}; //corner - center of mass
      typ distance_to_farthest_corner = 0.0;
      typ distance;

      /******** Computing the distance between the center of mass and the most distant corner ********/
      for (i = 0; i < 8; i ++){
            distance = sqrt(corner[i][0]*corner[i][0] + corner[i][1]*corner[i][1] + corner[i][2]*corner[i][2]);
            if (distance > distance_to_farthest_corner){
                  distance_to_farthest_corner = distance;
            }
      }
      
      typ * com_child; //Center of mass of a child
      typ rmax_child;  //Convergence radius of a child
      typ com_difference[3];
      typ max_ri = 0.;

      /******** Eq. (9) of Dehnen (2002) ********/
      for (i = idFirstChild; i < idLastChild; i++){
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
            (FlatTree + a) -> r_max = max_ri;
      }
      else{
            (FlatTree + a) -> r_max = distance_to_farthest_corner;
      }
}


#if expansion_order >= 3 && mutual_bool
void get_Mn_from_children(struct node * FlatTree, int a){

      /******** Computes the multipole moments of node a of FlatTree from that of its children ********/
      /******** Initializes the corresponding fields of the FlatTree                           ********/
      /******** This function is called only if the expansion order is at least 3              ********/

      int idFirstChild = (FlatTree + a) -> idFirstChild;
      int idLastChild  = idFirstChild + (FlatTree + a) -> how_many_children;
      typ * com = (FlatTree + a) -> com;
      typ * com_child;
      typ M0_child;
      typ Y[4]; //com - com_child. Index 0 is unused
      int p;
      
      #if expansion_order >= 3
      typ M2[6]  = {0., 0., 0., 0., 0., 0.};                                                             //The 6  distinct components of the quadrupole
      typ * M2_child;
      #endif
      #if expansion_order >= 4
      typ M3[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};                                             //The 10 distinct components of the octupole
      typ * M3_child;
      int * array_of_ijklm;
      int q;
      int i, j, k;
      #endif
      #if expansion_order >= 5
      typ M4[15] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};                         //The 15 distinct components of the fourth multipole moment
      typ * M4_child;
      int l;
      #endif
      #if expansion_order >= 6
      typ M5[21] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 21 distinct components of the fifth multipole moment
      typ * M5_child;
      int m;
      #endif
      #if expansion_order >= 7
      typ * M6_child;
      int n;
      typ M6[28] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      #endif                                                                                             //The 28 distinct components of the sixth multipole moment
      #if expansion_order >= 8
      typ * M7_child;
      int o;
      typ M7[36] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      #endif                                                                                             //The 36 distinct components of the seventh multipole moment
      
      
      for (p = idFirstChild; p < idLastChild; p ++){
            /******** Retrieving the multipole moments and the center of mass of the child ********/
            M0_child  = (FlatTree + p) -> M0;
            com_child = (FlatTree + p) -> com;
            #if expansion_order >= 3
                  M2_child = (FlatTree + p) -> M2;
                  #if expansion_order >= 4
                        M3_child = (FlatTree + p) -> M3;
                        #if expansion_order >= 5
                              M4_child = (FlatTree + p) -> M4;
                              #if expansion_order >= 6
                                    M5_child = (FlatTree + p) -> M5;
                                    #if expansion_order >= 7
                                          M6_child = (FlatTree + p) -> M6;
                                          #if expansion_order >= 8
                                                M7_child = (FlatTree + p) -> M7;
                                          #endif
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
            
            /******** Retrieving the difference Y between the expansion centers ********/
            Y[1] = com[0] - com_child[0];  Y[2] = com[1] - com_child[1];  Y[3] = com[2] - com_child[2];  
            
            /******** Computing the multipole moments ********/
            #if expansion_order >= 3
                  M2[0] += M2_child[0] + M0_child*Y[1]*Y[1];
                  M2[1] += M2_child[1] + M0_child*Y[1]*Y[2];
                  M2[2] += M2_child[2] + M0_child*Y[1]*Y[3];
                  M2[3] += M2_child[3] + M0_child*Y[2]*Y[2];
                  M2[4] += M2_child[4] + M0_child*Y[2]*Y[3];
                  M2[5] += M2_child[5] + M0_child*Y[3]*Y[3];
            #if expansion_order >= 4 //M3 is computed
                  for (q = 0; q < 10; q ++){
                        array_of_ijklm = ijklmn_from_kn[q][3];
                        i = array_of_ijklm[0]; j = array_of_ijklm[1]; k = array_of_ijklm[2]; 
                        M3[q] += M3_child[q]; //Term in M^3
                        M3[q] -= M0_child*Y[i]*Y[j]*Y[k]; //Term in M^0
                        M3[q] -= Y[k]*M2_child[k_from_ijklmn[i][j][0][0][0][0]]; //Terms in M^2
                        M3[q] -= Y[i]*M2_child[k_from_ijklmn[j][k][0][0][0][0]];
                        M3[q] -= Y[j]*M2_child[k_from_ijklmn[i][k][0][0][0][0]];
                  }
            #if expansion_order >= 5 //M4 is computed
                  for (q = 0; q < 15; q ++){
                        array_of_ijklm = ijklmn_from_kn[q][4];
                        i = array_of_ijklm[0]; j = array_of_ijklm[1]; k = array_of_ijklm[2]; l = array_of_ijklm[3];
                        M4[q] += M4_child[q]; //Term in M^4
                        M4[q] += M0_child*Y[i]*Y[j]*Y[k]*Y[l]; //Term in M^0
                        M4[q] -= Y[l]*M3_child[k_from_ijklmn[i][j][k][0][0][0]]; //Terms in M^3
                        M4[q] -= Y[k]*M3_child[k_from_ijklmn[i][j][l][0][0][0]];
                        M4[q] -= Y[j]*M3_child[k_from_ijklmn[i][k][l][0][0][0]];
                        M4[q] -= Y[i]*M3_child[k_from_ijklmn[j][k][l][0][0][0]];
                        M4[q] += Y[k]*Y[l]*M2_child[k_from_ijklmn[i][j][0][0][0][0]]; //Terms in M^2
                        M4[q] += Y[j]*Y[l]*M2_child[k_from_ijklmn[i][k][0][0][0][0]];
                        M4[q] += Y[j]*Y[k]*M2_child[k_from_ijklmn[i][l][0][0][0][0]];
                        M4[q] += Y[i]*Y[l]*M2_child[k_from_ijklmn[j][k][0][0][0][0]];
                        M4[q] += Y[i]*Y[k]*M2_child[k_from_ijklmn[j][l][0][0][0][0]];
                        M4[q] += Y[i]*Y[j]*M2_child[k_from_ijklmn[k][l][0][0][0][0]];
                  }
            #if expansion_order >= 6 //M5 is computed
                  for (q = 0; q < 21; q ++){
                        array_of_ijklm = ijklmn_from_kn[q][5];
                        i = array_of_ijklm[0]; j = array_of_ijklm[1]; k = array_of_ijklm[2]; l = array_of_ijklm[3]; m = array_of_ijklm[4];
                        M5[q] += M5_child[q]; //Term in M^5
                        M5[q] -= M0_child*Y[i]*Y[j]*Y[k]*Y[l]*Y[m]; //Term in M^0
                        M5[q] -= Y[m]*M4_child[k_from_ijklmn[i][j][k][l][0][0]]; //Terms in M^4
                        M5[q] -= Y[l]*M4_child[k_from_ijklmn[i][j][k][m][0][0]];
                        M5[q] -= Y[k]*M4_child[k_from_ijklmn[i][j][l][m][0][0]];
                        M5[q] -= Y[j]*M4_child[k_from_ijklmn[i][k][l][m][0][0]];
                        M5[q] -= Y[i]*M4_child[k_from_ijklmn[j][k][l][m][0][0]];
                        M5[q] += Y[l]*Y[m]*M3_child[k_from_ijklmn[i][j][k][0][0][0]]; //Terms in M^3
                        M5[q] += Y[k]*Y[m]*M3_child[k_from_ijklmn[i][j][l][0][0][0]];
                        M5[q] += Y[k]*Y[l]*M3_child[k_from_ijklmn[i][j][m][0][0][0]];
                        M5[q] += Y[j]*Y[m]*M3_child[k_from_ijklmn[i][k][l][0][0][0]];
                        M5[q] += Y[j]*Y[l]*M3_child[k_from_ijklmn[i][k][m][0][0][0]];
                        M5[q] += Y[j]*Y[k]*M3_child[k_from_ijklmn[i][l][m][0][0][0]];
                        M5[q] += Y[i]*Y[m]*M3_child[k_from_ijklmn[j][k][l][0][0][0]];
                        M5[q] += Y[i]*Y[l]*M3_child[k_from_ijklmn[j][k][m][0][0][0]];
                        M5[q] += Y[i]*Y[k]*M3_child[k_from_ijklmn[j][l][m][0][0][0]];
                        M5[q] += Y[i]*Y[j]*M3_child[k_from_ijklmn[k][l][m][0][0][0]];
                        M5[q] -= Y[i]*Y[j]*Y[k]*M2_child[k_from_ijklmn[l][m][0][0][0][0]]; //Terms in M^2
                        M5[q] -= Y[i]*Y[j]*Y[l]*M2_child[k_from_ijklmn[k][m][0][0][0][0]];
                        M5[q] -= Y[i]*Y[j]*Y[m]*M2_child[k_from_ijklmn[k][l][0][0][0][0]];
                        M5[q] -= Y[i]*Y[k]*Y[l]*M2_child[k_from_ijklmn[j][m][0][0][0][0]];
                        M5[q] -= Y[i]*Y[k]*Y[m]*M2_child[k_from_ijklmn[j][l][0][0][0][0]];
                        M5[q] -= Y[i]*Y[l]*Y[m]*M2_child[k_from_ijklmn[j][k][0][0][0][0]];
                        M5[q] -= Y[j]*Y[k]*Y[l]*M2_child[k_from_ijklmn[i][m][0][0][0][0]];
                        M5[q] -= Y[j]*Y[k]*Y[m]*M2_child[k_from_ijklmn[i][l][0][0][0][0]];
                        M5[q] -= Y[j]*Y[l]*Y[m]*M2_child[k_from_ijklmn[i][k][0][0][0][0]];
                        M5[q] -= Y[k]*Y[l]*Y[m]*M2_child[k_from_ijklmn[i][j][0][0][0][0]];
                  }
            #if expansion_order >= 7 //M6 is computed
                  for (q = 0; q < 28; q ++){
                        array_of_ijklm = ijklmn_from_kn[q][6];
                        i = array_of_ijklm[0]; j = array_of_ijklm[1]; k = array_of_ijklm[2];
                        l = array_of_ijklm[3]; m = array_of_ijklm[4]; n = array_of_ijklm[5];
                        M6[q] += M6_child[q]; //Term in M^6
                        M6[q] -= M5_child[k_from_ijklmn[i][j][k][l][m][0]]*Y[n]; //Terms in M^5
                        M6[q] -= M5_child[k_from_ijklmn[i][j][k][l][n][0]]*Y[m];
                        M6[q] -= M5_child[k_from_ijklmn[i][j][k][m][n][0]]*Y[l];
                        M6[q] -= M5_child[k_from_ijklmn[i][j][l][m][n][0]]*Y[k];
                        M6[q] -= M5_child[k_from_ijklmn[i][k][l][m][n][0]]*Y[j];
                        M6[q] -= M5_child[k_from_ijklmn[j][k][l][m][n][0]]*Y[i];
                        M6[q] += M4_child[k_from_ijklmn[i][j][k][l][0][0]]*Y[m]*Y[n]; //Terms in M^4
                        M6[q] += M4_child[k_from_ijklmn[i][j][k][m][0][0]]*Y[l]*Y[n];
                        M6[q] += M4_child[k_from_ijklmn[i][j][k][n][0][0]]*Y[l]*Y[m];
                        M6[q] += M4_child[k_from_ijklmn[i][j][l][m][0][0]]*Y[k]*Y[n];
                        M6[q] += M4_child[k_from_ijklmn[i][j][l][n][0][0]]*Y[k]*Y[m];
                        M6[q] += M4_child[k_from_ijklmn[i][j][m][n][0][0]]*Y[k]*Y[l];
                        M6[q] += M4_child[k_from_ijklmn[i][k][l][m][0][0]]*Y[j]*Y[n];
                        M6[q] += M4_child[k_from_ijklmn[i][k][l][n][0][0]]*Y[j]*Y[m];
                        M6[q] += M4_child[k_from_ijklmn[i][k][m][n][0][0]]*Y[j]*Y[l];
                        M6[q] += M4_child[k_from_ijklmn[i][l][m][n][0][0]]*Y[j]*Y[k];
                        M6[q] += M4_child[k_from_ijklmn[j][k][l][m][0][0]]*Y[i]*Y[n];
                        M6[q] += M4_child[k_from_ijklmn[j][k][l][n][0][0]]*Y[i]*Y[m];
                        M6[q] += M4_child[k_from_ijklmn[j][k][m][n][0][0]]*Y[i]*Y[l];
                        M6[q] += M4_child[k_from_ijklmn[j][l][m][n][0][0]]*Y[i]*Y[k];
                        M6[q] += M4_child[k_from_ijklmn[k][l][m][n][0][0]]*Y[i]*Y[j];
                        M6[q] -= M3_child[k_from_ijklmn[i][j][k][0][0][0]]*Y[l]*Y[m]*Y[n]; //Terms in M^3
                        M6[q] -= M3_child[k_from_ijklmn[i][j][l][0][0][0]]*Y[k]*Y[m]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[i][j][m][0][0][0]]*Y[k]*Y[l]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[i][j][n][0][0][0]]*Y[k]*Y[l]*Y[m];
                        M6[q] -= M3_child[k_from_ijklmn[i][k][l][0][0][0]]*Y[j]*Y[m]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[i][k][m][0][0][0]]*Y[j]*Y[l]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[i][k][n][0][0][0]]*Y[j]*Y[l]*Y[m];
                        M6[q] -= M3_child[k_from_ijklmn[i][l][m][0][0][0]]*Y[j]*Y[k]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[i][l][n][0][0][0]]*Y[j]*Y[k]*Y[m];
                        M6[q] -= M3_child[k_from_ijklmn[i][m][n][0][0][0]]*Y[j]*Y[k]*Y[l];
                        M6[q] -= M3_child[k_from_ijklmn[j][k][l][0][0][0]]*Y[i]*Y[m]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[j][k][m][0][0][0]]*Y[i]*Y[l]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[j][k][n][0][0][0]]*Y[i]*Y[l]*Y[m];
                        M6[q] -= M3_child[k_from_ijklmn[j][l][m][0][0][0]]*Y[i]*Y[k]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[j][l][n][0][0][0]]*Y[i]*Y[k]*Y[m];
                        M6[q] -= M3_child[k_from_ijklmn[j][m][n][0][0][0]]*Y[i]*Y[k]*Y[l];
                        M6[q] -= M3_child[k_from_ijklmn[k][l][m][0][0][0]]*Y[i]*Y[j]*Y[n];
                        M6[q] -= M3_child[k_from_ijklmn[k][l][n][0][0][0]]*Y[i]*Y[j]*Y[m];
                        M6[q] -= M3_child[k_from_ijklmn[k][m][n][0][0][0]]*Y[i]*Y[j]*Y[l];
                        M6[q] -= M3_child[k_from_ijklmn[l][m][n][0][0][0]]*Y[i]*Y[j]*Y[k];
                        M6[q] += M2_child[k_from_ijklmn[i][j][0][0][0][0]]*Y[k]*Y[l]*Y[m]*Y[n]; //Terms in M^2
                        M6[q] += M2_child[k_from_ijklmn[i][k][0][0][0][0]]*Y[j]*Y[l]*Y[m]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[i][l][0][0][0][0]]*Y[j]*Y[k]*Y[m]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[i][m][0][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[i][n][0][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[m];
                        M6[q] += M2_child[k_from_ijklmn[j][k][0][0][0][0]]*Y[i]*Y[l]*Y[m]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[j][l][0][0][0][0]]*Y[i]*Y[k]*Y[m]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[j][m][0][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[j][n][0][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[m];
                        M6[q] += M2_child[k_from_ijklmn[k][l][0][0][0][0]]*Y[i]*Y[j]*Y[m]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[k][m][0][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[k][n][0][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[m];
                        M6[q] += M2_child[k_from_ijklmn[l][m][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[n];
                        M6[q] += M2_child[k_from_ijklmn[l][n][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[m];
                        M6[q] += M2_child[k_from_ijklmn[m][n][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[l];
                        M6[q] += M0_child*Y[i]*Y[j]*Y[k]*Y[l]*Y[m]*Y[n]; //Term in M^0
                  }
            #if expansion_order >= 8 //M7 is computed
                  for (q = 0; q < 36; q ++){
                        array_of_ijklm = ijklmn_from_kn[q][7];
                        i = array_of_ijklm[0]; j = array_of_ijklm[1]; k = array_of_ijklm[2];
                        l = array_of_ijklm[3]; m = array_of_ijklm[4]; n = array_of_ijklm[5]; o = array_of_ijklm[6];
                        M7[q] += M7_child[q]; //Term in M^7
                        M7[q] -= M6_child[k_from_ijklmn[i][j][k][l][m][n]]*Y[o]; //Terms in M^6
                        M7[q] -= M6_child[k_from_ijklmn[i][j][k][l][m][o]]*Y[n];
                        M7[q] -= M6_child[k_from_ijklmn[i][j][k][l][n][o]]*Y[m];
                        M7[q] -= M6_child[k_from_ijklmn[i][j][k][m][n][o]]*Y[l];
                        M7[q] -= M6_child[k_from_ijklmn[i][j][l][m][n][o]]*Y[k];
                        M7[q] -= M6_child[k_from_ijklmn[i][k][l][m][n][o]]*Y[j];
                        M7[q] -= M6_child[k_from_ijklmn[j][k][l][m][n][o]]*Y[i];
                        M7[q] += M5_child[k_from_ijklmn[i][j][k][l][m][0]]*Y[n]*Y[o]; //Terms in M^5
                        M7[q] += M5_child[k_from_ijklmn[i][j][k][l][n][0]]*Y[m]*Y[o];
                        M7[q] += M5_child[k_from_ijklmn[i][j][k][l][o][0]]*Y[m]*Y[n];
                        M7[q] += M5_child[k_from_ijklmn[i][j][k][m][n][0]]*Y[l]*Y[o];
                        M7[q] += M5_child[k_from_ijklmn[i][j][k][m][o][0]]*Y[l]*Y[n];
                        M7[q] += M5_child[k_from_ijklmn[i][j][k][n][o][0]]*Y[l]*Y[m];
                        M7[q] += M5_child[k_from_ijklmn[i][j][l][m][n][0]]*Y[k]*Y[o];
                        M7[q] += M5_child[k_from_ijklmn[i][j][l][m][o][0]]*Y[k]*Y[n];
                        M7[q] += M5_child[k_from_ijklmn[i][j][l][n][o][0]]*Y[k]*Y[m];
                        M7[q] += M5_child[k_from_ijklmn[i][j][m][n][o][0]]*Y[k]*Y[l];
                        M7[q] += M5_child[k_from_ijklmn[i][k][l][m][n][0]]*Y[j]*Y[o];
                        M7[q] += M5_child[k_from_ijklmn[i][k][l][m][o][0]]*Y[j]*Y[n];
                        M7[q] += M5_child[k_from_ijklmn[i][k][l][n][o][0]]*Y[j]*Y[m];
                        M7[q] += M5_child[k_from_ijklmn[i][k][m][n][o][0]]*Y[j]*Y[l];
                        M7[q] += M5_child[k_from_ijklmn[i][l][m][n][o][0]]*Y[j]*Y[k];
                        M7[q] += M5_child[k_from_ijklmn[j][k][l][m][n][0]]*Y[i]*Y[o];
                        M7[q] += M5_child[k_from_ijklmn[j][k][l][m][o][0]]*Y[i]*Y[n];
                        M7[q] += M5_child[k_from_ijklmn[j][k][l][n][o][0]]*Y[i]*Y[m];
                        M7[q] += M5_child[k_from_ijklmn[j][k][m][n][o][0]]*Y[i]*Y[l];
                        M7[q] += M5_child[k_from_ijklmn[j][l][m][n][o][0]]*Y[i]*Y[k];
                        M7[q] += M5_child[k_from_ijklmn[k][l][m][n][o][0]]*Y[i]*Y[j];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][k][l][0][0]]*Y[m]*Y[n]*Y[o]; //Terms in M^4
                        M7[q] -= M4_child[k_from_ijklmn[i][j][k][m][0][0]]*Y[l]*Y[n]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][k][n][0][0]]*Y[l]*Y[m]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][k][o][0][0]]*Y[l]*Y[m]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][l][m][0][0]]*Y[k]*Y[n]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][l][n][0][0]]*Y[k]*Y[m]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][l][o][0][0]]*Y[k]*Y[m]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][m][n][0][0]]*Y[k]*Y[l]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][m][o][0][0]]*Y[k]*Y[l]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[i][j][n][o][0][0]]*Y[k]*Y[l]*Y[m];
                        M7[q] -= M4_child[k_from_ijklmn[i][k][l][m][0][0]]*Y[j]*Y[n]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][k][l][n][0][0]]*Y[j]*Y[m]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][k][l][o][0][0]]*Y[j]*Y[m]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[i][k][m][n][0][0]]*Y[j]*Y[l]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][k][m][o][0][0]]*Y[j]*Y[l]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[i][k][n][o][0][0]]*Y[j]*Y[l]*Y[m];
                        M7[q] -= M4_child[k_from_ijklmn[i][l][m][n][0][0]]*Y[j]*Y[k]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[i][l][m][o][0][0]]*Y[j]*Y[k]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[i][l][n][o][0][0]]*Y[j]*Y[k]*Y[m];
                        M7[q] -= M4_child[k_from_ijklmn[i][m][n][o][0][0]]*Y[j]*Y[k]*Y[l];
                        M7[q] -= M4_child[k_from_ijklmn[j][k][l][m][0][0]]*Y[i]*Y[n]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[j][k][l][n][0][0]]*Y[i]*Y[m]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[j][k][l][o][0][0]]*Y[i]*Y[m]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[j][k][m][n][0][0]]*Y[i]*Y[l]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[j][k][m][o][0][0]]*Y[i]*Y[l]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[j][k][n][o][0][0]]*Y[i]*Y[l]*Y[m];
                        M7[q] -= M4_child[k_from_ijklmn[j][l][m][n][0][0]]*Y[i]*Y[k]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[j][l][m][o][0][0]]*Y[i]*Y[k]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[j][l][n][o][0][0]]*Y[i]*Y[k]*Y[m];
                        M7[q] -= M4_child[k_from_ijklmn[j][m][n][o][0][0]]*Y[i]*Y[k]*Y[l];
                        M7[q] -= M4_child[k_from_ijklmn[k][l][m][n][0][0]]*Y[i]*Y[j]*Y[o];
                        M7[q] -= M4_child[k_from_ijklmn[k][l][m][o][0][0]]*Y[i]*Y[j]*Y[n];
                        M7[q] -= M4_child[k_from_ijklmn[k][l][n][o][0][0]]*Y[i]*Y[j]*Y[m];
                        M7[q] -= M4_child[k_from_ijklmn[k][m][n][o][0][0]]*Y[i]*Y[j]*Y[l];
                        M7[q] -= M4_child[k_from_ijklmn[l][m][n][o][0][0]]*Y[i]*Y[j]*Y[k];
                        M7[q] += M3_child[k_from_ijklmn[m][n][o][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[l]; //Terms in M^3
                        M7[q] += M3_child[k_from_ijklmn[l][n][o][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[m];
                        M7[q] += M3_child[k_from_ijklmn[l][m][o][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[l][m][n][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[k][n][o][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[m];
                        M7[q] += M3_child[k_from_ijklmn[k][m][o][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[k][m][n][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[o];                       
                        M7[q] += M3_child[k_from_ijklmn[k][l][o][0][0][0]]*Y[i]*Y[j]*Y[m]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[k][l][n][0][0][0]]*Y[i]*Y[j]*Y[m]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[k][l][m][0][0][0]]*Y[i]*Y[j]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[j][n][o][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[m];
                        M7[q] += M3_child[k_from_ijklmn[j][m][o][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[j][m][n][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[j][l][o][0][0][0]]*Y[i]*Y[k]*Y[m]*Y[n];                       
                        M7[q] += M3_child[k_from_ijklmn[j][l][n][0][0][0]]*Y[i]*Y[k]*Y[m]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[j][l][m][0][0][0]]*Y[i]*Y[k]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[j][k][o][0][0][0]]*Y[i]*Y[l]*Y[m]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[j][k][n][0][0][0]]*Y[i]*Y[l]*Y[m]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[j][k][m][0][0][0]]*Y[i]*Y[l]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[j][k][l][0][0][0]]*Y[i]*Y[m]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][n][o][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[m];                      
                        M7[q] += M3_child[k_from_ijklmn[i][m][o][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[i][m][n][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][l][o][0][0][0]]*Y[j]*Y[k]*Y[m]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[i][l][n][0][0][0]]*Y[j]*Y[k]*Y[m]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][l][m][0][0][0]]*Y[j]*Y[k]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][k][o][0][0][0]]*Y[j]*Y[l]*Y[m]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[i][k][n][0][0][0]]*Y[j]*Y[l]*Y[m]*Y[o];                        
                        M7[q] += M3_child[k_from_ijklmn[i][k][m][0][0][0]]*Y[j]*Y[l]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][k][l][0][0][0]]*Y[j]*Y[m]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][j][o][0][0][0]]*Y[k]*Y[l]*Y[m]*Y[n];
                        M7[q] += M3_child[k_from_ijklmn[i][j][n][0][0][0]]*Y[k]*Y[l]*Y[m]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][j][m][0][0][0]]*Y[k]*Y[l]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][j][l][0][0][0]]*Y[k]*Y[m]*Y[n]*Y[o];
                        M7[q] += M3_child[k_from_ijklmn[i][j][k][0][0][0]]*Y[l]*Y[m]*Y[n]*Y[o];                       
                        M7[q] -= M2_child[k_from_ijklmn[n][o][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[l]*Y[m]; //Terms in M^2
                        M7[q] -= M2_child[k_from_ijklmn[m][o][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[l]*Y[n];
                        M7[q] -= M2_child[k_from_ijklmn[m][n][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[l]*Y[o];                        
                        M7[q] -= M2_child[k_from_ijklmn[l][o][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[m]*Y[n];
                        M7[q] -= M2_child[k_from_ijklmn[l][n][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[m]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[l][m][0][0][0][0]]*Y[i]*Y[j]*Y[k]*Y[n]*Y[o];                        
                        M7[q] -= M2_child[k_from_ijklmn[k][o][0][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[m]*Y[n];
                        M7[q] -= M2_child[k_from_ijklmn[k][n][0][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[m]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[k][m][0][0][0][0]]*Y[i]*Y[j]*Y[l]*Y[n]*Y[o];                    
                        M7[q] -= M2_child[k_from_ijklmn[k][l][0][0][0][0]]*Y[i]*Y[j]*Y[m]*Y[n]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[j][o][0][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[m]*Y[n];
                        M7[q] -= M2_child[k_from_ijklmn[j][n][0][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[m]*Y[o];                      
                        M7[q] -= M2_child[k_from_ijklmn[j][m][0][0][0][0]]*Y[i]*Y[k]*Y[l]*Y[n]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[j][l][0][0][0][0]]*Y[i]*Y[k]*Y[m]*Y[n]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[j][k][0][0][0][0]]*Y[i]*Y[l]*Y[m]*Y[n]*Y[o];                       
                        M7[q] -= M2_child[k_from_ijklmn[i][o][0][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[m]*Y[n];
                        M7[q] -= M2_child[k_from_ijklmn[i][n][0][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[m]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[i][m][0][0][0][0]]*Y[j]*Y[k]*Y[l]*Y[n]*Y[o];                        
                        M7[q] -= M2_child[k_from_ijklmn[i][l][0][0][0][0]]*Y[j]*Y[k]*Y[m]*Y[n]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[i][k][0][0][0][0]]*Y[j]*Y[l]*Y[m]*Y[n]*Y[o];
                        M7[q] -= M2_child[k_from_ijklmn[i][j][0][0][0][0]]*Y[k]*Y[l]*Y[m]*Y[n]*Y[o];
                        M7[q] -= M0_child*Y[i]*Y[j]*Y[k]*Y[l]*Y[m]*Y[n]*Y[o]; //Term in M^0
                  }
            #endif
            #endif
            #endif
            #endif
            #endif
            #endif
      }
      
      /******** Initializing the fields M2 to M5 of the FlatTree ********/
      #if expansion_order >= 3
            for (p = 0; p < 6; p ++){
                  ((FlatTree + a) -> M2)[p] = M2[p];
            }
            #if expansion_order >= 4
                  for (p = 0; p < 10; p ++){
                        ((FlatTree + a) -> M3)[p] = M3[p];
                  }
                  #if expansion_order >= 5
                        for (p = 0; p < 15; p ++){
                              ((FlatTree + a) -> M4)[p] = M4[p];
                        }
                        #if expansion_order >= 6
                              for (p = 0; p < 21; p ++){
                                    ((FlatTree + a) -> M5)[p] = M5[p];
                              }
                              #if expansion_order >= 7
                                    for (p = 0; p < 28; p ++){
                                          ((FlatTree + a) -> M6)[p] = M6[p];
                                    }
                                    #if expansion_order >= 8
                                          for (p = 0; p < 36; p ++){
                                                ((FlatTree + a) -> M7)[p] = M7[p];
                                          }
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      #endif
}
#endif


void com_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Computes and initializes the masses and center of masses of all the nodes of FlatTree ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that could not be treated due to their child not having been treated yet
      int i; //Id of the current node
      int j = 0; //Index of where to put a node in the stack
      int how_many_children;
      
      /******** I travel the flattree. If a node is childless, I compute its mass and center of mass ********/
      /******** Otherwise, I store it in the stack for future treatment                              ********/
      for (i = 0; i < cell_id; i ++){
            how_many_children = (FlatTree + i) -> how_many_children;
            if (how_many_children == 0){ //If the node has no children
                  get_com(FlatTree, moonlets, i);                  
            }
            else{ //If the node has children
                  stack[j] = i;
                  j ++;
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function com_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** I now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j --;
            i = stack[j];
            get_com_from_children(FlatTree, i);
      }
      
      free(stack);
      stack = NULL;
}


void rmax_flattree(struct node * FlatTree, struct moonlet * moonlets){

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
                  get_rmax(FlatTree, moonlets, i);
            }
            else{ //If the node has children
                  stack[j] = i;
                  j ++;
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function rmax_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** I now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j --;
            i = stack[j];
            get_rmax_from_children(FlatTree, i);
      }
      
      free(stack);
      stack = NULL;
}


void rcrit_flattree(struct node * FlatTree){

      /******** Computes and initializes r_crit for all the nodes of FlatTree ********/
      
      int i;
      
      /******** I travel the flattree and compute the tolerance parameter theta and r_crit for each node ********/
      for (i = 0; i < cell_id; i ++){
            get_tolerance_parameter(FlatTree, i, 1.0e-6);
      }
}


#if expansion_order >= 3 && mutual_bool
void multipole_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Computes and initializes the multipole moments of all the nodes of FlatTree ********/
      /******** This function is called only if expansion_order is at least 3               ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that could not be treated due to their child not having been treated yet
      int i; //Id of the current node
      int j = 0; //Index of where to put a node in the stack
      int how_many_children;
      int how_many_dots;
      int child_multipole_threshold = expansion_order <= 4 ? 1 : 2*expansion_order; //Probably not optimal
      
      /******** I travel the flattree. If a node is childless or contains few bodies, I compute ********/
      /******** its multipole moments. Otherwise, I store it in the stack for future treatment  ********/
      for (i = 0; i < cell_id; i ++){
            how_many_children = (FlatTree + i) -> how_many_children;
            if (how_many_children == 0){ //If the node has no children
                  get_Mn(FlatTree, moonlets, i);
            }
            else{ //If the node has children
                  how_many_dots = (FlatTree + i) -> how_many_dots;
                  if (how_many_dots / how_many_children <= child_multipole_threshold){
                        get_Mn(FlatTree, moonlets, i);
                  }
                  else{
                        stack[j] = i;
                        j ++;
                  }
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function multipole_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** I now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j --;
            i = stack[j];
            get_Mn_from_children(FlatTree, i);
      } 
      free(stack);
      stack = NULL;
}
#endif


/*************************************************************************************/
/******** I now write functions to compute the gradient of the Green function ********/
/******** and the inner product between two tensors                           ********/
/*************************************************************************************/


void gradR(typ * R, typ * grad){

      /******** This function computes the n^th gradient of G/R, for 1 <= n <= expansion_order, and stores them   ********/
      /******** into grad. The vector R has three components and is of the form R = {Rx, Ry, Rz}. It is assumed   ********/
      /******** that the array grad is indexed from 0 to (p + 1)(p(2p + 10)/12 + 1)-2, where p = expansion_order, ********/
      /******** so it can hold all the independant components of the p first gradients of G/R                     ********/

      typ R1 = R[0], R2 = R[1], R3 = R[2];
      typ normR2 = R1*R1 + R2*R2 + R3*R3;
      typ normR  = sqrt(normR2);
      typ D1     = -G/(normR2*normR);
      #if expansion_order >= 2
            typ D2 = -3.*D1/normR2;
            #if expansion_order >= 3
                  typ D3  = -5.*D2/normR2;
                  typ R12 = R1*R1, R22 = R2*R2, R32 = R3*R3;
                  #if expansion_order >= 4
                        typ D4    = -7.*D3/normR2;
                        typ R1_R2 = R1*R2,  R1_R3 = R1*R3,  R2_R3 = R2*R3;
                        typ R13   = R12*R1, R23   = R22*R2, R33   = R32*R3;
                        #if expansion_order >= 5
                              typ D5  = -9.*D4/normR2;
                              typ R14 = R12*R12, R24 = R22*R22, R34 = R32*R32;
                              #if expansion_order >= 6
                                    typ D6  = -11.*D5/normR2;
                                    typ R15 = R13*R12, R25 = R23*R22, R35 = R33*R32;
                                    #if expansion_order >= 7
                                          typ D7  = -13.*D6/normR2;
                                          typ R16 = R13*R13, R26 = R23*R23, R36 = R33*R33;
                                          #if expansion_order >= 8
                                                typ D8  = -15.*D7/normR2;
                                                typ R17 = R14*R13, R27 = R24*R23, R37 = R34*R33;
                                          #endif
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      #endif
      
      grad[0]   = D1*R1;       //Computing the first gradient
      grad[1]   = D1*R2;
      grad[2]   = D1*R3;
      #if expansion_order >= 2 //Computing the second gradient
      grad[3]   = D2*R1*R1 + D1;
      grad[4]   = D2*R1*R2;
      grad[5]   = D2*R1*R3;
      grad[6]   = D2*R2*R2 + D1;
      grad[7]   = D2*R2*R3;
      grad[8]   = D2*R3*R3 + D1;
      #if expansion_order >= 3 //Computing the third gradient
      grad[9]   = D3*R12*R1 + 3.0*D2*R1;
      grad[10]  = D3*R12*R2 +     D2*R2;
      grad[11]  = D3*R12*R3 +     D2*R3;
      grad[12]  = D3*R22*R1 +     D2*R1;
      grad[13]  = D3*R1*R2*R3;
      grad[14]  = D3*R32*R1 +     D2*R1;
      grad[15]  = D3*R22*R2 + 3.0*D2*R2;
      grad[16]  = D3*R22*R3 +     D2*R3;
      grad[17]  = D3*R32*R2 +     D2*R2;
      grad[18]  = D3*R32*R3 + 3.0*D2*R3;
      #if expansion_order >= 4 //Computing the fourth gradient
      grad[19]  = D4*R12*R12 + 6.0*D3*R12 +   3.0*D2;
      grad[20]  = D4*R13*R2  + 3.0*D3*R1_R2;
      grad[21]  = D4*R13*R3  + 3.0*D3*R1_R3;
      grad[22]  = D4*R12*R22 +     D3*(R12+R22) + D2;
      grad[23]  = D4*R12*R2_R3 +   D3*R2_R3;
      grad[24]  = D4*R12*R32 +     D3*(R12+R32) + D2;
      grad[25]  = D4*R23*R1  + 3.0*D3*R1_R2;
      grad[26]  = D4*R22*R1_R3 +   D3*R1_R3;
      grad[27]  = D4*R1_R2*R32 +   D3*R1_R2;
      grad[28]  = D4*R33*R1  + 3.0*D3*R1_R3;
      grad[29]  = D4*R22*R22 + 6.0*D3*R22 +   3.0*D2;
      grad[30]  = D4*R23*R3  + 3.0*D3*R2_R3;
      grad[31]  = D4*R22*R32 +     D3*(R22+R32) + D2;
      grad[32]  = D4*R33*R2  + 3.0*D3*R2_R3;
      grad[33]  = D4*R32*R32 + 6.0*D3*R32 +   3.0*D2;
      #if expansion_order >= 5 //Computing the fifth gradient
      grad[34]  = D5*R13*R12 +  10.0*D4*R13 +             15.0*D3*R1;
      grad[35]  = D5*R14*R2  +   6.0*D4*R12*R2 +           3.0*D3*R2;
      grad[36]  = D5*R14*R3  +   6.0*D4*R12*R3 +           3.0*D3*R3;
      grad[37]  = D5*R13*R22 +       D4*R1*(R12+3.0*R22) + 3.0*D3*R1;
      grad[38]  = D5*R13*R2_R3 + 3.0*D4*R1_R2*R3;
      grad[39]  = D5*R13*R32 +       D4*R1*(R12+3.0*R32) + 3.0*D3*R1;
      grad[40]  = D5*R23*R12 +       D4*R2*(R22+3.0*R12) + 3.0*D3*R2;
      grad[41]  = D5*R12*R22*R3 +    D4*R3*(R12+R22) +         D3*R3;
      grad[42]  = D5*R12*R2*R32 +    D4*R2*(R12+R32) +         D3*R2;
      grad[43]  = D5*R33*R12 +       D4*R3*(R32+3.0*R12) + 3.0*D3*R3;
      grad[44]  = D5*R24*R1  +   6.0*D4*R22*R1 +           3.0*D3*R1;
      grad[45]  = D5*R23*R1_R3 + 3.0*D4*R1_R2*R3;
      grad[46]  = D5*R1*R22*R32 +    D4*R1*(R32+R22) +         D3*R1;
      grad[47]  = D5*R1_R2*R33 + 3.0*D4*R1_R2*R3;
      grad[48]  = D5*R34*R1  +   6.0*D4*R32*R1 +           3.0*D3*R1;
      grad[49]  = D5*R23*R22 +  10.0*D4*R23 +             15.0*D3*R2;
      grad[50]  = D5*R24*R3  +   6.0*D4*R22*R3 +           3.0*D3*R3;
      grad[51]  = D5*R23*R32 +       D4*R2*(R22+3.0*R32) + 3.0*D3*R2;
      grad[52]  = D5*R33*R22 +       D4*R3*(R32+3.0*R22) + 3.0*D3*R3;
      grad[53]  = D5*R34*R2  +   6.0*D4*R32*R2 +           3.0*D3*R2;
      grad[54]  = D5*R33*R32 +  10.0*D4*R33 +             15.0*D3*R3;
      #if expansion_order >= 6 //Computing the sixth gradient
      grad[55]  = D6*R13*R13 +  15.0*D5*R14 +                  45.0*D4*R12 +          15.0*D3;
      grad[56]  = D6*R15*R2 +   10.0*D5*R13*R2 +               15.0*D4*R1_R2;
      grad[57]  = D6*R15*R3 +   10.0*D5*R13*R3 +               15.0*D4*R1_R3;
      grad[58]  = D6*R14*R22 +       D5*R12*(R12+6.0*R22) +     3.0*D4*(R22+2.0*R12) + 3.0*D3;
      grad[59]  = D6*R14*R2_R3 + 6.0*D5*R12*R2_R3 +             3.0*D4*R2_R3;
      grad[60]  = D6*R14*R32 +       D5*R12*(R12+6.0*R32) +     3.0*D4*(R32+2.0*R12) + 3.0*D3;
      grad[61]  = D6*R13*R23 +   3.0*D5*R1_R2*(R12+R22) +       9.0*D4*R1_R2;
      grad[62]  = D6*R13*R22*R3 +    D5*R1_R3*(R12+3.0*R22) +   3.0*D4*R1_R3;
      grad[63]  = D6*R13*R32*R2 +    D5*R1_R2*(R12+3.0*R32) +   3.0*D4*R1_R2;
      grad[64]  = D6*R13*R33 +   3.0*D5*R1_R3*(R12+R32) +       9.0*D4*R1_R3;
      grad[65]  = D6*R24*R12 +       D5*R22*(R22+6.0*R12) +     3.0*D4*(R12+2.0*R22) + 3.0*D3;
      grad[66]  = D6*R23*R12*R3 +    D5*R2_R3*(R22+3.0*R12) +   3.0*D4*R2_R3;
      grad[67]  = D6*R12*R22*R32 +   D5*(R12*R22+R12*R32+R22*R32) + D4*(R12+R22+R32) +     D3;
      grad[68]  = D6*R33*R12*R2 +    D5*R2_R3*(R32+3.0*R12) +   3.0*D4*R2_R3;
      grad[69]  = D6*R34*R12 +       D5*R32*(R32+6.0*R12) +     3.0*D4*(R12+2.0*R32) + 3.0*D3;
      grad[70]  = D6*R25*R1 +   10.0*D5*R23*R1 +               15.0*D4*R1_R2;
      grad[71]  = D6*R24*R1_R3 + 6.0*D5*R22*R1_R3 +             3.0*D4*R1_R3; 
      grad[72]  = D6*R23*R32*R1 +    D5*R1_R2*(R22+3.0*R32) +   3.0*D4*R1_R2;
      grad[73]  = D6*R33*R22*R1 +    D5*R1_R3*(R32+3.0*R22) +   3.0*D4*R1_R3;
      grad[74]  = D6*R34*R1_R2 + 6.0*D5*R32*R1_R2 +             3.0*D4*R1_R2;
      grad[75]  = D6*R35*R1 +   10.0*D5*R33*R1 +               15.0*D4*R1_R3;
      grad[76]  = D6*R23*R23 +  15.0*D5*R24 +                  45.0*D4*R22 +          15.0*D3;
      grad[77]  = D6*R25*R3 +   10.0*D5*R23*R3 +               15.0*D4*R2_R3;
      grad[78]  = D6*R24*R32 +       D5*R22*(R22+6.0*R32) +     3.0*D4*(R32+2.0*R22) + 3.0*D3;
      grad[79]  = D6*R33*R23 +   3.0*D5*R2_R3*(R32+R22) +       9.0*D4*R2_R3;
      grad[80]  = D6*R34*R22 +       D5*R32*(R32+6.0*R22) +     3.0*D4*(R22+2.0*R32) + 3.0*D3;
      grad[81]  = D6*R35*R2 +   10.0*D5*R33*R2 +               15.0*D4*R2_R3;
      grad[82]  = D6*R33*R33 +  15.0*D5*R34 +                  45.0*D4*R32 +          15.0*D3;
      #if expansion_order >= 7 //Computing the seventh gradient
      grad[83]  = D7*R14*R13   + 21.*D6*R15                        + 105.*D5*R13                     + 105.*D4*R1;
      grad[84]  = D7*R16*R2    + 15.*D6*R14*R2                      + 45.*D5*R12*R2                   + 15.*D4*R2;
      grad[85]  = D7*R16*R3    + 15.*D6*R14*R3                      + 45.*D5*R12*R3                   + 15.*D4*R3;
      grad[86]  = D7*R15*R22       + D6*(10.*R13*R22+R15)               + D5*(15.*R1*R22+10.*R13)     + 15.*D4*R1;
      grad[87]  = D7*R15*R2_R3 + 10.*D6*R13*R2_R3                   + 15.*D5*R1_R2*R3;
      grad[88]  = D7*R15*R32       + D6*(10.*R13*R32+R15)               + D5*(15.*R1*R32+10.*R13)     + 15.*D4*R1;
      grad[89]  = D7*R14*R23       + D6*(6.*R12*R23+3.*R14*R2)          + D5*(3.*R23+18.*R12*R2)       + 9.*D4*R2;
      grad[90]  = D7*R14*R22*R3    + D6*(6.*R12*R22*R3+R14*R3)          + D5*(3.*R22*R3+6.*R12*R3)     + 3.*D4*R3;
      grad[91]  = D7*R14*R2*R32    + D6*(6.*R12*R2*R32+R14*R2)          + D5*(3.*R2*R32+6.*R12*R2)     + 3.*D4*R2;
      grad[92]  = D7*R14*R33       + D6*(6.*R12*R33+3.*R14*R3)          + D5*(3.*R33+18.*R12*R3)       + 9.*D4*R3;
      grad[93]  = D7*R13*R24       + D6*(3.*R1*R24+6.*R13*R22)          + D5*(18.*R1*R22+3.*R13)       + 9.*D4*R1;
      grad[94]  = D7*R13*R23*R3 + 3.*D6*(R1_R3*R23+R13*R2_R3)        + 9.*D5*R1_R2*R3;
      grad[95]  = D7*R13*R22*R32   + D6*(3.*R1*R22*R32+R13*R32+R13*R22) + D5*(3.*R1*R32+3.*R1*R22+R13) + 3.*D4*R1;
      grad[96]  = D7*R13*R2*R33 + 3.*D6*(R1_R2*R33+R13*R2_R3)        + 9.*D5*R1_R2*R3;
      grad[97]  = D7*R13*R34       + D6*(3.*R1*R34+6.*R13*R32)          + D5*(18.*R1*R32+3.*R13)       + 9.*D4*R1;
      grad[98]  = D7*R12*R25       + D6*(R25+10.*R12*R23)               + D5*(10.*R23+15.*R12*R2)     + 15.*D4*R2;
      grad[99]  = D7*R12*R24*R3    + D6*(R24*R3+6.*R12*R22*R3)          + D5*(6.*R22*R3+3.*R12*R3)     + 3.*D4*R3;
      grad[100] = D7*R12*R23*R32   + D6*(R23*R32+3.*R12*R2*R32+R12*R23) + D5*(3.*R2*R32+R23+3.*R12*R2) + 3.*D4*R2;
      grad[101] = D7*R12*R22*R33   + D6*(R22*R33+R12*R33+3.*R12*R22*R3) + D5*(R33+3.*R22*R3+3.*R12*R3) + 3.*D4*R3;
      grad[102] = D7*R12*R2*R34    + D6*(R2*R34+6.*R12*R2*R32)          + D5*(6.*R2*R32+3.*R12*R2)     + 3.*D4*R2;
      grad[103] = D7*R12*R35       + D6*(R35+10.*R12*R33)               + D5*(10.*R33+15.*R12*R3)     + 15.*D4*R3;
      grad[104] = D7*R1*R26    + 15.*D6*R1*R24                      + 45.*D5*R1*R22                   + 15.*D4*R1;
      grad[105] = D7*R1_R3*R25 + 10.*D6*R1_R3*R23                   + 15.*D5*R1_R2*R3;
      grad[106] = D7*R1*R24*R32    + D6*(6.*R1*R22*R32+R1*R24)          + D5*(3.*R1*R32+6.*R1*R22)     + 3.*D4*R1;
      grad[107] = D7*R1*R23*R33    + D6*(3.*R1_R2*R33+3.*R1_R3*R23)  + 9.*D5*R1_R2*R3;
      grad[108] = D7*R1*R22*R34    + D6*(R1*R34+6.*R1*R22*R32)          + D5*(6.*R1*R32+3.*R1*R22)     + 3.*D4*R1;
      grad[109] = D7*R1_R2*R35 + 10.*D6*R1_R2*R33                   + 15.*D5*R1_R2*R3;
      grad[110] = D7*R1*R36    + 15.*D6*R1*R34                      + 45.*D5*R1*R32                   + 15.*D4*R1;
      grad[111] = D7*R24*R23   + 21.*D6*R25                        + 105.*D5*R23                     + 105.*D4*R2;
      grad[112] = D7*R26*R3    + 15.*D6*R24*R3                      + 45.*D5*R22*R3                   + 15.*D4*R3;
      grad[113] = D7*R25*R32       + D6*(10.*R23*R32+R25)               + D5*(15.*R2*R32+10.*R23)     + 15.*D4*R2;
      grad[114] = D7*R24*R33       + D6*(6.*R22*R33+3.*R24*R3)          + D5*(3.*R33+18.*R22*R3)       + 9.*D4*R3;
      grad[115] = D7*R23*R34       + D6*(3.*R2*R34+6.*R23*R32)          + D5*(18.*R2*R32+3.*R23)       + 9.*D4*R2;
      grad[116] = D7*R22*R35       + D6*(R35+10.*R22*R33)               + D5*(10.*R33+15.*R22*R3)     + 15.*D4*R3;
      grad[117] = D7*R2*R36    + 15.*D6*R2*R34                      + 45.*D5*R2*R32                   + 15.*D4*R2;
      grad[118] = D7*R34*R33   + 21.*D6*R35                        + 105.*D5*R33                     + 105.*D4*R3;
      #if expansion_order >= 8 //Computing the eighth gradient
      grad[119] = D8*R14*R14   +28.*D7*R16                              +210.*D6*R14                               +420.*D5*R12                  +105.*D4;
      grad[120] = D8*R17*R2    +21.*D7*R15*R2                           +105.*D6*R13*R2                            +105.*D5*R1_R2;
      grad[121] = D8*R17*R3    +21.*D7*R15*R3                           +105.*D6*R13*R3                            +105.*D5*R1_R3;
      grad[122] = D8*R16*R22       +D7*(15.*R14*R22+R16)                     +D6*(45.*R12*R22+15.*R14)                  +D5*(15.*R22+45.*R12)     +15.*D4;
      grad[123] = D8*R16*R2_R3 +15.*D7*R14*R2_R3                         +45.*D6*R12*R2_R3                          +15.*D5*R2_R3;
      grad[124] = D8*R16*R32       +D7*(15.*R14*R32+R16)                     +D6*(45.*R12*R32+15.*R14)                  +D5*(15.*R32+45.*R12)     +15.*D4;
      grad[125] = D8*R15*R23       +D7*(10.*R13*R23+3.*R15*R2)               +D6*(15.*R1*R23+30.*R13*R2)            +45.*D5*R1_R2;
      grad[126] = D8*R15*R22*R3    +D7*(10.*R13*R22*R3+R15*R3)               +D6*(15.*R1*R22*R3+10.*R13*R3)         +15.*D5*R1_R3;
      grad[127] = D8*R15*R2*R32    +D7*(10.*R13*R2*R32+R15*R2)               +D6*(15.*R1*R2*R32+10.*R13*R2)         +15.*D5*R1_R2;
      grad[128] = D8*R15*R33       +D7*(10.*R13*R33+3.*R15*R3)               +D6*(15.*R1*R33+30.*R13*R3)            +45.*D5*R1_R3;
      grad[129] = D8*R14*R24    +6.*D7*(R12*R24+R14*R22)                  +3.*D6*(R24+12.*R12*R22+R14)              +18.*D5*(R22+R12)              +9.*D4;
      grad[130] = D8*R14*R23*R3    +D7*(6.*R12*R23*R3+3.*R14*R2_R3)          +D6*(3.*R23*R3+18.*R12*R2_R3)           +9.*D5*R2_R3;
      grad[131] = D8*R14*R22*R32   +D7*(6.*R12*R22*R32+R14*R32+R14*R22)      +D6*(3.*R22*R32+6.*R12*R32+6.*R12*R22+R14) +D5*(3.*R32+3.*R22+6.*R12) +3.*D4;
      grad[132] = D8*R14*R2*R33    +D7*(6.*R12*R2*R33+3.*R14*R2_R3)          +D6*(3.*R2*R33+18.*R12*R2_R3)           +9.*D5*R2_R3;
      grad[133] = D8*R14*R34       +D7*(6.*R12*R34+6.*R14*R32)               +D6*(3.*R34+36.*R12*R32+3.*R14)        +18.*D5*(R32+R12)              +9.*D4;
      grad[134] = D8*R13*R25       +D7*(3.*R1*R25+10.*R13*R23)               +D6*(30.*R1*R23+15.*R13*R2)            +45.*D5*R1_R2;
      grad[135] = D8*R13*R24*R3    +D7*(3.*R1_R3*R24+6.*R13*R22*R3)          +D6*(18.*R1_R3*R22+3.*R13*R3)           +9.*D5*R1_R3;
      grad[136] = D8*R13*R23*R32   +D7*(3.*R1*R23*R32+3.*R13*R2*R32+R13*R23) +D6*(9.*R1_R2*R32+3.*R1*R23+3.*R13*R2)  +9.*D5*R1_R2;
      grad[137] = D8*R13*R22*R33   +D7*(3.*R1*R22*R33+R13*R33+3.*R13*R22*R3) +D6*(3.*R1*R33+9.*R1_R3*R22+3.*R13*R3)  +9.*D5*R1_R3;
      grad[138] = D8*R13*R2*R34    +D7*(3.*R1_R2*R34+6.*R13*R2*R32)          +D6*(18.*R1_R2*R32+3.*R13*R2)           +9.*D5*R1_R2;
      grad[139] = D8*R13*R35       +D7*(3.*R1*R35+10.*R13*R33)               +D6*(30.*R1*R33+15.*R13*R3)            +45.*D5*R1_R3;
      grad[140] = D8*R12*R26       +D7*(R26+15.*R12*R24)                     +D6*(15.*R24+45.*R12*R22)                  +D5*(45.*R22+15.*R12)     +15.*D4;
      grad[141] = D8*R12*R25*R3    +D7*(R25*R3+10.*R12*R23*R3)               +D6*(10.*R23*R3+15.*R12*R2_R3)         +15.*D5*R2_R3;
      grad[142] = D8*R12*R24*R32   +D7*(R24*R32+6.*R12*R22*R32+R12*R24)      +D6*(6.*R22*R32+3.*R12*R32+R24+6.*R12*R22) +D5*(3.*R32+6.*R22+3.*R12) +3.*D4;
      grad[143] = D8*R12*R23*R33   +D7*(R23*R33+3.*R12*R2*R33+3.*R12*R23*R3) +D6*(3.*R2*R33+3.*R23*R3+9.*R12*R2_R3)  +9.*D5*R2_R3;
      grad[144] = D8*R12*R22*R34   +D7*(R22*R34+R12*R34+6.*R12*R22*R32)      +D6*(R34+6.*R22*R32+6.*R12*R32+3.*R12*R22) +D5*(6.*R32+3.*R22+3.*R12) +3.*D4;
      grad[145] = D8*R12*R2*R35    +D7*(R2*R35+10.*R12*R2*R33)               +D6*(10.*R2*R33+15.*R12*R2_R3)         +15.*D5*R2_R3;
      grad[146] = D8*R12*R36       +D7*(R36+15.*R12*R34)                     +D6*(15.*R34+45.*R12*R32)                  +D5*(45.*R32+15.*R12)     +15.*D4;
      grad[147] = D8*R1*R27    +21.*D7*R1*R25                           +105.*D6*R1*R23                            +105.*D5*R1_R2;
      grad[148] = D8*R1_R3*R26 +15.*D7*R1_R3*R24                         +45.*D6*R1*R22*R3                          +15.*D5*R1_R3;
      grad[149] = D8*R1*R25*R32    +D7*(10.*R1*R23*R32+R1*R25)               +D6*(15.*R1_R2*R32+10.*R1*R23)         +15.*D5*R1_R2;
      grad[150] = D8*R1*R24*R33    +D7*(6.*R1*R22*R33+3.*R1_R3*R24)          +D6*(3.*R1*R33+18.*R1_R3*R22)           +9.*D5*R1_R3;
      grad[151] = D8*R1*R23*R34    +D7*(3.*R1_R2*R34+6.*R1*R23*R32)          +D6*(18.*R1_R2*R32+3.*R1*R23)           +9.*D5*R1_R2;
      grad[152] = D8*R1*R22*R35    +D7*(R1*R35+10.*R1*R22*R33)               +D6*(10.*R1*R33+15.*R1_R3*R22)         +15.*D5*R1_R3;
      grad[153] = D8*R1_R2*R36 +15.*D7*R1_R2*R34                         +45.*D6*R1_R2*R32                          +15.*D5*R1_R2;
      grad[154] = D8*R1*R37    +21.*D7*R1*R35                           +105.*D6*R1*R33                            +105.*D5*R1_R3;
      grad[155] = D8*R24*R24   +28.*D7*R26                              +210.*D6*R24                               +420.*D5*R22                  +105.*D4;
      grad[156] = D8*R27*R3    +21.*D7*R25*R3                           +105.*D6*R23*R3                            +105.*D5*R2_R3;
      grad[157] = D8*R26*R32       +D7*(15.*R24*R32+R26)                     +D6*(45.*R22*R32+15.*R24)                  +D5*(15.*R32+45.*R22)     +15.*D4;
      grad[158] = D8*R25*R33       +D7*(10.*R23*R33+3.*R25*R3)               +D6*(15.*R2*R33+30.*R23*R3)            +45.*D5*R2_R3;
      grad[159] = D8*R24*R34    +6.*D7*(R22*R34+R24*R32)                  +3.*D6*(R34+12.*R22*R32+R24)              +18.*D5*(R32+R22)              +9.*D4;
      grad[160] = D8*R23*R35       +D7*(3.*R2*R35+10.*R23*R33)               +D6*(30.*R2*R33+15.*R23*R3)            +45.*D5*R2_R3;
      grad[161] = D8*R22*R36       +D7*(R36+15.*R22*R34)                     +D6*(15.*R34+45.*R22*R32)                  +D5*(45.*R32+15.*R22)     +15.*D4;
      grad[162] = D8*R2*R37    +21.*D7*R2*R35                           +105.*D6*R2*R33                            +105.*D5*R2_R3;
      grad[163] = D8*R34*R34   +28.*D7*R36                              +210.*D6*R34                               +420.*D5*R32                  +105.*D4;
      #endif
      #endif
      #endif
      #endif
      #endif
      #endif
      #endif
}


void inner_product(typ * T1, typ * T2, typ * T3, int p, int q, typ factor){

      /******** Computes the inner product between the tensor T1 of order p + q (max 6) and the tensor T2 ********/
      /******** of order p. The result is accumulated into the tensor T3 of order q with a factor factor. ********/
      /******** T1, T2 and T3 are assumed to be arrays of respective sizes (p+q+1)(p+q+2)/2, (p+1)(p+2)/2 ********/
      /******** and (q+1)(q+2)/2. The tensor T3 is not overwritten, but accumulated instead               ********/

      int q1, q2, q3; //Indexes into the arrays T1, T2 and T3
      int T2_size = ((p + 1)*(p + 2))/2;
      int T3_size = ((q + 1)*(q + 2))/2;
      int perm; //Number of tensor index permutations
      
      for (q3 = 0; q3 < T3_size; q3 ++){ //Looping over the components of the tensor to be computed
            for (q2 = 0; q2 < T2_size; q2 ++){   //Looping over the components of T2, in order to accumulate T3[q3]
                  perm    = perm_from_kn[q2][p]; //Number of permutations that can be done with the indexes of T2
                  q1      = q1fromq2q3[q2][q3];  //Retrieving the array index of T1
                  T3[q3] += T1[q1] * T2[q2] * (typ) perm * factor; //Accumulating into tensor T3 with a factor factor
            }
      }
}


void double_IP(typ * T1, typ * T2, typ * TT2, typ * T3, typ * TT3, int p, int q, typ factor1, typ factor2){

      /******** Same as above but computes two inner products at the same time ********/
      /******** Accumulates T1 X T2 into T3 and T1 X TT2 into TT3              ********/

      int q1, q2, q3; //Indexes into the arrays T1, T2 and T3
      int T2_size = ((p + 1)*(p + 2))/2;
      int T3_size = ((q + 1)*(q + 2))/2;
      int perm; //Number of tensor index permutations
      
      for (q3 = 0; q3 < T3_size; q3 ++){ //Looping over the components of the tensors to be computed
            for (q2 = 0; q2 < T2_size; q2 ++){    //Looping over the components of T2, in order to accumulate T3[q3]
                  perm     = perm_from_kn[q2][p]; //Number of permutations that can be done with the indexes of T2
                  q1       = q1fromq2q3[q2][q3];  //Retrieving the array index of T1
                  T3[q3]  += T1[q1] * T2[q2]  * (typ) perm * factor1; //Accumulating into tensor T3  with a factor factor1
                  TT3[q3] += T1[q1] * TT2[q2] * (typ) perm * factor2; //Accumulating into tensor TT3 with a factor factor2
            }
      }
}


/***************************************************************************/
/******** I now write functions relative two the second stage of    ********/
/******** Dehnen's algorithm, accumulation of the C^(m) in the tree ********/
/***************************************************************************/


#if mutual_bool
void Cm_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Performs the tree walk of Dehnen's algorithm, called the interaction phase      ********/
      /******** This function computes the tensors C^(m) between interacting cells, or the      ********/
      /******** acceleration C^(1) in case of interaction between bodies, but does not pass     ********/
      /******** the C^(m) down the tree. It is assumed that the array C1Moonlets, whose indexes ********/
      /******** 3*a to 3*a + 2 contain the acceleration of body a, is given initialized to 0.0  ********/
      /******** Similarly, it is assumed that the C^(m) are 0.0 upon calling this function      ********/
      
      /******** It is hard to tell in advance how many pairs will be treated by the tree walk   ********/
      /******** I expect it will be at most factor * cell_id, but that might have to be changed ********/
      int factor = (int) floor(350.0 * 0.5 / theta_min);
      struct pair * stack = (struct pair *)malloc(factor * cell_id * sizeof(struct pair)); //Stack of pairs of nodes' id that have to be treated
      if (stack == NULL){
            fprintf(stderr, "Error : Cannot allocate memory for the stack in function Cm_flattree.\n");
            abort();
      }
      
      int i = 0; //Index in the stack of the current pair of ids
      int j = 0; //Index of where to put a pair in the stack
      int a, b;  //Ids of the nodes of the current pair
      int p, q;  //Loop indexes
      int s, u;  //Body indexes
      int * dots_a; //All the bodies in cell a
      int * dots_b; //All the bodies in cell b
      int Na, Nb; //Number of bodies in cells a and b
      int how_many_children_a, how_many_children_b; //Number of children of cells a and b
      int idFirstChild; //Id of first child
      int idLastChild;  //Id of last child
      typ M0_a, M0_b;   //Nodes' masses
      typ * com_a, * com_b; //Nodes' centers of mass
      #if expansion_order   == 1
      typ Grad[3]; //Independant components of all the gradients of G/R from order 1 to expansion_order
      #elif expansion_order == 2
      typ Grad[9];
      #elif expansion_order == 3
      typ Grad[19];
      #elif expansion_order == 4
      typ Grad[34];
      #elif expansion_order == 5
      typ Grad[55];
      #elif expansion_order == 6
      typ Grad[83];
      #elif expansion_order == 7
      typ Grad[119];
      #elif expansion_order == 8
      typ Grad[164];
      #endif
      typ C1_a[3];
      typ C1_b[3];
      typ R[3];  //com_a - com_b
      
      #if expansion_order >= 2
            typ * grad2;
            typ C2_a[6];
            typ C2_b[6];
            #if expansion_order >= 3
                  typ * grad3;
                  typ C3_a[10];       //interaction tensor of node a
                  typ C3_b[10];       //interaction tensor of node b
                  typ * M2_a, * M2_b; //Node's multipole moments
                  #if expansion_order >= 4
                        typ * grad4;
                        typ C4_a[15];
                        typ C4_b[15];
                        typ * M3_a, * M3_b;
                        #if expansion_order >= 5
                              typ * grad5;
                              typ C5_a[21];
                              typ C5_b[21];
                              typ * M4_a, * M4_b;
                              #if expansion_order >= 6
                                    typ * grad6;
                                    typ C6_a[28];
                                    typ C6_b[28];
                                    typ * M5_a, * M5_b;
                                    #if expansion_order >= 7
                                          typ * grad7;
                                          typ C7_a[36];
                                          typ C7_b[36];
                                          typ * M6_a, * M6_b;
                                          #if expansion_order >= 8
                                                typ * grad8;
                                                typ C8_a[45];
                                                typ C8_b[45];
                                                typ * M7_a, * M7_b;
                                          #endif
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      #endif

      typ r_crit_a, r_crit_b; // Critical radii of node a and b
      typ ms, Rs, xs, ys, zs, mu, Ru, xu, yu, zu, dx, dy, dz, r, r3, softening; //Body coordinates
      typ omega2_x, omega2_y, omega2_z;
      
      /******** Putting the pair (root_cell, root_cell) in the stack ********/
      (stack + j) -> fst = 0;
      (stack + j) -> snd = 0;
      j ++;
      
      
      /******** I travel the stack of pairs of nodes. At each pair, if NaNb < N_cc_pre, I treat it brute-forcely,  ********/
      /******** otherwise, if the nodes are well-separated, I treat them by a multipole expansion, otherwise, if   ********/
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
                  if (Na < N_cs || how_many_children_a == 0){ //Direct interaction
                        dots_a = (FlatTree + a) -> dots;
                        for (p = 0; p < Na; p ++){
                              for (q = p + 1; q < Na; q ++){
                                    s = dots_a[p]; //Id of first  moonlet
                                    u = dots_a[q]; //Id of second monnlet 
                                    xs = (moonlets + s) -> x;
                                    ys = (moonlets + s) -> y;
                                    zs = (moonlets + s) -> z;
                                    ms = (moonlets + s) -> mass;
                                    Rs = (moonlets + s) -> radius;
                                    xu = (moonlets + u) -> x;
                                    yu = (moonlets + u) -> y;
                                    zu = (moonlets + u) -> z;
                                    mu = (moonlets + u) -> mass;
                                    Ru = (moonlets + u) -> radius;
                                    dx = xs - xu;
                                    dy = ys - yu;
                                    dz = zs - zu;
                                    softening = softening_parameter*(Rs + Ru);
                                    r  = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                                    r3 = r*r*r;
                                    omega2_x = G*dx/r3;
                                    omega2_y = G*dy/r3;
                                    omega2_z = G*dz/r3;
                                    C1Moonlets[3*s]     -= mu*omega2_x;
                                    C1Moonlets[3*s + 1] -= mu*omega2_y;
                                    C1Moonlets[3*s + 2] -= mu*omega2_z;
                                    C1Moonlets[3*u]     += ms*omega2_x;
                                    C1Moonlets[3*u + 1] += ms*omega2_y;
                                    C1Moonlets[3*u + 2] += ms*omega2_z;
                              }
                        }
                  }
                  else{ //Subdivision
                        idFirstChild = (FlatTree + a) -> idFirstChild;
                        idLastChild  = idFirstChild + how_many_children_a;
                        /******** Making sure that the stack is big enough. To be removed when the code is robust ********/
                        if (j + 36 >= factor * cell_id){
                              fprintf(stderr, "Error : The stack is not big enough in function Cm_flattree. Try increasing the value of factor.\n");
                              abort();
                        }
                        for (p = idFirstChild; p < idLastChild; p ++){
                              for (q = p; q < idLastChild; q ++){
                                    (stack + j) -> fst = p;
                                    (stack + j) -> snd = q;
                                    j ++;
                              }
                        }
                  }
            }
            else{ //Interaction between two different cells
                  if (!(Na > N_cc_pre && Nb > N_cc_pre) && Na*Nb < N_cc_pre){ //Direct interaction
                        dots_a = (FlatTree + a) -> dots;
                        dots_b = (FlatTree + b) -> dots;
                        for (p = 0; p < Na; p ++){
                              for (q = 0; q < Nb; q ++){
                                    s = dots_a[p]; //Id of first  body
                                    u = dots_b[q]; //Id of second body
                                    xs = (moonlets + s) -> x;
                                    ys = (moonlets + s) -> y;
                                    zs = (moonlets + s) -> z;
                                    ms = (moonlets + s) -> mass;
                                    Rs = (moonlets + s) -> radius;
                                    xu = (moonlets + u) -> x;
                                    yu = (moonlets + u) -> y;
                                    zu = (moonlets + u) -> z;
                                    mu = (moonlets + u) -> mass;
                                    Ru = (moonlets + u) -> radius;
                                    dx = xs - xu;
                                    dy = ys - yu;
                                    dz = zs - zu;
                                    softening = softening_parameter*(Rs + Ru);
                                    r  = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                                    r3 = r*r*r;
                                    omega2_x = G*dx/r3;
                                    omega2_y = G*dy/r3;
                                    omega2_z = G*dz/r3;
                                    C1Moonlets[3*s]     -= mu*omega2_x;
                                    C1Moonlets[3*s + 1] -= mu*omega2_y;
                                    C1Moonlets[3*s + 2] -= mu*omega2_z;
                                    C1Moonlets[3*u]     += ms*omega2_x;
                                    C1Moonlets[3*u + 1] += ms*omega2_y;
                                    C1Moonlets[3*u + 2] += ms*omega2_z;
                              }
                        }  
                  }
                  else{
                        r_crit_a = (FlatTree + a) -> r_crit;
                        com_a    = (FlatTree + a) -> com;
                        r_crit_b = (FlatTree + b) -> r_crit;
                        com_b    = (FlatTree + b) -> com;
                        R[0]     = com_a[0] - com_b[0];
                        R[1]     = com_a[1] - com_b[1];
                        R[2]     = com_a[2] - com_b[2];
                        r        = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
                        if (r > r_crit_a + r_crit_b){ //The nodes are well-separated -> Multipole expansion
                              /******** Computing the tensors of interaction C^(m) ********/
                              M0_a = (FlatTree + a) -> M0;
                              M0_b = (FlatTree + b) -> M0;
                              gradR(R, Grad);  //Computing all the gradients of G/R
                              C1_a[0] =  M0_b * Grad[0];  C1_a[1] =  M0_b * Grad[1];  C1_a[2] =  M0_b * Grad[2];
                              C1_b[0] = -M0_a * Grad[0];  C1_b[1] = -M0_a * Grad[1];  C1_b[2] = -M0_a * Grad[2];
                              
                              #if expansion_order >= 2
                                    grad2 = Grad + 3;
                                    for (p = 0; p < 6; p ++){
                                          C2_a[p] = M0_b * grad2[p];
                                          C2_b[p] = M0_a * grad2[p];
                                    }
                                    
                                    #if expansion_order >= 3
                                          M2_a  = (FlatTree + a) -> M2;
                                          M2_b  = (FlatTree + b) -> M2;
                                          grad3 = Grad + 9;
                                          for (p = 0; p < 10; p ++){
                                                C3_a[p] =  M0_b * grad3[p];
                                                C3_b[p] = -M0_a * grad3[p];
                                          }
                                          C1_a[0] += 0.5*(M2_b[0]*grad3[0]+2.0*M2_b[1]*grad3[1]+2.0*M2_b[2]*grad3[2]+M2_b[3]*grad3[3]+2.0*M2_b[4]*grad3[4]+M2_b[5]*grad3[5]);
                                          C1_a[1] += 0.5*(M2_b[0]*grad3[1]+2.0*M2_b[1]*grad3[3]+2.0*M2_b[2]*grad3[4]+M2_b[3]*grad3[6]+2.0*M2_b[4]*grad3[7]+M2_b[5]*grad3[8]);
                                          C1_a[2] += 0.5*(M2_b[0]*grad3[2]+2.0*M2_b[1]*grad3[4]+2.0*M2_b[2]*grad3[5]+M2_b[3]*grad3[7]+2.0*M2_b[4]*grad3[8]+M2_b[5]*grad3[9]);
                                          C1_b[0] -= 0.5*(M2_a[0]*grad3[0]+2.0*M2_a[1]*grad3[1]+2.0*M2_a[2]*grad3[2]+M2_a[3]*grad3[3]+2.0*M2_a[4]*grad3[4]+M2_a[5]*grad3[5]);
                                          C1_b[1] -= 0.5*(M2_a[0]*grad3[1]+2.0*M2_a[1]*grad3[3]+2.0*M2_a[2]*grad3[4]+M2_a[3]*grad3[6]+2.0*M2_a[4]*grad3[7]+M2_a[5]*grad3[8]);
                                          C1_b[2] -= 0.5*(M2_a[0]*grad3[2]+2.0*M2_a[1]*grad3[4]+2.0*M2_a[2]*grad3[5]+M2_a[3]*grad3[7]+2.0*M2_a[4]*grad3[8]+M2_a[5]*grad3[9]);
                                          
                                          #if expansion_order >= 4
                                                M3_a  = (FlatTree + a) -> M3;
                                                M3_b  = (FlatTree + b) -> M3;
                                                grad4 = Grad + 19;
                                                for (p = 0; p < 15; p ++){
                                                      C4_a[p] = M0_b * grad4[p];
                                                      C4_b[p] = M0_a * grad4[p];
                                                }
                                                double_IP(grad4, M3_b, M3_a, C1_a, C1_b, 3, 1, -0.166666666666666666666666667, -0.166666666666666666666666667);
                                                double_IP(grad4, M2_b, M2_a, C2_a, C2_b, 2, 2, 0.5, 0.5);
                                                
                                                #if expansion_order >= 5
                                                      M4_a  = (FlatTree + a) -> M4;
                                                      M4_b  = (FlatTree + b) -> M4;
                                                      grad5 = Grad + 34;
                                                      for (p = 0; p < 21; p ++){
                                                            C5_a[p] =  M0_b * grad5[p];
                                                            C5_b[p] = -M0_a * grad5[p];
                                                      }
                                                      double_IP(grad5, M4_b, M4_a, C1_a, C1_b, 4, 1, 0.041666666666666666666666667, -0.041666666666666666666666667);
                                                      double_IP(grad5, M3_b, M3_a, C2_a, C2_b, 3, 2, -0.166666666666666666666666667, 0.166666666666666666666666667);
                                                      double_IP(grad5, M2_b, M2_a, C3_a, C3_b, 2, 3, 0.5, -0.5);
                                                      
                                                      #if expansion_order >= 6
                                                            M5_a  = (FlatTree + a) -> M5;
                                                            M5_b  = (FlatTree + b) -> M5;
                                                            grad6 = Grad + 55;
                                                            for (p = 0; p < 28; p ++){
                                                                  C6_a[p] = M0_b * grad6[p];
                                                                  C6_b[p] = M0_a * grad6[p];
                                                            }
                                                            double_IP(grad6, M5_b, M5_a, C1_a, C1_b, 5, 1, -0.008333333333333333333333333, -0.008333333333333333333333333);
                                                            double_IP(grad6, M4_b, M4_a, C2_a, C2_b, 4, 2, 0.041666666666666666666666667, 0.041666666666666666666666667);
                                                            double_IP(grad6, M3_b, M3_a, C3_a, C3_b, 3, 3, -0.166666666666666666666666667, -0.166666666666666666666666667);
                                                            double_IP(grad6, M2_b, M2_a, C4_a, C4_b, 2, 4, 0.5, 0.5);
                                                            #if expansion_order >= 7
                                                                  M6_a  = (FlatTree + a) -> M6;
                                                                  M6_b  = (FlatTree + b) -> M6;
                                                                  grad7 = Grad + 83;
                                                                  for (p = 0; p < 36; p ++){
                                                                        C7_a[p] =  M0_b * grad7[p];
                                                                        C7_b[p] = -M0_a * grad7[p];
                                                                  }
                                                                  double_IP(grad7, M6_b, M6_a, C1_a, C1_b, 6, 1, 0.0013888888888888888888889, -0.0013888888888888888888889);
                                                                  double_IP(grad7, M5_b, M5_a, C2_a, C2_b, 5, 2, -0.0083333333333333333333333, 0.0083333333333333333333333);
                                                                  double_IP(grad7, M4_b, M4_a, C3_a, C3_b, 4, 3, 0.0416666666666666666666667, -0.0416666666666666666666667);
                                                                  double_IP(grad7, M3_b, M3_a, C4_a, C4_b, 3, 4, -0.1666666666666666666666667, 0.1666666666666666666666667);
                                                                  double_IP(grad7, M2_b, M2_a, C5_a, C5_b, 2, 5, 0.5, -0.5);
                                                                  #if expansion_order >= 8
                                                                        M7_a  = (FlatTree + a) -> M7;
                                                                        M7_b  = (FlatTree + b) -> M7;
                                                                        grad8 = Grad + 119;
                                                                        for (p = 0; p < 45; p ++){
                                                                              C8_a[p] = M0_b * grad8[p];
                                                                              C8_b[p] = M0_a * grad8[p];
                                                                        }
                                                                        double_IP(grad8, M7_b, M7_a, C1_a, C1_b, 7, 1, -0.0001984126984126984126984, -0.0001984126984126984126984);
                                                                        double_IP(grad8, M6_b, M6_a, C2_a, C2_b, 6, 2, 0.0013888888888888888888889, 0.0013888888888888888888889);
                                                                        double_IP(grad8, M5_b, M5_a, C3_a, C3_b, 5, 3, -0.0083333333333333333333333, -0.0083333333333333333333333);
                                                                        double_IP(grad8, M4_b, M4_a, C4_a, C4_b, 4, 4, 0.0416666666666666666666667, 0.0416666666666666666666667);
                                                                        double_IP(grad8, M3_b, M3_a, C5_a, C5_b, 3, 5, -0.1666666666666666666666667, -0.1666666666666666666666667);
                                                                        double_IP(grad8, M2_b, M2_a, C6_a, C6_b, 2, 6, 0.5, 0.5);
                                                                  #endif
                                                            #endif
                                                      #endif
                                                #endif
                                          #endif
                                    #endif
                              #endif
                              /******** Actualizing the FlatTrees C^(m) ********/
                              ((FlatTree + a) -> C1)[0] += C1_a[0];  ((FlatTree + a) -> C1)[1] += C1_a[1];  ((FlatTree + a) -> C1)[2] += C1_a[2];
                              ((FlatTree + b) -> C1)[0] += C1_b[0];  ((FlatTree + b) -> C1)[1] += C1_b[1];  ((FlatTree + b) -> C1)[2] += C1_b[2];
                              #if expansion_order >= 2
                                    for (p = 0; p < 6; p ++){
                                          ((FlatTree + a) -> C2)[p] += C2_a[p];  ((FlatTree + b) -> C2)[p] += C2_b[p];
                                    }
                                    #if expansion_order >= 3
                                          for (p = 0; p < 10; p ++){
                                                ((FlatTree + a) -> C3)[p] += C3_a[p];  ((FlatTree + b) -> C3)[p] += C3_b[p];
                                          }
                                          #if expansion_order >= 4
                                                for (p = 0; p < 15; p ++){
                                                      ((FlatTree + a) -> C4)[p] += C4_a[p];  ((FlatTree + b) -> C4)[p] += C4_b[p];
                                                }
                                                #if expansion_order >= 5
                                                      for (p = 0; p < 21; p ++){
                                                            ((FlatTree + a) -> C5)[p] += C5_a[p];  ((FlatTree + b) -> C5)[p] += C5_b[p];
                                                      }
                                                      #if expansion_order >= 6
                                                            for (p = 0; p < 28; p ++){
                                                                  ((FlatTree + a) -> C6)[p] += C6_a[p];  ((FlatTree + b) -> C6)[p] += C6_b[p];
                                                            }
                                                            #if expansion_order >= 7
                                                                  for (p = 0; p < 36; p ++){
                                                                        ((FlatTree + a) -> C7)[p] += C7_a[p];  ((FlatTree + b) -> C7)[p] += C7_b[p];
                                                                  }
                                                                  #if expansion_order >= 8
                                                                        for (p = 0; p < 45; p ++){
                                                                              ((FlatTree + a) -> C8)[p] += C8_a[p];  ((FlatTree + b) -> C8)[p] += C8_b[p];
                                                                        }
                                                                  #endif
                                                            #endif
                                                      #endif
                                                #endif
                                          #endif
                                    #endif
                              #endif
                        }
                        else{ //The nodes are not well-separated
                              how_many_children_a = (FlatTree + a) -> how_many_children;
                              how_many_children_b = (FlatTree + b) -> how_many_children;
                              if ((!(Na > N_cc_post && Nb > N_cc_post) && Na*Nb < N_cc_post) || (how_many_children_a == 0 && how_many_children_b == 0)){ //Direct interaction
                                    dots_a = (FlatTree + a) -> dots;
                                    dots_b = (FlatTree + b) -> dots;
                                    for (p = 0; p < Na; p ++){
                                          for (q = 0; q < Nb; q ++){
                                                s  = dots_a[p]; //Id of first  body
                                                u  = dots_b[q]; //Id of second body
                                                xs = (moonlets + s) -> x;
                                                ys = (moonlets + s) -> y;
                                                zs = (moonlets + s) -> z;
                                                ms = (moonlets + s) -> mass;
                                                Rs = (moonlets + s) -> radius;
                                                xu = (moonlets + u) -> x;
                                                yu = (moonlets + u) -> y;
                                                zu = (moonlets + u) -> z;
                                                mu = (moonlets + u) -> mass;
                                                Ru = (moonlets + u) -> radius;
                                                dx = xs - xu;
                                                dy = ys - yu;
                                                dz = zs - zu;
                                                softening = softening_parameter*(Rs + Ru);
                                                r  = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                                                r3 = r*r*r;
                                                omega2_x = G*dx/r3;
                                                omega2_y = G*dy/r3;
                                                omega2_z = G*dz/r3;
                                                C1Moonlets[3*s]     -= mu*omega2_x;
                                                C1Moonlets[3*s + 1] -= mu*omega2_y;
                                                C1Moonlets[3*s + 2] -= mu*omega2_z;
                                                C1Moonlets[3*u]     += ms*omega2_x;
                                                C1Moonlets[3*u + 1] += ms*omega2_y;
                                                C1Moonlets[3*u + 2] += ms*omega2_z;
                                          }
                                    }
                              }
                              else{ //Subdivision
                                    /******** Making sure that the stack is big enough. To be removed when the code is robust ********/
                                    if (j + 8 >= factor * cell_id){
                                          fprintf(stderr, "Error : The stack is not big enough in function Cm_flattree. Try increasing the value of factor.\n");
                                          abort();
                                    }
                                    if (how_many_children_a == 0){ //Subdivising b
                                          idFirstChild = (FlatTree + b) -> idFirstChild;
                                          idLastChild  = idFirstChild + how_many_children_b;
                                          for (p = idFirstChild; p < idLastChild; p ++){
                                                (stack + j) -> fst = a;
                                                (stack + j) -> snd = p;
                                                j ++;
                                          }
                                    }
                                    else if (how_many_children_b == 0){ //Subdivising a
                                          idFirstChild = (FlatTree + a) -> idFirstChild;
                                          idLastChild  = idFirstChild + how_many_children_a;
                                          for (p = idFirstChild; p < idLastChild; p ++){
                                                (stack + j) -> fst = b;
                                                (stack + j) -> snd = p;
                                                j ++;
                                          }
                                    }
                                    else{ //Both cells have children, subdivising the largest (in term of critical radius)
                                          if (r_crit_a < r_crit_b){ //Subdivising b
                                                idFirstChild = (FlatTree + b) -> idFirstChild;
                                                idLastChild  = idFirstChild + how_many_children_b;
                                                for (p = idFirstChild; p < idLastChild; p ++){
                                                      (stack + j) -> fst = a;
                                                      (stack + j) -> snd = p;
                                                      j ++;
                                                }
                                          }
                                          else{ //Subdivising a
                                                idFirstChild = (FlatTree + a) -> idFirstChild;
                                                idLastChild  = idFirstChild + how_many_children_a;
                                                for (p = idFirstChild; p < idLastChild; p ++){
                                                      (stack + j) -> fst = b;
                                                      (stack + j) -> snd = p;
                                                      j ++;
                                                }
                                          }
                                    }
                              }
                        }
                  }
            }
            i ++;
      }
      free(stack);
      stack = NULL;
}
#endif


/***********************************************************************/
/******** I now write functions relative two the third stage of ********/
/******** Dehnen's algorithm, passing the C^(n) down the tree   ********/
/***********************************************************************/


#if mutual_bool
void Cm_downtree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Performs the third stage of Dehnen's algorithm, called the evaluation phase, ********/
      /******** where the C^(m) are passed down the three and the accelerations C^(1) of     ********/
      /******** each body computed. It is assumed that this function follows Cm_flattree     ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of pairs of ids of nodes that have to be treated
      int i = 0; //Index in the stack of the current id
      int j = 0; //Index of where to put an id in the stack
      int a;     //Id of the current node
      int p;     //Loop index
      int b;     //Body index
      int * dots;//All the bodies in the current cell
      int Na;    //Number of bodies in cells a and b
      int how_many_children; //Number of children of the current node
      int idFirstChild; //Id of first child 
      int idLastChild;  //Id of last child
      int idParent;
      typ * C1, * C1p; //interaction tensor of the current node and of its parent
      
      #if expansion_order >= 2
            typ * com, * com_parent; //centers of mass of the cell and of its parent
            typ * C2, * C2p;
            typ R[3];                //com - com_parent
            #if expansion_order >= 3
                  typ * C3, * C3p;
                  typ X2[6];
                  #if expansion_order >= 4
                        int q;  //Loop index
                        typ * C4, * C4p;
                        typ X3[10];
                        #if expansion_order >= 5
                              typ * C5, * C5p;
                              typ X4[15];
                              #if expansion_order >= 6
                                    typ * C6, * C6p;
                                    typ X5[21];
                                    #if expansion_order >= 7
                                          typ * C7, * C7p;
                                          typ X6[28];
                                          #if expansion_order >= 8
                                                typ * C8, * C8p;
                                                typ X7[36];
                                          #endif
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      #endif

      /******** Putting the root cell in the stack ********/
      stack[j] = 0;
      j ++;

      while (j > i){
            a = stack[i];

            /******** If the node has a parent, I shift and accumulate its C^(n) ********/
            idParent = (FlatTree + a) -> idParent;
            if (idParent != -1){
                  C1  = (FlatTree + a)        -> C1;
                  C1p = (FlatTree + idParent) -> C1;
                  C1[0] += C1p[0]; C1[1] += C1p[1]; C1[2] += C1p[2];
                  
                  #if expansion_order >= 2
                        com        = (FlatTree + a) -> com;
                        com_parent = (FlatTree + idParent) -> com;
                        R[0] = com[0] - com_parent[0];
                        R[1] = com[1] - com_parent[1];
                        R[2] = com[2] - com_parent[2];
                        C2   = (FlatTree + a)        -> C2;
                        C2p  = (FlatTree + idParent) -> C2;
                        C2[0] += C2p[0]; C2[1] += C2p[1]; C2[2] += C2p[2]; C2[3] += C2p[3]; C2[4] += C2p[4]; C2[5] += C2p[5];
                        C1[0] += R[0]*C2p[0] + R[1]*C2p[1] + R[2]*C2p[2];
                        C1[1] += R[0]*C2p[1] + R[1]*C2p[3] + R[2]*C2p[4];
                        C1[2] += R[0]*C2p[2] + R[1]*C2p[4] + R[2]*C2p[5];

                        #if expansion_order >= 3
                              C3  = (FlatTree + a)        -> C3;
                              C3p = (FlatTree + idParent) -> C3;
                              for (p = 0; p < 10; p ++){
                                    C3[p] += C3p[p];
                              }
                              X2[0]  = R[0]*R[0]; X2[1] = R[0]*R[1]; X2[2] = R[0]*R[2]; X2[3] = R[1]*R[1]; X2[4] = R[1]*R[2]; X2[5] = R[2]*R[2];
                              C1[0] += 0.5*(X2[0]*C3p[0] + 2.0*X2[1]*C3p[1] + 2.0*X2[2]*C3p[2] + X2[3]*C3p[3] + 2.0*X2[4]*C3p[4] + X2[5]*C3p[5]);
                              C1[1] += 0.5*(X2[0]*C3p[1] + 2.0*X2[1]*C3p[3] + 2.0*X2[2]*C3p[4] + X2[3]*C3p[6] + 2.0*X2[4]*C3p[7] + X2[5]*C3p[8]);
                              C1[2] += 0.5*(X2[0]*C3p[2] + 2.0*X2[1]*C3p[4] + 2.0*X2[2]*C3p[5] + X2[3]*C3p[7] + 2.0*X2[4]*C3p[8] + X2[5]*C3p[9]);
                              C2[0] += R[0]*C3p[0] + R[1]*C3p[1] + R[2]*C3p[2];
                              C2[1] += R[0]*C3p[1] + R[1]*C3p[3] + R[2]*C3p[4];
                              C2[2] += R[0]*C3p[2] + R[1]*C3p[4] + R[2]*C3p[5];
                              C2[3] += R[0]*C3p[3] + R[1]*C3p[6] + R[2]*C3p[7];
                              C2[4] += R[0]*C3p[4] + R[1]*C3p[7] + R[2]*C3p[8];
                              C2[5] += R[0]*C3p[5] + R[1]*C3p[8] + R[2]*C3p[9];

                              #if expansion_order >= 4
                                    C4  = (FlatTree + a)        -> C4;
                                    C4p = (FlatTree + idParent) -> C4;
                                    for (p = 0; p < 15; p ++){
                                          C4[p] += C4p[p];
                                    }
                                    for (p = 0; p < 10; p ++){
                                          get_Xn_overwrite(p, 3, R, 1.0, X3);
                                    }
                                    inner_product(C4p, X3, C1, 3, 1, 0.1666666666666666666666666667);
                                    inner_product(C4p, X2, C2, 2, 2, 0.5);
                                    inner_product(C4p,  R, C3, 1, 3, 1.0);

                                    #if expansion_order >= 5
                                          C5  = (FlatTree + a)        -> C5;
                                          C5p = (FlatTree + idParent) -> C5;
                                          for (p = 0; p < 21; p ++){
                                                C5[p] += C5p[p];
                                          }
                                          for (p = 0; p < 15; p ++){
                                                get_Xn_overwrite(p, 4, R, 1.0, X4);
                                          }
                                          inner_product(C5p, X4, C1, 4, 1, 0.0416666666666666666666666667);
                                          inner_product(C5p, X3, C2, 3, 2, 0.1666666666666666666666666667);
                                          inner_product(C5p, X2, C3, 2, 3, 0.5);
                                          inner_product(C5p,  R, C4, 1, 4, 1.0);

                                          #if expansion_order >= 6
                                                C6  = (FlatTree + a)        -> C6;
                                                C6p = (FlatTree + idParent) -> C6;
                                                for (p = 0; p < 28; p ++){
                                                      C6[p] += C6p[p];
                                                }
                                                for (p = 0; p < 21; p ++){
                                                      get_Xn_overwrite(p, 5, R, 1.0, X5);
                                                }
                                                inner_product(C6p, X5, C1, 5, 1, 0.0083333333333333333333333333);
                                                inner_product(C6p, X4, C2, 4, 2, 0.0416666666666666666666666667);
                                                inner_product(C6p, X3, C3, 3, 3, 0.1666666666666666666666666667);
                                                inner_product(C6p, X2, C4, 2, 4, 0.5);
                                                inner_product(C6p,  R, C5, 1, 5, 1.0);
                                                #if expansion_order >= 7
                                                      C7  = (FlatTree + a)        -> C7;
                                                      C7p = (FlatTree + idParent) -> C7;
                                                      for (p = 0; p < 36; p ++){
                                                            C7[p] += C7p[p];
                                                      }
                                                      for (p = 0; p < 28; p ++){
                                                            get_Xn_overwrite(p, 6, R, 1.0, X6);
                                                      }
                                                      inner_product(C7p, X6, C1, 6, 1, 0.0013888888888888888888889);
                                                      inner_product(C7p, X5, C2, 5, 2, 0.0083333333333333333333333);
                                                      inner_product(C7p, X4, C3, 4, 3, 0.0416666666666666666666667);
                                                      inner_product(C7p, X3, C4, 3, 4, 0.1666666666666666666666667);
                                                      inner_product(C7p, X2, C5, 2, 5, 0.5);
                                                      inner_product(C7p,  R, C6, 1, 6, 1.0);
                                                      #if expansion_order >= 8
                                                            C8  = (FlatTree + a)        -> C8;
                                                            C8p = (FlatTree + idParent) -> C8;
                                                            for (p = 0; p < 45; p ++){
                                                                  C8[p] += C8p[p];
                                                            }
                                                            for (p = 0; p < 36; p ++){
                                                                  get_Xn_overwrite(p, 7, R, 1.0, X7);
                                                            }
                                                            inner_product(C8p, X7, C1, 7, 1, 0.0001984126984126984126984);
                                                            inner_product(C8p, X6, C2, 6, 2, 0.0013888888888888888888889);
                                                            inner_product(C8p, X5, C3, 5, 3, 0.0083333333333333333333333);
                                                            inner_product(C8p, X4, C4, 4, 4, 0.0416666666666666666666667);
                                                            inner_product(C8p, X3, C5, 3, 5, 0.1666666666666666666666667);
                                                            inner_product(C8p, X2, C6, 2, 6, 0.5);
                                                            inner_product(C8p,  R, C7, 1, 7, 1.0);
                                                      #endif
                                                #endif
                                          #endif
                                    #endif
                              #endif
                        #endif
                  #endif
            }
            
            /******** If the node has children, I put them in the stack ********/
            how_many_children = (FlatTree + a) -> how_many_children;
            if (how_many_children > 0){
                  idFirstChild = (FlatTree + a) -> idFirstChild;
                  idLastChild  = idFirstChild + how_many_children;
                  for (p = idFirstChild; p < idLastChild; p ++){
                        stack[j] = p;
                        j ++;
                  }
            }
            
            /******** If the node has no children, I shift and accumulate its C^(n) to its body ********/
            else{
                  Na   = (FlatTree + a) -> how_many_dots;
                  dots = (FlatTree + a) -> dots;
                  for (p = 0; p < Na; p ++){
                        b  = dots[p];
                        C1 = (FlatTree + a) -> C1;
                        C1Moonlets[3*b] += C1[0];  C1Moonlets[3*b + 1] += C1[1];  C1Moonlets[3*b + 2] += C1[2];
                        
                        #if expansion_order >= 2
                              com  = (FlatTree + a) -> com;
                              C2   = (FlatTree + a) -> C2;
                              R[0] = (moonlets + b) -> x - com[0];
                              R[1] = (moonlets + b) -> y - com[1];
                              R[2] = (moonlets + b) -> z - com[2];
                              C1Moonlets[3*b]     += R[0]*C2[0] + R[1]*C2[1] + R[2]*C2[2];
                              C1Moonlets[3*b + 1] += R[0]*C2[1] + R[1]*C2[3] + R[2]*C2[4];
                              C1Moonlets[3*b + 2] += R[0]*C2[2] + R[1]*C2[4] + R[2]*C2[5];
                              
                              #if expansion_order >= 3
                                    C3 = (FlatTree + a) -> C3;
                                    X2[0]  = R[0]*R[0]; X2[1] = R[0]*R[1]; X2[2] = R[0]*R[2]; X2[3] = R[1]*R[1]; X2[4] = R[1]*R[2]; X2[5] = R[2]*R[2];
                                    C1Moonlets[3*b]     += 0.5*(X2[0]*C3[0] + 2.0*X2[1]*C3[1] + 2.0*X2[2]*C3[2] + X2[3]*C3[3] + 2.0*X2[4]*C3[4] + X2[5]*C3[5]);
                                    C1Moonlets[3*b + 1] += 0.5*(X2[0]*C3[1] + 2.0*X2[1]*C3[3] + 2.0*X2[2]*C3[4] + X2[3]*C3[6] + 2.0*X2[4]*C3[7] + X2[5]*C3[8]);
                                    C1Moonlets[3*b + 2] += 0.5*(X2[0]*C3[2] + 2.0*X2[1]*C3[4] + 2.0*X2[2]*C3[5] + X2[3]*C3[7] + 2.0*X2[4]*C3[8] + X2[5]*C3[9]);
                                    #if expansion_order >= 4
                                          C4 = (FlatTree + a) -> C4;
                                          for (q = 0; q < 10; q ++){
                                                get_Xn_overwrite(q, 3, R, 1.0, X3);
                                          }
                                          inner_product(C4, X3, C1Moonlets + 3*b, 3, 1, 0.1666666666666666666666666667);
                                          
                                          #if expansion_order >= 5
                                                C5 = (FlatTree + a) -> C5;
                                                for (q = 0; q < 15; q ++){
                                                      get_Xn_overwrite(q, 4, R, 1.0, X4);
                                                }
                                                inner_product(C5, X4, C1Moonlets + 3*b, 4, 1, 0.0416666666666666666666666667);
                                                
                                                #if expansion_order >= 6
                                                      C6 = (FlatTree + a) -> C6;
                                                      for (q = 0; q < 21; q ++){
                                                            get_Xn_overwrite(q, 5, R, 1.0, X5);
                                                      }
                                                      inner_product(C6, X5, C1Moonlets + 3*b, 5, 1, 0.0083333333333333333333333333);
                                                      #if expansion_order >= 7
                                                            C7 = (FlatTree + a) -> C7;
                                                            for (q = 0; q < 28; q ++){
                                                                  get_Xn_overwrite(q, 6, R, 1.0, X6);
                                                            }
                                                            inner_product(C7, X6, C1Moonlets + 3*b, 6, 1, 0.0013888888888888888888889);
                                                            #if expansion_order >= 8
                                                                  C8 = (FlatTree + a) -> C8;
                                                                  for (q = 0; q < 36; q ++){
                                                                        get_Xn_overwrite(q, 7, R, 1.0, X7);
                                                                  }
                                                                  inner_product(C8, X7, C1Moonlets + 3*b, 7, 1, 0.0001984126984126984126984);
                                                            #endif
                                                      #endif
                                                #endif
                                          #endif
                                    #endif
                              #endif
                        #endif
                  }
            }
            i ++;
      }
      free(stack);
      stack = NULL;
}
#endif


/******************************************************************************************/
/******** I now write functions allowing to compute gravity with the tree-code of  ********/
/******** Barnes & Hut. The original algorithm of Barnes & Hut performs a first    ********/
/******** order Taylor expansion, but here, the order is given by expansion_order. ********/
/******************************************************************************************/


#if mutual_bool
void standard_tree_acceleration(struct node * FlatTree, struct moonlet * moonlets, int b){

      /******** Computes the acceleration of body b from the multipole moments of the          ********/
      /******** flattree FlatTree. Stores the result into C1Moonlets[3*b] to C1Moonlets[3*b+2] ********/
      /******** This is the standard tree code at expansion order expansion_order              ********/
      
      /******** Retrieving the bodies's coordinates ********/
      typ X, Y, Z, Rb;
      X  = (moonlets + b) -> x;
      Y  = (moonlets + b) -> y;
      Z  = (moonlets + b) -> z;
      Rb = (moonlets + b) -> radius;
      
      /******** Initializing the acceleration ********/
      typ acc[3] = {0.0, 0.0, 0.0};
      
      typ M0;    //Node's mass
      typ * com; //Node's center of mass
      typ R[3];  //bodies's position - center of mass
      typ distance; // |R|
      typ r_crit; // Critical distance for well-separation
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that have to be considered
      int i = 0; //Index in the stack of the id of the current node
      int j = 0; //Index of where to put a node in the stack
      int a;     //Id of the current node
      int k;
      int how_many_children; 
      int how_many_dots;
      int idFirstChild;
      int idLastChild;
      int * dots; //All the bodies in that cell
      int index; //Index of a body in dots
      typ m, Rad, x, y, z, dx, dy, dz, r, r3, softening;
      #if expansion_order >= 3
            typ * M2;  //Node's multipole moments
            #if expansion_order >= 4
                  typ * M3;
                  #if expansion_order >= 5
                        typ * M4;
                        #if expansion_order >= 6
                              typ * M5;
                              #if expansion_order >= 7
                                    typ * M6;
                                    #if expansion_order >= 8
                                          typ * M7;
                                    #endif
                              #endif
                        #endif
                  #endif
            #endif
      #endif
      #if expansion_order   == 1
      typ Grad[3]; //Independant components of all the gradients of G/R from order 1 to expansion_order
      #elif expansion_order == 2
      typ Grad[9];
      #elif expansion_order == 3
      typ Grad[19];
      #elif expansion_order == 4
      typ Grad[34];
      #elif expansion_order == 5
      typ Grad[55];
      #elif expansion_order == 6
      typ Grad[83];
      #elif expansion_order == 7
      typ Grad[119];
      #elif expansion_order == 8
      typ Grad[164];
      #endif
      
      stack[j] = 0;
      j++;
      
      /******** At each node, if it has at most N_cb_pre bodies, I treat it directly, otherwise,     ********/
      /******** if it is well-separated, I treat it by multipole expansion, otherwise, if it has at  ********/
      /******** most N_cb_post bodies or no children, I treat it directly, otherwise, I subdivise it ********/
      while (j > i){
            a = stack[i];
            how_many_dots = (FlatTree + a) -> how_many_dots;
            if (how_many_dots < N_cb_pre){ //Direct summation
                  dots = (FlatTree + a) -> dots;
                  for (index = 0; index < how_many_dots; index++){
                        k = dots[index];
                        if (k != b){ //If the body is different from body b, I accumulate its contribution to the acceleration
                              m   = (moonlets + k) -> mass;
                              Rad = (moonlets + k) -> radius;
                              x   = (moonlets + k) -> x;
                              y   = (moonlets + k) -> y;
                              z   = (moonlets + k) -> z;
                              dx  = x - X;
                              dy  = y - Y;
                              dz  = z - Z;
                              softening = softening_parameter * (Rb + Rad);
                              r   = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                              r3  = r*r*r;
                              acc[0] += G*m*dx/r3;
                              acc[1] += G*m*dy/r3;
                              acc[2] += G*m*dz/r3;
                        }
                  }   
            }
            else{
                  r_crit = (FlatTree + a) -> r_crit;
                  com    = (FlatTree + a) -> com;
                  R[0] = X - com[0];
                  R[1] = Y - com[1];
                  R[2] = Z - com[2];
                  distance = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
                  if (distance >= r_crit){ //Multipole expansion
                        gradR(R, Grad);
                        M0 = (FlatTree + a) -> M0;
                        inner_product(Grad, &M0, acc, 0, 1, 1.0);
                        #if expansion_order >= 3
                              M2 = (FlatTree + a) -> M2;
                              inner_product(Grad + 9, M2, acc, 2, 1, 0.5);
                              #if expansion_order >= 4
                                    M3 = (FlatTree + a) -> M3;
                                    inner_product(Grad + 19, M3, acc, 3, 1, -0.1666666666666666666666667);
                                    #if expansion_order >= 5
                                          M4 = (FlatTree + a) -> M4;
                                          inner_product(Grad + 34, M4, acc, 4, 1, 0.0416666666666666666666667);
                                          #if expansion_order >= 6
                                                M5 = (FlatTree + a) -> M5;
                                                inner_product(Grad + 55, M5, acc, 5, 1, -0.0083333333333333333333333);
                                                #if expansion_order >= 7
                                                      M6 = (FlatTree + a) -> M6;
                                                      inner_product(Grad + 83, M6, acc, 6, 1, 0.0013888888888888888888889);
                                                      #if expansion_order >= 8
                                                            M7 = (FlatTree + a) -> M7;
                                                            inner_product(Grad + 119, M7, acc, 7, 1, -0.0001984126984126984126984);
                                                      #endif
                                                #endif
                                          #endif
                                    #endif
                              #endif
                        #endif
                  }
                  else{
                        how_many_children = (FlatTree + a) -> how_many_children;
                        if (how_many_dots < N_cb_post || how_many_children == 0){ //Direct summation
                              dots = (FlatTree + a) -> dots;
                              for (index = 0; index < how_many_dots; index ++){
                                    k = dots[index];
                                    if (k != b){ //If the body is different from body b, I accumulate its contribution to the acceleration
                                          m   = (moonlets + k) -> mass;
                                          Rad = (moonlets + k) -> radius;
                                          x   = (moonlets + k) -> x;
                                          y   = (moonlets + k) -> y;
                                          z   = (moonlets + k) -> z;
                                          dx  = x - X;
                                          dy  = y - Y;
                                          dz  = z - Z;
                                          softening = softening_parameter*(Rb + Rad);
                                          r   = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                                          r3  = r*r*r;
                                          acc[0] += G*m*dx/r3;
                                          acc[1] += G*m*dy/r3;
                                          acc[2] += G*m*dz/r3;
                                    }
                              }   
                        }
                        else{ //Subdivision
                              idFirstChild = (FlatTree + a) -> idFirstChild;
                              idLastChild = idFirstChild + how_many_children;
                              for (k = idFirstChild; k < idLastChild; k++){
                                    stack[j] = k;
                                    j++; 
                              }
                        }
                  }
            }
            i++;
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function standard_tree_acceleration. Aborting before segmentation fault.\n");
            abort();
      }
      free(stack);
      stack = NULL;
      
      /******** Storing the acceleration into C1Moonlets ********/
      C1Moonlets[3*b]     = acc[0];
      C1Moonlets[3*b + 1] = acc[1];
      C1Moonlets[3*b + 2] = acc[2];
}
#endif


/*********************************************************/
/******** Functions relative to tree manipulation ********/
/*********************************************************/



void create_boxdot(struct boxdot ** BoxDot, typ * corner_coordinates, typ D){

      /******** Allocates memory for a boxdot at address *Boxdot     ********/
      /******** Initializes the fields dots and oct[i] to NULL       ********/
      /******** D is the sidelength of BoxDot and corner_coordinates ********/
      /******** will initialize the field corner                     ********/

      *BoxDot = (struct boxdot *)malloc(sizeof(struct boxdot));
      (*BoxDot) -> dots = NULL;
      struct boxdot ** octants = (*BoxDot) -> oct;
      int i;
      for (i = 0; i < 8; i ++){
            *(octants + i) = NULL;
      }
      (*BoxDot)  -> sidelength = D;
      (*BoxDot)  -> how_many   = 0;
      ((*BoxDot) -> corner)[0] = * corner_coordinates;
      ((*BoxDot) -> corner)[1] = *(corner_coordinates + 1);
      ((*BoxDot) -> corner)[2] = *(corner_coordinates + 2);
}


void clear_boxdot(struct boxdot ** BoxDot){

      /******** Frees the memory allocated to the boxdot Boxdot and to its descendants ********/
      /******** A stack is used to avoid recursivity                                   ********/
      
      /******** To be removed when the code is robust ********/
      if (*BoxDot == NULL){
            fprintf(stderr, "Error : clear_boxdot has been called on a NULL boxdot.");
            abort();
      }
      
      int p;
      
      /******** Stack of nodes still to be deleted ********/
      struct boxdot ** stack = (struct boxdot **)malloc(how_many_cells*sizeof(struct boxdot *));
      
      int i = 0, j = 0; // i is the index of the next node to be treated. j is the index of where to store a node to be treated
      
      stack[j] = *BoxDot;
      j ++;
      
      struct boxdot * to_be_deleted;
      struct boxdot * octant;
      
      while (j > i){ //While there are still nodes to be deleted
            
            to_be_deleted = stack[i];
            clear_chain(&(to_be_deleted -> dots)); //Clearing the chain of dots
            for (p = 0; p < 8; p ++){ //Adding to the stack the children of the current node
                  octant = (to_be_deleted -> oct)[p];
                  if (octant != NULL){
                        stack[j] = octant;
                        j ++;
                  }
            }
            free(to_be_deleted);
            to_be_deleted = NULL;
            i ++;
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > how_many_cells){
            fprintf(stderr, "Error : The stack is not big enough in function clear_boxdot. Aborting before segmentation fault.\n");
            abort();
      }
      free(stack);
      stack = NULL;
}


int get_octant(typ x, typ y, typ z, typ xa, typ ya, typ za, typ D){

      /******** Returns the octant in which a body with coordinates (xa,ya,za) ********/
      /******** should go, given the corner coordinates xyz of the parent cell ********/
      /******** D is half the parent sidelength                                ********/

      
      if (xa < x + D){             //Body is in octant 0, 2, 4 or 6
            if (ya < y + D){       //Body is in octant 0 or 2
                  if (za > z - D){ //Body is in octant 0
                        return 0;
                  }
                  else{            //Body is in octant 2
                        return 2;
                  }
            }
            else{                  //Body is in octant 4 or 6
                  if (za > z - D){ //Body is in octant 4
                        return 4;
                  }
                  else{            //Body is in octant 6
                        return 6;
                  }
            }
      }
      else {                       //Body is in octant 1, 3, 5 or 7
            if (ya < y + D){       //Body is in octant 1 or 3
                  if (za > z - D){ //Body is in octant 1
                        return 1;
                  }
                  else{            //Body is in octant 3
                        return 3;
                  }
            }
            else{                  //Body is in octant 5 or 7
                  if (za > z - D){ //Body is in octant 5
                        return 5;
                  }
                  else{            //Body is in octant 7
                        return 7;
                  }
            }
      }
}


void get_corner_coordinates(typ X, typ Y, typ Z, typ D, int i, typ * corner){

      /******** Returns the corner coordinates of the i^th child, given the corner ********/
      /******** coordinates of the parent. D is half the parent's sidelength       ********/
      
      typ x,y,z; //Child's coordinates
      
      /******** Retrieving the child's coordinates ********/
      if (i == 0){
            x = X;  y = Y;  z = Z;
      }
      else if (i == 1){
            x = X + D;  y = Y;  z = Z;
      }
      else if (i == 2){
            x = X;  y = Y;  z = Z - D;
      }
      else if (i == 3){
            x = X + D;  y = Y;  z = Z - D;
      }
      else if (i == 4){
            x = X;  y = Y + D;  z = Z;
      }
      else if (i == 5){
            x = X + D;  y = Y + D;  z = Z;
      }
      else if (i == 6){
            x = X;  y = Y + D;  z = Z - D;
      }
      else{
            x = X + D;  y = Y + D;  z = Z - D;
      }
      
      /******** Returning the child's coordinates ********/
      *corner       = x;
      *(corner + 1) = y;
      *(corner + 2) = z;
}

void print_boxdot(struct boxdot * BoxDot){

      /******** Prints in terminal the boxdot BoxDot  ********/
      /******** This function is not used for the FFM ********/
      /******** To be removed eventually              ********/

      int p;
      typ * corner;
      typ D;
      int octant;
      int rotation;
      typ distance;
      typ previous_center[3] = {0.0, 0.0, 0.0};
      typ center[3]          = {0.0, 0.0, 0.0};
      
      /******** Stack of nodes still to be printed ********/
      struct boxdot ** stack = (struct boxdot **)malloc(how_many_cells*sizeof(struct boxdot *));

      int i = 0, j = 0; // i is the index of the next node to be treated. j is the index of where to store a node to be treated
      
      stack[j] = BoxDot;
      j++;
      
      struct boxdot * to_be_printed;
      struct boxdot * child;
      while (j>i){ //While there are still nodes to be printed
            to_be_printed = stack[i];
            corner = to_be_printed -> corner;
            D = to_be_printed -> sidelength;
            rotation = to_be_printed -> rotation;
            center[0] = corner[0]+D/2.0;  center[1] = corner[1]+D/2.0;  center[2] = corner[2]-D/2.0;
            distance  = sqrt((center[0]-previous_center[0])*(center[0]-previous_center[0]) + (center[1]-previous_center[1])*(center[1]-previous_center[1]) + 
            (center[2]-previous_center[2])*(center[2]-previous_center[2]));
            previous_center[0] = center[0];  previous_center[1] = center[1];  previous_center[2] = center[2];
            
            printf("Node n° %d : N = %d, id = %d, level = %d, distance from previous node = %.6lf, center = {x,y,z} = {%.6lf, %.6lf, %.6lf}\n", i, 
            to_be_printed -> how_many, to_be_printed -> id, to_be_printed -> level, distance, center[0], center[1], center[2]);
            
            for (p = 0; p < 8; p ++){ // p is the Hilbert-Peano digit
                  octant = OctantFromDigit[rotation][p];
                  child = (to_be_printed -> oct)[octant];
                  if (child != NULL){
                        stack[j] = child;
                        j++;
                  }
            }
            i++;
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > how_many_cells){
            fprintf(stderr, "Error : The stack is not big enough in function print_boxdot. Aborting before segmentation fault.\n");
            abort();
      }
      free(stack);
      stack = NULL;
}
