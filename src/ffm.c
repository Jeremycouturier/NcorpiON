#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "structure.h"
#include "ffm.h"
#include <errno.h>
#include <math.h>

#include "parameters.h"

/*******************************************************************************************************************/
/******** In this file, we implement the fast multipole algorithm falcON for O(N) mutual gravity evaluation ********/
/*******************************************************************************************************************/


typ Mtot;
int how_many_cells;
int cell_id;
int k_from_s2s3[7][7];
int s1s2s3_from_kn[28][7][3];
int s2s3_from_k[28][2];
int ijklmn_from_kn[28][7][6];
int k_from_ijklmn[4][4][4][4][4][4];
int factorial[7] = {1, 1, 2, 6, 24, 120, 720};
int perm_from_kn[28][7];
int q1fromq2q3[28][28];
int * PeanoHilbertOrder;
int IndexPeanoHilbertOrder;
typ * C1Moonlets;
typ * C2FlatTree;
typ * C3FlatTree;
typ * C4FlatTree;
typ * C5FlatTree;
typ * C6FlatTree;
typ * M2FlatTree;
typ * M3FlatTree;
typ * M4FlatTree;
typ * M5FlatTree;

/******** Arrays relative to Hilbert-Peano order ********/
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

int DigitFromOctant[48][8] = {
      { 0, 3, 1, 2, 7, 4, 6, 5},
      { 7, 0, 6, 1, 4, 3, 5, 2},
      { 4, 7, 5, 6, 3, 0, 2, 1},
      { 3, 4, 2, 5, 0, 7, 1, 6},
      { 1, 2, 6, 5, 0, 3, 7, 4},
      { 0, 1, 7, 6, 3, 2, 4, 5},
      { 3, 0, 4, 7, 2, 1, 5, 6},
      { 2, 3, 5, 4, 1, 0, 6, 7},
      { 6, 5, 7, 4, 1, 2, 0, 3},
      { 1, 6, 0, 7, 2, 5, 3, 4},
      { 2, 1, 3, 0, 5, 6, 4, 7},
      { 5, 2, 4, 3, 6, 1, 7, 0},
      { 7, 4, 0, 3, 6, 5, 1, 2},
      { 6, 7, 1, 0, 5, 4, 2, 3},
      { 5, 6, 2, 1, 4, 7, 3, 0},
      { 4, 5, 3, 2, 7, 6, 0, 1},
      { 6, 1, 5, 2, 7, 0, 4, 3},
      { 7, 6, 4, 5, 0, 1, 3, 2},
      { 0, 7, 3, 4, 1, 6, 2, 5},
      { 1, 0, 2, 3, 6, 7, 5, 4},
      { 2, 5, 1, 6, 3, 4, 0, 7},
      { 3, 2, 0, 1, 4, 5, 7, 6},
      { 4, 3, 7, 0, 5, 2, 6, 1},
      { 5, 4, 6, 7, 2, 3, 1, 0},
      { 7, 4, 6, 5, 0, 3, 1, 2},
      { 0, 7, 1, 6, 3, 4, 2, 5},
      { 3, 0, 2, 1, 4, 7, 5, 6},
      { 4, 3, 5, 2, 7, 0, 6, 1},
      { 6, 5, 1, 2, 7, 4, 0, 3},
      { 7, 6, 0, 1, 4, 5, 3, 2},
      { 4, 7, 3, 0, 5, 6, 2, 1},
      { 5, 4, 2, 3, 6, 7, 1, 0},
      { 1, 2, 0, 3, 6, 5, 7, 4},
      { 6, 1, 7, 0, 5, 2, 4, 3},
      { 5, 6, 4, 7, 2, 1, 3, 0},
      { 2, 5, 3, 4, 1, 6, 0, 7},
      { 0, 3, 7, 4, 1, 2, 6, 5},
      { 1, 0, 6, 7, 2, 3, 5, 4},
      { 2, 1, 5, 6, 3, 0, 4, 7},
      { 3, 2, 4, 5, 0, 1, 7, 6},
      { 1, 6, 2, 5, 0, 7, 3, 4},
      { 0, 1, 3, 2, 7, 6, 4, 5},
      { 7, 0, 4, 3, 6, 1, 5, 2},
      { 6, 7, 5, 4, 1, 0, 2, 3},
      { 5, 2, 6, 1, 4, 3, 7, 0},
      { 4, 5, 7, 6, 3, 2, 0, 1},
      { 3, 4, 0, 7, 2, 5, 1, 6},
      { 2, 3, 1, 0, 5, 4, 6, 7}
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

/*********************************************************************************************************************/
/******** First, we implement functions allowing to put moonlets in a tree-like structure called boxdot       ********/
/*********************************************************************************************************************/

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
                  current_power = (current_power-1)/2;
                  to_be_multiplied_by *= to_be_returned;
            }
            to_be_returned *= to_be_returned;
      }
      return to_be_multiplied_by*to_be_returned;
}


void subdivision(struct boxdot * BoxDot, struct moonlet * moonlets){

      /******** Auxiliary function to add_boxdot.            ********/
      /******** Adds the moonlets of BoxDot to its children. ********/

      int a; //The current moonlet
      typ x, y, z, xa, ya, za, D; //Coordinates
      int octant; //The octant where the moonlet must go
      struct chain * dots = BoxDot -> dots; //Chain of ids of moonlets inside BoxDot
      int index = dots -> how_many - 1; //The index of a in BoxDot -> dots
      
      /******** Getting the box coordinates ********/
      x = (BoxDot -> corner)[0];
      y = (BoxDot -> corner)[1];
      z = (BoxDot -> corner)[2];
      D = BoxDot -> sidelength / 2.0;

      while (dots != NULL){
            a = (dots -> ids)[index]; //The moonlet's id
            
            /******** Getting the moonlet coordinates ********/
            xa = (moonlets + a) -> x;
            ya = (moonlets + a) -> y;
            za = (moonlets + a) -> z;
            
            octant = get_octant(x, y, z, xa, ya, za, D); //Retrieving the octant where the moonlet belongs
            
            /******** Initializing the child if necessary ********/
            if ((BoxDot -> oct)[octant] == NULL){
                  typ corner[3];
                  get_corner_coordinates(x, y, z, D, octant, corner);
                  create_boxdot(&((BoxDot -> oct)[octant]), corner, D);
                  fill_boxdot_int((BoxDot -> oct)[octant], BoxDot, octant);
                  how_many_cells++;
            }
            
            add_boxdot((BoxDot -> oct)[octant], moonlets, a); //The while loop of add_boxdot will never be reached.
                                                              //subdivision will be called again by add_boxdot only if all moonlets of BoxDot go into the same octant
            /******** Switching to the next moonlet ********/                                                  
            index--;
            if (index < 0){
                  dots = dots -> queue;
                  index = max_ids_per_node - 1;
            }
      }
}


void add_boxdot(struct boxdot * BoxDot, struct moonlet * moonlets, int a){

      /******** Adds moonlet a to the boxdot BoxDot, and to its descendants.  ********/
      /******** It is assumed that it has been checked, prior to calling this ********/
      /******** function, that the moonlet a belongs to the boxdot BoxDot.    ********/
      /******** It is also assumed that BoxDot was previously initialized.    ********/
      
      int n = BoxDot -> how_many; //Number of moonlets currently in BoxDot
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
      D = BoxDot -> sidelength;
      x = (BoxDot -> corner)[0];
      y = (BoxDot -> corner)[1];
      z = (BoxDot -> corner)[2];
      
      /******** Retrieving the moonlet coordinates ********/
      typ xa, ya, za;
      xa = (moonlets + a) -> x;
      ya = (moonlets + a) -> y;
      za = (moonlets + a) -> z;
      
      int octant;
      
      while (1){ //While the box has children, we keep getting deeper into the tree
            D /= 2.0;
            octant = get_octant(x, y, z, xa, ya, za, D); //Retrieving the corresponding octant
            
            if ((BoxDot -> oct)[octant] == NULL){ //If the child is still NULL, we initialize it 
                  typ corner[3];
                  get_corner_coordinates(x, y, z, D, octant, corner);
                  create_boxdot(&((BoxDot -> oct)[octant]), corner, D);
                  fill_boxdot_int((BoxDot -> oct)[octant], BoxDot, octant);
                  how_many_cells++;
            }
            
            BoxDot = (BoxDot -> oct)[octant]; //The new considered box is now the child
            n = BoxDot -> how_many;           //Retrieving the number of moonlets it contains
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
      /******** Adds all moonlets to it               ********/
      
      int i, a;
      typ x, y, z; //Moonlet coordinates
      for (i = 0; i <= largest_id; i++){
            already_in_tree[i] = 0;
      }
      
      /******** Creating and initializing the root cell ********/
      struct boxdot * root = NULL;
      typ corner[3] = {-root_sidelength/2.0, -root_sidelength/2.0, root_sidelength/2.0};
      create_boxdot(&root, corner, root_sidelength);
      root -> rotation = 0;
      root -> level = 0; //The root cell is at level 0
      how_many_cells++;
      
      /******** Adding all moonlets to it ********/
      for (i = 0; i < IndexPeanoHilbertOrder; i++){ //Adding the moonlets in Hilbert order
            a = PeanoHilbertOrder[i];
            if (*(exists+a)){
                  x = (moonlets + a) -> x;
                  y = (moonlets + a) -> y;
                  z = (moonlets + a) -> z;
                  if (absolute(x) <= root_sidelength/2.0 && absolute(y) <= root_sidelength/2.0 && absolute(z) <= root_sidelength/2.0){ //If moonlet a belongs to the root cell
                        add_boxdot(root, moonlets, a);
                  }
                  already_in_tree[a] = 1;
            }
      }
      for (i = 0; i <= largest_id; i++){ //Adding the few remaining moonlets in random order
            if (exists[i] && !already_in_tree[i]){
                  x = (moonlets + i) -> x;
                  y = (moonlets + i) -> y;
                  z = (moonlets + i) -> z;
                  if (absolute(x) <= root_sidelength/2.0 && absolute(y) <= root_sidelength/2.0 && absolute(z) <= root_sidelength/2.0){ //If moonlet i belongs to the root cell
                        add_boxdot(root, moonlets, i);
                  }
            }
      }
      
      return root;
}


/***********************************************************************************************************************/
/******** We now link the boxdot tree to a simple array containing the multipole moments and field tensor C^(n) ********/
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
      BoxDot -> level = levelChild;
}


struct node * flattree_init(struct boxdot * BoxDot){

      /******** This function, when called on the root cell, initializes its unique id and that of all its descendants            ********/
      /******** It is assumed that the global variable cell_id is 0 prior to calling this function and that the global variable   ********/
      /******** how_many_cells is the total number of cells. This function also returns an array FlatTree on which the three      ********/
      /******** stages of Dehnen's FalcON algorithm are performed, and initializes the fields idFirstChild, how_many_children,    ********/
      /******** idParent, how_many_dots, dots, sidelength and corner of each node of FlatTree                                     ********/
      
      int p;
      int octant;
      int rotation;
      int how_many_child = 0;
      int isFirstChild = 1;
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
      (FlatTree + cell_id) -> idParent = -1; //Arbitrarily setting to -1 the unique id of the parent of the root cell
      
      /******** Stack of nodes still to be treated ********/
      struct boxdot ** stack = (struct boxdot **)malloc(how_many_cells * sizeof(struct boxdot *));

      int j = 0; // j is the index of where to store a node to be treated
      
      stack[j] = BoxDot;
      j++;
      
      struct boxdot * to_be_treated;
      struct boxdot * child;
      
      while (j > cell_id){ //While there are still nodes to be treated
            to_be_treated = stack[cell_id];
            to_be_treated -> id = cell_id; //Initializing the unique id of that cell
            rotation = to_be_treated -> rotation;
            D = to_be_treated -> sidelength;
            
            for (p = 0; p < 8; p++){ // p is the Hilbert-Peano digit
                  octant = OctantFromDigit[rotation][p];
                  child = (to_be_treated -> oct)[octant];
                  if (child != NULL){
                        stack[j] = child;
                        (FlatTree + j) -> idParent = cell_id; //Initializing the field idParent of all the children
                        if (isFirstChild){ //Initializing the field idFirstChild 
                              (FlatTree + cell_id) -> idFirstChild = j;
                              isFirstChild = 0;
                        }
                        j++;
                        how_many_child++;
                  }
            }
            (FlatTree + cell_id) -> how_many_children = how_many_child;
            how_many_dots = to_be_treated -> how_many;
            (FlatTree + cell_id) -> how_many_dots     = how_many_dots;
            
            /******** Initializing the field array dots ********/
            dots = (int *)malloc(how_many_dots * sizeof(int));
            (FlatTree + cell_id) -> dots              = dots;
            if (dots == NULL){
                  fprintf(stderr, "Error : Could not allocate memory for %d dots in function flattree_init.\n", how_many_dots);
                  abort();
            }
            ch = to_be_treated -> dots;
            index = ch -> how_many - 1;
            for (p = 0; p < how_many_dots; p++){
                  *(dots + p) = (ch -> ids)[index];
                  if (how_many_child == 0){
                        PeanoHilbertOrder[IndexPeanoHilbertOrder] = (ch -> ids)[index];
                        IndexPeanoHilbertOrder ++;
                  }
                  /******** Switching to the next moonlet ********/                                                  
                  index--;
                  if (index < 0){
                        ch = ch -> queue;
                        index = max_ids_per_node - 1;
                  }
            }
            
            /******** Initializing the C^(1) ********/
            ((FlatTree + cell_id) -> C1)[0] = 0.0;  ((FlatTree + cell_id) -> C1)[1] = 0.0;  ((FlatTree + cell_id) -> C1)[2] = 0.0;
            
            (FlatTree + cell_id) -> sidelength        = D;
            ((FlatTree + cell_id) -> center)[0]       = (to_be_treated -> corner)[0] + D/2.0;
            ((FlatTree + cell_id) -> center)[1]       = (to_be_treated -> corner)[1] + D/2.0;
            ((FlatTree + cell_id) -> center)[2]       = (to_be_treated -> corner)[2] - D/2.0;
            if (how_many_child == 0){ //If to_be_treated has no children, then we arbitrarily set the unique id of its first child to -1 
                  (FlatTree + cell_id) -> idFirstChild = -1;
            }
            how_many_child = 0;
            isFirstChild = 1;
            cell_id++;
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


void tensor_initialization(){

      int p;

      /******** Allocating and initializing the global FlatTree C^(m) and M^(n) ********/
      for (p = 0; p <= 3*largest_id+2; p++){
            *(C1Moonlets+p) = 0.0;
      }
      if (expansion_order >= 2){
            C2FlatTree = (typ *)malloc(how_many_cells * 6 * sizeof(typ));
            for (p = 0; p < 6 * how_many_cells; p++){
                  *(C2FlatTree+p) = 0.0;
            }
            if (expansion_order >= 3){
                  C3FlatTree = (typ *)malloc(how_many_cells * 10 * sizeof(typ));
                  M2FlatTree = (typ *)malloc(how_many_cells * 6  * sizeof(typ));
                  for (p = 0; p < 10 * how_many_cells; p++){
                        *(C3FlatTree+p) = 0.0;
                  }
                  if (expansion_order >= 4){
                        C4FlatTree = (typ *)malloc(how_many_cells * 15 * sizeof(typ));
                        M3FlatTree = (typ *)malloc(how_many_cells * 10 * sizeof(typ));
                        for (p = 0; p < 15 * how_many_cells; p++){
                              *(C4FlatTree+p) = 0.0;
                        }
                        if (expansion_order >= 5){
                              C5FlatTree = (typ *)malloc(how_many_cells * 21 * sizeof(typ));
                              M4FlatTree = (typ *)malloc(how_many_cells * 15 * sizeof(typ));
                              for (p = 0; p < 21 * how_many_cells; p++){
                                    *(C5FlatTree+p) = 0.0;
                              }
                              if (expansion_order >= 6){
                                    C6FlatTree = (typ *)malloc(how_many_cells * 28 * sizeof(typ));
                                    M5FlatTree = (typ *)malloc(how_many_cells * 21 * sizeof(typ));
                                    for (p = 0; p < 28 * how_many_cells; p++){
                                          *(C6FlatTree+p) = 0.0;
                                    }
                              }
                        }
                  }
            }
      }
}


void tensor_free(){

      /******** Frees the global arrays C^(n) and M^(n) ********/

      if (expansion_order >= 2){
            free(C2FlatTree);
            C2FlatTree = NULL;
            if (expansion_order >= 3){
                  free(C3FlatTree);
                  C3FlatTree = NULL;
                  free(M2FlatTree);
                  M2FlatTree = NULL;
                  if (expansion_order >= 4){
                        free(C4FlatTree);
                        C4FlatTree = NULL;
                        free(M3FlatTree);
                        M3FlatTree = NULL;
                        if (expansion_order >= 5){
                              free(C5FlatTree);
                              C5FlatTree = NULL;
                              free(M4FlatTree);
                              M4FlatTree = NULL;
                              if (expansion_order >= 6){
                                    free(C6FlatTree);
                                    C6FlatTree = NULL;
                                    free(M5FlatTree);
                                    M5FlatTree = NULL;
                              }
                        }
                  }
            }
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
      
      int n1 = (int) integral(0.5*(sqrt(1.0+8.0*(typ) k)-1.0));
      int n2 = k - (n1*(n1 + 1))/2;
      *s1 = n - n1;
      *s2 = n1 - n2;
      *s3 = n2;

}


void get_s2_s3(int k, int * s2, int * s3){

      /******** Same as above but depends only on k and computes s2 and s3 only ********/
      
      int n1 = (int) integral(0.5*(sqrt(1.0+8.0*(typ) k)-1.0));
      int n2 = k - (n1*(n1 + 1))/2;
      *s2 = n1 - n2;
      *s3 = n2;

}


void s1s2s3_from_kn_init(){

      /******** In order not to call function get_s1_s2_s3 too many times, we store its return values in a table s1s2s3_from_kn ********/
      
      int k, n, s1, s2, s3;
      for (k = 0; k < 28; k++){
            for (n = 0; n < 7; n++){
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
      for (k = 0; k < 28; k++){
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

      /******** In order not to call function get_k too many times, we store its return values in a table k_from_s2s3 ********/
      int k;
      int i,j;
      
      for (i = 0; i <= 6; i++){ // i is the number of 2's
            for (j = 0; j <= 6; j++){ // j is the number of 3's
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
            p++;
            s1--;
      }
      while (s2 > 0){
            ijklmn[p] = 2;
            p++;
            s2--;
      }
      while (s3 > 0){
            ijklmn[p] = 3;
            p++;
            s3--;
      }
      while (p < 6){
            ijklmn[p] = 0;
            p++;
      }
}


void ijklmn_from_kn_init(){

      /******** Stores the return value of indexes_from_kn into ijklmn_from_kn ********/
      
      int ijklmn[6];
      int k, n;
      
      for (k = 0; k < 28; k++){
            for (n = 0; n < 7; n++){
                  indexes_from_kn(k, n, ijklmn);
                  ijklmn_from_kn[k][n][0] = ijklmn[0];
                  ijklmn_from_kn[k][n][1] = ijklmn[1];
                  ijklmn_from_kn[k][n][2] = ijklmn[2];
                  ijklmn_from_kn[k][n][3] = ijklmn[3];
                  ijklmn_from_kn[k][n][4] = ijklmn[4];
                  ijklmn_from_kn[k][n][5] = ijklmn[5];
            }
      }
}


void k_from_indexes(int * k, int * ijklmn){

      /******** Inverse of function indexes_from_kn. Retrieves k from the indexes ********/
      
      int s2 = 0, s3 = 0;
      for (int i = 0; i < 6; i++){
            if (ijklmn[i] == 2){
                  s2++;
            }
            else if (ijklmn[i] == 3){
                  s3++;
            }
      }
      *k = k_from_s2s3[s2][s3];
}


void k_from_ijklmn_init(){

      /******** Stores the return value of k_from_indexes into k_from_ijklmn ********/
      
      int i, j, k, l, m, n;
      int p;
      int ijklmn[6];
      
      for (i = 0; i < 4; i++){
            for (j = 0; j < 4; j++){
                  for (k = 0; k < 4; k++){
                        for (l = 0; l < 4; l++){
                              for (m = 0; m < 4; m++){
                                    for (n = 0; n < 4; n++){
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
      
      int s1, s2, s3;
      s1 = s1s2s3_from_kn[k][n][0];
      s2 = s1s2s3_from_kn[k][n][1];
      s3 = s1s2s3_from_kn[k][n][2];
      if (s1 >= 0 && s2 >= 0 && s3 >= 0){
            *perm = factorial[n]/(factorial[s1]*factorial[s2]*factorial[s3]);
      }
}


void perm_from_kn_init(){

      /******** Stores the return value of perm_from_kn into the array perm_from_kn ********/

      int perm = 0;
      int k, n;
      
      for (k = 0; k < 28; k++){
            for (n = 0; n < 7; n++){
                  permutation_from_kn(k, n, &perm);
                  perm_from_kn[k][n] = perm;
            }
      }
}

void q1_from_q2q3(int * q1, int q2, int q3){

      /******** Function auxiliary to inner_product. Gives q1 as a function of q2 and q3 ********/

      int n1_2, n1_3, n2_2, n2_3, n1, n2, s2_2, s2_3, s3_2, s3_3, s2, s3;
      n1_2 = (int) integral(0.5*(sqrt(1.0+8.0*(typ) q2)-1.0));
      n1_3 = (int) integral(0.5*(sqrt(1.0+8.0*(typ) q3)-1.0));
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
      *q1 = (n1*(n1 + 1))/2 + n2;
}


void q1fromq2q3_init(){

      /******** Stores the return value of q1_from_q2q3 into the array q1fromq2q3 ********/

      int q1, q2, q3;
      for (q2 = 0; q2 < 28; q2++){
            for (q3 = 0; q3 < 28; q3++){
                  q1_from_q2q3(&q1, q2, q3);
                  q1fromq2q3[q2][q3] = q1;
            }
      }
}


/******************************************************************************************/
/******** We now write functions relative to the first stage of Dehnen's algorithm ********/
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
      
      typ x, y, z; //The current moonlet's coordinates
      typ m; //The current moonlet's mass
      
      int * dots = (FlatTree + a) -> dots; //All the moonlets in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int j; //Id of the current moonlet
      
      for (i = 0; i < how_many_dots; i++){
            j = dots[i];
            m = (moonlets + j) -> mass;
            x = (moonlets + j) -> x;
            y = (moonlets + j) -> y;
            z = (moonlets + j) -> z;
            
            M0 += m;
            com[0] += m*x;
            com[1] += m*y;
            com[2] += m*z;
      }
      
      /******** Initializing the relevant fields ********/
      (FlatTree + a) -> M0 = M0;
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
      
      int * dots = (FlatTree + a) -> dots; //All the moonlets in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int j; //Id of the moonlet whose distance to the center of mass is to be computed
      
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
      
      typ M0 = (FlatTree + a) -> M0;
      typ one_theta_min = 1.0-theta_min;
      typ power = ((typ) expansion_order)+2.0;
      typ rhs = fast_pow(theta_min,power)/one_theta_min/one_theta_min*pow(M0/Mtot,-1.0/3.0);

      typ current_precision = 1.0;
      typ theta = 0.7;
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
      typ fX0,dfX0,one_theta,one_theta_2,one_theta_3,theta_power,dtheta;
      int nstep=0;
      while(current_precision > precision){
            one_theta = 1.0-theta;
            one_theta_2 = one_theta*one_theta;
            one_theta_3 = one_theta_2*one_theta;
            theta_power = fast_pow(theta,power);
            fX0 = theta_power/one_theta_2-rhs;
            dfX0 = (2.0*theta_power+power*one_theta*theta_power/theta)/one_theta_3;
            dtheta = -fX0/dfX0;
            theta += dtheta;
            current_precision = absolute(dtheta/theta);
            nstep++;
            
            /******** If the method gets lost, we try to put it back on a right path ********/
            if (theta >= 1.0 || theta <= theta_min){
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


void get_Mn(struct node * FlatTree, struct moonlet * moonlets, int a){

      /******** Computes the multipole moments of node a of FlatTree, assuming that it has ********/
      /******** no children. Initializes the corresponding fields of (FlatTree + a)        ********/
      /******** This function is called only if the expansion order is at least 3          ********/
      
      /******** Retrieving the mass and center of mass of the cell ********/
      typ com[3];
      com[0] = ((FlatTree + a) -> com)[0];  com[1] = ((FlatTree + a) -> com)[1];  com[2] = ((FlatTree + a) -> com)[2];  
      
      typ M2[6] =  {0., 0., 0., 0., 0., 0.}; //The 6 distinct components of the quadrupole
      typ M3[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 10 distinct components of the octupole
      typ M4[15] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 15 distinct components of the fourth multipole moment
      typ M5[21] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 21 distinct components of the fifth multipole moment
      
                                
      typ m; //The moonlet's mass
      
      int * dots = (FlatTree + a) -> dots; //All the moonlets in that cell
      int how_many_dots = (FlatTree + a) -> how_many_dots;
      int i;
      int p; //Id of the current moonlet
      int k; //Array index
      typ dX[3]; // x_i - com
      
      for (i = 0; i < how_many_dots; i++){
            p = dots[i];
            m = (moonlets + p) -> mass;
            dX[0] = (moonlets + p) -> x - com[0];
            dX[1] = (moonlets + p) -> y - com[1];
            dX[2] = (moonlets + p) -> z - com[2];
            
            /******** Accumulating the multipole moments ********/
            M2[0] += m * dX[0] * dX[0];
            M2[1] += m * dX[0] * dX[1];
            M2[2] += m * dX[0] * dX[2];
            M2[3] += m * dX[1] * dX[1];
            M2[4] += m * dX[1] * dX[2];
            M2[5] += m * dX[2] * dX[2];
            if (expansion_order >= 4){ // M3 need to be computed
                  for (k = 0; k < 10; k++){
                        get_Xn(k, 3, dX, m, M3);
                  }
                  if (expansion_order >= 5){ // M4 need to be computed
                        for (k = 0; k < 15; k++){
                              get_Xn(k, 4, dX, m, M4);
                        }
                        if (expansion_order >= 6){ // M5 need to be computed
                              for (k = 0; k < 21; k++){
                                    get_Xn(k, 5, dX, m, M5);
                              }
                        }
                  }
            }
      }
      
      /******** Initializing the fields M2 to M5 of the FlatTree ********/
      for (k = 0; k < 6; k++){
            (M2FlatTree + 6*a)[k] = M2[k];
      }
      if (expansion_order >= 4){
            for (k = 0; k < 10; k++){
                  (M3FlatTree + 10*a)[k] = M3[k];
            }
            if (expansion_order >= 5){
                  for (k = 0; k < 15; k++){
                        (M4FlatTree + 15*a)[k] = M4[k];
                  }
                  if (expansion_order >= 6){
                        for (k = 0; k < 21; k++){
                              (M5FlatTree + 21*a)[k] = M5[k];
                        }
                  }
            }
      }
}


void get_com_from_children(struct node * FlatTree, int a){

      /******** Computes the mass and center of mass of node a of FlatTree from that of ********/
      /******** its children. Initializes the corresponding fields of (FlatTree + a)    ********/

      int idFirstChild = (FlatTree + a) -> idFirstChild;
      int idLastChild  = idFirstChild + (FlatTree + a) -> how_many_children;
      
      typ M0 = 0.0;
      typ com[3] = {0.0, 0.0, 0.0};
      typ M0_child;
      typ * com_child; //Center of mass of a child
      
      for (int i = idFirstChild; i < idLastChild; i++){
            M0_child = (FlatTree + i) -> M0;
            com_child = (FlatTree + i) -> com;
            
            M0 += M0_child;
            com[0] += M0_child * com_child[0];
            com[1] += M0_child * com_child[1];
            com[2] += M0_child * com_child[2];
            
      }

      /******** Initializing the relevant fields ********/
      (FlatTree + a) -> M0 = M0;
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

      typ * center = (FlatTree + a) -> center; //The center of the node
      typ * com = (FlatTree + a) -> com;       //The center of mass of the node
      typ D = (FlatTree + a) -> sidelength;    //The sidelength of the node
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
      for (i = 0; i < 8; i++){
            distance = sqrt(corner[i][0]*corner[i][0] + corner[i][1]*corner[i][1] + corner[i][2]*corner[i][2]);
            if (distance > distance_to_farthest_corner){
                  distance_to_farthest_corner = distance;
            }
      }
      
      typ * com_child; //Center of mass of a child
      typ rmax_child;  //Convergence radius of a child
      typ com_difference[3];
      typ max_ri = 0.0;

      /******** Eq. (9) of Dehnen (2002) ********/
      for (int i = idFirstChild; i < idLastChild; i++){
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


void get_Mn_from_children(struct node * FlatTree, int a){

      /******** Computes the multipole moments of node a of FlatTree from that of its children        ********/
      /******** Initializes the corresponding fields of (FlatTree + a)                                ********/
      /******** This function is called only if the expansion order is at least 3 and if              ********/
      /******** FlatTree + a has at least child_multipole_threshold times more children than moonlets ********/

      int idFirstChild = (FlatTree + a) -> idFirstChild;
      int idLastChild  = idFirstChild + (FlatTree + a) -> how_many_children;
      typ * com = (FlatTree + a) -> com;
      typ * com_child;
      typ * M5_child;
      typ * M4_child;
      typ * M3_child;
      typ * M2_child;
      typ M0_child;
      typ Y[4]; //com - com_child. Index 0 is unused
      int p, q;
      int i, j, k, l, m;
      int * array_of_ijklm;
      
      typ M2[6] =  {0., 0., 0., 0., 0., 0.}; //The 6 distinct components of the quadrupole
      typ M3[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 10 distinct components of the octupole
      typ M4[15] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 15 distinct components of the fourth multipole moment
      typ M5[21] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; //The 21 distinct components of the fifth multipole moment
      
      for (p = idFirstChild; p < idLastChild; p++){
            /******** Retrieving the multipole moments and the center of mass of the child ********/
            M0_child = (FlatTree + p) -> M0;
            com_child = (FlatTree + p) -> com;
            M2_child = M2FlatTree + 6*p;
            if (expansion_order >= 4){
                  M3_child = M3FlatTree + 10*p;
                  if (expansion_order >= 5){
                        M4_child = M4FlatTree + 15*p;
                        if (expansion_order >= 6){
                              M5_child = M5FlatTree + 21*p;
                        }
                  }
            }
            
            /******** Retrieving the difference Y between the expansion centers ********/
            Y[1] = com[0] - com_child[0];  Y[2] = com[1] - com_child[1];  Y[3] = com[2] - com_child[2];  
            
            /******** Computing the multipole moments ********/
            M2[0] += M2_child[0] + M0_child*Y[1]*Y[1];
            M2[1] += M2_child[1] + M0_child*Y[1]*Y[2];
            M2[2] += M2_child[2] + M0_child*Y[1]*Y[3];
            M2[3] += M2_child[3] + M0_child*Y[2]*Y[2];
            M2[4] += M2_child[4] + M0_child*Y[2]*Y[3];
            M2[5] += M2_child[5] + M0_child*Y[3]*Y[3];
            if (expansion_order >= 4){ //M3 is computed
                  for (q = 0; q < 10; q++){
                        array_of_ijklm = ijklmn_from_kn[q][3];
                        i = array_of_ijklm[0]; j = array_of_ijklm[1]; k = array_of_ijklm[2]; 
                        M3[q] += M3_child[q]; //Term in M^3
                        M3[q] -= M0_child*Y[i]*Y[j]*Y[k]; //Term in M^0
                        M3[q] -= Y[k]*M2_child[k_from_ijklmn[i][j][0][0][0][0]]; //Terms in M^2
                        M3[q] -= Y[i]*M2_child[k_from_ijklmn[j][k][0][0][0][0]];
                        M3[q] -= Y[j]*M2_child[k_from_ijklmn[i][k][0][0][0][0]];
                  }
                  if (expansion_order >= 5){ //M4 is computed
                        for (q = 0; q < 15; q++){
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
                        if (expansion_order >= 6){ //M5 is computed
                              for (q = 0; q < 21; q++){
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
                        }
                  }
            }
      }
      
      /******** Initializing the fields M2 to M5 of the FlatTree ********/
      for (k = 0; k < 6; k++){
            (M2FlatTree + 6*a)[k] = M2[k];
      }
      if (expansion_order >= 4){
            for (k = 0; k < 10; k++){
                  (M3FlatTree + 10*a)[k] = M3[k];
            }
            if (expansion_order >= 5){
                  for (k = 0; k < 15; k++){
                        (M4FlatTree + 15*a)[k] = M4[k];
                  }
                  if (expansion_order >= 6){
                        for (k = 0; k < 21; k++){
                              (M5FlatTree + 21*a)[k] = M5[k];
                        }
                  }
            }
      }
}


void com_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Computes and initializes the masses and center of masses of all the nodes of FlatTree ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that could not be treated due to their child not having been treated yet
      int i; //Id of the current node
      int j = 0; //Index of where to put a node in the stack
      int how_many_children;
      
      /******** We travel the flattree. If a node is childless, we compute its mass and center of mass ********/
      /******** Otherwise, we store it in the stack for future treatment                               ********/
      for (i = 0; i < cell_id; i++){
            how_many_children = (FlatTree + i) -> how_many_children;
            if (how_many_children == 0){ //If the node has no children
                  get_com(FlatTree, moonlets, i);
            }
            else{ //If the node has children
                  stack[j] = i;
                  j++;
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function com_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** We now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j--;
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
      
      /******** We travel the flattree. If a node is childless, we compute its convergence radius ********/
      /******** Otherwise, we store it in the stack for future treatment                          ********/
      for (i = 0; i < cell_id; i++){
            how_many_children = (FlatTree + i) -> how_many_children;
            if (how_many_children == 0){ //If the node has no children
                  get_rmax(FlatTree, moonlets, i);
            }
            else{ //If the node has children
                  stack[j] = i;
                  j++;
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function rmax_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** We now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j--;
            i = stack[j];
            get_rmax_from_children(FlatTree, i);
      }
      
      free(stack);
      stack = NULL;
}


void rcrit_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Computes and initializes r_crit for all the nodes of FlatTree ********/
      
      int i;
      
      
      /******** We travel the flattree and compute the tolerance parameter theta and r_crit for each node ********/
      for (i = 0; i < cell_id; i++){
            get_tolerance_parameter(FlatTree, i, 1.0e-6);
      }
}


void multipole_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Computes and initializes the multipole moments of all the nodes of FlatTree ********/
      /******** This function is called only if expansion_order is at least 3               ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of ids of nodes that could not be treated due to their child not having been treated yet
      int i; //Id of the current node
      int j = 0; //Index of where to put a node in the stack
      int how_many_children;
      int how_many_dots;
      
      /******** We travel the flattree. If a node is childless or contains few moonlets, we compute ********/
      /******** its multipole moments. Otherwise, we store it in the stack for future treatment     ********/
      for (i = 0; i < cell_id; i++){
            how_many_children = (FlatTree + i) -> how_many_children;
            if (how_many_children == 0){ //If the node has no children
                  get_Mn(FlatTree, moonlets, i);
            }
            else{ //If the node has children
                  how_many_dots = (FlatTree + i) -> how_many_dots;
                  if (how_many_dots / how_many_children < child_multipole_threshold){
                        get_Mn(FlatTree, moonlets, i);
                  }
                  else{
                        stack[j] = i;
                        j++;
                  }
            }
      }
      
      /******** Making sure that the stack was big enough. To be removed when the code is robust ********/
      if (j > cell_id){
            fprintf(stderr, "Error : The stack is not big enough in function multipole_flattree. Aborting before segmentation fault.\n");
            abort();
      }
      
      /******** We now travel the stack from the end to treat nodes that were not treated previously ********/
      while(j > 0){
            j--;
            i = stack[j];
            get_Mn_from_children(FlatTree, i);
      } 
      free(stack);
      stack = NULL;
}


/**************************************************************************************/
/******** We now write functions to compute the gradient of the Green function ********/
/******** and the inner product between two tensors                            ********/
/**************************************************************************************/


void gradR(typ * R, typ * grad, int p){

      /******** This function computes the p^th gradient of G/R, for p >= 0, and stores it into grad ********/
      /******** The vector R must be given under the form R = {Rx, Ry, Rz}                           ********/
      /******** It is assumed that the array grad is indexed from 0 to (p+1)(p+2)/2 - 1,             ********/
      /******** corresponding to the (p+1)(p+2)/2 independent components of the tensor grad^(p)(1/R) ********/

      typ R1 = R[0], R2 = R[1], R3 = R[2];
      typ R12, R13, R14, R15, R22, R23, R24, R25, R32, R33, R34, R35;
      typ R1_R2, R1_R3, R2_R3;
      typ normR  = sqrt(R1*R1 + R2*R2 + R3*R3); // |R|
      typ normR3; // |R|^3
      typ normR5; // |R|^5
      typ normR7; // |R|^7
      typ normR9; // |R|^9
      typ normR11;// |R|^11
      typ normR13;// |R|^13
      typ D1;
      typ D2;
      typ D3;
      typ D4;
      typ D5;
      typ D6;

      if (p == 0){
            *grad = G/normR;
      }
      else if (p == 1){
            normR3 = normR*normR*normR;
            D1 = -1.0*G/normR3;
            grad[0] = D1*R1;
            grad[1] = D1*R2;
            grad[2] = D1*R3;
      }
      else if (p == 2){
            normR3 = normR*normR*normR;
            normR5 = normR3*normR*normR;
            D1 = -1.0*G/normR3;
            D2 = 3.0*G/normR5;
            grad[0] = D2*R1*R1 + D1;
            grad[1] = D2*R1*R2;
            grad[2] = D2*R1*R3;
            grad[3] = D2*R2*R2 + D1;
            grad[4] = D2*R2*R3;
            grad[5] = D2*R3*R3 + D1;
      }
      else if (p == 3){
            normR5 = normR*normR*normR*normR*normR;
            normR7 = normR5*normR*normR;
            D2 = 3.0*G/normR5;
            D3 = -15.0*G/normR7;
            R12 = R1*R1;  R22 = R2*R2;  R32 = R3*R3;
            grad[0] = D3*R12*R1 + 3.0*D2*R1;
            grad[1] = D3*R12*R2 +     D2*R2;
            grad[2] = D3*R12*R3 +     D2*R3;
            grad[3] = D3*R22*R1 +     D2*R1;
            grad[4] = D3*R1*R2*R3;
            grad[5] = D3*R32*R1 +     D2*R1;
            grad[6] = D3*R22*R2 + 3.0*D2*R2;
            grad[7] = D3*R22*R3 +     D2*R3;
            grad[8] = D3*R32*R2 +     D2*R2;
            grad[9] = D3*R32*R3 + 3.0*D2*R3;
      }
      else if (p == 4){
            normR5 = normR*normR*normR*normR*normR;
            normR7 = normR5*normR*normR;
            normR9 = normR7*normR*normR;
            D2 = 3.0*G/normR5;
            D3 = -15.0*G/normR7;
            D4 = 105.0*G/normR9;
            R12 = R1*R1;   R22 = R2*R2;   R32 = R3*R3;
            R13 = R12*R1;  R23 = R22*R2;  R33 = R32*R3;
            R1_R2 = R1*R2; R1_R3 = R1*R3; R2_R3 = R2*R3;
            grad[0]  = D4*R12*R12 + 6.0*D3*R12 +   3.0*D2;
            grad[1]  = D4*R13*R2  + 3.0*D3*R1_R2;
            grad[2]  = D4*R13*R3  + 3.0*D3*R1_R3;
            grad[3]  = D4*R12*R22 +     D3*(R12+R22) + D2;
            grad[4]  = D4*R12*R2_R3 +   D3*R2_R3;
            grad[5]  = D4*R12*R32 +     D3*(R12+R32) + D2;
            grad[6]  = D4*R23*R1  + 3.0*D3*R1_R2;
            grad[7]  = D4*R22*R1_R3 +   D3*R1_R3;
            grad[8]  = D4*R1_R2*R32 +   D3*R1_R2;
            grad[9]  = D4*R33*R1  + 3.0*D3*R1_R3;
            grad[10] = D4*R22*R22 + 6.0*D3*R22 +   3.0*D2;
            grad[11] = D4*R23*R3  + 3.0*D3*R2_R3;
            grad[12] = D4*R22*R32 +     D3*(R22+R32) + D2;
            grad[13] = D4*R33*R2  + 3.0*D3*R2_R3;
            grad[14] = D4*R32*R32 + 6.0*D3*R32 +   3.0*D2;
      }
      else if (p == 5){
            normR7 = normR*normR*normR*normR*normR*normR*normR;
            normR9 = normR7*normR*normR;
            normR11 = normR9*normR*normR;
            D3 = -15.0*G/normR7;
            D4 = 105.0*G/normR9;
            D5 = -945.0*G/normR11;
            R12 = R1*R1;    R22 = R2*R2;    R32 = R3*R3;
            R13 = R12*R1;   R23 = R22*R2;   R33 = R32*R3;
            R14 = R12*R12;  R24 = R22*R22;  R34 = R32*R32;
            R1_R2 = R1*R2; R1_R3 = R1*R3; R2_R3 = R2*R3;
            grad[0]  = D5*R13*R12 +  10.0*D4*R13 +             15.0*D3*R1;
            grad[1]  = D5*R14*R2  +   6.0*D4*R12*R2 +           3.0*D3*R2;
            grad[2]  = D5*R14*R3  +   6.0*D4*R12*R3 +           3.0*D3*R3;
            grad[3]  = D5*R13*R22 +       D4*R1*(R12+3.0*R22) + 3.0*D3*R1;
            grad[4]  = D5*R13*R2_R3 + 3.0*D4*R1_R2*R3;
            grad[5]  = D5*R13*R32 +       D4*R1*(R12+3.0*R32) + 3.0*D3*R1;
            grad[6]  = D5*R23*R12 +       D4*R2*(R22+3.0*R12) + 3.0*D3*R2;
            grad[7]  = D5*R12*R22*R3 +    D4*R3*(R12+R22) +         D3*R3;
            grad[8]  = D5*R12*R2*R32 +    D4*R2*(R12+R32) +         D3*R2;
            grad[9]  = D5*R33*R12 +       D4*R3*(R32+3.0*R12) + 3.0*D3*R3;
            grad[10] = D5*R24*R1  +   6.0*D4*R22*R1 +           3.0*D3*R1;
            grad[11] = D5*R23*R1_R3 + 3.0*D4*R1_R2*R3;
            grad[12] = D5*R1*R22*R32 +    D4*R1*(R32+R22) +         D3*R1;
            grad[13] = D5*R1_R2*R33 + 3.0*D4*R1_R2*R3;
            grad[14] = D5*R34*R1  +   6.0*D4*R32*R1 +           3.0*D3*R1;
            grad[15] = D5*R23*R22 +  10.0*D4*R23 +             15.0*D3*R2;
            grad[16] = D5*R24*R3  +   6.0*D4*R22*R3 +           3.0*D3*R3;
            grad[17] = D5*R23*R32 +       D4*R2*(R22+3.0*R32) + 3.0*D3*R2;
            grad[18] = D5*R33*R22 +       D4*R3*(R32+3.0*R22) + 3.0*D3*R3;
            grad[19] = D5*R34*R2  +   6.0*D4*R32*R2 +           3.0*D3*R2;
            grad[20] = D5*R33*R32 +  10.0*D4*R33 +             15.0*D3*R3;
      }
      else if (p == 6){
            normR7 = normR*normR*normR*normR*normR*normR*normR;
            normR9 = normR7*normR*normR;
            normR11 = normR9*normR*normR;
            normR13 = normR11*normR*normR;
            D3 = -15.0*G/normR7;
            D4 = 105.0*G/normR9;
            D5 = -945.0*G/normR11;
            D6 = 10395.0*G/normR13;
            R12 = R1*R1;    R22 = R2*R2;    R32 = R3*R3;
            R13 = R12*R1;   R23 = R22*R2;   R33 = R32*R3;
            R14 = R12*R12;  R24 = R22*R22;  R34 = R32*R32;
            R15 = R13*R12;  R25 = R23*R22;  R35 = R33*R32;
            R1_R2 = R1*R2; R1_R3 = R1*R3; R2_R3 = R2*R3;
            grad[0]  = D6*R13*R13 +  15.0*D5*R14 +                  45.0*D4*R12 +          15.0*D3;
            grad[1]  = D6*R15*R2 +   10.0*D5*R13*R2 +               15.0*D4*R1_R2;
            grad[2]  = D6*R15*R3 +   10.0*D5*R13*R3 +               15.0*D4*R1_R3;
            grad[3]  = D6*R14*R22 +       D5*R12*(R12+6.0*R22) +     3.0*D4*(R22+2.0*R12) + 3.0*D3;
            grad[4]  = D6*R14*R2_R3 + 6.0*D5*R12*R2_R3 +             3.0*D4*R2_R3;
            grad[5]  = D6*R14*R32 +       D5*R12*(R12+6.0*R32) +     3.0*D4*(R32+2.0*R12) + 3.0*D3;
            grad[6]  = D6*R13*R23 +   3.0*D5*R1_R2*(R12+R22) +       9.0*D4*R1_R2;
            grad[7]  = D6*R13*R22*R3 +    D5*R1_R3*(R12+3.0*R22) +   3.0*D4*R1_R3;
            grad[8]  = D6*R13*R32*R2 +    D5*R1_R2*(R12+3.0*R32) +   3.0*D4*R1_R2;
            grad[9]  = D6*R13*R33 +   3.0*D5*R1_R3*(R12+R32) +       9.0*D4*R1_R3;
            grad[10] = D6*R24*R12 +       D5*R22*(R22+6.0*R12) +     3.0*D4*(R12+2.0*R22) + 3.0*D3;
            grad[11] = D6*R23*R12*R3 +    D5*R2_R3*(R22+3.0*R12) +   3.0*D4*R2_R3;
            grad[12] = D6*R12*R22*R32 +   D5*(R12*R22+R12*R32+R22*R32) + D4*(R12+R22+R32) +     D3;
            grad[13] = D6*R33*R12*R2 +    D5*R2_R3*(R32+3.0*R12) +   3.0*D4*R2_R3;
            grad[14] = D6*R34*R12 +       D5*R32*(R32+6.0*R12) +     3.0*D4*(R12+2.0*R32) + 3.0*D3;
            grad[15] = D6*R25*R1 +   10.0*D5*R23*R1 +               15.0*D4*R1_R2;
            grad[16] = D6*R24*R1_R3 + 6.0*D5*R22*R1_R3 +             3.0*D4*R1_R3; 
            grad[17] = D6*R23*R32*R1 +    D5*R1_R2*(R22+3.0*R32) +   3.0*D4*R1_R2;
            grad[18] = D6*R33*R22*R1 +    D5*R1_R3*(R32+3.0*R22) +   3.0*D4*R1_R3;
            grad[19] = D6*R34*R1_R2 + 6.0*D5*R32*R1_R2 +             3.0*D4*R1_R2;
            grad[20] = D6*R35*R1 +   10.0*D5*R33*R1 +               15.0*D4*R1_R3;
            grad[21] = D6*R23*R23 +  15.0*D5*R24 +                  45.0*D4*R22 +          15.0*D3;
            grad[22] = D6*R25*R3 +   10.0*D5*R23*R3 +               15.0*D4*R2_R3;
            grad[23] = D6*R24*R32 +       D5*R22*(R22+6.0*R32) +     3.0*D4*(R32+2.0*R22) + 3.0*D3;
            grad[24] = D6*R33*R23 +   3.0*D5*R2_R3*(R32+R22) +       9.0*D4*R2_R3;
            grad[25] = D6*R34*R22 +       D5*R32*(R32+6.0*R22) +     3.0*D4*(R22+2.0*R32) + 3.0*D3;
            grad[26] = D6*R35*R2 +   10.0*D5*R33*R2 +               15.0*D4*R2_R3;
            grad[27] = D6*R33*R33 +  15.0*D5*R34 +                  45.0*D4*R32 +          15.0*D3;
      }
}


void inner_product(typ * T1, typ * T2, typ * T3, int p, int q, typ factor){

      /******** Computes the inner product between the tensor T1 of order p+q (max 6) and the tensor T2          ********/
      /******** of order p. The result is accumulated into the tensor T3 of order q with a factor factor.        ********/
      /******** T1, T2 and T3 are assumed to be arrays of sizes (p+q+1)(p+q+2)/2, (p+1)(p+2)/2 and (q+1)(q+2)/2. ********/
      /******** The tensor T3 is not overwritten                                                                 ********/

      int q1, q2, q3; //Indexes into the arrays T1, T2 and T3
      int T2_size = ((p+1)*(p+2))/2;
      int T3_size = ((q+1)*(q+2))/2;
      int perm; //Number of tensor index permutations
      
      for (q3 = 0; q3 < T3_size; q3++){ //Looping over the components of the tensor to be computed
            for (q2 = 0; q2 < T2_size; q2++){  //Looping over the components of T2, in order to accumulate T3[q3]
                  perm = perm_from_kn[q2][p];  //number of permutations that can be done with the indexes of T2
                  q1 = q1fromq2q3[q2][q3];     //Retrieving the array index of T1
                  T3[q3] += T1[q1] * T2[q2] * (typ) perm * factor; //Accumulating into tensor T3 with a factor factor
            }
      }
}


/***************************************************************************/
/******** We now write functions relative two the second stage of   ********/
/******** Dehnen's algorithm, accumulation of the C^(m) in the tree ********/
/***************************************************************************/


void Cm_flattree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Performs the tree walk of Dehnen's algorithm, called the interaction phase      ********/
      /******** This function computes the tensors C^(m) between interacting cells, or the      ********/
      /******** acceleration C^(1) in case of interaction between moonlets, but does not pass   ********/
      /******** the C^(m) down the tree. It is assumed that the array C1Moonlets, whose indexes ********/
      /******** 3*a to 3*a+2 contain the acceleration of moonlet a, is given initialized to 0.0 ********/
      /******** Similarly, it is assumed that the C^(m) are 0.0 upon calling this function      ********/
      
      
      /******** It is hard to tell in advance how many pairs will be treated by the tree walk    ********/
      /******** We expect it will be at most factor * cell_id, but that might have to be changed ********/
      int factor = (int) integral(250.0 * 0.5 / theta_min);
      struct pair * stack = (struct pair *)malloc(factor * cell_id * sizeof(struct pair)); //Stack of pairs of ids of nodes that have to be treated
      if (stack == NULL){
            fprintf(stderr, "Error : Cannot allocate memory for the stack in function Cm_flattree.\n");
            abort();
      }
      
      int i = 0; //Index in the stack of the current pair of ids
      int j = 0; //Index of where to put a pair in the stack
      int a, b;  //Ids of the nodes of the current pair
      int p, q;  //Loop indexes
      int s, u;  //Moonlet indexes
      int * dots_a; //All the moonlets in cell a
      int * dots_b; //All the moonlets in cell b
      int Na, Nb; //Number of moonlets in cells a and b
      int how_many_children_a, how_many_children_b; //Number of children of cells a and b
      int idFirstChild; //Id of first child
      int idLastChild;  //Id of last child
      typ M0_a, M0_b;   //Nodes' masses
      typ * com_a, * com_b; //Nodes' centers of mass
      typ * M2_a, * M3_a, * M4_a, * M5_a; //Node a's multipole moments
      typ * M2_b, * M3_b, * M4_b, * M5_b; //Node b's multipole moments
      typ R[3];  //com_a - com_b
      typ grad[3], grad2[6], grad3[10], grad4[15], grad5[21], grad6[28]; //p^th gradient of G/R
      typ C1_a[3], C2_a[6], C3_a[10], C4_a[15], C5_a[21], C6_a[28]; //p^th order interaction tensor of node a
      typ C1_b[3], C2_b[6], C3_b[10], C4_b[15], C5_b[21], C6_b[28]; //p^th order interaction tensor of node b
      typ r_crit_a, r_crit_b; // Critical radii of node a and b
      typ ms, Rs, xs, ys, zs, mu, Ru, xu, yu, zu, dx, dy, dz, r, r3, softening; //moonlet coordinates
      typ omega2_x, omega2_y, omega2_z;
      
      /******** Putting the pair (root_cell, root_cell) in the stack ********/
      (stack + j) -> fst = 0;
      (stack + j) -> snd = 0;
      j++;
      
      
      /******** We travel the stack of pairs of nodes. At each pair, if NaNb < N_cc_pre, we treat it brute-forcely,  ********/
      /******** otherwise, if the nodes are well-separated, we treat them by a multipole expansion, otherwise, if    ********/
      /******** NaNb < N_cc_post or the pair has no children, we treat it brute-forcely, otherwise, we subdivise the ********/
      /******** largest node of the pair (or the only one that has children). If it is a pair of the same node, we   ********/
      /******** treat it brute-forcely if Na < N_cs or if it has no children, and we subdivise it else               ********/
      while (j > i){
            a = (stack + i) -> fst; //Id of first  node
            b = (stack + i) -> snd; //Id of second node
            Na = (FlatTree + a) -> how_many_dots;
            Nb = (FlatTree + b) -> how_many_dots;
            if (a == b){ //Cell self-interation
                  how_many_children_a = (FlatTree + a) -> how_many_children;
                  if (Na < N_cs || how_many_children_a == 0){ //Direct interaction
                        dots_a = (FlatTree + a) -> dots;
                        for (p = 0; p < Na; p++){
                              for (q = p + 1; q < Na; q++){
                                    s = dots_a[p]; //Id of first  moonlet
                                    u = dots_a[q]; //Id of second monnlet
                                    ms = (moonlets + s) -> mass;
                                    Rs = (moonlets + s) -> radius;
                                    xs = (moonlets + s) -> x;
                                    ys = (moonlets + s) -> y;
                                    zs = (moonlets + s) -> z;
                                    mu = (moonlets + u) -> mass;
                                    Ru = (moonlets + u) -> radius;
                                    xu = (moonlets + u) -> x;
                                    yu = (moonlets + u) -> y;
                                    zu = (moonlets + u) -> z;
                                    dx = xs - xu;
                                    dy = ys - yu;
                                    dz = zs - zu;
                                    softening = softening_parameter*(Rs + Ru);
                                    r  = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                                    r3 = r*r*r;
                                    omega2_x = G*dx/r3;
                                    omega2_y = G*dy/r3;
                                    omega2_z = G*dz/r3;
                                    C1Moonlets[3*s]   -= mu*omega2_x;
                                    C1Moonlets[3*s+1] -= mu*omega2_y;
                                    C1Moonlets[3*s+2] -= mu*omega2_z;
                                    C1Moonlets[3*u]   += ms*omega2_x;
                                    C1Moonlets[3*u+1] += ms*omega2_y;
                                    C1Moonlets[3*u+2] += ms*omega2_z;
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
                        for (p = idFirstChild; p < idLastChild; p++){
                              for (q = p; q < idLastChild; q++){
                                    (stack + j) -> fst = p;
                                    (stack + j) -> snd = q;
                                    j++;
                              }
                        }
                  }
            }
            else{ //Interaction between two different cells
                  if (!(Na > N_cc_pre && Nb > N_cc_pre) && Na*Nb < N_cc_pre){ //Direct interaction
                        dots_a = (FlatTree + a) -> dots;
                        dots_b = (FlatTree + b) -> dots;
                        for (p = 0; p < Na; p++){
                              for (q = 0; q < Nb; q++){
                                    s = dots_a[p]; //Id of first  moonlet
                                    u = dots_b[q]; //Id of second moonlet
                                    ms = (moonlets + s) -> mass;
                                    Rs = (moonlets + s) -> radius;
                                    xs = (moonlets + s) -> x;
                                    ys = (moonlets + s) -> y;
                                    zs = (moonlets + s) -> z;
                                    mu = (moonlets + u) -> mass;
                                    Ru = (moonlets + u) -> radius;
                                    xu = (moonlets + u) -> x;
                                    yu = (moonlets + u) -> y;
                                    zu = (moonlets + u) -> z;
                                    dx = xs - xu;
                                    dy = ys - yu;
                                    dz = zs - zu;
                                    softening = softening_parameter*(Rs + Ru);
                                    r  = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                                    r3 = r*r*r;
                                    omega2_x = G*dx/r3;
                                    omega2_y = G*dy/r3;
                                    omega2_z = G*dz/r3;
                                    C1Moonlets[3*s]   -= mu*omega2_x;
                                    C1Moonlets[3*s+1] -= mu*omega2_y;
                                    C1Moonlets[3*s+2] -= mu*omega2_z;
                                    C1Moonlets[3*u]   += ms*omega2_x;
                                    C1Moonlets[3*u+1] += ms*omega2_y;
                                    C1Moonlets[3*u+2] += ms*omega2_z;
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
                              gradR(R, grad, 1);  //Computing the first gradient of G/R
                              C1_a[0] =  M0_b * grad[0];  C1_a[1] =  M0_b * grad[1];  C1_a[2] =  M0_b * grad[2];
                              C1_b[0] = -M0_a * grad[0];  C1_b[1] = -M0_a * grad[1];  C1_b[2] = -M0_a * grad[2];
                              
                              if (expansion_order >= 2){
                                    gradR(R, grad2, 2);  //Computing the second gradient
                                    for (p = 0; p < 6; p++){
                                          C2_a[p] = M0_b * grad2[p];
                                          C2_b[p] = M0_a * grad2[p];
                                    }
                                    
                                    if (expansion_order >= 3){
                                          M2_a = M2FlatTree + 6*a;
                                          M2_b = M2FlatTree + 6*b;
                                          gradR(R, grad3, 3);  //Computing the third gradient
                                          for (p = 0; p < 10; p++){
                                                C3_a[p] =  M0_b * grad3[p];
                                                C3_b[p] = -M0_a * grad3[p];
                                          }
                                          C1_a[0] += 0.5*(M2_b[0]*grad3[0]+2.0*M2_b[1]*grad3[1]+2.0*M2_b[2]*grad3[2]+M2_b[3]*grad3[3]+2.0*M2_b[4]*grad3[4]+M2_b[5]*grad3[5]);
                                          C1_a[1] += 0.5*(M2_b[0]*grad3[1]+2.0*M2_b[1]*grad3[3]+2.0*M2_b[2]*grad3[4]+M2_b[3]*grad3[6]+2.0*M2_b[4]*grad3[7]+M2_b[5]*grad3[8]);
                                          C1_a[2] += 0.5*(M2_b[0]*grad3[2]+2.0*M2_b[1]*grad3[4]+2.0*M2_b[2]*grad3[5]+M2_b[3]*grad3[7]+2.0*M2_b[4]*grad3[8]+M2_b[5]*grad3[9]);
                                          C1_b[0] -= 0.5*(M2_a[0]*grad3[0]+2.0*M2_a[1]*grad3[1]+2.0*M2_a[2]*grad3[2]+M2_a[3]*grad3[3]+2.0*M2_a[4]*grad3[4]+M2_a[5]*grad3[5]);
                                          C1_b[1] -= 0.5*(M2_a[0]*grad3[1]+2.0*M2_a[1]*grad3[3]+2.0*M2_a[2]*grad3[4]+M2_a[3]*grad3[6]+2.0*M2_a[4]*grad3[7]+M2_a[5]*grad3[8]);
                                          C1_b[2] -= 0.5*(M2_a[0]*grad3[2]+2.0*M2_a[1]*grad3[4]+2.0*M2_a[2]*grad3[5]+M2_a[3]*grad3[7]+2.0*M2_a[4]*grad3[8]+M2_a[5]*grad3[9]);
                                          
                                          if (expansion_order >= 4){
                                                M3_a = M3FlatTree + 10*a;
                                                M3_b = M3FlatTree + 10*b;
                                                gradR(R, grad4, 4);  //Computing the fourth gradient
                                                for (p = 0; p < 15; p++){
                                                      C4_a[p] = M0_b * grad4[p];
                                                      C4_b[p] = M0_a * grad4[p];
                                                }
                                                inner_product(grad4, M3_b, C1_a, 3, 1, -0.166666666666666666666666667);
                                                inner_product(grad4, M3_a, C1_b, 3, 1, -0.166666666666666666666666667);
                                                inner_product(grad4, M2_b, C2_a, 2, 2, 0.5);
                                                inner_product(grad4, M2_a, C2_b, 2, 2, 0.5);
                                                
                                                if (expansion_order >= 5){
                                                      M4_a = M4FlatTree + 15*a;
                                                      M4_b = M4FlatTree + 15*b;
                                                      gradR(R, grad5, 5);  //Computing the fifth gradient
                                                      for (p = 0; p < 21; p++){
                                                            C5_a[p] =  M0_b * grad5[p];
                                                            C5_b[p] = -M0_a * grad5[p];
                                                      }
                                                      inner_product(grad5, M4_b, C1_a, 4, 1,  0.041666666666666666666666667);
                                                      inner_product(grad5, M4_a, C1_b, 4, 1, -0.041666666666666666666666667);
                                                      inner_product(grad5, M3_b, C2_a, 3, 2, -0.166666666666666666666666667);
                                                      inner_product(grad5, M3_a, C2_b, 3, 2,  0.166666666666666666666666667);
                                                      inner_product(grad5, M2_b, C3_a, 2, 3,  0.5);
                                                      inner_product(grad5, M2_a, C3_b, 2, 3, -0.5);
                                                      
                                                      if (expansion_order >= 6){
                                                            M5_a = M5FlatTree + 21*a;
                                                            M5_b = M5FlatTree + 21*b;
                                                            gradR(R, grad6, 6);  //Computing the sixth gradient
                                                            for (p = 0; p < 28; p++){
                                                                  C6_a[p] = M0_b * grad6[p];
                                                                  C6_b[p] = M0_a * grad6[p];
                                                            }
                                                            inner_product(grad6, M5_b, C1_a, 5, 1, -0.008333333333333333333333333);
                                                            inner_product(grad6, M5_a, C1_b, 5, 1, -0.008333333333333333333333333);
                                                            inner_product(grad6, M4_b, C2_a, 4, 2,  0.041666666666666666666666667);
                                                            inner_product(grad6, M4_a, C2_b, 4, 2,  0.041666666666666666666666667);
                                                            inner_product(grad6, M3_b, C3_a, 3, 3, -0.166666666666666666666666667);
                                                            inner_product(grad6, M3_a, C3_b, 3, 3, -0.166666666666666666666666667);
                                                            inner_product(grad6, M2_b, C4_a, 2, 4,  0.5);
                                                            inner_product(grad6, M2_a, C4_b, 2, 4,  0.5);
                                                      }
                                                }
                                          }
                                    }
                              }
                              /******** Actualizing the FlatTrees C^(m) ********/
                              ((FlatTree + a) -> C1)[0] += C1_a[0];  ((FlatTree + a) -> C1)[1] += C1_a[1];  ((FlatTree + a) -> C1)[2] += C1_a[2];
                              ((FlatTree + b) -> C1)[0] += C1_b[0];  ((FlatTree + b) -> C1)[1] += C1_b[1];  ((FlatTree + b) -> C1)[2] += C1_b[2];
                              if(expansion_order >= 2){
                                    for (p = 0; p < 6; p++){
                                          C2FlatTree[6*a+p] += C2_a[p];  C2FlatTree[6*b+p] += C2_b[p];
                                    }
                                    if(expansion_order >= 3){
                                          for (p = 0; p < 10; p++){
                                                C3FlatTree[10*a+p] += C3_a[p];  C3FlatTree[10*b+p] += C3_b[p];
                                          }
                                          if(expansion_order >= 4){
                                                for (p = 0; p < 15; p++){
                                                      C4FlatTree[15*a+p] += C4_a[p];  C4FlatTree[15*b+p] += C4_b[p];
                                                }
                                                if(expansion_order >= 5){
                                                      for (p = 0; p < 21; p++){
                                                            C5FlatTree[21*a+p] += C5_a[p];  C5FlatTree[21*b+p] += C5_b[p];
                                                      }
                                                      if(expansion_order >= 6){
                                                            for (p = 0; p < 28; p++){
                                                                  C6FlatTree[28*a+p] += C6_a[p];  C6FlatTree[28*b+p] += C6_b[p];
                                                            }
                                                      }
                                                }
                                          }
                                    }
                              }
                        }
                        else{ //The nodes are not well-separated
                              how_many_children_a = (FlatTree + a) -> how_many_children;
                              how_many_children_b = (FlatTree + b) -> how_many_children;
                              if ((!(Na > N_cc_post && Nb > N_cc_post) && Na*Nb < N_cc_post) || (how_many_children_a == 0 && how_many_children_b == 0)){ //Direct interaction
                                    dots_a = (FlatTree + a) -> dots;
                                    dots_b = (FlatTree + b) -> dots;
                                    for (p = 0; p < Na; p++){
                                          for (q = 0; q < Nb; q++){
                                                s = dots_a[p]; //Id of first  moonlet
                                                u = dots_b[q]; //Id of second moonlet
                                                ms = (moonlets + s) -> mass;
                                                Rs = (moonlets + s) -> radius;
                                                xs = (moonlets + s) -> x;
                                                ys = (moonlets + s) -> y;
                                                zs = (moonlets + s) -> z;
                                                mu = (moonlets + u) -> mass;
                                                Ru = (moonlets + u) -> radius;
                                                xu = (moonlets + u) -> x;
                                                yu = (moonlets + u) -> y;
                                                zu = (moonlets + u) -> z;
                                                dx = xs - xu;
                                                dy = ys - yu;
                                                dz = zs - zu;
                                                softening = softening_parameter*(Rs + Ru);
                                                r  = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
                                                r3 = r*r*r;
                                                omega2_x = G*dx/r3;
                                                omega2_y = G*dy/r3;
                                                omega2_z = G*dz/r3;
                                                C1Moonlets[3*s]   -= mu*omega2_x;
                                                C1Moonlets[3*s+1] -= mu*omega2_y;
                                                C1Moonlets[3*s+2] -= mu*omega2_z;
                                                C1Moonlets[3*u]   += ms*omega2_x;
                                                C1Moonlets[3*u+1] += ms*omega2_y;
                                                C1Moonlets[3*u+2] += ms*omega2_z;
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
                                          for (p = idFirstChild; p < idLastChild; p++){
                                                (stack + j) -> fst = a;
                                                (stack + j) -> snd = p;
                                                j++;
                                          }
                                    }
                                    else if (how_many_children_b == 0){ //Subdivising a
                                          idFirstChild = (FlatTree + a) -> idFirstChild;
                                          idLastChild  = idFirstChild + how_many_children_a;
                                          for (p = idFirstChild; p < idLastChild; p++){
                                                (stack + j) -> fst = b;
                                                (stack + j) -> snd = p;
                                                j++;
                                          }
                                    }
                                    else{ //Both cells have children, subdivising the largest (in term of critical radius)
                                          if (r_crit_a < r_crit_b){ //Subdivising b
                                                idFirstChild = (FlatTree + b) -> idFirstChild;
                                                idLastChild  = idFirstChild + how_many_children_b;
                                                for (p = idFirstChild; p < idLastChild; p++){
                                                      (stack + j) -> fst = a;
                                                      (stack + j) -> snd = p;
                                                      j++;
                                                }
                                          }
                                          else{ //Subdivising a
                                                idFirstChild = (FlatTree + a) -> idFirstChild;
                                                idLastChild  = idFirstChild + how_many_children_a;
                                                for (p = idFirstChild; p < idLastChild; p++){
                                                      (stack + j) -> fst = b;
                                                      (stack + j) -> snd = p;
                                                      j++;
                                                }
                                          }
                                    }
                              }
                        }
                  }
            }
            i++;
      }
      free(stack);
      stack = NULL;
}


/************************************************************************/
/******** We now write functions relative two the third stage of ********/
/******** Dehnen's algorithm, passing the C^(n) down the tree    ********/
/************************************************************************/


void Cm_downtree(struct node * FlatTree, struct moonlet * moonlets){

      /******** Performs the third stage of Dehnen's algorithm, called the evaluation phase, ********/
      /******** where the C^(m) are passed down the three and the accelerations C^(1) of     ********/
      /******** each moonlet computed. It is assumed that this function follows Cm_flattree  ********/
      
      int * stack = (int *)malloc(cell_id * sizeof(int)); //Stack of pairs of ids of nodes that have to be treated
      int i = 0; //Index in the stack of the current id
      int j = 0; //Index of where to put an id in the stack
      int a;     //Id of the current node
      int p, q;     //Loop indexes
      int b;     //Moonlet index
      int * dots;//All the moonlets in the current cell
      int Na;    //Number of moonlets in cells a and b
      int how_many_children; //Number of children of the current node
      int idFirstChild; //Id of first child
      int idLastChild;  //Id of last child
      int idParent;
      typ M0;     //Node' mass
      typ * com, * com_parent; //centers of mass of the cell and of its parent
      typ R[3];  //com - com_parent
      typ X2[6], X3[10], X4[15], X5[21];
      typ * C1, * C2, * C3, * C4, * C5, * C6; //p^th order interaction tensor of the current node
      typ * C1p, * C2p, * C3p, * C4p, * C5p, * C6p; //p^th order interaction tensor of the parent of the current node

      /******** Putting the root cell in the stack ********/
      stack[j] = 0;
      j++;

      while (j > i){
            a = stack[i];

            /******** If the node has a parent, we shift and accumulate its C^(n) ********/
            idParent = (FlatTree + a) -> idParent;
            if (idParent != -1){
                  C1 =  (FlatTree + a)        -> C1;
                  C1p = (FlatTree + idParent) -> C1;
                  C1[0] += C1p[0]; C1[1] += C1p[1]; C1[2] += C1p[2];
                  if (expansion_order >= 2){
                        com = (FlatTree + a) -> com;
                        com_parent = (FlatTree + idParent) -> com;
                        R[0] = com[0] - com_parent[0];
                        R[1] = com[1] - com_parent[1];
                        R[2] = com[2] - com_parent[2];
                        C2  = C2FlatTree + 6*a;
                        C2p = C2FlatTree + 6*idParent;
                        C2[0] += C2p[0]; C2[1] += C2p[1]; C2[2] += C2p[2]; C2[3] += C2p[3]; C2[4] += C2p[4]; C2[5] += C2p[5];
                        C1[0] += R[0]*C2p[0] + R[1]*C2p[1] + R[2]*C2p[2];
                        C1[1] += R[0]*C2p[1] + R[1]*C2p[3] + R[2]*C2p[4];
                        C1[2] += R[0]*C2p[2] + R[1]*C2p[4] + R[2]*C2p[5];

                        if (expansion_order >= 3){
                              C3  = C3FlatTree + 10*a;
                              C3p = C3FlatTree + 10*idParent;
                              for (p = 0; p < 10; p++){
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

                              if (expansion_order >= 4){
                                    C4  = C4FlatTree + 15*a;
                                    C4p = C4FlatTree + 15*idParent;
                                    for (p = 0; p < 15; p++){
                                          C4[p] += C4p[p];
                                    }
                                    for (p = 0; p < 10; p++){
                                          get_Xn_overwrite(p, 3, R, 1.0, X3);
                                    }
                                    inner_product(C4p, X3, C1, 3, 1, 0.1666666666666666666666666667);
                                    inner_product(C4p, X2, C2, 2, 2, 0.5);
                                    inner_product(C4p,  R, C3, 1, 3, 1.0);

                                    if (expansion_order >= 5){
                                          C5  = C5FlatTree + 21*a;
                                          C5p = C5FlatTree + 21*idParent;
                                          for (p = 0; p < 21; p++){
                                                C5[p] += C5p[p];
                                          }
                                          for (p = 0; p < 15; p++){
                                                get_Xn_overwrite(p, 4, R, 1.0, X4);
                                          }
                                          inner_product(C5p, X4, C1, 4, 1, 0.0416666666666666666666666667);
                                          inner_product(C5p, X3, C2, 3, 2, 0.1666666666666666666666666667);
                                          inner_product(C5p, X2, C3, 2, 3, 0.5);
                                          inner_product(C5p,  R, C4, 1, 4, 1.0);

                                          if (expansion_order >= 6){
                                                C6  = C6FlatTree + 28*a;
                                                C6p = C6FlatTree + 28*idParent;
                                                for (p = 0; p < 28; p++){
                                                      C6[p] += C6p[p];
                                                }
                                                for (p = 0; p < 21; p++){
                                                      get_Xn_overwrite(p, 5, R, 1.0, X5);
                                                }
                                                inner_product(C6p, X5, C1, 5, 1, 0.0083333333333333333333333333);
                                                inner_product(C6p, X4, C2, 4, 2, 0.0416666666666666666666666667);
                                                inner_product(C6p, X3, C3, 3, 3, 0.1666666666666666666666666667);
                                                inner_product(C6p, X2, C4, 2, 4, 0.5);
                                                inner_product(C6p,  R, C5, 1, 5, 1.0);
                                          }
                                    }
                              }
                        }
                  }
            }
            
            /******** If the node has children, we put them in the stack ********/
            how_many_children = (FlatTree + a) -> how_many_children;
            if (how_many_children > 0){
                  idFirstChild = (FlatTree + a) -> idFirstChild;
                  idLastChild  = idFirstChild + how_many_children;
                  for (p = idFirstChild; p < idLastChild; p++){
                        stack[j] = p;
                        j++;
                  }
            }
            
            /******** If the node has no children, we shift and accumulate its C^(n) to its moonlets ********/
            else{
                  Na = (FlatTree + a) -> how_many_dots;
                  dots = (FlatTree + a) -> dots;
                  com = (FlatTree + a) -> com;
                  for (p = 0; p < Na; p++){
                        b = dots[p];
                        C1 = (FlatTree + a) -> C1;
                        C1Moonlets[3*b] += C1[0];  C1Moonlets[3*b+1] += C1[1];  C1Moonlets[3*b+2] += C1[2];
                        
                        if (expansion_order >= 2){
                              C2 = C2FlatTree + 6*a;
                              R[0] = (moonlets + b) -> x - com[0];
                              R[1] = (moonlets + b) -> y - com[1];
                              R[2] = (moonlets + b) -> z - com[2];
                              C1Moonlets[3*b]   += R[0]*C2[0] + R[1]*C2[1] + R[2]*C2[2];
                              C1Moonlets[3*b+1] += R[0]*C2[1] + R[1]*C2[3] + R[2]*C2[4];
                              C1Moonlets[3*b+2] += R[0]*C2[2] + R[1]*C2[4] + R[2]*C2[5];
                              
                              if (expansion_order >= 3){
                                    C3 = C3FlatTree + 10*a;
                                    X2[0]  = R[0]*R[0]; X2[1] = R[0]*R[1]; X2[2] = R[0]*R[2]; X2[3] = R[1]*R[1]; X2[4] = R[1]*R[2]; X2[5] = R[2]*R[2];
                                    C1Moonlets[3*b]   += 0.5*(X2[0]*C3[0] + 2.0*X2[1]*C3[1] + 2.0*X2[2]*C3[2] + X2[3]*C3[3] + 2.0*X2[4]*C3[4] + X2[5]*C3[5]);
                                    C1Moonlets[3*b+1] += 0.5*(X2[0]*C3[1] + 2.0*X2[1]*C3[3] + 2.0*X2[2]*C3[4] + X2[3]*C3[6] + 2.0*X2[4]*C3[7] + X2[5]*C3[8]);
                                    C1Moonlets[3*b+2] += 0.5*(X2[0]*C3[2] + 2.0*X2[1]*C3[4] + 2.0*X2[2]*C3[5] + X2[3]*C3[7] + 2.0*X2[4]*C3[8] + X2[5]*C3[9]);
                                    if (expansion_order >= 4){
                                          C4 = C4FlatTree + 15*a;
                                          for (q = 0; q < 10; q++){
                                                get_Xn_overwrite(q, 3, R, 1.0, X3);
                                          }
                                          inner_product(C4, X3, C1Moonlets + 3*b, 3, 1, 0.1666666666666666666666666667);
                                          
                                          if (expansion_order >= 5){
                                                C5 = C5FlatTree + 21*a;
                                                for (q = 0; q < 15; q++){
                                                      get_Xn_overwrite(q, 4, R, 1.0, X4);
                                                }
                                                inner_product(C5, X4, C1Moonlets + 3*b, 4, 1, 0.0416666666666666666666666667);
                                                
                                                if (expansion_order >= 6){
                                                      C6 = C6FlatTree + 28*a;
                                                      for (q = 0; q < 21; q++){
                                                            get_Xn_overwrite(q, 5, R, 1.0, X5);
                                                      }
                                                      inner_product(C6, X5, C1Moonlets + 3*b, 5, 1, 0.0083333333333333333333333333);
                                                }
                                          }
                                    }
                              }
                        }
                  }
            }
            i++;
      }
      free(stack);
      stack = NULL;
}


/******************************************************************************************/
/******** We now write functions allowing to compute gravity with the tree-code of ********/
/******** Barnes & Hut. The original algorithm of Barnes & Hut performs a first    ********/
/******** order Taylor expansion, but here, the order is given by expansion_order. ********/
/******************************************************************************************/



void standard_tree_acceleration(struct node * FlatTree, struct moonlet * moonlets, int b){

      /******** Computes the acceleration of moonlet b from the multipole moments of the       ********/
      /******** flattree FlatTree. Stores the result into C1Moonlets[3*b] to C1Moonlets[3*b+2] ********/
      /******** This is the standard tree code at expansion order expansion_order              ********/
      
      /******** Retrieving the moonlet's coordinates ********/
      typ X, Y, Z, Rb;
      X  = (moonlets + b) -> x;
      Y  = (moonlets + b) -> y;
      Z  = (moonlets + b) -> z;
      Rb = (moonlets + b) -> radius;
      
      /******** Initializing the acceleration ********/
      typ acc[3] = {0.0, 0.0, 0.0};
      
      typ M0;    //Node's mass
      typ * com; //Node's center of mass
      typ * M2;  //Node's multipole moments
      typ * M3;
      typ * M4;
      typ * M5;
      typ R[3];  //moonlet's position - center of mass
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
      int * dots; //All the moonlets in that cell
      int index; //Index of a moonlet in dots
      typ m, Rad, x, y, z, dx, dy, dz, r, r3, softening;
      
      stack[j] = 0;
      j++;
      
      /******** At each node, if it has at most N_cb_pre moonlets, we treat it directly, otherwise,      ********/
      /******** if it is well-separated, we treat it by multipole expansion, otherwise, if it has at     ********/
      /******** most N_cb_post moonlets or no children, we treat it directly, otherwise, we subdivise it ********/
      while (j > i){
            a = stack[i];
            how_many_dots = (FlatTree + a) -> how_many_dots;
            if (how_many_dots < N_cb_pre){ //Direct summation
                  dots = (FlatTree + a) -> dots;
                  for (index = 0; index < how_many_dots; index++){
                        k = dots[index];
                        if (k != b){ //If the moonlet is different from moonlet b, we accumulate its contribution to the acceleration
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
                        typ grad[3];
                        gradR(R, grad, 1);
                        M0 = (FlatTree + a) -> M0;
                        inner_product(grad, &M0, acc, 0, 1, 1.0);
                        if (expansion_order >= 3){
                              typ grad3[10];
                              gradR(R, grad3, 3);
                              M2 = M2FlatTree + 6*a;
                              inner_product(grad3, M2, acc, 2, 1, 0.5);
                              if (expansion_order >= 4){
                                    typ grad4[15];
                                    gradR(R, grad4, 4);
                                    M3 = M3FlatTree + 10*a;
                                    inner_product(grad4, M3, acc, 3, 1, -0.166666666666666666666666666666666666666);
                                    if (expansion_order >= 5){
                                          typ grad5[21];
                                          gradR(R, grad5, 5);
                                          M4 = M4FlatTree + 15*a;
                                          inner_product(grad5, M4, acc, 4, 1, 0.0416666666666666666666666666666666666666);
                                          if (expansion_order >= 6){
                                                typ grad6[28];
                                                gradR(R, grad6, 6);
                                                M5 = M5FlatTree + 21*a;
                                                inner_product(grad6, M5, acc, 5, 1, -0.008333333333333333333333333333333333333);
                                          }
                                    }
                              }
                        }
                  }
                  else{
                        how_many_children = (FlatTree + a) -> how_many_children;
                        if (how_many_dots < N_cb_post || how_many_children == 0){ //Direct summation
                              dots = (FlatTree + a) -> dots;
                              for (index = 0; index < how_many_dots; index++){
                                    k = dots[index];
                                    if (k != b){ //If the moonlet is different from moonlet b, we accumulate its contribution to the acceleration
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
      for (i=0; i<8; i++){
            *(octants+i) = NULL;
      }
      (*BoxDot)  -> sidelength = D;
      (*BoxDot)  -> how_many = 0;
      ((*BoxDot) -> corner)[0] = *corner_coordinates;
      ((*BoxDot) -> corner)[1] = *(corner_coordinates+1);
      ((*BoxDot) -> corner)[2] = *(corner_coordinates+2);
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
      j++;
      
      struct boxdot * to_be_deleted;
      struct boxdot * octant;
      
      while (j>i){ //While there are still nodes to be deleted
            
            to_be_deleted = stack[i];
            clear_chain(&(to_be_deleted -> dots)); //Clearing the chain of dots
            for (p=0; p<8; p++){ //Adding to the stack the children of the current node
                  octant = (to_be_deleted -> oct)[p];
                  if (octant != NULL){
                        stack[j] = octant;
                        j++;
                  }
            }
            free(to_be_deleted);
            to_be_deleted = NULL;
            i++;
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

      /******** Returns the octant in which a moonlet with coordinates xyz_a   ********/
      /******** should go, given the corner coordinates xyz of the parent cell ********/
      /******** D is half the parent sidelength                                ********/
      
      int octant;
      
      if (xa < x+D){ //Moonlet is in octant 0, 2, 4 or 6
            if (ya < y+D){ //Moonlet is in octant 0 or 2
                  if (za > z-D){ //Moonlet is in octant 0
                        octant = 0;
                  }
                  else{ //Moonlet is in octant 2
                        octant = 2;
                  }
            }
            else{ //Moonlet is in octant 4 or 6
                  if (za > z-D){ //Moonlet is in octant 4
                        octant = 4;
                  }
                  else{ //Moonlet is in octant 6
                        octant = 6;
                  }
            }
      }
      else { //Moonlet is in octant 1, 3, 5 or 7
            if (ya < y+D){ //Moonlet is in octant 1 or 3
                  if (za > z-D){ //Moonlet is in octant 1
                        octant = 1;
                  }
                  else{ //Moonlet is in octant 3
                        octant = 3;
                  }
            }
            else{ //Moonlet is in octant 5 or 7
                  if (za > z-D){ //Moonlet is in octant 5
                        octant = 5;
                  }
                  else{ //Moonlet is in octant 7
                        octant = 7;
                  }
            }
      }
      return octant;
}


void get_corner_coordinates(typ X, typ Y, typ Z, typ D, int i, typ * corner){

      /******** Returns the corner coordinates of the i^th child, given the corner ********/
      /******** coordinates of the parent. D is half the parent's sidelength       ********/
      
      typ x,y,z; //Child's coordinates
      
      /******** Retrieving the child's coordinates ********/
      if (i == 0){
            x=X;  y=Y;  z=Z;
      }
      else if (i==1){
            x=X+D;  y=Y;  z=Z;
      }
      else if (i==2){
            x=X;  y=Y;  z=Z-D;
      }
      else if (i==3){
            x=X+D;  y=Y;  z=Z-D;
      }
      else if (i==4){
            x=X;  y=Y+D;  z=Z;
      }
      else if (i==5){
            x=X+D;  y=Y+D;  z=Z;
      }
      else if (i==6){
            x=X;  y=Y+D;  z=Z-D;
      }
      else{
            x=X+D;  y=Y+D;  z=Z-D;
      }
      
      /******** Returning the child's coordinates ********/
      *corner     = x;
      *(corner+1) = y;
      *(corner+2) = z;
}

void print_boxdot(struct boxdot * BoxDot){

      /******** Prints in terminal the boxdot BoxDot  ********/
      /******** This function is not used for the FFM ********/

      int p;
      typ * corner;
      typ D;
      int octant;
      int rotation;
      typ distance;
      typ previous_center[3] = {0.0, 0.0, 0.0};
      typ center[3] = {0.0, 0.0, 0.0};
      
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
            distance = sqrt((center[0]-previous_center[0])*(center[0]-previous_center[0]) + (center[1]-previous_center[1])*(center[1]-previous_center[1]) + 
            (center[2]-previous_center[2])*(center[2]-previous_center[2]));
            previous_center[0] = center[0];  previous_center[1] = center[1];  previous_center[2] = center[2];
            
            printf("Node n° %d : N = %d, id = %d, level = %d, distance from previous node = %.6lf, center = {x,y,z} = {%.6lf, %.6lf, %.6lf}\n", i, 
            to_be_printed -> how_many, to_be_printed -> id, to_be_printed -> level, distance, center[0], center[1], center[2]);
            
            for (p=0; p<8; p++){ // p is the Hilbert-Peano digit
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




