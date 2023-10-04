#ifndef _COLLISION_H_
#define _COLLISION_H_

#include "parameters.h"
#include "structure.h"
#include "ffm.h"


typ * closest_approach(struct moonlet * moonlets, int a, int b);


void hash_table_cell_index(struct moonlet * moonlets, int a);


void hash_table_cell_index_z_order(struct moonlet * moonlets, int a);


void neighbours(struct moonlet * moonlets, int a);


void mesh(struct moonlet * moonlets);


void brute_force(struct moonlet * moonlets);


void get_center_and_maxR(struct node * FlatTree, struct moonlet * moonlets, int a);


void get_rmax_and_rcrit(struct node * FlatTree, struct moonlet * moonlets, int a);


void get_center_and_maxR_from_children(struct node * FlatTree, int a);


void get_rmax_and_rcrit_from_children(struct node * FlatTree, int a);


void center_and_maxR_flattree(struct node * FlatTree, struct moonlet * moonlets);


void rmax_and_rcrit_flattree(struct node * FlatTree, struct moonlet * moonlets);


void collision_flattree(struct node * FlatTree, struct moonlet * moonlets);


void standard_tree_collision(struct node * FlatTree, struct moonlet * moonlets, int b);


#endif
