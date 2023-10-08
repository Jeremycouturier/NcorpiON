#ifndef _RK4_H_
#define _RK4_H_

#include "parameters.h"
#include "structure.h"
#include "ffm.h"


extern struct node * FlatTree;


void kick(struct moonlet * X, void (*F)(struct moonlet *));


void drift(struct moonlet * X);


FILE ** file_opening();


void file_closing();


void display(struct moonlet * moonlets, typ * aei);


void end_of_timestep(struct moonlet * moonlets, int progressed);


int integration(typ t);


int integration_brute_force(typ t);


int integration_tree(typ t);


int integration_mesh(typ t);


int integration_brute_force_SABA1(typ t);


#endif
