#ifndef _PHYSICS_H_
#define _PHYSICS_H_


#include "parameters.h"
#include "structure.h"


void vector_field(struct moonlet * moonlets);


void collision(struct moonlet * moonlets, int a, int b, typ f);


void merger(struct moonlet * moonlets, int a, int b);


void fragmentation(struct moonlet * moonlets, int a, int b);


void collision_treatment(struct moonlet * moonlets, int a, int b, int type_of_collision);


void get_neighbours_mesh(struct moonlet * moonlets);



#endif
