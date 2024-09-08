/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    spring.h                                                    ********/
/******** @brief   Header file to spring.c                                     ********/
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


#ifndef _SPRING_H_
#define _SPRING_H_

#include "parameters.h"


/******** Defining the connection structure between two particles (Kelvin-Voigt model) ********/
struct connection {
      struct pair Pair; //The pair of particles linked by the connection
      typ rest_length;  //The rest length of the spring
};


/******** Defining a quaternion structure to perform rotations from one vector to another ********/
struct quaternion {
      typ w;
      typ x;
      typ y;
      typ z;
};


extern struct connection * connections; //The array of connections between particles of the viscoelastic body
extern struct chain * first;            //This chain will contain the index of the first  particle of the connection
extern struct chain * second;           //This chain will contain the index of the second particle of the connection
extern int N_connections;               //The total number of connections in the viscoelastic body
extern typ shapeV;                      //The volume of the shape model


void precision(struct moonlet * viscoelastic); //To be removed


struct moonlet * generate_visco_elastic_body();


void generate_connections(struct moonlet * viscoelastic);


void make_rotate(struct moonlet * viscoelastic);


void point_angular_momentum(struct moonlet * viscoelastic);


void overlap(struct moonlet * viscoelastic);


void deOverlap(struct moonlet * viscoelastic, int a, int b);


void viscoelastic_flattree(struct node * FlatTree, struct moonlet * viscoelastic, int generating);


int connects(struct moonlet * viscoelastic, int a, int b);


struct connection make_connection(struct moonlet * viscoelastic, int a, int b);


typ get_perturbing_true_anomaly(typ time);


void get_pert_coordinates(typ time, typ * x, typ * y, typ * z);


void quaternion_norm(struct quaternion * q);


struct quaternion get_quaternion(typ ux, typ uy, typ uz, typ vx, typ vy, typ vz);


void rotate_with_quaternion(typ x, typ y, typ z, struct quaternion q, typ * xr, typ * yr, typ * zr);


void three_closest_nodes(struct moonlet * viscoelastic, int k, int * indexes);


#endif





