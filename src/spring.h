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
      struct pair Pair;       //The pair of particles linked by the connection
      typ rest_length;        //The rest length of the spring
      typ equilibrium_length; //The equilibrium length of the spring. Unknown until the viscous body has been integrated to equilibrium
};


extern struct connection * connections; //The array of connections between particles of the viscoelastic body
extern struct chain * first;            //This chain will contain the index of the first  particle of the connection
extern struct chain * second;           //This chain will contain the index of the second particle of the connection
extern int N_connections;               //The total number of connections in the viscoelastic body


struct moonlet * generate_visco_elastic_body();


void generate_connections(struct moonlet * viscoelastic);


void overlap(struct moonlet * viscoelastic);


void deOverlap(struct moonlet * viscoelastic, int a, int b);


void viscoelastic_flattree(struct node * FlatTree, struct moonlet * viscoelastic, int generating);


int connects(struct moonlet * viscoelastic, int a, int b);


struct connection make_connection(struct moonlet * viscoelastic, int a, int b);


#endif





