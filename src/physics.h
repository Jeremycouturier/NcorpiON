/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/******** @file    physics.h                                                   ********/
/******** @brief   Header file to physics.c                                    ********/
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


#ifndef _PHYSICS_H_
#define _PHYSICS_H_


#include "parameters.h"
#include "structure.h"


void vector_field(struct moonlet * moonlets);


void tides(struct moonlet * X);


void collision(struct moonlet * moonlets, int a, int b, typ f);


void merger(struct moonlet * moonlets, int a, int b);


void fragmentation(struct moonlet * moonlets, int a, int b);


void collision_treatment(struct moonlet * moonlets, int a, int b, int type_of_collision);


void get_neighbours_mesh(struct moonlet * moonlets);



#endif
